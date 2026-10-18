"""后处理：读场、可视化。

分两个层次，别混用：

* :func:`read_field` / :func:`plot` —— **只读已有结果**，只需要文件名，
  不需要 :class:`design.Case`。日常画图用这个。
* :func:`plot_mesh` —— 只读网格文件，按 physical group（label）给不同部分上色，
  用来看网格分布：端口最细、波导内部次之、包层最粗。
* :func:`post_process` / :func:`visualization` —— 需要 ``save_field=1`` 的场，
  而优化循环里是 ``save_field=0``，所以它们会**重跑一次求解器**，因此要 Case。

例::

    import postproc
    postproc.plot("out_cross_waveguide.msh", "result.h5")
    postproc.plot_mesh("cross_waveguide.msh")

    # 或者：从 design.h5 还原算例，重跑一次带场的正问题再画
    import design
    postproc.visualization(design.Case.load())
"""

import numpy as np
import scipy.io as spio

import dataio
import solver
import paths
from design import Case, make_case

DEFAULT_RESULT = "result.h5"


def read_field(result_file=DEFAULT_RESULT):
    """从已经算好的结果里读场，不重跑求解器。"""
    res = dataio.read_result(result_file, want_field=True)
    if "Ez" not in res:
        raise RuntimeError(
            "%s 里没有 /fields/Ez；求解时需要把 infile.pot 的 save_field 设为 1"
            % result_file)
    return res["Ez"]


def plot(mesh_file, result_file=DEFAULT_RESULT, scalars="Ez", Ez=None):
    """把场画出来（只读文件，不跑求解器）。"""
    import pyvista

    if Ez is None:
        Ez = read_field(result_file)
    mesh = pyvista.read(mesh_file)
    mesh["E_norm"] = np.abs(Ez)
    mesh["Ez"] = Ez.real
    return mesh.plot(scalars=scalars)


# ------------------------------------------------------------------ 网格可视化

# 网格里没写 $PhysicalNames 时（例如 libMesh 写回的 out_*.msh）按标准 id 回落
DEFAULT_LABELS = {
    1: "waveguide", 2: "air",
    3: "input", 4: "output1", 5: "output21", 6: "output22",
    7: "ex_opt", 8: "in_opt",
}

# 默认配色：子域浅色填充，端口/优化边界用饱和色描边
DEFAULT_COLORS = {
    "waveguide": "#e8cf7a",
    "air": "#fbf8ea",
    "input": "#d62728",
    "output1": "#1f77b4",
    "output21": "#2ca02c",
    "output22": "#9467bd",
    "ex_opt": "#ff7f0e",
    "in_opt": "#17becf",
}

# 线单元的 VTK 类型（用来判断一组单元是"边界曲线"还是"子域面"）
_LINE_CELL_TYPES = {3, 4, 21}

# 图例顺序：先子域，再端口，最后两条优化边界
LABEL_ORDER = ["waveguide", "air", "input", "output1", "output21", "output22",
               "ex_opt", "in_opt"]


def _label_map(mesh):
    """physical tag -> label。优先用网格里的 ``$PhysicalNames``，没有就用标准 id 表。"""
    labels = dict(DEFAULT_LABELS)
    for name, value in (mesh.field_data or {}).items():
        tag = np.asarray(value).ravel()[0]
        labels[int(tag)] = name
    return labels


def _line_polydata(points, cells):
    """把线单元子集转成 ``pyvista.PolyData`` 的 polyline。

    VTK 的扁平连接格式是 ``[每条线的点数, i0, i1, ..., 点数, ...]``。
    二阶线单元（LINE3）的节点顺序是 ``(v0, v1, mid)``，中间那个节点要插在两端
    之间，曲边的弯曲才画得出来。

    只保留这一组真正用到的节点：否则每个边界 actor 的包围盒都是整个域，
    图例/相机取景会算错。
    """
    import pyvista as pv

    ordered = cells[:, [0, 2, 1]] if cells.shape[1] >= 3 else cells[:, [0, 1]]
    used, remap = np.unique(ordered, return_inverse=True)
    remap = remap.reshape(ordered.shape)
    conn = np.hstack([np.full((len(ordered), 1), ordered.shape[1]), remap])
    return pv.PolyData(np.asarray(points)[used], lines=conn.astype(np.int64).ravel())


def _order_of(name):
    return LABEL_ORDER.index(name) if name in LABEL_ORDER else len(LABEL_ORDER)


def mesh_label_groups(mesh_file, colors=None):
    """按 physical group（label）把网格拆成几块，返回 ``{label: dict}``。

    每条记录的键：

    ``part``
        可直接交给 ``Plotter.add_mesh`` 的对象：子域是 ``UnstructuredGrid``
        （保留二阶三角形，曲边能画出来），边界是 ``PolyData`` 的 polyline。
    ``color``
        颜色（``DEFAULT_COLORS``，可用 ``colors`` 参数覆盖）。
    ``kind``
        ``"subdomain"`` 或 ``"boundary"``，决定怎么画（填充面 + 单元边 vs 粗线）。

    这一步跟渲染无关，所以没有显示环境时也能单独检查分组是否正确。
    """
    import pyvista as pv

    mesh = pv.read(mesh_file)
    if "gmsh:physical" not in mesh.cell_data:
        raise RuntimeError(
            "%s 里没有 gmsh:physical（网格没写 physical group），无法按 label 上色"
            % mesh_file)

    tags = np.abs(np.asarray(mesh.cell_data["gmsh:physical"])).astype(int)
    labels = _label_map(mesh)
    palette = dict(DEFAULT_COLORS)
    palette.update(colors or {})

    # 每个单元的节点索引：cells 是扁平的连接表，offset 是它的分段边界
    conn, offset = np.asarray(mesh.cell_connectivity), np.asarray(mesh.offset)

    groups = {}
    for tag in np.unique(tags):
        idx = np.flatnonzero(tags == tag)
        width = offset[idx + 1] - offset[idx]
        if not np.all(width == width[0]):
            raise RuntimeError("%s 里同一 physical group 混了不同的单元类型" % mesh_file)
        rows = offset[idx][:, None] + np.arange(width[0])[None, :]
        cells = conn[rows]

        types = set(np.unique(mesh.celltypes[idx]).tolist())
        is_boundary = bool(types) and types <= _LINE_CELL_TYPES
        if is_boundary:
            part = _line_polydata(mesh.points, cells)
        else:
            # 这里**不能**用 extract_surface()：它会把二阶单元三角化成 4 个线性
            # 三角形（3320 个 TRI6 → 13280 个 TRI3），单元边就多出对角线，看起来
            # 比实际网格密。直接把单元子集交给渲染器，VTK 原生支持二阶单元。
            part = mesh.extract_cells(idx)

        name = labels.get(int(tag), str(int(tag)))
        groups[name] = {
            "part": part,
            "color": palette.get(name, "#bbbbbb"),
            "kind": "boundary" if is_boundary else "subdomain",
        }
    return dict(sorted(groups.items(), key=lambda kv: _order_of(kv[0])))


def plot_mesh(mesh_file, colors=None, annotate=None, centerlines=False,
              show_edges=True, edge_color="#8f8f8f", edge_width=0.6,
              boundary_width=3.0, opacity=1.0, legend=True, title=None,
              window_size=(900, 800), show=True):
    """按 label 给网格的不同部分上色（pyvista 交互窗口），用来看剖分分布。

    ``label`` 就是 gmsh 的 physical group 名字：``waveguide`` / ``air`` 两个子域
    填充上色，``input`` / ``output1`` / ``output21`` / ``output22`` 四个端口和
    ``ex_opt`` / ``in_opt`` 两条优化边界描边。这样一眼能看出"端口最细、波导内部
    次之、包层最粗"对不对。

    参数
    ----
    mesh_file : str
        网格文件。``cross_waveguide.msh``（gmsh 写，带名字）和
        ``out_cross_waveguide.msh``（libMesh 写，没名字，按标准 id 回落）都行。
    colors : dict, optional
        覆盖默认配色，键是 label。
    annotate : list of (x, y, text), optional
        在图上加点文字标注，例如 ``[(0.25, 0.30, "optimization region")]``。
        注意 VTK 的默认字体没有中文字形，中文可能显示成方框。
    centerlines : bool
        是否画穿过原点的坐标轴线（默认不画）。
    show : bool
        是否立刻 ``show()``（默认是）。``show=False`` 时返回的 Plotter 可以自己
        再 ``show()`` 或叠加其它图元。

    返回 ``pyvista.Plotter``。
    """
    import pyvista as pv

    groups = mesh_label_groups(mesh_file, colors=colors)
    plotter = pv.Plotter(window_size=window_size)
    plotter.set_background("white")

    if centerlines:
        span = 1.1 * max(np.abs(np.asarray(g["part"].points)[:, :2]).max()
                         for g in groups.values() if g["part"].n_points)
        plotter.add_mesh(pv.Line((-span, 0, 0), (span, 0, 0)), color="#d62728",
                         line_width=2, lighting=False)
        plotter.add_mesh(pv.Line((0, -span, 0), (0, span, 0)), color="#d62728",
                         line_width=2, lighting=False)

    for name, group in groups.items():
        if group["kind"] == "subdomain":
            plotter.add_mesh(group["part"], color=group["color"],
                             show_edges=show_edges, edge_color=edge_color,
                             line_width=edge_width, opacity=opacity,
                             lighting=False, label=name)
        else:
            plotter.add_mesh(group["part"], color=group["color"],
                             line_width=boundary_width, lighting=False,
                             label=name)

    if annotate:
        points = np.array([[float(a[0]), float(a[1]), 0.0] for a in annotate])
        plotter.add_point_labels(points, [a[2] for a in annotate],
                                 font_size=14, shape=None, fill_shape=False,
                                 text_color="black", always_visible=True,
                                 show_points=False)

    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.camera.zoom(1.05)
    if legend:
        plotter.add_legend(bcolor="white", border=False, size=(0.22, 0.30))
    if title:
        plotter.add_title(title, font_size=11)
    if show:
        plotter.show()
    return plotter


def post_process(case, out_mat=None):
    """重跑一次正问题（``save_field=1``）并把场导出成 .mat。

    只有这个（和下面的 :func:`visualization`）需要 ``case``，因为要真的
    再解一遍才能拿到场。
    """
    res = solver.run(case, save_field=True)
    if "Ez" not in res:
        raise RuntimeError("结果里没有场；确认求解器成功写出了 /fields")
    spio.savemat(out_mat or case.name + "_post.mat", {"Ez": res["Ez"]})
    return res["Ez"]


def visualization(case, scalars="Ez", out_mat=None):
    """重跑一次并作图（旧 ``visualization()`` 的行为）。"""
    Ez = post_process(case, out_mat=out_mat)
    return plot(case.out_mesh_file, case.result_file, scalars=scalars, Ez=Ez)


def get_case(idx: int =  -1):
    dof_data = np.loadtxt("DOF_list.csv", delimiter=',')
    return make_case(dof_data[idx, :])


if __name__ == "__main__":
    paths.chdir_data()
    cs = get_case(-3)
    print(cs)
    # 网格按 label 上色（快，只读文件、不跑求解器）
    plot_mesh(cs.mesh_file, centerlines=True,
              title="cross_waveguide: coloured by physical group")
    visualization(cs)