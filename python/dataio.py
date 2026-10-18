"""C++ 后端 <-> Python 前端的数据契约。

所有跨语言数据都集中在这里，避免每个脚本各自解释文件名和变量名。

格式约定
--------
``design.h5``  (Python -> C++，几何与参数)

    /design/in_Bezier_P   (n,2) float64   Bezier 控制点
    /design/in_Bezier_V   (n,2) float64   Bezier 切向量
    /design/ex_Bezier_P   (n,2) float64
    /design/ex_Bezier_V   (n,2) float64
    /design/dofs          (n,)  float64   设计变量本身（v2 起）
    /design/node_t        (n_nodes,) float64   每个网格节点在优化边界曲线上的
                          参数 t，非边界节点为 NaN（v3 起）。**下标 == ``.msh``
                          文件里的节点顺序**，也就是 libMesh 读入后
                          （``allow_renumbering(false)``）的节点 id。
    /ports/<name>         (1,m) int32 / (k,m) float64   端口表，名字由
                          infile.pot 里的 ``*_varname`` 决定
    /meta/*               属性：format_version / name / freq / epsilon_* /
                          mu_* / w_wg / L_sim / L_wg_fixed / 网格尺寸参数

``design.h5`` 是"算例的权威描述"：只靠它就能还原出完整的
:class:`design.Case`（见 ``Case.load``）。

``node_t`` 是"网格前移"的产物：以前 C++ 用 ``mesh_snap_to_Bezier`` 把节点
吸附到曲线上，同时顺手记下参数 t；现在网格由 gmsh 直接在真曲线上生成，
t 就由 Python 一次算好。C++ 侧会逐节点校验 ``‖node − P(t)‖``，所以
"node id == 文件里的节点顺序"这个约定一旦被破坏会立刻报错，不会静默错算。

``result.h5``  (C++ -> Python，求解结果)

    /scalars/S_params         (4,)
    /scalars/ex_opt_sens      (n,2)
    /scalars/in_opt_sens      (n,2)
    /fields/Ez                (n_dofs,2)   仅当 save_field 打开
    /meta/*                   属性：format_version / freq / n_process ...

复数统一存成 ``(n,2)`` 的实数数据集并带属性 ``complex = 1``（实部第 0 列、
虚部第 1 列）。这与 PETSc 的 HDF5 viewer 写复数 Vec 的约定一致。

各 rank 的局部量（S 参数的单元积分、按节点拼出来的场）现在由 C++ 侧用 MPI
归约后再写出，所以这里读到的一律是**全局量**，不需要求和。
"""

import os

import h5py
import numpy as np

FORMAT_VERSION = 2          # result.h5 的布局版本
DESIGN_FORMAT_VERSION = 3   # design.h5 的布局版本（3 起含 /design/node_t）

DESIGN_KEYS = ("in_Bezier_P", "in_Bezier_V", "ex_Bezier_P", "ex_Bezier_V")


def write_design(path, *, in_P, in_V, ex_P, ex_V, port_id, port_is_analytic,
                 port_coord, port_x_axis, incident_on, port_other_params,
                 dofs=None, node_t=None, port_names=None, extra_meta=None):
    """把几何（Bezier 控制点）、端口表和节点参数 t 写成 design.h5。

    ``node_t`` 见模块开头的说明；不给就跳过（旧行为，只对不需要网格的
    检查有用——求解器要求它存在）。

    ``port_names`` 用于覆盖端口表在文件里的名字，键是 ``infile.pot`` 里
    对应的 ``*_varname``；默认与 ``original_mesh.py`` 写出的一致。
    """
    names = {
        "port_id": "port_id",
        "port_is_analytic": "port_is_analytic",
        "port_coord": "port_coord",
        "port_x_axis": "port_x_axis",
        "incident_on": "incident_on",
        "port_other_params": "port_other_params",
    }
    if port_names:
        names.update(port_names)

    _remove_quietly(path)
    with h5py.File(path, "w") as f:
        d = f.create_group("design")
        for name, arr in (("in_Bezier_P", in_P), ("in_Bezier_V", in_V),
                          ("ex_Bezier_P", ex_P), ("ex_Bezier_V", ex_V)):
            d.create_dataset(name, data=np.asarray(arr, dtype=np.float64))
        if dofs is not None:
            d.create_dataset("dofs", data=np.asarray(dofs, dtype=np.float64))
        if node_t is not None:
            # 1 维数据集：C++ 侧的 read_matrix 会当成 (1,n) 读回
            d.create_dataset("node_t", data=np.asarray(node_t, dtype=np.float64))

        p = f.create_group("ports")
        for key, arr in (("port_id", port_id),
                         ("port_is_analytic", port_is_analytic),
                         ("port_coord", port_coord),
                         ("port_x_axis", port_x_axis),
                         ("incident_on", incident_on),
                         ("port_other_params", port_other_params)):
            a = np.atleast_2d(np.asarray(arr))
            # 整数表（端口 id / 开关）用 int32，实数表保持 float64
            if key in ("port_id", "port_is_analytic", "incident_on"):
                a = a.astype(np.int32)
            else:
                a = a.astype(np.float64)
            p.create_dataset(names[key], data=a)

        m = f.create_group("meta")
        m.attrs["format_version"] = DESIGN_FORMAT_VERSION
        for k, v in (extra_meta or {}).items():
            m.attrs[k] = v
    return path


def read_design(path):
    """读回 design.h5（主要给测试和检查用）。"""
    with h5py.File(path, "r") as f:
        out = {k: f["design/" + k][...] for k in DESIGN_KEYS}
        # dofs / node_t 是后来才有的；旧文件没有
        out["dofs"] = f["design/dofs"][...] if "design/dofs" in f else None
        out["node_t"] = f["design/node_t"][...] if "design/node_t" in f else None
        out["ports"] = {k: f["ports/" + k][...] for k in f["ports"]}
        out["meta"] = dict(f["meta"].attrs)
    return out


def _complex(arr):
    """(n,2) -> complex；也兼容本来就存成 complex 的情况。"""
    a = np.asarray(arr)
    if np.iscomplexobj(a):
        return a.ravel()
    if a.ndim == 2 and a.shape[1] == 2:
        return a[:, 0] + 1j * a[:, 1]
    if a.ndim == 1:
        return a.astype(np.complex128)
    raise ValueError("不认识的复数布局: %r" % (a.shape,))


def read_result(path, want_field=False):
    """读 result.h5，返回统一结构的 dict。

    返回的键：``S_params``、``sens_ex``、``sens_in``、``Ez``（可选）、``meta``。
    文件里的量已经是全局的（各 rank 的部分量在 C++ 侧归约过了）。
    """
    if not os.path.exists(path):
        raise FileNotFoundError("找不到结果文件: %s" % path)

    out = {}
    with h5py.File(path, "r") as f:
        meta = dict(f["meta"].attrs) if "meta" in f else {}
        version = int(meta.get("format_version", 0))
        if version != FORMAT_VERSION:
            raise ValueError(
                "result.h5 的格式版本是 %d，本代码需要 %d。若是旧文件（每 rank 一份 "
                "S_params_<rank>），请重新跑一次求解器。" % (version, FORMAT_VERSION))

        out["S_params"] = _complex(f["scalars/S_params"][...])
        out["sens_ex"] = f["scalars/ex_opt_sens"][...]
        out["sens_in"] = f["scalars/in_opt_sens"][...]

        if want_field and "fields/Ez" in f:
            out["Ez"] = _complex(f["fields/Ez"][...])

        out["meta"] = meta
        out["n_process"] = int(meta.get("n_process", 1))
    return out


def mat_to_design(mat_path, h5_path, port_names=None):
    """迁移辅助：把旧的 cross_waveguide.mat 转成 design.h5。

    有了这个就不必为了换格式而重新跑一遍 gmsh。
    """
    import scipy.io as spio

    m = spio.loadmat(mat_path)

    def get(name):
        return np.asarray(m[name])

    return write_design(
        h5_path,
        in_P=get("in_Bezier_P"), in_V=get("in_Bezier_V"),
        ex_P=get("ex_Bezier_P"), ex_V=get("ex_Bezier_V"),
        port_id=get("port_id"), port_is_analytic=get("port_is_analytic"),
        port_coord=get("port_coord"), port_x_axis=get("port_x_axis"),
        incident_on=get("incident_on"), port_other_params=get("port_other_params"),
        port_names=port_names,
    )


def _remove_quietly(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass
