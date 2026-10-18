"""几何与参数设置：设计变量布局、Bézier 控制点、端口表、gmsh 建模。

这个模块是"网格产生"的唯一入口，对外只暴露 :func:`make_case`。

设计变量布局（``Layout``，只有 4N+4 这一套）::

    [ x_in_bd, theta_in, (x,y) x N_ctrl, y_ex_bd, theta_ex, (x,y) x N_ctrl ]

``x_in_bd``/``y_ex_bd`` 是两个端口的边界点位置，``theta_*`` 是端点切向角，
中间 2*N_ctrl 个是两条 Bézier 曲线内部的控制点坐标。
"""

import os
import sys
from dataclasses import dataclass, field
from math import cos, pi, sin

import numpy as np

# ``Bezier`` 是编在 python/ 下的 pybind11 扩展模块。如果调用方一开始就在
# python/ 里跑（``python3 optimize.py``），它天然可导入；但如果脚本先
# ``chdir`` 到 data/ 再用相对路径导入，``sys.path`` 里的 '' 会跟着变，
# 就找不到它了。这里把本文件所在目录固定加进去，避免这类脆弱性。
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import dataio
import physics as physics_mod

# gmsh 剖分参数（沿用旧值，改动会影响网格）
MESH_SIZE_FACTOR = 0.5   # 单元目标尺寸 = MESH_SIZE_FACTOR * w_wg
PORT_REFINE = 3          # 端口处再加密的倍数
FAR_FIELD_FACTOR = 10    # 远场（空气外圈）尺寸倍数
N_BPOINTS = 50           # 每条边界曲线的采样点数

# 端口编号（对应 gmsh 里 add_physical 的顺序，见 build_mesh 末尾）
PORT_IDS = np.array([3, 4, 5, 6], dtype=int)


# --------------------------------------------------------------------- 参数


@dataclass(frozen=True)
class Geometry:
    """几何尺寸（不随优化变化的部分）。"""

    w_wg: float = 0.08
    L_sim: float = 0.5 / cos(pi / 6)
    L_wg_fixed: float = 0.3

    @property
    def e1(self):
        return np.array([cos(pi / 6), sin(pi / 6)])

    @property
    def e2(self):
        return np.array([-sin(pi / 6), cos(pi / 6)])

    @property
    def e3(self):
        return np.array([cos(pi / 3), sin(pi / 3)])

    @property
    def wg_input_cen(self):
        return self.L_sim * cos(pi / 6) * self.e1

    @property
    def wg_input_bd_ex(self):
        return self.wg_input_cen + self.w_wg * self.e2 / 2

    @property
    def wg_input_bd_in(self):
        return self.wg_input_cen - self.w_wg * self.e2 / 2

    @property
    def wg_fixed_bd_ex(self):
        return self.wg_input_bd_ex - self.L_wg_fixed * self.e1

    @property
    def wg_fixed_bd_in(self):
        return self.wg_input_bd_in - self.L_wg_fixed * self.e1


@dataclass(frozen=True)
class Layout:
    """设计变量布局：只需要控制点个数。"""

    n_ctrl: int

    @classmethod
    def from_dofs(cls, dofs):
        return cls(int((len(dofs) - 4) / 4 + 0.5))

    @property
    def size(self):
        return 4 * self.n_ctrl + 4


@dataclass
class Case:
    """①（几何 + 参数）交给 ②（求解器）的完整描述。"""

    dofs: np.ndarray
    physics: physics_mod.Physics = field(default_factory=physics_mod.Physics)
    geometry: Geometry = field(default_factory=Geometry)
    name: str = "cross_waveguide"
    # gmsh 里实际的 physical id；只有 build_mesh 跑过之后才有值，
    # 拿不到时求解器会回落到按名字查。
    ids: dict = None

    @property
    def layout(self):
        return Layout.from_dofs(self.dofs)

    @property
    def mesh_file(self):
        return self.name + ".msh"

    @property
    def out_mesh_file(self):
        return "out_" + self.name + ".msh"

    @property
    def design_file(self):
        return "design.h5"

    @property
    def result_file(self):
        return "result.h5"

    @classmethod
    def load(cls, design_file="design.h5", name=None):
        """从磁盘上的 ``design.h5`` 还原一个算例。

        ``design.h5`` 是自描述的（控制点 + 设计变量 + 物理/几何参数），所以
        后处理或重算时不需要调用方手里还留着 DOF 向量。
        """
        data = dataio.read_design(design_file)
        meta = data["meta"]
        version = int(meta.get("format_version", 0))
        if version != dataio.DESIGN_FORMAT_VERSION:
            raise ValueError(
                "%s 的格式版本是 %d，本代码需要 %d。旧文件缺少 /design/dofs，"
                "请用 make_mesh.py 重新生成。" % (design_file, version,
                                                  dataio.DESIGN_FORMAT_VERSION))
        if data["dofs"] is None:
            raise ValueError("%s 里没有 /design/dofs，无法还原算例。" % design_file)

        return cls(
            dofs=np.asarray(data["dofs"], dtype=float),
            physics=physics_mod.Physics(
                freq=float(meta["freq"]),
                epsilon_wg=float(meta["epsilon_wg"]),
                epsilon_amb=float(meta["epsilon_amb"]),
                mu_wg=float(meta["mu_wg"]),
                mu_amb=float(meta["mu_amb"]),
            ),
            geometry=Geometry(
                w_wg=float(meta["w_wg"]),
                L_sim=float(meta["L_sim"]),
                L_wg_fixed=float(meta["L_wg_fixed"]),
            ),
            name=name or str(meta.get("name", "cross_waveguide")),
        )


# --------------------------------------------------------- 设计变量 <-> 控制点


def initial_dofs(geometry, layout):
    """未优化的初始设计：两条边界都是直线（旧 ``use_original`` 分支）。"""
    g = geometry
    x0 = g.w_wg / 2 / sin(pi / 6)
    y0 = g.w_wg / 2 / cos(pi / 6)
    ex_origin = np.linspace(g.wg_fixed_bd_ex, [0, y0], layout.n_ctrl + 2)
    in_origin = np.linspace(g.wg_fixed_bd_in, [x0, 0], layout.n_ctrl + 2)
    return np.hstack([
        [x0, 7 * pi / 6], in_origin[1:-1].flatten(),
        [y0, 7 * pi / 6], ex_origin[1:-1].flatten(),
    ])


def control_points(dofs, geometry, layout):
    """DOF -> 两条边界的控制点与端点切向。"""
    g, n = geometry, layout.n_ctrl
    in_points_bd = np.array([[dofs[0], 0]])
    in_end_vec = np.array([[cos(dofs[1]), sin(dofs[1])]])
    in_points = np.reshape(dofs[2:2 * n + 2], (n, 2))
    in_points = np.vstack([g.wg_fixed_bd_in, in_points, in_points_bd])

    ex_points_bd = np.array([[0, dofs[2 * n + 2]]])
    ex_end_vec = np.array([[cos(dofs[2 * n + 3]), sin(dofs[2 * n + 3])]])
    ex_points = np.reshape(dofs[2 * n + 4:], (n, 2))
    ex_points = np.vstack([g.wg_fixed_bd_ex, ex_points, ex_points_bd])
    return in_points, in_end_vec, ex_points, ex_end_vec


def _points_to_dofs(in_points, in_theta, ex_points, ex_theta):
    """控制点 -> DOF（``control_points`` 的逆）。"""
    n = in_points.shape[0] - 1
    return np.hstack([
        [in_points[-1, 0], in_theta], in_points[:-1, :].flatten(),
        [ex_points[-1, 1], ex_theta], ex_points[:-1, :].flatten(),
    ])


def sensitivity_to_dofs(ex_sens, in_sens):
    """求解器给出的边界灵敏度 -> 按 DOF 顺序排好的梯度向量。

    ``ex_sens``/``in_sens`` 的形状是 ``(n_ctrl+3, 2)``：第 0 行对应固定端
    （不参与设计变量），最后一行是两个端点角度，中间是控制点。
    """
    return _points_to_dofs(in_sens[1:-1, :], in_sens[-1, 1],
                           ex_sens[1:-1, :], ex_sens[-1, 1])


def bezier_curves(dofs, geometry, layout):
    """算出供 gmsh 和 ``design.h5`` 用的 Bézier 控制点与切向。"""
    import Bezier

    in_points, in_end_vec, ex_points, ex_end_vec = control_points(dofs, geometry, layout)
    e1 = geometry.e1
    ex_tangent = Bezier.Bezier_tangent_vector(ex_points, 0, 0, -e1, ex_end_vec.flatten(), True)
    in_tangent = Bezier.Bezier_tangent_vector(in_points, 0, 0, -e1, in_end_vec.flatten(), True)
    return in_points, in_tangent, ex_points, ex_tangent


# ------------------------------------------------------------------ 端口表


def port_tables(case):
    """构造 C++ 侧需要的端口参数表。"""
    g, ph = case.geometry, case.physics
    beta = physics_mod.beta(ph, g.w_wg)

    mx = np.array([1, -1])
    my = np.array([-1, 1])
    mxy = -1
    cen = g.wg_input_cen

    port_coords = np.array([cen, cen * mxy, cen * mx, cen * my])
    port_coords = np.vstack([port_coords.transpose(), np.zeros((1, 4))])

    port_x_axis = np.array([g.e2, -g.e2, g.e3, -g.e3])
    port_x_axis = np.vstack([port_x_axis.transpose(), np.zeros((1, 4))])

    port_other_params = np.array([
        [1, 1, 1, 1],                          # amplitude
        [g.w_wg, g.w_wg, g.w_wg, g.w_wg],      # width
        [0, 0, 0, 0],                          # mode index
        [beta, beta, beta, beta],
        [ph.epsilon_wg] * 4,
        [ph.epsilon_amb] * 4,
        [ph.mu_wg] * 4,
        [ph.mu_amb] * 4,
    ])
    return {
        "port_id": PORT_IDS,
        "port_is_analytic": np.array([1, 1, 1, 1], dtype=int),
        "port_coord": port_coords,
        "port_x_axis": port_x_axis,
        "incident_on": np.array([1, 0, 0, 0], dtype=int),
        "port_other_params": port_other_params,
    }


# ------------------------------------------------------------------ gmsh 建模


def build_mesh(case):
    """用 gmsh 生成 ``case.mesh_file``，返回各 physical 的数字 id。"""
    # gmsh / pygmsh 只在真正建模时才导入，方便单独使用本模块的纯计算函数
    import gmsh
    import pygmsh

    g, ph = case.geometry, case.physics
    layout = case.layout
    mesh_size = MESH_SIZE_FACTOR * g.w_wg
    K = FAR_FIELD_FACTOR
    port_refine = PORT_REFINE
    N_BPoints = N_BPOINTS

    e1 = g.e1
    wg_input_bd_ex = g.wg_input_bd_ex
    wg_input_bd_in = g.wg_input_bd_in

    in_points, in_tangent, ex_points, ex_tangent = bezier_curves(case.dofs, g, layout)

    import Bezier
    t = np.linspace(0, layout.n_ctrl + 1, N_BPoints)
    ex_line_pts = Bezier.Bezier_interp(ex_points, ex_tangent, t, 0)
    in_line_pts = Bezier.Bezier_interp(in_points, in_tangent, t, 0)

    geom = pygmsh.geo.Geometry()
    model = geom.__enter__()

    mx = np.array([1, -1])
    my = np.array([-1, 1])
    mxy = -1

    origin = model.add_point((0, 0, 0), mesh_size=mesh_size)

    input_ex_pt = model.add_point(wg_input_bd_ex, mesh_size=mesh_size / port_refine)
    input_in_pt = model.add_point(wg_input_bd_in, mesh_size=mesh_size / port_refine)
    ex_opt_pts = [model.add_point(ex_line_pts[i, :], mesh_size=mesh_size) for i in range(N_BPoints)]
    in_opt_pts = [model.add_point(in_line_pts[i, :], mesh_size=mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mx = model.add_point(wg_input_bd_ex * mx, mesh_size=mesh_size / port_refine)
    input_in_pt_mx = model.add_point(wg_input_bd_in * mx, mesh_size=mesh_size / port_refine)
    ex_opt_pts_mx = [model.add_point(ex_line_pts[i, :] * mx, mesh_size=mesh_size) for i in range(N_BPoints)]
    in_opt_pts_mx = [model.add_point(in_line_pts[i, :] * mx, mesh_size=mesh_size) for i in range(N_BPoints - 1)] + [in_opt_pts[-1]]

    input_ex_pt_my = model.add_point(wg_input_bd_ex * my, mesh_size=mesh_size / port_refine)
    input_in_pt_my = model.add_point(wg_input_bd_in * my, mesh_size=mesh_size / port_refine)
    ex_opt_pts_my = [model.add_point(ex_line_pts[i, :] * my, mesh_size=mesh_size) for i in range(N_BPoints - 1)] + [ex_opt_pts[-1]]
    in_opt_pts_my = [model.add_point(in_line_pts[i, :] * my, mesh_size=mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mxy = model.add_point(wg_input_bd_ex * mxy, mesh_size=mesh_size / port_refine)
    input_in_pt_mxy = model.add_point(wg_input_bd_in * mxy, mesh_size=mesh_size / port_refine)
    ex_opt_pts_mxy = [model.add_point(ex_line_pts[i, :] * mxy, mesh_size=mesh_size) for i in range(N_BPoints - 1)] + [ex_opt_pts_mx[-1]]
    in_opt_pts_mxy = [model.add_point(in_line_pts[i, :] * mxy, mesh_size=mesh_size) for i in range(N_BPoints - 1)] + [in_opt_pts_my[-1]]

    # ---- 第 0 象限（m0）----
    curves = []
    ex_opt_m0 = [model.add_line(ex_opt_pts[i], ex_opt_pts[i - 1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_m0
    ex_fixed_m0 = [model.add_line(ex_opt_pts[0], input_ex_pt)]
    curves += ex_fixed_m0
    input_bd_m0 = [model.add_line(input_ex_pt, input_in_pt)]
    curves += input_bd_m0
    in_fixed_m0 = [model.add_line(input_in_pt, in_opt_pts[0])]
    curves += in_fixed_m0
    in_opt_m0 = [model.add_line(in_opt_pts[i], in_opt_pts[i + 1]) for i in range(N_BPoints - 1)]
    curves += in_opt_m0
    curves += [model.add_line(in_opt_pts[-1], origin), model.add_line(origin, ex_opt_pts[-1])]
    m0_curves = curves
    wg1 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第三象限（mxy）----
    curves = []
    ex_opt_mxy = [model.add_line(ex_opt_pts_mxy[i], ex_opt_pts_mxy[i - 1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_mxy
    ex_fixed_mxy = [model.add_line(ex_opt_pts_mxy[0], input_ex_pt_mxy)]
    curves += ex_fixed_mxy
    input_bd_mxy = [model.add_line(input_ex_pt_mxy, input_in_pt_mxy)]
    curves += input_bd_mxy
    in_fixed_mxy = [model.add_line(input_in_pt_mxy, in_opt_pts_mxy[0])]
    curves += in_fixed_mxy
    in_opt_mxy = [model.add_line(in_opt_pts_mxy[i], in_opt_pts_mxy[i + 1]) for i in range(N_BPoints - 1)]
    curves += in_opt_mxy
    curves += [model.add_line(in_opt_pts_mxy[-1], origin), model.add_line(origin, ex_opt_pts_mxy[-1])]
    mxy_curves = curves
    wg3 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第二象限（mx）----
    curves = []
    in_opt_mx = [model.add_line(in_opt_pts_mx[i], in_opt_pts_mx[i - 1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_mx
    in_fixed_mx = [model.add_line(in_opt_pts_mx[0], input_in_pt_mx)]
    curves += in_fixed_mx
    input_bd_mx = [model.add_line(input_in_pt_mx, input_ex_pt_mx)]
    curves += input_bd_mx
    ex_fixed_mx = [model.add_line(input_ex_pt_mx, ex_opt_pts_mx[0])]
    curves += ex_fixed_mx
    ex_opt_mx = [model.add_line(ex_opt_pts_mx[i], ex_opt_pts_mx[i + 1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_mx
    curves += [-mxy_curves[-1], -m0_curves[-2]]
    wg2 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第四象限（my）----
    curves = []
    in_opt_my = [model.add_line(in_opt_pts_my[i], in_opt_pts_my[i - 1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_my
    in_fixed_my = [model.add_line(in_opt_pts_my[0], input_in_pt_my)]
    curves += in_fixed_my
    input_bd_my = [model.add_line(input_in_pt_my, input_ex_pt_my)]
    curves += input_bd_my
    ex_fixed_my = [model.add_line(input_ex_pt_my, ex_opt_pts_my[0])]
    curves += ex_fixed_my
    ex_opt_my = [model.add_line(ex_opt_pts_my[i], ex_opt_pts_my[i + 1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_my
    curves += [-m0_curves[-1], -mxy_curves[-2]]
    wg4 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 外圈空气 ----
    sim_region_pts = [model.add_point((g.L_sim * cos(2 * pi * n / 6), g.L_sim * sin(2 * pi * n / 6)),
                                      mesh_size=K * mesh_size) for n in range(6)]

    curves = []
    port_m0_air1 = [model.add_line(sim_region_pts[0], input_in_pt)]
    curves += port_m0_air1
    curves += in_fixed_m0
    curves += in_opt_m0
    curves += in_opt_mx
    curves += in_fixed_mx
    port_mx_air1 = [model.add_line(input_in_pt_mx, sim_region_pts[0])]
    curves += port_mx_air1
    air1 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[1], sim_region_pts[2])]
    port_my_air2 = [model.add_line(sim_region_pts[2], input_ex_pt_my)]
    curves += port_my_air2
    curves += ex_fixed_my
    curves += ex_opt_my
    curves += ex_opt_m0
    curves += ex_fixed_m0
    port_m0_air2 = [model.add_line(input_ex_pt, sim_region_pts[1])]
    curves += port_m0_air2
    air2 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    port_mxy_air3 = [model.add_line(sim_region_pts[3], input_in_pt_mxy)]
    curves += port_mxy_air3
    curves += in_fixed_mxy
    curves += in_opt_mxy
    curves += in_opt_my
    curves += in_fixed_my
    port_my_air3 = [model.add_line(input_in_pt_my, sim_region_pts[3])]
    curves += port_my_air3
    air3 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[4], sim_region_pts[5])]
    port_mx_air4 = [model.add_line(sim_region_pts[5], input_ex_pt_mx)]
    curves += port_mx_air4
    curves += ex_fixed_mx
    curves += ex_opt_mx
    curves += ex_opt_mxy
    curves += ex_fixed_mxy
    port_mxy_air4 = [model.add_line(input_ex_pt_mxy, sim_region_pts[4])]
    curves += port_mxy_air4
    air4 = model.add_plane_surface(model.add_curve_loop(curves))

    model.synchronize()
    model.add_physical([wg1, wg2, wg3, wg4], "waveguide")
    model.add_physical([air1, air2, air3, air4], "air")
    model.add_physical(input_bd_m0 + port_m0_air1 + port_m0_air2, "input")
    model.add_physical(input_bd_mxy + port_mxy_air3 + port_mxy_air4, "output1")
    model.add_physical(input_bd_mx + port_mx_air1 + port_mx_air4, "output21")
    model.add_physical(input_bd_my + port_my_air2 + port_my_air3, "output22")
    model.add_physical(ex_opt_m0 + ex_opt_mx + ex_opt_my + ex_opt_mxy, "ex_opt")
    model.add_physical(in_opt_m0 + in_opt_mx + in_opt_my + in_opt_mxy, "in_opt")

    model.generate_mesh(dim=2, order=2)
    gmsh.write(case.mesh_file)

    # 记下数字 id：求解器优先用它们，就不再依赖 physical name
    by_name = {gmsh.model.getPhysicalName(dim, tag): tag
               for dim, tag in gmsh.model.getPhysicalGroups()}
    ids = {}
    for key, label in (("wg", "waveguide"), ("ex_opt", "ex_opt"), ("in_opt", "in_opt")):
        if label in by_name:
            ids[key] = int(by_name[label])

    gmsh.clear()
    geom.__exit__()
    return ids


# ------------------------------------------------------------------ 入口


def write_design(case):
    """把设计变量与端口表写成 ``design.h5``。"""
    in_points, in_tangent, ex_points, ex_tangent = bezier_curves(
        case.dofs, case.geometry, case.layout)
    ports = port_tables(case)
    ph, g = case.physics, case.geometry
    return dataio.write_design(
        case.design_file,
        in_P=in_points, in_V=in_tangent,
        ex_P=ex_points, ex_V=ex_tangent,
        dofs=case.dofs,
        extra_meta={
            # 有了这些，design.h5 就能独立还原出完整的 Case（见 Case.load）
            "name": case.name,
            "freq": float(ph.freq),
            "epsilon_wg": float(ph.epsilon_wg),
            "epsilon_amb": float(ph.epsilon_amb),
            "mu_wg": float(ph.mu_wg),
            "mu_amb": float(ph.mu_amb),
            "w_wg": float(g.w_wg),
            "L_sim": float(g.L_sim),
            "L_wg_fixed": float(g.L_wg_fixed),
        },
        **ports,
    )


def make_case(dofs, physics=None, geometry=None, name="cross_waveguide", build=True):
    """① 的完整动作：生成网格 + 写 design.h5，返回 :class:`Case`。

    ``build=False`` 只构造对象（用于查询或复用已有网格）。
    """
    case = Case(dofs=np.asarray(dofs, dtype=float),
                physics=physics or physics_mod.Physics(),
                geometry=geometry or Geometry(),
                name=name)
    if build:
        case.ids = build_mesh(case)
        write_design(case)
    return case
