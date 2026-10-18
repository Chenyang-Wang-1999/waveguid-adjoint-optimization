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
from math import cos, pi, sin, sqrt

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

# 网格尺寸参数。迁移之后，剖分完全由这几个数决定（以前是 C++ 里的自适应加密，
# 见已删除的 EM_mesh_flagging.cpp）。
PORT_SIZE = 0.0037       # 端口曲线（input/output*）上的目标边长。S 参数是沿端口
                         # 积分出来的，这里的密度直接影响积分精度，所以刻意给得
                         # 比内部细得多。
WL_PER_SIZE = 16         # 介质内部：目标边长 ≤ 介质内波长 / WL_PER_SIZE
GRADE_LEN = 0.025        # 从端口尺寸过渡到内部尺寸的距离
DISTANCE_SAMPLING = 200  # Distance 场每条曲线的采样点数

# 允许优化边界越出对称轴多少（相对内部目标边长）。见 geometry_margins。
GEOMETRY_TOL_FACTOR = 0.05

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


def bezier_polygons(dofs, geometry, layout):
    """两条优化边界的分段三次 Bézier 控制多边形，形状 ``(n_seg, 4, 2)``。

    ``Bezier_interp`` 在第 i 段上的取值就是控制多边形
    ``[P_i, P_i + V_i/3, P_{i+1} - V_{i+1}/3, P_{i+1}]`` 的三次 Bernstein 形式，
    所以把每段建成一条 gmsh 三次 Bézier 之后，几何与 ``Bezier_interp``
    逐位相同（已验证）。端点切向 ``θ_in``/``θ_ex`` 由 ``V`` 的两端决定、
    段间 C² 连续，都随控制多边形自动保留，不需要额外约束。
    """
    in_P, in_V, ex_P, ex_V = bezier_curves(dofs, geometry, layout)
    return _cubic_polygons(in_P, in_V), _cubic_polygons(ex_P, ex_V)


def geometry_margins(dofs, geometry, layout, n_samples=2000):
    """两条优化边界"离对称轴还有多远"：正数表示还在自己象限里，负数表示越界。

    四个象限的曲线是把第一象限的那条**镜像**出来的，所以曲线一旦越过对称轴，
    相邻象限的镜像就会真的相交。这在物理上意味着波导壁越过中心线（端口几何
    已经不成立），在数值上意味着 gmsh 报 "intersections in the 1D mesh" 然后
    **直接 abort 整个进程**（gmsh 4.9.5 的致命错误是 C++ throw，Python 侧拦不住）。
    所以 :func:`build_mesh` 会先用它拦一道。
    """
    import Bezier

    in_P, in_V, ex_P, ex_V = bezier_curves(dofs, geometry, layout)
    t = np.linspace(0.0, layout.n_ctrl + 1, n_samples)
    ex = Bezier.Bezier_interp(ex_P, ex_V, t.copy(), 0)
    inn = Bezier.Bezier_interp(in_P, in_V, t.copy(), 0)
    return {
        "ex_min_x": float(ex[:, 0].min()),
        "ex_min_y": float(ex[:, 1].min()),
        "in_min_x": float(inn[:, 0].min()),
        "in_min_y": float(inn[:, 1].min()),
    }


def _cubic_polygons(P, V):
    """控制点 + 切向 → 逐段的三次 Bézier 控制多边形。"""
    return np.stack([
        np.stack([P[i], P[i] + V[i] / 3.0, P[i + 1] - V[i + 1] / 3.0, P[i + 1]])
        for i in range(P.shape[0] - 1)
    ])



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


def _guidance_wavelength(freq, epsilon):
    """判据里用的介质内波长。

    沿用旧实现（``EM_mesh_flagging.cpp``）的写法 ``1/(freq*max(1/freq, √ε))``：
    当 ``freq < 1`` 时把"折射率"抬到 ``1/freq``，等价于把波长上限卡在 1。
    当前算例（``freq=1.55``）下这个保护项不起作用。
    """
    index = max(1.0 / freq, sqrt(epsilon))
    return 1.0 / (freq * index)


def mesh_sizes(case):
    """剖分用的目标边长（含单位说明，全部是"边长"而不是"面积"）。

    旧实现是让 libMesh 反复二分：单元体积 > (λ/16)^dim 就再加密一层，直到
    不再有单元超标。现在直接把等价的**目标边长**交给 gmsh 的尺寸场，一轮
    生成的网格就满足同一个判据，而且各向同性、没有二分产生的层级痕迹。
    """
    ph = case.physics
    return {
        "h_port": PORT_SIZE,
        "h_wg": _guidance_wavelength(ph.freq, ph.epsilon_wg) / WL_PER_SIZE,
        "h_air": _guidance_wavelength(ph.freq, ph.epsilon_amb) / WL_PER_SIZE,
        "grade_len": GRADE_LEN,
    }


class _PointCache:
    """按坐标去重的 gmsh 点工厂。

    gmsh 要求相邻曲线共享**同一个点实体**（pygmsh 的 ``CurveLoop`` 直接用
    ``points[-1] is points[0]`` 断言端点是同一个对象），而四个象限是镜像
    出来的、共享端点很多，靠手工索引最容易出错。这里的镜像就是乘以 ±1，
    不引入舍入误差，所以"坐标相同"严格等价于"端点相同"，用缓存最省事。
    """

    def __init__(self, model):
        self._model = model
        self._cache = {}

    def __call__(self, xy):
        key = (float(xy[0]), float(xy[1]))
        if key not in self._cache:
            self._cache[key] = self._model.add_point(key + (0.0,))
        return self._cache[key]


def _backward(curves):
    """把一串正向曲线按相反方向使用（顺序反转 + 每条取负）。"""
    return [-c for c in reversed(curves)]


def _apply_size_field(sizes, surfaces_wg, surfaces_air, port_curves):
    """把目标边长交给 gmsh 的背景尺寸场。

    端口处 ``PORT_SIZE``，沿距离 ``GRADE_LEN`` 过渡到介质内部的目标尺寸
    （波导和空气各一套，因为两者的波长不同）。

    两个坑（都是实测出来的）：

    * ``Constant`` 字段有 ``VIn/VOut`` 语义，不指定区域时会静默给出巨大
      尺寸，所以宁可用 ``MathEval``；
    * ``Restrict`` 里的字段按"点所属的实体"取值，所以再单独给端口曲线
      一条 ``Restrict(MathEval(h_port), CurvesList=…)`` 兜底，保证一维网格
      在端口上一定是细的（否则 ``Min`` 里两条 ``Threshold`` 在共享边界上
      可能都退化成 1e22）。
    """
    import gmsh

    F = gmsh.model.mesh.field

    dist = F.add("Distance")
    F.setNumbers(dist, "CurvesList", [int(c._id) for c in port_curves])
    F.setNumber(dist, "Sampling", DISTANCE_SAMPLING)

    def threshold(size_max, surfaces):
        tag = F.add("Threshold")
        F.setNumber(tag, "InField", dist)
        F.setNumber(tag, "SizeMin", sizes["h_port"])
        F.setNumber(tag, "SizeMax", size_max)
        F.setNumber(tag, "DistMin", 0.0)
        F.setNumber(tag, "DistMax", sizes["grade_len"])
        restricted = F.add("Restrict")
        F.setNumber(restricted, "InField", tag)
        F.setNumbers(restricted, "SurfacesList", [int(s._id) for s in surfaces])
        return restricted

    port_value = F.add("MathEval")
    F.setString(port_value, "F", repr(float(sizes["h_port"])))
    port_field = F.add("Restrict")
    F.setNumber(port_field, "InField", port_value)
    F.setNumbers(port_field, "CurvesList", [int(c._id) for c in port_curves])

    combined = F.add("Min")
    F.setNumbers(combined, "FieldsList", [
        threshold(sizes["h_wg"], surfaces_wg),
        threshold(sizes["h_air"], surfaces_air),
        port_field,
    ])
    F.setAsBackgroundMesh(combined)

    # 只让尺寸场决定单元大小：不要点尺寸、不要边界外插、不要曲率加密
    # （后者在 gmsh 4.9.5 里默认关着，这里显式写出来，避免换版本后行为漂移）
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def build_mesh(case):
    """用 gmsh 生成 ``case.mesh_file``，返回各 physical 的数字 id。

    两条优化边界建成**分段三次 Bézier**（每个控制点区间一段），不是折线：
    段的控制多边形 ``[P_i, P_i+V_i/3, P_{i+1}-V_{i+1}/3, P_{i+1}]`` 与
    ``Bezier_interp`` 逐位等价，所以网格节点天生就在真曲线上，不再需要
    "先生成折线网格、再吸附到曲线"那一步（那正是并行下崩溃的来源）。
    """
    # gmsh / pygmsh 只在真正建模时才导入，方便单独使用本模块的纯计算函数
    import gmsh
    import pygmsh

    g = case.geometry
    sizes = mesh_sizes(case)

    # 越界检查必须放在 gmsh 之前：越界的几何会让 gmsh 直接 abort 掉进程
    margins = geometry_margins(case.dofs, g, case.layout)
    tol = GEOMETRY_TOL_FACTOR * sizes["h_wg"]
    worst = min(margins.values())
    if worst < -tol:
        raise ValueError(
            "优化边界越出了对称轴（%s，最深处 %.2e，允许 %.2e）。四个象限是镜像"
            "出来的，越界会让相邻象限的曲线真的相交；gmsh 遇到这种几何会直接"
            " abort 整个进程，所以这里先拦下来。"
            % ({k: round(v, 6) for k, v in margins.items()}, worst, tol))

    in_polys, ex_polys = bezier_polygons(case.dofs, g, case.layout)

    geom = pygmsh.geo.Geometry()
    model = geom.__enter__()
    pts = _PointCache(model)

    # 4 个象限：坐标乘以 scale。±1 乘法是精确的，所以镜像出来的共享端点
    # 会被 _PointCache 合并成同一个 gmsh 点。
    SCALES = {"m0": (1, 1), "mx": (1, -1), "my": (-1, 1), "mxy": (-1, -1)}

    def scale_of(polys, scale):
        return np.asarray(polys, dtype=float) * np.asarray(scale, dtype=float)

    def curves_of(polys, scale):
        return [model.add_bezier([pts(xy) for xy in poly])
                for poly in scale_of(polys, scale)]

    def end_of(polys, scale, first):
        p = scale_of(polys, scale)
        return pts(p[0, 0] if first else p[-1, 3])

    def at(xy, scale):
        return pts(np.asarray(xy, dtype=float) * np.asarray(scale, dtype=float))

    ex_curves = {k: curves_of(ex_polys, s) for k, s in SCALES.items()}
    in_curves = {k: curves_of(in_polys, s) for k, s in SCALES.items()}
    # 每条曲线两个端点：固定端（t=0）与端口端（t=n_ctrl+1）
    ex_fixed_pt = {k: end_of(ex_polys, s, True) for k, s in SCALES.items()}
    ex_opt_pt = {k: end_of(ex_polys, s, False) for k, s in SCALES.items()}
    in_fixed_pt = {k: end_of(in_polys, s, True) for k, s in SCALES.items()}
    in_opt_pt = {k: end_of(in_polys, s, False) for k, s in SCALES.items()}

    origin = pts((0.0, 0.0))

    # ---- 第 0 象限（m0）----
    s = SCALES["m0"]
    curves = []
    # ex_opt 由端口端走向固定端，in_opt 由固定端走向端口端（闭合环的走向要求）
    ex_opt_m0 = _backward(ex_curves["m0"])
    curves += ex_opt_m0
    ex_fixed_m0 = [model.add_line(ex_fixed_pt["m0"], at(g.wg_input_bd_ex, s))]
    curves += ex_fixed_m0
    input_bd_m0 = [model.add_line(at(g.wg_input_bd_ex, s), at(g.wg_input_bd_in, s))]
    curves += input_bd_m0
    in_fixed_m0 = [model.add_line(at(g.wg_input_bd_in, s), in_fixed_pt["m0"])]
    curves += in_fixed_m0
    in_opt_m0 = list(in_curves["m0"])
    curves += in_opt_m0
    curves += [model.add_line(in_opt_pt["m0"], origin), model.add_line(origin, ex_opt_pt["m0"])]
    m0_curves = curves
    wg1 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第三象限（mxy）----
    s = SCALES["mxy"]
    curves = []
    ex_opt_mxy = _backward(ex_curves["mxy"])
    curves += ex_opt_mxy
    ex_fixed_mxy = [model.add_line(ex_fixed_pt["mxy"], at(g.wg_input_bd_ex, s))]
    curves += ex_fixed_mxy
    input_bd_mxy = [model.add_line(at(g.wg_input_bd_ex, s), at(g.wg_input_bd_in, s))]
    curves += input_bd_mxy
    in_fixed_mxy = [model.add_line(at(g.wg_input_bd_in, s), in_fixed_pt["mxy"])]
    curves += in_fixed_mxy
    in_opt_mxy = list(in_curves["mxy"])
    curves += in_opt_mxy
    curves += [model.add_line(in_opt_pt["mxy"], origin), model.add_line(origin, ex_opt_pt["mxy"])]
    mxy_curves = curves
    wg3 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第二象限（mx，y<0）----
    s = SCALES["mx"]
    curves = []
    in_opt_mx = _backward(in_curves["mx"])
    curves += in_opt_mx
    in_fixed_mx = [model.add_line(in_fixed_pt["mx"], at(g.wg_input_bd_in, s))]
    curves += in_fixed_mx
    input_bd_mx = [model.add_line(at(g.wg_input_bd_in, s), at(g.wg_input_bd_ex, s))]
    curves += input_bd_mx
    ex_fixed_mx = [model.add_line(at(g.wg_input_bd_ex, s), ex_fixed_pt["mx"])]
    curves += ex_fixed_mx
    ex_opt_mx = list(ex_curves["mx"])
    curves += ex_opt_mx
    curves += [-mxy_curves[-1], -m0_curves[-2]]
    wg2 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 第四象限（my，x<0）----
    s = SCALES["my"]
    curves = []
    in_opt_my = _backward(in_curves["my"])
    curves += in_opt_my
    in_fixed_my = [model.add_line(in_fixed_pt["my"], at(g.wg_input_bd_in, s))]
    curves += in_fixed_my
    input_bd_my = [model.add_line(at(g.wg_input_bd_in, s), at(g.wg_input_bd_ex, s))]
    curves += input_bd_my
    ex_fixed_my = [model.add_line(at(g.wg_input_bd_ex, s), ex_fixed_pt["my"])]
    curves += ex_fixed_my
    ex_opt_my = list(ex_curves["my"])
    curves += ex_opt_my
    curves += [-m0_curves[-1], -mxy_curves[-2]]
    wg4 = model.add_plane_surface(model.add_curve_loop(curves))

    # ---- 外圈空气 ----
    sim_region_pts = [pts((g.L_sim * cos(2 * pi * n / 6), g.L_sim * sin(2 * pi * n / 6)))
                      for n in range(6)]

    curves = []
    port_m0_air1 = [model.add_line(sim_region_pts[0], at(g.wg_input_bd_in, SCALES["m0"]))]
    curves += port_m0_air1
    curves += in_fixed_m0
    curves += in_opt_m0
    curves += in_opt_mx
    curves += in_fixed_mx
    port_mx_air1 = [model.add_line(at(g.wg_input_bd_in, SCALES["mx"]), sim_region_pts[0])]
    curves += port_mx_air1
    air1 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[1], sim_region_pts[2])]
    port_my_air2 = [model.add_line(sim_region_pts[2], at(g.wg_input_bd_ex, SCALES["my"]))]
    curves += port_my_air2
    curves += ex_fixed_my
    curves += ex_opt_my
    curves += ex_opt_m0
    curves += ex_fixed_m0
    port_m0_air2 = [model.add_line(at(g.wg_input_bd_ex, SCALES["m0"]), sim_region_pts[1])]
    curves += port_m0_air2
    air2 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    port_mxy_air3 = [model.add_line(sim_region_pts[3], at(g.wg_input_bd_in, SCALES["mxy"]))]
    curves += port_mxy_air3
    curves += in_fixed_mxy
    curves += in_opt_mxy
    curves += in_opt_my
    curves += in_fixed_my
    port_my_air3 = [model.add_line(at(g.wg_input_bd_in, SCALES["my"]), sim_region_pts[3])]
    curves += port_my_air3
    air3 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[4], sim_region_pts[5])]
    port_mx_air4 = [model.add_line(sim_region_pts[5], at(g.wg_input_bd_ex, SCALES["mx"]))]
    curves += port_mx_air4
    curves += ex_fixed_mx
    curves += ex_opt_mx
    curves += ex_opt_mxy
    curves += ex_fixed_mxy
    port_mxy_air4 = [model.add_line(at(g.wg_input_bd_ex, SCALES["mxy"]), sim_region_pts[4])]
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

    # 端口曲线（4 个端口 × 口径 + 两条径线）：尺寸场在这里加密到 PORT_SIZE
    port_curves = (input_bd_m0 + port_m0_air1 + port_m0_air2
                   + input_bd_mxy + port_mxy_air3 + port_mxy_air4
                   + input_bd_mx + port_mx_air1 + port_mx_air4
                   + input_bd_my + port_my_air2 + port_my_air3)
    _apply_size_field(sizes, [wg1, wg2, wg3, wg4], [air1, air2, air3, air4],
                      port_curves)

    model.generate_mesh(dim=2, order=2)
    # 写 MSH 2.2 而不是默认的 4.1：libMesh 1.9 读 4.x 时会校验 $Entities 里
    # 每条 curve 的包围盒，而 gmsh 4.9.5 给**曲边**算的包围盒比真实曲线小
    # 一点点（实测 5.5e-6，而 libMesh 的相对容差只有 1e-8 量级），于是读网格
    # 直接报 "outside entity physical bounding box"。2.2 没有 $Entities 段，
    # physical tag / $PhysicalNames 一样保留。
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
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


def _physical_tag(ids, mesh, label):
    """优化边界的 physical id：优先用建网格时记下的数字，其次查网格里的名字。"""
    if ids and label in ids:
        return int(ids[label])
    if label in mesh.field_data:
        return int(mesh.field_data[label][0])
    return None


def boundary_param_t(case):
    """每个网格节点在优化边界曲线上的参数 t，非边界节点记 ``NaN``。

    返回数组的**下标 == ``.msh`` 文件里的节点顺序**，也就是 libMesh 读入
    （``allow_renumbering(false)``）之后的节点 id。libMesh 的 gmsh reader
    明确不支持 parametric 坐标（``libmesh_error_msg_if(parametric, …)``），
    而 ``prepare_for_use()`` 又会按单元遍历顺序重排节点 id，所以这条路上
    唯一干净的传递方式就是"文件里的顺序"。

    C++ 侧会逐节点校验 ``‖node − P(t)‖``，所以这个约定一旦错位会立刻报错，
    不会静默算错灵敏度。
    """
    import Bezier
    import meshio

    mesh = meshio.read(case.mesh_file)
    ids = case.ids or {}
    in_P, in_V, ex_P, ex_V = bezier_curves(case.dofs, case.geometry, case.layout)
    curves = {
        _physical_tag(ids, mesh, "ex_opt"): (ex_P, ex_V),
        _physical_tag(ids, mesh, "in_opt"): (in_P, in_V),
    }
    curves.pop(None, None)

    t = np.full(len(mesh.points), np.nan)
    for block, tags in zip(mesh.cells, mesh.cell_data.get("gmsh:physical", [])):
        if block.type not in ("line", "line3"):
            continue
        # 逐元素看 physical tag：MSH 4.x 里 gmsh 按**几何实体**分块（一块一个
        # physical），2.2 里按**单元类型**分块（一个 line3 块里混着所有端口）。
        # meshio 用带符号的 tag 表示单元方向，物理组本身取绝对值。
        tags = np.abs(np.asarray(tags).ravel())
        for pid, (P, V) in curves.items():
            for idx in np.unique(block.data[tags == pid].ravel()):
                # 四个象限是镜像出来的：把点折回第一象限再求参数（与 C++ 一致）。
                # 容差取 1e-8 而不是 C++ 里的 1e-5：Newton 二阶收敛，多走两步
                # 不花钱，却能把参数误差从 ~5e-3（t 的量纲）压到 1e-12 量级。
                pt = np.abs(mesh.points[int(idx)][:2])
                t[int(idx)] = float(Bezier.Nearest_t(P, V, pt, 1e-8))
    return t


def write_design(case):
    """把设计变量、端口表和节点参数 t 写成 ``design.h5``。

    ``node_t`` 是当场从 ``case.mesh_file`` 读回来算的，所以 ``design.h5``
    永远和磁盘上那张网格一致——这也是求解器能直接按节点 id 查表的前提。
    """
    in_points, in_tangent, ex_points, ex_tangent = bezier_curves(
        case.dofs, case.geometry, case.layout)
    ports = port_tables(case)
    sizes = mesh_sizes(case)
    ph, g = case.physics, case.geometry
    return dataio.write_design(
        case.design_file,
        in_P=in_points, in_V=in_tangent,
        ex_P=ex_points, ex_V=ex_tangent,
        dofs=case.dofs,
        node_t=boundary_param_t(case),
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
            # 网格尺寸参数：没有它们就复现不出同一张网格
            "h_port": float(sizes["h_port"]),
            "h_wg": float(sizes["h_wg"]),
            "h_air": float(sizes["h_air"]),
            "grade_len": float(sizes["grade_len"]),
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
