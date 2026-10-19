"""把一版优化结果导出成 COMSOL 能直接导入的 DXF，并生成导入后的校核表。

用法::

    cd data && python3 ../python/export_comsol.py
    cd data && python3 ../python/export_comsol.py --design design.h5 --out comsol

产物（默认 ``data/comsol/``）::

    geometry.dxf   2D 平面几何，按 label 分图层
    CHECK.md       导入后逐项核对的数字 + COMSOL 侧的物理设置

为什么自己写 DXF：gmsh 4.9.5 的**内置** geometry kernel 根本导不出 DXF
（``gmsh.write("x.dxf")`` 报 ``Unknown output file format``），而 STEP/IGES 只认
OCC kernel，我们的几何又是用内置 kernel 建的。好在这条优化壁曲线本就等价于
**一条**三次 B-spline（见 :func:`wall_spline`），DXF 的 SPLINE 实体能逐位承载它
（控制点 + 节点向量原样写进去，不插值、不近似），所以自己写这一百来行比中间隔
一层 CAD 转换更可靠。

坐标系与 ``design.h5`` **完全相同，不做任何缩放**。COMSOL 里把长度单位设成
``m``、频率设成 ``1.55*c_const``，就与本项目的归一化逐位等价；想要 λ₀ = 1.55 µm
的物理器件，把所有长度乘 2.4025、频率换 193.4 THz 即可（物理不变）。
"""

import argparse
import os
from dataclasses import dataclass
from math import cos, pi, sin, sqrt

import numpy as np

import design
import paths

# ---------------------------------------------------------------- 图层约定

#: DXF 图层 -> (颜色号, 说明)。COMSOL 的 DXF 导入不一定保留图层，所以真正的
#: 依据是 ``CHECK.md`` 里的坐标；图层主要用于自查和"按层框选"。
LAYERS = {
    "port_input":    (1, "入射端口：六边形边 0-60°，在波导口处切成 3 段"),
    "port_output1":  (3, "端口：边 180-240°"),
    "port_output21": (4, "端口：边 300-360°"),
    "port_output22": (6, "端口：边 120-180°"),
    "outer":         (8, "非端口的 2 条外边界（y = +-0.5）→ 散射边界条件"),
    "stub":          (7, "固定直壁段：波导口 → 壁的固定端（长 L_wg_fixed）"),
    "wall_in":       (5, "in 优化壁：一条三次 B-spline（14 控制点 / 18 节点）"),
    "wall_ex":       (2, "ex 优化壁：同上"),
}

#: 象限 (符号 x, 符号 y) -> (端口名, 六边形边号)。边 k 从顶点 k 到顶点 k+1，
#: 顶点 k 在 60k 度、半径 L_sim 处。
QUADRANTS = (
    ((1, 1), "input", 0),
    ((-1, 1), "output22", 2),
    ((-1, -1), "output1", 3),
    ((1, -1), "output21", 5),
)

#: 两条非端口的外边界（90° 与 270°）：求解器在它们上面加一阶散射边界条件。
SCATTER_SIDES = (1, 4)

#: ``--polyline`` 退路里每条壁的采样点数。当前算例下弦-弧偏差 ~3e-6，
#: 而网格尺寸是 λ_guided/16 ≈ 0.015，所以这点误差完全看不出来。
POLYLINE_POINTS = 400


@dataclass
class Entity:
    """DXF 里的一条曲线。"""

    layer: str
    kind: str                     # "line" | "spline"
    points: np.ndarray            # LINE: (2, 2)；SPLINE: (n_cp, 2)
    knots: np.ndarray = None      # 仅 SPLINE
    degree: int = 0               # 仅 SPLINE


# ---------------------------------------------------------------- 曲线表示


def wall_spline(P, V):
    """一条优化壁（n 段复合三次 Bézier）→ 一条**完全等价**的三次 B-spline。

    第 i 段的控制多边形是 ``[P_i, P_i+V_i/3, P_{i+1}-V_{i+1}/3, P_{i+1}]``
    （Bernstein 形式，与 ``Bezier_interp`` 逐位相同）；而整条曲线在段间 C² 连续
    （三对角系统的性质），所以它同时落在"两端四重节点、中间全为单重节点"的三次
    B-spline 空间里。用 Greville 横标做配置法（14×14 小线性系统，条件数 3.4）
    解出控制点即可，实测残差 ~1e-16、弧长吻合到 1e-9。

    这样每条壁在 DXF 里就是一个 SPLINE 实体，而不是 11 段首尾相接的曲线——
    COMSOL 里选边界、看面积都干净得多。
    """
    from scipy.interpolate import BSpline

    import Bezier

    n = len(P) - 1
    knots = np.concatenate([np.zeros(4), np.arange(1, n, dtype=float),
                            np.full(4, float(n))])
    n_cp = len(knots) - 4
    greville = np.array([knots[i + 1:i + 4].sum() / 3.0 for i in range(n_cp)])
    basis = np.column_stack([BSpline(knots, np.eye(n_cp)[i], 3)(greville)
                             for i in range(n_cp)])
    target = Bezier.Bezier_interp(P, V, greville.copy(), 0)
    return np.linalg.solve(basis, target), knots


def spline_error(control, knots, P, V, n_samples=8001):
    """B-spline 与原复合 Bézier 的最大偏差（导出前自检用）。"""
    from scipy.interpolate import BSpline

    import Bezier

    t = np.linspace(knots[3], knots[-4], n_samples)
    mine = np.column_stack([BSpline(knots, control[:, 0], 3)(t),
                            BSpline(knots, control[:, 1], 3)(t)])
    return float(np.abs(mine - Bezier.Bezier_interp(P, V, t.copy(), 0)).max())


def hexagon_vertices(geometry):
    """外框（正六边形）的 6 个顶点，第 k 个在 60k 度、半径 L_sim 处。"""
    L = geometry.L_sim
    return np.array([[L * cos(pi / 3 * k), L * sin(pi / 3 * k)] for k in range(6)])


def _ordered_on_segment(a, b, points):
    """把若干点按它们在线段 a→b 上的投影位置排序。"""
    d = b - a
    t = [np.dot(p - a, d) / np.dot(d, d) for p in points]
    return [points[i] for i in np.argsort(t)]


# ---------------------------------------------------------------- 实体表


def build_entities(case):
    """算出整个 2D 几何的曲线表（不含任何网格信息）。

    集合与 :func:`design.build_mesh` 建的曲线一一对应，只有一处不同：原点附近
    连接两片波导的 4 条轴向线段只写**一次**（gmsh 里它们被相邻象限共用）。
    """
    geometry, layout = case.geometry, case.layout
    n_ctrl = layout.n_ctrl

    in_P, in_V, ex_P, ex_V = design.bezier_curves(case.dofs, geometry, layout)
    in_cp, in_knots = wall_spline(in_P, in_V)
    ex_cp, ex_knots = wall_spline(ex_P, ex_V)

    # 自检：写进 DXF 的样条必须就是求解器用的那条曲线
    err = max(spline_error(in_cp, in_knots, in_P, in_V),
              spline_error(ex_cp, ex_knots, ex_P, ex_V))
    if err > 1e-12:
        raise RuntimeError(
            "B-spline 与原复合 Bézier 的偏差 %.2e 过大，导出的几何会和求解器用的"
            "不一致（正常应在 1e-15 量级）。" % err)

    # 样条是夹紧的，首末控制点就是曲线端点，正好是壁的两个端
    tip_in = np.array([case.dofs[0], 0.0])
    tip_ex = np.array([0.0, case.dofs[2 * n_ctrl + 2]])
    if not (np.allclose(in_cp[0], geometry.wg_fixed_bd_in)
            and np.allclose(in_cp[-1], tip_in)
            and np.allclose(ex_cp[0], geometry.wg_fixed_bd_ex)
            and np.allclose(ex_cp[-1], tip_ex)):
        raise RuntimeError("样条端点与几何定义不符：%s %s %s %s"
                           % (in_cp[0], in_cp[-1], ex_cp[0], ex_cp[-1]))

    verts = hexagon_vertices(geometry)
    entities = []

    for sign, port, side in QUADRANTS:
        s = np.asarray(sign, dtype=float)
        # ±1 乘法是精确的，所以镜像出来的共享端点逐位相同
        mouth_in = s * geometry.wg_input_bd_in
        mouth_ex = s * geometry.wg_input_bd_ex
        fixed_in = s * geometry.wg_fixed_bd_in
        fixed_ex = s * geometry.wg_fixed_bd_ex

        # 端口边：整条六边形边被切成 5 段 —— 中间 3 段是端口（共 6w），两端各
        # 留 2w 的端部（不打端口标签）。参考 cross_waveguide_test.java 的 ls1：
        # 端口画到顶点时端点落在 120° 折角上，COMSOL 的一维模式解算会出错。
        a, b = verts[side], verts[(side + 1) % 6]
        chain = _ordered_on_segment(
            a, b, [a, s * geometry.port_edge_bd_in, mouth_in,
                   mouth_ex, s * geometry.port_edge_bd_ex, b])
        for k, (p, q) in enumerate(zip(chain[:-1], chain[1:])):
            layer = "outer" if k in (0, 4) else "port_" + port
            entities.append(Entity(layer, "line", np.array([p, q])))

        # 固定直壁段 + 两条优化壁
        entities.append(Entity("stub", "line", np.array([fixed_in, mouth_in])))
        entities.append(Entity("stub", "line", np.array([fixed_ex, mouth_ex])))
        entities.append(Entity("wall_in", "spline", s * in_cp, in_knots, 3))
        entities.append(Entity("wall_ex", "spline", s * ex_cp, ex_knots, 3))

    # 原点附近连接各象限的同材料内部线段**不放**：它们是波导域内部的分割线，
    # 留着只会把波导切成 4 块（材料要挂 4 次），对物理毫无影响。
    # 不放的话 Form Union 得到 5 个域：1 个波导 + 4 个空气楔。

    # 两条非端口外边界（y = ±0.5）
    for side in SCATTER_SIDES:
        a, b = verts[side], verts[(side + 1) % 6]
        entities.append(Entity("outer", "line", np.array([a, b])))

    return entities


# ---------------------------------------------------------------- DXF 写出


class _DxfWriter:
    """最小的 AC1015（R2000）ASCII DXF 写出器。

    R12（AC1009）没有 SPLINE 实体，所以版本号必须到 R2000；相应地每个实体都要
    带一个唯一句柄（组码 5），否则严格的读取器会拒绝。这里句柄从 0x30 递增。
    """

    def __init__(self, fp):
        self.fp = fp
        self._next = 0x30

    def pair(self, code, value):
        self.fp.write("%d\n%s\n" % (code, value))

    def real(self, code, value):
        # repr(float) 是最短往返表示：读回来与写进去逐位相同
        self.fp.write("%d\n%s\n" % (code, repr(float(value))))

    def handle(self):
        self._next += 1
        return "%X" % self._next


def write_dxf(path, entities, polyline=False, scale=1.0):
    """写 DXF。``polyline=True`` 时把样条按密集折线写（SPLINE 读不动时的退路）。

    ``scale`` 把一个几何单位换算成多少米，必须和 ``model_java()`` 的
    ``length_scale`` 一致，否则 DXF 和 .java 两个产物尺寸对不上。
    """
    if scale != 1.0:
        entities = [Entity(e.layer, e.kind, e.points * scale, e.knots, e.degree)
                    for e in entities]
    layers = sorted({e.layer for e in entities})
    # 包围盒：样条的控制多边形包住曲线本身，用控制点算出来是安全的超集
    allpts = np.vstack([e.points for e in entities])
    lo, hi = allpts.min(axis=0), allpts.max(axis=0)
    with open(path, "w", newline="\n") as fp:
        w = _DxfWriter(fp)

        # ---- HEADER
        w.pair(0, "SECTION")
        w.pair(2, "HEADER")
        w.pair(9, "$ACADVER")
        w.pair(1, "AC1015")
        w.pair(9, "$HANDSEED")
        w.pair(5, "FFFF")
        w.pair(9, "$INSUNITS")
        w.pair(70, 6)               # 6 = 米
        for name, (x, y) in (("$EXTMIN", lo), ("$EXTMAX", hi)):
            w.pair(9, name)
            w.real(10, x)
            w.real(20, y)
            w.real(30, 0.0)
        w.pair(0, "ENDSEC")

        # ---- TABLES（线型 + 图层 + 两个几乎所有写出器都会给的注册表）
        w.pair(0, "SECTION")
        w.pair(2, "TABLES")
        w.pair(0, "TABLE")
        w.pair(2, "LTYPE")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTable")
        w.pair(70, 1)
        w.pair(0, "LTYPE")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTableRecord")
        w.pair(100, "AcDbLinetypeTableRecord")
        w.pair(2, "CONTINUOUS")
        w.pair(70, 0)
        w.pair(3, "Solid line")
        w.pair(72, 65)
        w.pair(73, 0)
        w.real(40, 0.0)
        w.pair(0, "ENDTAB")

        w.pair(0, "TABLE")
        w.pair(2, "LAYER")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTable")
        w.pair(70, len(layers))
        for name in layers:
            w.pair(0, "LAYER")
            w.pair(5, w.handle())
            w.pair(100, "AcDbSymbolTableRecord")
            w.pair(100, "AcDbLayerTableRecord")
            w.pair(2, name)
            w.pair(70, 0)
            w.pair(62, LAYERS[name][0])
            w.pair(6, "CONTINUOUS")
        w.pair(0, "ENDTAB")

        w.pair(0, "TABLE")
        w.pair(2, "STYLE")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTable")
        w.pair(70, 1)
        w.pair(0, "STYLE")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTableRecord")
        w.pair(100, "AcDbTextStyleTableRecord")
        w.pair(2, "STANDARD")
        w.pair(70, 0)
        w.real(40, 0.0)
        w.real(41, 1.0)
        w.real(50, 0.0)
        w.pair(71, 0)
        w.real(42, 2.5)
        w.pair(3, "txt")
        w.pair(4, "")
        w.pair(0, "ENDTAB")

        w.pair(0, "TABLE")
        w.pair(2, "APPID")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTable")
        w.pair(70, 1)
        w.pair(0, "APPID")
        w.pair(5, w.handle())
        w.pair(100, "AcDbSymbolTableRecord")
        w.pair(100, "AcDbRegAppTableRecord")
        w.pair(2, "ACAD")
        w.pair(70, 0)
        w.pair(0, "ENDTAB")
        w.pair(0, "ENDSEC")

        # ---- ENTITIES
        w.pair(0, "SECTION")
        w.pair(2, "ENTITIES")
        for e in entities:
            if e.kind == "line":
                p, q = e.points
                _write_line(w, e.layer, p, q)
            elif polyline:
                # 退路：沿曲线采样成折线（不是控制多边形）
                u = np.linspace(e.knots[3], e.knots[-4], POLYLINE_POINTS)
                samples = evaluate(e, u)
                for p, q in zip(samples[:-1], samples[1:]):
                    _write_line(w, e.layer, p, q)
            else:
                _write_spline(w, e)
        w.pair(0, "ENDSEC")
        w.pair(0, "EOF")


def _write_line(w, layer, p, q):
    w.pair(0, "LINE")
    w.pair(5, w.handle())
    w.pair(100, "AcDbEntity")
    w.pair(8, layer)
    w.pair(100, "AcDbLine")
    w.real(10, p[0])
    w.real(20, p[1])
    w.real(30, 0.0)
    w.real(11, q[0])
    w.real(21, q[1])
    w.real(31, 0.0)


def _write_spline(w, e):
    w.pair(0, "SPLINE")
    w.pair(5, w.handle())
    w.pair(100, "AcDbEntity")
    w.pair(8, e.layer)
    w.pair(100, "AcDbSpline")
    w.real(210, 0.0)                # 法向：平面样条
    w.real(220, 0.0)
    w.real(230, 1.0)
    w.pair(70, 8)                   # 位 3：平面样条
    w.pair(71, int(e.degree))
    w.pair(72, len(e.knots))
    w.pair(73, len(e.points))
    w.pair(74, 0)                   # 不给拟合点
    for k in e.knots:
        w.real(40, k)
    # 权重不写：全为 1 时按规范应当省略（70 里的 Rational 位也没置）
    for p in e.points:
        w.real(10, p[0])
        w.real(20, p[1])
        w.real(30, 0.0)


# ---------------------------------------------------------------- DXF 读回


def read_dxf(path):
    """把 :func:`write_dxf` 写的 DXF 读回来（只为自检，不追求通用）。"""
    with open(path) as fp:
        text = [line.rstrip("\n") for line in fp]
    pairs = [(int(text[i]), text[i + 1]) for i in range(0, len(text) - 1, 2)]

    entities, cur = [], None
    for code, value in pairs:
        if code == 0:
            if cur is not None:
                entities.append(cur)
            cur = {"type": value, "codes": []} if value in ("LINE", "SPLINE") else None
        elif cur is not None:
            cur["codes"].append((code, value))

    out = []
    for raw in entities:
        if raw["type"] == "LINE":
            coords = {code: float(v) for code, v in raw["codes"]
                      if code in (10, 20, 11, 21)}
            layer = [v for c, v in raw["codes"] if c == 8][0]
            out.append(Entity(layer, "line",
                              np.array([[coords[10], coords[20]],
                                        [coords[11], coords[21]]])))
        else:
            code8 = [v for c, v in raw["codes"] if c == 8][0]
            degree = int([v for c, v in raw["codes"] if c == 71][0])
            knots = np.array([float(v) for c, v in raw["codes"] if c == 40])
            xs = [float(v) for c, v in raw["codes"] if c == 10]
            ys = [float(v) for c, v in raw["codes"] if c == 20]
            out.append(Entity(code8, "spline",
                              np.column_stack([xs, ys]), knots, degree))
    return out


def evaluate(entity, u):
    """在参数 ``u`` 处求值（折线用首末点做线性插值，仅用于粗检）。"""
    if entity.kind == "line":
        u = np.atleast_1d(np.asarray(u, dtype=float))
        p, q = entity.points
        return p[None, :] + np.reshape(u, (-1, 1)) * (q - p)[None, :]
    from scipy.interpolate import BSpline

    u = np.atleast_1d(np.asarray(u, dtype=float))
    return np.column_stack([BSpline(entity.knots, entity.points[:, 0], entity.degree)(u),
                            BSpline(entity.knots, entity.points[:, 1], entity.degree)(u)])


# ---------------------------------------------------------------- 校核数字


def geometry_numbers(case, entities):
    """导入 COMSOL 后应当核对到的几何量（全部由解析几何给出，与网格无关）。"""
    geometry = case.geometry
    L = geometry.L_sim
    verts = hexagon_vertices(geometry)

    def arclen(entity, n=20001):
        u = np.linspace(entity.knots[3], entity.knots[-4], n)
        return float(np.sum(np.sqrt(np.sum(np.diff(evaluate(entity, u), axis=0) ** 2,
                                           axis=1))))

    arcs = {}
    for name in ("in", "ex"):
        walls = [e for e in entities if e.layer == "wall_" + name]
        lengths = [arclen(e) for e in walls]
        # 4 个象限是同一个形状的镜像，弧长必须逐位相等
        if len(walls) != 4 or max(lengths) - min(lengths) > 1e-12:
            raise RuntimeError("%s 壁的 4 份镜像弧长不一致：%s" % (name, lengths))
        arcs[name] = lengths[0]

    mouth = float(np.linalg.norm(geometry.wg_input_bd_ex - geometry.wg_input_bd_in))
    stub = float(np.linalg.norm(geometry.wg_fixed_bd_in - geometry.wg_input_bd_in))
    port_len = float(np.linalg.norm(geometry.port_edge_bd_ex - geometry.port_edge_bd_in))
    apothem = float(np.linalg.norm(geometry.wg_input_cen))

    return {
        "L": L,
        "total_area": 1.5 * sqrt(3.0) * L * L,
        "mouth": mouth,
        "stub": stub,
        "port_len": port_len,
        "apothem": apothem,
        "margin": 0.5 * (L - port_len),
        "arc_in": arcs["in"],
        "arc_ex": arcs["ex"],
        "tips": [(case.dofs[0], 0.0), (0.0, case.dofs[2 * case.layout.n_ctrl + 2])],
        "verts": verts,
        "scatter_len": float(np.linalg.norm(verts[1] - verts[2])),
    }


def mesh_numbers(case, arc_in, arc_ex):
    """从 ``data`` 里的网格量出来的数字；网格不存在或与 design.h5 不同步时返回 None。"""
    import meshio

    path = case.mesh_file
    if not os.path.exists(path):
        return None

    mesh = meshio.read(path)
    pts = mesh.points[:, :2]
    tri_idx = next((i for i, c in enumerate(mesh.cells)
                    if c.type in ("triangle6", "triangle")), None)
    line_idx = next((i for i, c in enumerate(mesh.cells)
                     if c.type in ("line3", "line")), None)
    if tri_idx is None or line_idx is None:
        return None
    tri, line = mesh.cells[tri_idx].data, mesh.cells[line_idx].data
    tri_phys = mesh.cell_data["gmsh:physical"][tri_idx]
    line_phys = mesh.cell_data["gmsh:physical"][line_idx]

    def tri_area(t):
        a, b, c = pts[t[:, 0]], pts[t[:, 1]], pts[t[:, 2]]
        return 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                            - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))

    names = {tag: name for name, (tag, _dim) in mesh.field_data.items()}
    areas = {names.get(tag, tag): float(tri_area(tri[tri_phys == tag]).sum())
             for tag in np.unique(tri_phys)}
    lengths = {names.get(tag, tag): float(np.sum(np.sqrt(np.sum(
        (pts[line[line_phys == tag, 1]] - pts[line[line_phys == tag, 0]]) ** 2, axis=1))))
        for tag in np.unique(line_phys)}

    # 弦长比弧长短一点点（0.4% 量级）；差太多说明网格不是按这份 design.h5 建的
    stale = None
    if "in_opt" in lengths:
        stale = abs(lengths["in_opt"] - 4 * arc_in) / (4 * arc_in)
        if stale > 0.02:
            stale = "in_opt 弦长 %.6f vs 4x 弧长 %.6f" % (lengths["in_opt"], 4 * arc_in)
    return {"areas": areas, "lengths": lengths, "stale": stale,
            "n_tri": len(tri), "n_node": len(pts)}


def physics_numbers(case):
    import physics as physics_mod

    ph = case.physics
    beta = physics_mod.beta(ph, case.geometry.w_wg)
    k0 = ph.omega
    return {"beta": beta, "k0": k0, "n_eff": beta / k0,
            "lambda0": 1.0 / ph.freq, "lambda_guided": 2 * pi / beta}


def design_provenance():
    """这份 design.h5 是不是 ``DOF_list.csv`` 的最后一轮（不是的话多半正在跑）。"""
    try:
        arr = np.atleast_2d(np.loadtxt("DOF_list.csv", delimiter=","))
    except (OSError, ValueError):
        return None
    if arr.size == 0:
        return None
    return arr.shape[0], arr[-1, :]


# ---------------------------------------------------------------- 报告


def build_report(case, entities, design_file, dxf_name,
                 java_name="cross_waveguide_optimized.java", length_scale=1.0):
    """生成 ``CHECK.md`` 的内容。"""
    geo = geometry_numbers(case, entities)
    phys = physics_numbers(case)
    mesh = mesh_numbers(case, geo["arc_in"], geo["arc_ex"])

    counts = {}
    for e in entities:
        counts[e.layer] = counts.get(e.layer, 0) + 1

    lines = []
    add = lines.append
    add("# COMSOL 导入校核表\n")
    add("由 `python/export_comsol.py` 从 `%s` 生成。\n" % design_file)
    add("| 文件 | 用途 |")
    add("|---|---|")
    add("| `%s` | **主产物**：COMSOL 模型，几何直接建，跑起来就是完整算例 |" % java_name)
    add("| `%s` | 同一几何的 DXF（想手工搭、或者只想看图时用） |" % dxf_name)
    add("")
    add("两条路用的是同一个几何：DXF 里每条壁写成一个三次 B-spline，与求解器的曲线"
        "偏差 $10^{-16}$ 量级；`.java` 里用 `CubicBezier` 直接填控制多边形，"
        "连这一步转换都没有。\n")

    # ---- 单位
    add("## 1. 单位与频率\n")
    add("跟 `cross_waveguide_test.java` 一样，把量纲收进两个参数：\n")
    add("```text")
    add("L_unit = %s[m]      # 原始无量纲坐标里 “1” 对应多少米" % _num(length_scale))
    add("T_unit = L_unit/c_const")
    add("freq   = f/T_unit          # f = %.6f" % case.physics.freq)
    add("```\n")
    add("默认 `--length-scale 1` → `L_unit = 1[m]`，几何数字原样按米理解："
        "$\\lambda_0$ = %.6f m、$k_0$ = %.6f /m，与 `design.h5` 里的数字逐位相同。"
        "要 $\\lambda_0$ = 1.55 µm 的物理器件就用 `--length-scale 2.4025e-6` 重新"
        "生成（= 1.55/0.6452），坐标和频率一起缩放，比值不变。\n"
        % (phys["lambda0"], phys["k0"]))

    # ---- 几何
    add("## 2. 几何校核（建完立刻核对这几项）\n")
    add("下面的长度都以 `L_unit` 为单位；默认 `L_unit = 1[m]`，所以数字与"
        " `design.h5` 里的一模一样。\n")
    add("| 量 | 期望值 | 说明 |")
    add("|---|---|---|")
    add("| 六边形边长 | `%.9f` | = %.2f w（参考文件的 L/w） |"
        % (geo["L"], geo["L"] / case.geometry.w_wg))
    add("| **每条端口长** | `%.9f` | = %.2f w。端口只占中间 3 段，两端各留 "
        "`%.9f` = %.1f w 的端部 |"
        % (geo["port_len"], geo["port_len"] / case.geometry.w_wg,
           geo["margin"], geo["margin"] / case.geometry.w_wg))
    add("| 端口平面到中心 | `%.9f` | = %.2f w |"
        % (geo["apothem"], geo["apothem"] / case.geometry.w_wg))
    add("| 总面积 | `%.9f` | 必须严格等于 $(3\\sqrt3/2)L^2$；对不上就是单位错了 |"
        % geo["total_area"])
    add("| 域个数 | 5 个：1 个波导 + 4 个空气楔 | 原点附近的同材料分割线没有画，"
        "所以波导是一整块 |")
    add("| 波导口宽度 | `%.9f` | `w_wg`（在端口里居中） |" % geo["mouth"])
    add("| 固定直壁段长 | `%.9f` | `L_wg_fixed`：波导口 → 优化壁的固定端 |" % geo["stub"])
    add("| `in` 壁弧长（每条 / 4 条） | `%.9f` / `%.9f` |"
        % (geo["arc_in"], 4 * geo["arc_in"]))
    add("| `ex` 壁弧长（每条 / 4 条） | `%.9f` / `%.9f` |"
        % (geo["arc_ex"], 4 * geo["arc_ex"]))
    add("| 散射边界 | 2 条整边 + 8 段端部 | 共 10 条，与参考文件 `sctr1` 选中的"
        " 10 条边一一对应 |")
    add("| 4 个尖端 | `(±%.9f, 0)`、`(0, ±%.9f)` | 两条壁在对称轴上汇合处 |"
        % (geo["tips"][0][0], geo["tips"][1][1]))
    add("")
    add("端口边中点（`.java` 里的 Box 选择就是按这些坐标框的）：\n")
    add("| 端点 | 名称 | 中点 x | 中点 y | 对应 |")
    add("|---|---|---|---|---|")
    for k, port in enumerate(PORT_ORDER):
        mid = 0.5 * (geo["verts"][PORT_SIDE[port]] + geo["verts"][(PORT_SIDE[port] + 1) % 6])
        add("| %d | `%s` | `%.9f` | `%.9f` | `S_params[%d]` |"
            % (k + 1, port, mid[0], mid[1], k))
    add("")
    add("六边形顶点（`k` = 0…5，角度 60k°，半径 `%.9f`）：\n" % geo["L"])
    add("| k | x | y |")
    add("|---|---|---|")
    for k, v in enumerate(geo["verts"]):
        add("| %d | `%.9f` | `%.9f` |" % (k, v[0], v[1]))
    add("")

    if mesh is not None:
        add("### 我们网格里的对应值（供交叉检查）\n")
        add("| 量 | 值 |")
        add("|---|---|")
        for name in ("waveguide", "air"):
            if name in mesh["areas"]:
                add("| %s 域面积 | `%.9f` |" % (name, mesh["areas"][name]))
        add("| 单元 / 节点 | %d / %d |" % (mesh["n_tri"], mesh["n_node"]))
        for name in ("input", "output1", "output21", "output22",
                     "in_opt", "ex_opt"):
            if name in mesh["lengths"]:
                add("| boundary `%s` 弦长 | `%.6f` |" % (name, mesh["lengths"][name]))
        add("")
        if isinstance(mesh["stale"], str):
            add("> ⚠️ 网格与这份 `design.h5` **不同步**（%s）。上面的网格数字请忽略，"
                "重新 `make_mesh.py` 再看。\n" % mesh["stale"])

    # ---- 怎么建
    sizes = design.mesh_sizes(case)
    add("## 3. 怎么在 COMSOL 里建\n")
    add("### A. 直接跑生成的模型文件（推荐）\n")
    add("`%s` 与 `cross_waveguide_test.java` 是同一套东西，只把几何从"
        "“两条直矩形交叉”换成了优化后的形状。几何用 `CubicBezier` 直接建，"
        "控制点就是 `Bezier_interp` 用的控制多边形（%d 段 × 8 条壁），所以既没有"
        " DXF 导入的容差，也没有插值曲线的近似。\n" % (java_name, case.layout.n_ctrl + 1))
    add("| 项 | 设置 | 为什么 |")
    add("|---|---|---|")
    add("| 物理场 | `ewfd`，`components = outofplane` | 我们只解 $E_z$ |")
    add("| 端口 | 4 个 `Port`，`PortType = Numeric`，1 个模式 | 我们用的是解析板模，"
        "两边应该在 $n_{eff} = %.6f$ 上对上 |" % phys["n_eff"])
    add("| 非端口外边界 | 2 个 `Scattering`（$y = \\pm 0.5$） | 求解器在这里加的是"
        " $\\gamma = -i k_0 n_{amb}$，**不是 PEC**；`ewfd` 的默认外边界是 PEC"
        "（$E_z = 0$），不显式设就是另一个物理 |")
    add("| 材料 | `RefractiveIndex` 里 `n = sqrt(eps)` | 与参考文件一致"
        "（$\\varepsilon_r$ = %.2f / %.1f） |"
        % (case.physics.epsilon_wg, case.physics.epsilon_amb))
    add("| 网格 | 全局 `hmax` = `%.4f`（空气 λ/16）<br>波导域 `%.4f`<br>端口边 "
        "`%.4f`<br>4 个尖端 `1e-3` | 端口/波导两档是参考文件的做法；尖端那一档"
        "见下面的说明 |"
        % (sizes["h_air"], sizes["h_wg"], sizes["h_port"]))
    add("| 研究 | 4 步 `Boundary Mode Analysis`（每端口一步）+ 1 步频域 | COMSOL 5.5"
        " 的 `Numeric` 端口必须先跑 BMA 才有模式 |")
    add("")
    add("**第一个要核对的数是 Port 节点显示的有效折射率**：应该是 $%.6f$"
        "（我们的解析板模 $\\beta/k_0$）。对上了，端口模式就等价了。\n" % phys["n_eff"])
    add("第二个是功率守恒：`ewfd.Rtotal + ewfd.Ttotal + ewfd.Atotal = 1`。"
        "我们这边没有损耗材料，$A$ 应该几乎全是从两条散射边界漏掉的部分——"
        "它正好能回答下面 §5 里 $\\sum|S|^2 \\ne 1$ 的疑问。\n")
    add("> 尖端那一档网格不是可选项：端点被 $\\theta$ 惩罚推到几乎垂直于对称轴，"
        "尖端曲率半径只有 `3e-3`…`4e-3`（`in` / `ex`），全局 `%.4f` 的网格在那里"
        "会出很差的单元。参考文件没有这个坑，它的几何是两条直矩形。\n"
        % sizes["h_air"])

    add("### B. 纯 GUI（用 DXF）\n")
    add("1. 新建 2D，长度单位 `m`。")
    add("2. Geometry → Import → `%s` → **Convert to Solid** → Form Union。"
        "内部壁的边界要留着（材料分区靠它们）。" % dxf_name)
    add("3. 材料、物理场、网格、研究按上表设。")
    add("4. 端口边界按 §2 的中点坐标选（每条端口边是 3 段共线线段，一个 Port 选 3 条）。")
    add("")

    # ---- 物理参数
    add("## 4. 物理参数\n")
    add("| 量 | 值 | 说明 |")
    add("|---|---|---|")
    add("| `f`（归一化） | `%.6f` | $k_0 = 2\\pi f$ |" % case.physics.freq)
    add("| `eps_wg` / `eps_amb` | `%.6f` / `%.6f` | εr，μr 全为 1 |"
        % (case.physics.epsilon_wg, case.physics.epsilon_amb))
    add("| `mu_wg` / `mu_amb` | `%.6f` / `%.6f` | 非磁性 |"
        % (case.physics.mu_wg, case.physics.mu_amb))
    add("| `w_wg` / `L_sim` / `L_wg_fixed` | `%.6f` / `%.9f` / `%.6f` | "
        "波导宽度 / 六边形外接半径 / 固定直段 |"
        % (case.geometry.w_wg, case.geometry.L_sim, case.geometry.L_wg_fixed))
    add("| $k_0 = 2\\pi f$ | `%.9f` | 单位 1/长度 |" % phys["k0"])
    add("| $\\beta$（解析板模） | `%.9f` | 端口处的传播常数 |" % phys["beta"])
    add("| $n_{eff} = \\beta/k_0$ | `%.9f` | **Port 显示的有效折射率应该等于它** |"
        % phys["n_eff"])
    add("| $\\lambda_{guided} = 2\\pi/\\beta$ | `%.9f` | 网格尺寸的依据 |"
        % phys["lambda_guided"])
    add("")

    # ---- 我们的结果
    result = os.path.join(os.path.dirname(design_file) or ".", case.result_file)
    if os.path.exists(result):
        import h5py

        with h5py.File(result, "r") as fp:
            S = fp["scalars/S_params"][()] if "scalars/S_params" in fp else None
        if S is not None and S.shape == (4, 2):
            add("## 5. 我们这一版算出来的 S 参数\n")
            if os.path.getmtime(result) < os.path.getmtime(design_file):
                add("> ⚠️ `%s` 比 `%s` 旧，说明它算的是**上一轮**的几何。这份数字只能"
                    "当参考，要严格对比就等这一轮跑完再导一次。\n"
                    % (result, design_file))
            add("| 端口 | Re | Im | \\|S\\| |")
            add("|---|---|---|---|")
            for k, (re, im) in enumerate(S):
                add("| %d | `%.9f` | `%.9f` | `%.9f` |"
                    % (k + 1, re, im, sqrt(re * re + im * im)))
            total = float(np.sum(S[:, 0] ** 2 + S[:, 1] ** 2))
            add("")
            add("Σ|S|² = `%.6f`。" % total)
            if abs(total - 1.0) > 0.02:
                add("")
                add("> 无耗结构 + 两条吸收边，这个和必须 ≤ 1，现在偏了 %+.1f%%。"
                    "COMSOL 一搭起来最先能定位的就是这个（它会把 4 个 S 和 β 都独立"
                    "算一遍）。\n" % (100 * (total - 1.0)))
            add("")

    # ---- 图层
    add("## 6. DXF 图层\n")
    add("| 图层 | 实体数 | 颜色 | 说明 |")
    add("|---|---|---|---|")
    for name in sorted(counts):
        color, note = LAYERS[name]
        add("| `%s` | %d | %d | %s |" % (name, counts[name], color, note))
    add("")

    provenance = design_provenance()
    if provenance is not None:
        n_rows, last = provenance
        same = last.shape == case.dofs.shape and np.allclose(last, case.dofs)
        add("## 7. 来源\n")
        add("`DOF_list.csv` 共 %d 轮，最后一轮与这份 `design.h5` %s。"
            % (n_rows, "**一致**" if same else "**不一致**（很可能优化还在跑，"
               "这份是中间快照）"))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- COMSOL 模型

#: COMSOL 里的端口编号（1..4）对应的我们这边的边界名。``ewfd.S<k>1`` 就是
#: ``result.h5`` 里 ``S_params[k-1]``。
PORT_ORDER = ("input", "output1", "output21", "output22")

#: 端口名 -> 六边形边号
PORT_SIDE = {"input": 0, "output22": 2, "output1": 3, "output21": 5}


def _num(v):
    """浮点数写成 17 位往返字符串（COMSOL 读回来与这里逐位相同）。

    ``+ 0.0`` 只是把 ``-0.0`` 归一成 ``0.0``（IEEE 里 (-0)+(+0)=+0），
    其他值加 0 都是精确的。
    """
    return repr(float(v) + 0.0)


def scatter_groups(entities):
    """散射边界覆盖的实体下标：2 条整的非端口边 + 4 条端口边两端的 8 段端部。

    参考 ``cross_waveguide_test.java`` 里 ``sctr1`` 选中的正是 10 条边 —— 2 条
    整边加 8 段端部，与我们这里的 10 条一一对应。整边按"最长"认（整边长 L_sim，
    端部只有 2w），不依赖任何硬编码坐标。
    """
    outer = [i for i, e in enumerate(entities) if e.layer == "outer"]
    if len(outer) != 10:
        raise RuntimeError("散射边界应当有 2 条整边 + 8 段端部，实际 %d 条" % len(outer))
    length = {i: float(np.linalg.norm(entities[i].points[1] - entities[i].points[0]))
              for i in outer}
    full = sorted(outer, key=lambda i: -length[i])[:2]
    if not all(length[i] > 4 * max(length[j] for j in outer if j not in full)
               for i in full):
        raise RuntimeError("分不清哪两条是整边：%s" % {i: length[i] for i in outer})

    groups = {}
    for i in full:
        name = "top" if entities[i].points[:, 1].mean() > 0 else "bottom"
        corners = {tuple(p) for p in entities[i].points}
        groups[name] = sorted([i] + [j for j in outer if j != i
                                     and any(tuple(p) in corners
                                             for p in entities[j].points)])
    return groups


def selection_boxes(case, entities):
    """所有按坐标定义的 Box 选择 ``[(tag, entitydim, box, note)]``。

    用坐标而不是实体编号：编号是 COMSOL 内部的对象顺序，几何一变就全错；坐标是
    几何本身的性质。边界和点用 ``condition = "inside"``（实体整体落在框内），
    所以相邻壁蹭到框边的端点不会被误选；域用 ``"intersects"``（框落在域内部）。
    """
    g = case.geometry
    L_sim = g.L_sim
    pad = 0.02 * L_sim
    boxes = []

    for port in PORT_ORDER:
        pts = np.vstack([e.points for e in entities if e.layer == "port_" + port])
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        boxes.append(("sel_" + port, 1,
                      (lo[0] - pad, hi[0] + pad, lo[1] - pad, hi[1] + pad),
                      "端口 %s：3 条共线边，中点 (%+.4f, %+.4f)"
                      % (port, (lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2)))

    # 散射边界 = 2 条整的非端口边 + 4 条端口边两端的 8 段端部（参考文件里正好
    # 也是 10 条边）。一条整边的框顺带把贴在它两个端点上的端部一起框进去，
    # 所以只要 2 个框。
    y_top = L_sim * np.sin(np.pi / 3)
    for name, group in scatter_groups(entities).items():
        pts = np.vstack([entities[i].points for i in group])
        boxes.append(("sel_scatter_" + name, 1,
                      (pts[:, 0].min() - pad, pts[:, 0].max() + pad,
                       pts[:, 1].min() - pad, pts[:, 1].max() + pad),
                      "散射边界：y = %s%.4f 那条整边 + 它两端相邻的 2 段端部"
                      % ("+" if name == "top" else "-", y_top)))

    # 波导域：两条固定直段中点连线的中点是臂的正中心，一定在域内部
    center = 0.5 * (g.wg_fixed_bd_in + g.wg_fixed_bd_ex)
    half = 0.25 * g.w_wg
    boxes.append(("sel_wg", 2,
                  (center[0] - half, center[0] + half,
                   center[1] - half, center[1] + half),
                  "波导域（框落在直段内部）"))

    # 4 个空气楔：分别朝 0° / 90° / 180° / 270°
    y_top = L_sim * np.sin(np.pi / 3)
    half = 0.25 * g.w_wg
    for k, (x, y, where) in enumerate((
            (0.6 * L_sim, 0.0, "0°（两臂之间，六边形顶点侧）"),
            (0.0, 0.8 * y_top, "90°（贴上外边界）"),
            (-0.6 * L_sim, 0.0, "180°"),
            (0.0, -0.8 * y_top, "270°（贴下外边界）"))):
        boxes.append(("sel_air%d" % (k + 1), 2,
                      (x - half, x + half, y - half, y + half),
                      "空气楔 %s" % where))

    tip_in, tip_ex = case.dofs[0], case.dofs[2 * case.layout.n_ctrl + 2]
    half = 0.05 * min(tip_in, tip_ex)
    for k, (x, y) in enumerate(((tip_in, 0.0), (-tip_in, 0.0),
                                (0.0, tip_ex), (0.0, -tip_ex))):
        boxes.append(("sel_tip%d" % (k + 1), 0,
                      (x - half, x + half, y - half, y + half),
                      "尖端 %d：(%+.6f, %+.6f)" % (k + 1, x, y)))
    return boxes


def _vertices(entities):
    """几何里所有顶点（曲线端点去重）。"""
    pts = np.vstack([e.points for e in entities])
    return np.unique(np.round(pts, 12), axis=0)


def _dense_points(entity, n=401):
    """实体上足够密的采样点（用来判断一条曲线有没有穿过某个框）。"""
    return evaluate(entity, np.linspace(0.0, 1.0, n))


def check_selections(case, entities, boxes=None):
    """核对每个 Box 选择框住了什么（离线可做，不需要 COMSOL）。

    COMSOL 的 ``condition = "inside"`` 判的是"实体整体落在框内"；对轴对齐的框，
    这等价于"实体的包围盒落在框内"（包围盒是包含该实体的最小轴对齐矩形，两个
    方向都成立），所以这里用包围盒模拟是与 COMSOL 完全一致的判断。
    返回 ``{tag: 命中项}``，曲线命中的是下标、点命中的是坐标。
    """
    boxes = selection_boxes(case, entities) if boxes is None else boxes
    spans = [(e.points.min(axis=0), e.points.max(axis=0)) for e in entities]
    verts = _vertices(entities)
    hits = {}
    for tag, dim, (x0, x1, y0, y1), _note in boxes:
        if dim == 0:
            inside = ((verts[:, 0] >= x0) & (verts[:, 0] <= x1)
                      & (verts[:, 1] >= y0) & (verts[:, 1] <= y1))
            hits[tag] = [(float(v[0]), float(v[1])) for v in verts[inside]]
        else:
            hits[tag] = [i for i, (lo, hi) in enumerate(spans)
                         if lo[0] >= x0 and hi[0] <= x1
                         and lo[1] >= y0 and hi[1] <= y1]
    return hits


def verify_selections(case, entities):
    """选择框是否恰好框住该框的东西；返回问题列表（空 = 全部正确）。"""
    hits = check_selections(case, entities)
    layers = [e.layer for e in entities]
    problems = []

    for port in PORT_ORDER:
        want = [i for i, name in enumerate(layers) if name == "port_" + port]
        got = hits["sel_" + port]
        if got != want:
            problems.append("sel_%s 命中 %s，应当恰好是端口边的 %s" % (port, got, want))

    for name, want in scatter_groups(entities).items():
        got = sorted(hits["sel_scatter_" + name])
        if got != want:
            problems.append("sel_scatter_%s 命中 %s，应当是 %s" % (name, got, want))

    # 域用的 "intersects"：框必须完全落在某一个域内部 —— 既不能框住整条曲线
    # （上面已查），也不能被任何曲线穿过（否则会同时选中两个域）
    samples = np.vstack([_dense_points(e) for e in entities])
    for tag in ("sel_wg", "sel_air1", "sel_air2", "sel_air3", "sel_air4"):
        x0, x1, y0, y1 = [b[2] for b in selection_boxes(case, entities)
                          if b[0] == tag][0]
        if hits[tag]:
            problems.append("%s 本应在某个域内部，却整条框住了曲线 %s"
                            % (tag, hits[tag]))
        crossed = (((samples[:, 0] >= x0) & (samples[:, 0] <= x1)
                    & (samples[:, 1] >= y0) & (samples[:, 1] <= y1)).any())
        if crossed:
            problems.append("%s 被某条曲线穿过，会同时选中两个域" % tag)

    for k in range(1, 5):
        got = hits["sel_tip%d" % k]
        if len(got) != 1:
            problems.append("sel_tip%d 命中 %d 个点（应当恰好 1 个）：%s"
                            % (k, len(got), got))
    return problems


def model_java(case, entities, src_label, length_scale=1.0):
    """生成 COMSOL 5.5 的 ``.java`` 模型文件。

    几何直接建：优化壁用 ``CubicBezier``，它的 ``p`` 属性就是四个控制点，与
    ``design.bezier_polygons()`` 逐位相同，所以既没有 DXF 导入的容差问题，也没有
    插值曲线的近似。物理 / 网格 / 研究 / 结果照搬 ``cross_waveguide_test.java``。
    """
    S = float(length_scale)
    g, layout = case.geometry, case.layout
    in_polys, ex_polys = design.bezier_polygons(case.dofs, g, layout)
    n_seg = len(in_polys)
    phys = physics_numbers(case)
    sizes = design.mesh_sizes(case)
    verts = hexagon_vertices(g)

    m = 'model.component("comp1").geom("geom1")'
    L = []
    ln = L.append

    # ---- 文件头
    ln("/* ============================================================================")
    ln(" * cross_waveguide_optimized.java")
    ln(" *")
    ln(" *   ⚠ 由 python/export_comsol.py 从 %s 拼字符串生成，不是 COMSOL 的" % src_label)
    ln(" *     “Save as Model File for Java” 导出件（那份带一行 Model exported on ...")
    ln(" *     by COMSOL x.y 的戳，像 cross_waveguide_test.java）。语句照抄了参考文件")
    ln(" *     能抄的部分；其余（CubicBezier / Box 选择 / named 选择 / study().run）")
    ln(" *     是按 COMSOL 文档和 API 约定推的，没编译运行过。换设计请重新生成。")
    ln(" *")
    ln(" * 与 cross_waveguide_test.java 的关系：物理、网格、研究、结果那一套照搬")
    ln(" * （ewfd + outofplane；4 个 Numeric 端口 + 4 步 Boundary Mode Analysis；")
    ln(" * 非端口外边界用 Scattering 吸收；材料用 RefractiveIndex 填 n = sqrt(eps)），")
    ln(" * 只把几何从“两条直矩形交叉”换成优化后的形状：")
    ln(" *")
    ln(" *   正六边形外框                    1 个 Polygon")
    ln(" *   优化壁 in/ex x 4 象限 x %-2d 段  %d 个 CubicBezier" % (n_seg, 8 * n_seg))
    ln(" *   固定直壁段                      8 个 LineSegment")
    ln(" *")
    ln(" * 控制点就是求解器里 Bezier_interp 用的控制多边形，逐位相同，不做插值。")
    ln(" *")
    ln(" * 单位：几何数字的单位是模型长度单位（默认 m）。L_unit 是“原始无量纲坐标里")
    ln(" * 1 对应多少米”，只出现在频率表达式里；要换物理尺度用")
    ln(" *     python3 ../python/export_comsol.py --length-scale <s>")
    ln(" * 重新生成（坐标乘 s、L_unit 同时变成 s[m]），这样几何和频率一起缩放。")
    ln(" *")
    ln(" * 选择：端口 / 散射边界 / 材料域 / 尖端全部按坐标的 Box 选择，不写死实体编号")
    ln(" * （编号随几何变化，坐标不会）。哪个框没框住，点开对应的 Selection 就看得见。")
    ln(" *")
    ln(" * 端口编号 -> 我们的 S 参数（ewfd.S<k>1 = S_params[k-1]）：")
    for k, port in enumerate(PORT_ORDER):
        side = PORT_SIDE[port]
        mid = 0.5 * (verts[side] + verts[(side + 1) % 6])
        ln(" *   port%d %-9s 中点 (%+.4f, %+.4f)  <-> S_params[%d]"
           % (k + 1, port, mid[0], mid[1], k))
    ln(" * ============================================================================ */")
    ln("")
    ln("import com.comsol.model.*;")
    ln("import com.comsol.model.util.*;")
    ln("")
    ln("public class cross_waveguide_optimized {")
    ln("")
    ln("  public static Model run() {")
    ln('    Model model = ModelUtil.create("Model");')
    ln('    model.label("cross_waveguide_optimized.mph");')
    ln("")
    ln("    // ---------------- 参数 ----------------")
    ln('    model.param().set("L_unit", "%s[m]");' % _num(S))
    ln('    model.param().set("T_unit", "L_unit/c_const");')
    ln('    model.param().set("P_unit", "1[W]");')
    ln('    model.param().set("E_unit", "P_unit*T_unit");')
    ln('    model.param().set("S_unit", "P_unit/L_unit^2");')
    ln('    model.param().set("f", "%s");' % _num(case.physics.freq))
    ln('    model.param().set("epsilon_wg", "%s");' % _num(case.physics.epsilon_wg))
    ln('    model.param().set("epsilon_amb", "%s");' % _num(case.physics.epsilon_amb))
    ln("")
    ln("    // ---------------- 元件与几何 ----------------")
    ln('    model.component().create("comp1", true);')
    ln('    model.component("comp1").geom().create("geom1", 2);')
    ln("")
    ln("    // 外框：正六边形，顶点在半径 %s 处、角度 60k 度" % _num(g.L_sim * S))
    ln('    %s.create("pol1", "Polygon");' % m)
    ln('    %s.feature("pol1").set("x", "%s*cos(2*pi/6*range(0,1,6))");'
       % (m, _num(g.L_sim * S)))
    ln('    %s.feature("pol1").set("y", "%s*sin(2*pi/6*range(0,1,6))");'
       % (m, _num(g.L_sim * S)))
    ln("")

    signs = [np.asarray(q[0], dtype=float) for q in QUADRANTS]
    ports = [q[1] for q in QUADRANTS]
    tag = 0
    for name, polys in (("in", in_polys), ("ex", ex_polys)):
        ln("    // 优化壁 %s：4 个象限 x %d 段三次 Bézier（控制点 = 控制多边形）"
           % (name, n_seg))
        for sign, port in zip(signs, ports):
            for i, poly in enumerate(polys):
                tag += 1
                pts = poly * sign * S
                ln('    %s.create("cb%d", "CubicBezier");   // %s / %s / 第 %d 段'
                   % (m, tag, name, port, i + 1))
                for j, p in enumerate(pts):
                    ln('    %s.feature("cb%d").setIndex("p", %s, 0, %d);'
                       % (m, tag, _num(p[0]), j))
                    ln('    %s.feature("cb%d").setIndex("p", %s, 1, %d);'
                       % (m, tag, _num(p[1]), j))
        ln("")

    ltag = 0
    ln("    // 固定直壁段（长 %s）：波导口 -> 壁的固定端" % _num(g.L_wg_fixed * S))
    for sign in signs:
        for fixed, mouth in ((g.wg_fixed_bd_in, g.wg_input_bd_in),
                             (g.wg_fixed_bd_ex, g.wg_input_bd_ex)):
            ltag += 1
            a, b = fixed * sign * S, mouth * sign * S
            _java_line(m, ltag, a, b, ln)
    ln("")
    ln("    // 端口边界的印记线：与六边形边共线，把每条边切成 5 段（只有中间 3 段"
       "是端口）。")
    ln("    // 参考文件的 ls1 干的正是这件事。不加的话六边形边只会被波导口切成 3 段，")
    ln("    // Box 就框不出“中间 6w”这个集合了。")
    for sign, port, side in QUADRANTS:
        ltag += 1
        a = np.asarray(sign, dtype=float) * g.port_edge_bd_in * S
        b = np.asarray(sign, dtype=float) * g.port_edge_bd_ex * S
        _java_line(m, ltag, a, b, ln, note=port)
    ln("")
    ln("    // Form Union：这些曲线把六边形切成 5 个域（1 个波导 + 4 个空气楔）。")
    ln("    // 原点附近那 4 条同材料分割线故意不放 —— 它们会把波导切成 4 块，"
       "对物理没影响。")
    ln('    %s.run();' % m)
    ln("")

    ln("    // ---------------- 按坐标的 Box 选择 ----------------")
    for tag_, dim, (x0, x1, y0, y1), note in selection_boxes(case, entities):
        cond = "intersects" if dim == 2 else "inside"
        ln("    // %s" % note)
        ln('    model.component("comp1").selection().create("%s", "Box");' % tag_)
        ln('    model.component("comp1").selection("%s").set("entitydim", %d);'
           % (tag_, dim))
        ln('    model.component("comp1").selection("%s").set("xmin", %s);'
           % (tag_, _num(x0 * S)))
        ln('    model.component("comp1").selection("%s").set("xmax", %s);'
           % (tag_, _num(x1 * S)))
        ln('    model.component("comp1").selection("%s").set("ymin", %s);'
           % (tag_, _num(y0 * S)))
        ln('    model.component("comp1").selection("%s").set("ymax", %s);'
           % (tag_, _num(y1 * S)))
        ln('    model.component("comp1").selection("%s").set("condition", "%s");'
           % (tag_, cond))
    ln("")

    ln("    // ---------------- 材料（n = sqrt(eps)，与参考文件同一写法）"
       " ----------------")
    for mat, sel, eps in ([("mat_air%d" % k, "sel_air%d" % k, "epsilon_amb")
                           for k in range(1, 5)]
                          + [("mat_wg", "sel_wg", "epsilon_wg")]):
        _java_material(mat, sel, eps, ln)
    ln("")

    ln("    // ---------------- 物理场：ewfd + out-of-plane Ez ----------------")
    ln('    model.component("comp1").physics().create("ewfd", '
       '"ElectromagneticWavesFrequencyDomain", "geom1");')
    ln('    model.component("comp1").physics("ewfd").prop("components")'
       '.set("components", "outofplane");')
    for k, port in enumerate(PORT_ORDER):
        ln('    model.component("comp1").physics("ewfd").create("port%d", "Port", 1);'
           % (k + 1))
        ln('    model.component("comp1").physics("ewfd").feature("port%d")'
           '.selection().named("sel_%s");' % (k + 1, port))
        ln('    model.component("comp1").physics("ewfd").feature("port%d")'
           '.set("PortType", "Numeric");' % (k + 1))
        ln('    model.component("comp1").physics("ewfd").feature("port%d")'
           '.set("PortFreq", "f/T_unit");' % (k + 1))
    ln("    // 两条非端口外边界：散射边界条件。绝不能留默认 —— ewfd 的默认外边界")
    ln("    // 是 PEC（Ez = 0），那是完全不同的物理。")
    for name in ("top", "bottom"):
        ln('    model.component("comp1").physics("ewfd").create("sctr_%s", '
           '"Scattering", 1);' % name)
        ln('    model.component("comp1").physics("ewfd").feature("sctr_%s")'
           '.selection().named("sel_scatter_%s");' % (name, name))
    ln("")

    ln("    // ---------------- 网格 ----------------")
    ln("    // 全局 hmax = 空气里 λ/16，波导域再细一档，端口边按 h_port，"
       "4 个尖端单独压到 1e-3")
    ln('    model.component("comp1").mesh().create("mesh1");')
    ln('    model.component("comp1").mesh("mesh1").feature("size")'
       '.set("custom", "on");')
    ln('    model.component("comp1").mesh("mesh1").feature("size")'
       '.set("hmax", "%s");' % _num(sizes["h_air"] * S))
    ln('    model.component("comp1").mesh("mesh1").feature("size")'
       '.set("hmin", "%s");' % _num(0.005 * g.w_wg * S))
    mg = 'model.component("comp1").mesh("mesh1")'
    size_tag = 0
    for sel, hmax, note in (
            [("sel_wg", sizes["h_wg"], "波导域")]
            + [("sel_" + p, sizes["h_port"], "端口 %s" % p) for p in PORT_ORDER]
            + [("sel_tip%d" % k, 1e-3, "尖端 %d" % k) for k in range(1, 5)]):
        size_tag += 1
        ln("    // %s：hmax = %s" % (note, _num(hmax * S)))
        ln('    %s.create("size%d", "Size");' % (mg, size_tag))
        ln('    %s.feature("size%d").selection().named("%s");'
           % (mg, size_tag, sel))
        ln('    %s.feature("size%d").set("custom", "on");' % (mg, size_tag))
        ln('    %s.feature("size%d").set("hmax", "%s");'
           % (mg, size_tag, _num(hmax * S)))
        ln('    %s.feature("size%d").set("hmaxactive", true);' % (mg, size_tag))
    ln('    %s.create("ftri1", "FreeTri");' % mg)
    ln('    %s.run();' % mg)
    ln("")

    ln("    // ---------------- 研究：4 步边值模分析 + 1 步频域 ----------------")
    ln('    model.study().create("std1");')
    ln('    model.study("std1").create("bma", "BoundaryModeAnalysis");')
    for k in range(1, 4):
        ln('    model.study("std1").create("bma%d", "BoundaryModeAnalysis");' % k)
    ln('    model.study("std1").create("freq", "Frequency");')
    ln("    // 端口的数值模式：在 f 处求，shift 用解析板模的有效折射率当起点")
    ln('    model.study("std1").feature("bma").set("modeFreq", "f/T_unit");')
    ln('    model.study("std1").feature("bma").set("shift", "%s");'
       % _num(round(phys["n_eff"], 3)))
    for k in range(1, 4):
        ln('    model.study("std1").feature("bma%d").set("PortName", "%d");' % (k, k + 1))
        ln('    model.study("std1").feature("bma%d").set("modeFreq", "f/T_unit");' % k)
        ln('    model.study("std1").feature("bma%d").set("shift", "%s");'
           % (k, _num(round(phys["n_eff"], 3))))
    ln('    model.study("std1").feature("freq").set("plist", "f/T_unit");')
    ln('    model.study("std1").feature("freq").set("ngen", 1);')
    ln('    model.study("std1").run();')
    ln("")

    ln("    // ---------------- 结果 ----------------")
    ln("    // gev1：R / T / A 功率守恒；gev2：S11..S41（复数值）")
    ln('    model.result().table().create("tbl1", "Table");')
    ln('    model.result().table().create("tbl2", "Table");')
    ln('    model.result().numerical().create("gev1", "EvalGlobal");')
    ln('    model.result().numerical().create("gev2", "EvalGlobal");')
    ln('    model.result().numerical("gev1").set("probetag", "none");')
    ln('    model.result().numerical("gev1").set("table", "tbl1");')
    ln('    model.result().numerical("gev1").set("expr", new String[]{"ewfd.Rport_1", '
       '"ewfd.Rtotal", "ewfd.Tport_2", "ewfd.Tport_3", "ewfd.Tport_4", '
       '"ewfd.Ttotal", "ewfd.RTtotal", "ewfd.Atotal"});')
    ln('    model.result().numerical("gev1").set("descr", new String[]{"R (port 1)", '
       '"R total", "T (port 2)", "T (port 3)", "T (port 4)", "T total", '
       '"R+T", "A (= 从散射边界漏掉的)"});')
    ln('    model.result().numerical("gev2").set("probetag", "none");')
    ln('    model.result().numerical("gev2").set("table", "tbl2");')
    ln('    model.result().numerical("gev2").set("expr", new String[]{"ewfd.S11", '
       '"ewfd.S21", "ewfd.S31", "ewfd.S41"});')
    ln('    model.result().numerical("gev2").set("descr", new String[]{"S11", "S21", '
       '"S31", "S41"});')
    ln('    model.result().numerical("gev1").setResult();')
    ln('    model.result().numerical("gev2").setResult();')
    ln('    model.result().create("pg1", "PlotGroup2D");')
    ln('    model.result("pg1").create("surf1", "Surface");')
    ln('    model.result("pg1").feature("surf1").set("expr", "ewfd.normE");')
    ln('    model.result("pg1").feature("surf1").set("colortable", "RainbowLight");')
    ln("")
    ln("    return model;")
    ln("  }")
    ln("")
    ln("  public static void main(String[] args) {")
    ln("    run();")
    ln("  }")
    ln("}")
    return "\n".join(L) + "\n"


def _java_line(m, tag, a, b, ln, note=""):
    """一条 LineSegment（坐标用 String[]，和参考文件同一写法）。"""
    ln('    %s.create("ls%d", "LineSegment");%s'
       % (m, tag, ("   // 端口印记线 " + note) if note else ""))
    ln('    %s.feature("ls%d").set("specify1", "coord");' % (m, tag))
    ln('    %s.feature("ls%d").set("coord1", new String[]{"%s", "%s"});'
       % (m, tag, _num(a[0]), _num(a[1])))
    ln('    %s.feature("ls%d").set("specify2", "coord");' % (m, tag))
    ln('    %s.feature("ls%d").set("coord2", new String[]{"%s", "%s"});'
       % (m, tag, _num(b[0]), _num(b[1])))


def _java_material(tag, sel, eps, ln):
    """一个 Common 材料：折射率 n = sqrt(eps)（与参考文件同一写法）。"""
    c = 'model.component("comp1").material("%s")' % tag
    ln('    model.component("comp1").material().create("%s", "Common");' % tag)
    ln('    %s.selection().named("%s");' % (c, sel))
    ln('    %s.propertyGroup().create("RefractiveIndex", "Refractive index");' % c)
    ln('    %s.propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(%s)", '
       '"0", "0", "0", "sqrt(%s)", "0", "0", "0", "sqrt(%s)"});' % (c, eps, eps, eps))
    ln('    %s.propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", '
       '"0", "0", "0", "0", "0", "0", "0"});' % c)


# ---------------------------------------------------------------- 入口


def main(design_file="design.h5", out_dir="comsol", polyline=False,
         length_scale=1.0):
    paths.chdir_data()
    case = design.Case.load(design_file)
    entities = build_entities(case)

    # 选择框必须框住的正好是该框的东西 —— 这条离线就能严格验证（见
    # verify_selections 的注释），比到 COMSOL 里才发现材料挂错域强得多。
    problems = verify_selections(case, entities)
    if problems:
        raise RuntimeError("Box 选择框有问题：\n  " + "\n  ".join(problems))

    os.makedirs(out_dir, exist_ok=True)
    dxf_path = os.path.join(out_dir, "geometry.dxf")
    write_dxf(dxf_path, entities, polyline=polyline, scale=length_scale)

    java_path = os.path.join(out_dir, "cross_waveguide_optimized.java")
    with open(java_path, "w") as fp:
        fp.write(model_java(case, entities, design_file, length_scale))

    check_path = os.path.join(out_dir, "CHECK.md")
    report = build_report(case, entities, design_file,
                          os.path.basename(dxf_path),
                          os.path.basename(java_path), length_scale)
    with open(check_path, "w") as fp:
        fp.write(report)

    geo = geometry_numbers(case, entities)
    counts = {}
    for e in entities:
        counts[e.layer] = counts.get(e.layer, 0) + 1
    print("DXF:   %s（%d 个实体: %s）"
          % (dxf_path, len(entities),
             ", ".join("%s=%d" % kv for kv in sorted(counts.items()))))
    print("模型:  %s（%d 个几何特征 + %d 个 Box 选择）"
          % (java_path, 1 + 8 * (case.layout.n_ctrl + 1) + 12,
             len(selection_boxes(case, entities))))
    print("校核表: %s" % check_path)
    print("六边形边长 %.9f（= %.2f w）  端口长 %.9f（= %.2f w）  总面积 %.9f"
          % (geo["L"], geo["L"] / case.geometry.w_wg,
             geo["port_len"], geo["port_len"] / case.geometry.w_wg, geo["total_area"]))
    print("in/ex 壁弧长 %.9f / %.9f" % (geo["arc_in"], geo["arc_ex"]))
    return dxf_path, java_path, check_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--design", default="design.h5",
                        help="设计文件（默认 design.h5，相对 data/）")
    parser.add_argument("--out", default="comsol", help="输出目录（默认 data/comsol）")
    parser.add_argument("--polyline", action="store_true",
                        help="DXF 里的样条写成密集折线（SPLINE 读不动时的退路）")
    parser.add_argument("--length-scale", type=float, default=1.0,
                        help="把一个几何单位缩放成多少米（默认 1，即几何数字按米理解；"
                             "要 λ0 = 1.55 µm 用 1.55/0.645161 = 2.4025e-6）")
    args = parser.parse_args()
    main(design_file=args.design, out_dir=args.out, polyline=args.polyline,
         length_scale=args.length_scale)
