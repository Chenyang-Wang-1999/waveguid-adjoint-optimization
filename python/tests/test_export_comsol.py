"""COMSOL 导出的 DXF：样条等价性、往返、实体拓扑与几何不变量。

这里的用例都用**自己构造**的一版设计（确定性扰动，不读 ``data/``），所以
``data/design.h5`` 被优化器重写时测试也不会飘。
"""

import re
from collections import Counter

import numpy as np
import pytest

import Bezier
import design
import export_comsol as ec


@pytest.fixture(scope="module")
def case():
    """一版够弯、但不越轴的设计。"""
    geo, layout = design.Geometry(), design.Layout(10)
    dofs = design.initial_dofs(geo, layout)
    rng = np.random.default_rng(0)
    dofs[2:2 * layout.n_ctrl + 2] += rng.normal(scale=3e-3, size=2 * layout.n_ctrl)
    dofs[2 * layout.n_ctrl + 4:] += rng.normal(scale=3e-3, size=2 * layout.n_ctrl)
    crossings = design.axis_crossings(dofs, geo, layout)
    assert all(v is None for v in crossings.values()), crossings
    return design.make_case(dofs, build=False)


@pytest.fixture(scope="module")
def entities(case):
    return ec.build_entities(case)


def test_wall_spline_is_exact(case):
    """一条复合三次 Bézier 必须能被一条带节点向量的 B-spline 逐位表示。"""
    geo, layout = case.geometry, case.layout
    in_P, in_V, ex_P, ex_V = design.bezier_curves(case.dofs, geo, layout)
    n_seg = len(in_P) - 1

    for P, V in ((in_P, in_V), (ex_P, ex_V)):
        control, knots = ec.wall_spline(P, V)

        # 两端四重节点 + 中间 (n-1) 个单重节点，共 n+3 个控制点
        assert len(knots) == 4 + (n_seg - 1) + 4
        assert np.array_equal(knots[:4], np.zeros(4))
        assert np.array_equal(knots[-4:], np.full(4, float(n_seg)))
        assert np.array_equal(knots[4:-4], np.arange(1.0, float(n_seg)))
        assert len(control) == n_seg + 3

        # 夹紧样条：首末控制点就是曲线的两个端
        assert np.allclose(control[0], P[0])
        assert np.allclose(control[-1], P[-1])

        assert ec.spline_error(control, knots, P, V) < 1e-12


def test_entity_topology(case, entities):
    """实体集合必须和 `build_mesh` 建的曲线一一对应。"""
    counts = Counter(e.layer for e in entities)
    assert counts == {"port_input": 3, "port_output22": 3, "port_output1": 3,
                      "port_output21": 3, "outer": 10, "stub": 8,
                      "wall_in": 4, "wall_ex": 4}
    assert sum(1 for e in entities if e.kind == "spline") == 8

    geo = case.geometry
    for port in ("input", "output22", "output1", "output21"):
        segs = [e for e in entities if e.layer == "port_" + port]
        lens = [float(np.linalg.norm(s.points[1] - s.points[0])) for s in segs]
        # 端口 = 中间 3 段：口径 w 居中，两侧各 (6w-w)/2 = 2.5w
        assert sum(lens) == pytest.approx(6.0 * geo.w_wg, abs=1e-12)
        assert lens[1] == pytest.approx(geo.w_wg, abs=1e-12)
        assert lens[0] == pytest.approx(lens[2], abs=1e-12)
        assert lens[0] == pytest.approx(2.5 * geo.w_wg, abs=1e-12)

        unit = segs[1].points[1] - segs[1].points[0]
        unit /= np.linalg.norm(unit)
        for s in segs:
            d = s.points[1] - s.points[0]
            assert abs(np.cross(unit, d)) < 1e-12 * np.linalg.norm(d)

    # 散射边界的 10 条：2 条整边（长 L_sim）+ 8 段端部（各 2w）
    outer = [e for e in entities if e.layer == "outer"]
    lengths = [float(np.linalg.norm(e.points[1] - e.points[0])) for e in outer]
    assert sum(1 for x in lengths if x == pytest.approx(geo.L_sim)) == 2
    assert sum(1 for x in lengths
               if x == pytest.approx(2.0 * geo.w_wg)) == 8


def test_geometry_numbers(case, entities):
    numbers = ec.geometry_numbers(case, entities)
    geo = case.geometry
    assert numbers["L"] == pytest.approx(geo.L_sim)
    assert numbers["total_area"] == pytest.approx(1.5 * np.sqrt(3.0) * geo.L_sim ** 2)
    assert numbers["mouth"] == pytest.approx(geo.w_wg)
    assert numbers["stub"] == pytest.approx(geo.L_wg_fixed)
    # 尖端落在对称轴上
    assert numbers["tips"][0][1] == 0.0
    assert numbers["tips"][1][0] == 0.0


def test_dxf_roundtrip(case, entities, tmp_path):
    """写出去再读回来，8 条壁必须与求解器用的曲线逐点一致。"""
    path = str(tmp_path / "geometry.dxf")
    ec.write_dxf(path, entities)
    back = ec.read_dxf(path)
    assert len(back) == len(entities)
    assert Counter(e.layer for e in back) == Counter(e.layer for e in entities)

    geo, layout = case.geometry, case.layout
    in_P, in_V, ex_P, ex_V = design.bezier_curves(case.dofs, geo, layout)
    t = np.linspace(0.0, len(in_P) - 1, 2001)

    for name, (P, V) in (("in", (in_P, in_V)), ("ex", (ex_P, ex_V))):
        walls = [e for e in back if e.layer == "wall_" + name]
        assert len(walls) == 4
        ref0 = Bezier.Bezier_interp(P, V, t.copy(), 0)
        for (sign, _port, _side), wall in zip(ec.QUADRANTS, walls):
            ref = ref0 * np.asarray(sign, dtype=float)
            assert np.abs(ec.evaluate(wall, t) - ref).max() < 1e-12


def test_dxf_header(entities, tmp_path):
    path = str(tmp_path / "geometry.dxf")
    ec.write_dxf(path, entities)
    text = open(path).read().splitlines()

    assert "AC1015" in text
    splines = 0
    for i, line in enumerate(text):
        if line != "AcDbSpline":
            continue
        splines += 1
        codes = text[i + 1:i + 17]
        assert (codes[6], codes[7]) == ("70", "8")        # 平面样条
        assert (codes[9], codes[11], codes[13]) == ("3", "18", "14")
        assert (codes[14], codes[15]) == ("74", "0")       # 不给拟合点
    assert splines == 8
    # 权重全为 1：按规范应当整组省略 → SPLINE 实体里不该出现组码 41
    for i, line in enumerate(text):
        if line != "AcDbSpline":
            continue
        block = []
        j = i + 1
        while j + 1 < len(text) and text[j] != "0":
            block.append(text[j])
            j += 2
        assert "41" not in block
        assert block.count("40") == 18


def test_polyline_fallback(case, entities, tmp_path):
    """退路：沿曲线采样成折线，不能拿控制多边形冒充曲线。"""
    path = str(tmp_path / "poly.dxf")
    ec.write_dxf(path, entities, polyline=True)
    back = ec.read_dxf(path)
    assert all(e.kind == "line" for e in back)
    assert len(back) == 30 + 8 * (ec.POLYLINE_POINTS - 1)

    geo, layout = case.geometry, case.layout
    in_P, in_V, _, _ = design.bezier_curves(case.dofs, geo, layout)
    chain = [e for e in back if e.layer == "wall_in"][:ec.POLYLINE_POINTS - 1]
    pts = np.vstack([chain[0].points[0], np.array([e.points[1] for e in chain])])
    t = np.linspace(0.0, len(in_P) - 1, ec.POLYLINE_POINTS)
    ref = Bezier.Bezier_interp(in_P, in_V, t.copy(), 0)
    assert np.abs(pts - ref).max() < 1e-9


# ---------------------------------------------------------------- COMSOL 模型


def _parse_cubic_beziers(src):
    """从生成的 .java 里把 CubicBezier 的控制点读回来（只为自检）。"""
    out = {}
    for m in re.finditer(r'create\("(cb\d+)", "CubicBezier"\)', src):
        start = m.end()
        nxt = src.find('.create("cb', start)
        chunk = src[start:nxt if nxt > 0 else len(src)]
        pts = [[0.0, 0.0] for _ in range(4)]
        for mm in re.finditer(r'setIndex\("p", ([-+0-9.eE]+), (\d), (\d)\)', chunk):
            pts[int(mm.group(3))][int(mm.group(2))] = float(mm.group(1))
        out[m.group(1)] = np.array(pts)
    return out


def _parse_lines(src):
    """把 LineSegment 的端点读回来。"""
    out = {}
    for m in re.finditer(r'create\("(ls\d+)", "LineSegment"\)', src):
        start = m.end()
        nxt = src.find('.create("ls', start)
        chunk = src[start:nxt if nxt > 0 else len(src)]
        coords = re.findall(r'coord[12]", new String\[\]\{"([-+0-9.eE]+)", "([-+0-9.eE]+)"\}', chunk)
        out[m.group(1)] = np.array([[float(a), float(b)] for a, b in coords])
    return out


def test_model_java_control_points(case, entities, tmp_path):
    """生成的 .java 里 88 个 CubicBezier 必须就是求解器用的控制多边形。"""
    src = ec.model_java(case, entities, "design.h5")
    bez = _parse_cubic_beziers(src)
    n_seg = case.layout.n_ctrl + 1
    assert len(bez) == 8 * n_seg

    in_polys, ex_polys = design.bezier_polygons(case.dofs, case.geometry, case.layout)
    for k, (sign, _port, _side) in enumerate(ec.QUADRANTS):
        s = np.asarray(sign, dtype=float)
        for i in range(n_seg):
            for base, polys in ((1, in_polys), (1 + 4 * n_seg, ex_polys)):
                tag = "cb%d" % (base + k * n_seg + i)
                # 逐位相同：都是 poly * sign，浮点运算顺序也一致
                assert np.array_equal(bez[tag], polys[i] * s), tag
                if i:
                    prev = "cb%d" % (base + k * n_seg + i - 1)
                    assert np.array_equal(bez[prev][3], bez[tag][0])


def test_model_java_lines_and_selections(case, entities):
    src = ec.model_java(case, entities, "design.h5")
    lns = _parse_lines(src)
    assert len(lns) == 12

    # 后 4 条是端口印记线：与六边形边共线，用来把边切成 5 段
    for k, (sign, _port, _side) in enumerate(ec.QUADRANTS):
        want = (np.array([case.geometry.port_edge_bd_in, case.geometry.port_edge_bd_ex])
                * np.asarray(sign, float))
        assert np.allclose(lns["ls%d" % (9 + k)], want)

    # 前 8 条是固定直壁段：一端是壁的固定端，另一端是波导口
    g = case.geometry
    want = [g.wg_fixed_bd_in, g.wg_fixed_bd_ex] * 4
    for i, sign in enumerate([np.asarray(q[0], float) for q in ec.QUADRANTS]):
        for j in range(2):
            seg = lns["ls%d" % (1 + 2 * i + j)]
            assert np.allclose(seg[0], want[2 * i + j] * sign)
    # 15 个 Box 选择都必须恰好框住该框的东西
    assert ec.verify_selections(case, entities) == []
    assert src.count('selection().create("') == 15
    for tag, _dim, _box, _note in ec.selection_boxes(case, entities):
        assert 'selection("%s")' % tag in src


def test_model_java_uses_reference_idioms(case, entities):
    """与 cross_waveguide_test.java 对齐的几处关键设置，掉了就会静默出错。"""
    src = ec.model_java(case, entities, "design.h5")
    for needle in (
            'set("components", "outofplane")',
            'set("PortType", "Numeric")',
            'create("sctr_top", "Scattering", 1)',
            'create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")',
            '"RefractiveIndex"',
            '"BoundaryModeAnalysis"',
            'ewfd.S11',
            'ewfd.Atotal',
            'set("L_unit"'):
        assert needle in src, needle
    # 端口必须挂在按坐标的选择上，不能写死编号
    assert "selection().set(" not in src
