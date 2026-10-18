"""设计变量布局与 DOF <-> 控制点映射。"""

import numpy as np
import pytest

import design


def test_layout_from_dofs():
    for n in (5, 10, 17):
        assert design.Layout.from_dofs(np.zeros(4 * n + 4)).n_ctrl == n
        assert design.Layout(n).size == 4 * n + 4


def test_initial_dofs_size():
    geo = design.Geometry()
    layout = design.Layout(10)
    dofs = design.initial_dofs(geo, layout)
    assert dofs.shape == (layout.size,)


def test_control_points_inverse():
    """control_points 与 _points_to_dofs 必须互为逆映射。"""
    geo, layout = design.Geometry(), design.Layout(6)
    rng = np.random.default_rng(0)
    dofs = rng.random(layout.size)

    in_P, in_theta_vec, ex_P, ex_theta_vec = design.control_points(dofs, geo, layout)
    # _points_to_dofs 吃的是"内部点 + 端口边界点"（不含固定锚点），
    # 也就是 in_P[1:] 那一组，与求解器给出的灵敏度行布局一致
    back = design._points_to_dofs(in_P[1:], np.arctan2(in_theta_vec[0, 1], in_theta_vec[0, 0]),
                                 ex_P[1:], np.arctan2(ex_theta_vec[0, 1], ex_theta_vec[0, 0]))
    assert np.allclose(back, dofs)

    # 控制点矩阵的首行是固定锚点，末行是端口边界点
    assert np.allclose(in_P[0], geo.wg_fixed_bd_in)
    assert np.allclose(ex_P[0], geo.wg_fixed_bd_ex)
    assert np.allclose(in_P[-1], [dofs[0], 0.0])
    assert np.allclose(ex_P[-1], [0.0, dofs[2 * layout.n_ctrl + 2]])


def test_geometry_margins_detects_axis_crossing():
    """几何自检必须能发现"边界越过对称轴"——那是 gmsh 直接 abort 的入口。"""
    pytest.importorskip("Bezier")
    geo, layout = design.Geometry(), design.Layout(10)
    dofs = design.initial_dofs(geo, layout)
    # 直线设计完全落在象限内：曲线的端口端点按构造就在坐标轴上（== 0），
    # 所以判据是"不越界"而不是"严格大于 0"
    assert min(design.geometry_margins(dofs, geo, layout).values()) >= -1e-12

    bad = dofs.copy()
    bad[0] = -0.01          # 端点 x_in_bd 变负 -> in 曲线越过 y 轴
    assert design.geometry_margins(bad, geo, layout)["in_min_x"] < 0


def test_sensitivity_to_dofs_shape():
    """求解器给的 (n_ctrl+3, 2) 灵敏度应映射回 4*n_ctrl+4 长度。"""
    n = 6
    layout = design.Layout(n)
    ex_sens = np.arange((n + 3) * 2, dtype=float).reshape(n + 3, 2)
    in_sens = -ex_sens
    grad = design.sensitivity_to_dofs(ex_sens, in_sens)
    assert grad.shape == (layout.size,)


def test_case_load_roundtrip(tmp_path, monkeypatch):
    """design.h5 必须能独立还原出完整的 Case（不需要调用方手里有 DOF）。

    网格跟着一起建：``write_design`` 会从磁盘上的 ``.msh`` 读回节点、算出
    ``/design/node_t``，所以"没有网格的 design.h5"现在是不合法的输入。
    """
    pytest.importorskip("Bezier")
    monkeypatch.chdir(tmp_path)

    geo = design.Geometry(w_wg=0.05, L_sim=0.7, L_wg_fixed=0.25)
    ph = design.physics_mod.Physics(freq=1.3, epsilon_wg=9.0, epsilon_amb=1.1,
                                    mu_wg=1.2, mu_amb=1.3)
    layout = design.Layout(3)
    dofs = design.initial_dofs(geo, layout)
    design.make_case(dofs, physics=ph, geometry=geo, name="t")

    back = design.Case.load()
    assert np.allclose(back.dofs, dofs)
    assert back.physics == ph
    assert back.geometry == geo
    assert back.name == "t"

    # node_t 只应该在两条优化边界的节点上有值，而且 P(t) 必须还原出节点
    data = design.dataio.read_design("design.h5")
    node_t, mesh_file = data["node_t"], back.mesh_file
    assert node_t is not None and np.isfinite(node_t).sum() > 0
    case = design.Case(dofs=dofs, physics=ph, geometry=geo, name="t",
                       ids={"wg": 1, "ex_opt": 7, "in_opt": 8})
    import Bezier
    import meshio
    in_P, in_V, ex_P, ex_V = design.bezier_curves(dofs, geo, layout)
    mesh = meshio.read(mesh_file)
    worst = 0.0
    for idx in np.flatnonzero(np.isfinite(node_t)):
        target = np.abs(mesh.points[idx][:2])
        t = np.array([node_t[idx]])
        worst = max(worst, min(
            np.linalg.norm(Bezier.Bezier_interp(P, V, t, 0)[0] - target)
            for P, V in ((in_P, in_V), (ex_P, ex_V))))
    assert worst < 1e-3, worst


def test_bezier_polygons_match_interp():
    """分段三次 Bézier 的控制多边形必须和 Bezier_interp 逐位等价。

    这条是整个几何迁移的基石：gmsh 用这四条点建三次 Bézier，等价于
    Bezier_interp 在该段上的 Bernstein 形式，所以"网格节点天生在曲线上"
    成立，不需要再吸附。
    """
    Bezier = pytest.importorskip("Bezier")
    geo, layout = design.Geometry(), design.Layout(4)
    rng = np.random.default_rng(3)
    dofs = design.initial_dofs(geo, layout) + rng.normal(scale=0.005, size=layout.size)

    in_P, in_V, ex_P, ex_V = design.bezier_curves(dofs, geo, layout)
    in_polys, ex_polys = design.bezier_polygons(dofs, geo, layout)
    n_seg = layout.n_ctrl + 1

    for polys, P, V in ((in_polys, in_P, in_V), (ex_polys, ex_P, ex_V)):
        t = np.linspace(0.0, n_seg, 400)
        ref = Bezier.Bezier_interp(P, V, t.copy(), 0)
        out = np.zeros_like(ref)
        for k, tk in enumerate(t):
            tk = min(max(tk, 0.0), n_seg - 1e-15)
            i = min(int(tk), n_seg - 1)
            s, u = tk - i, 1.0 - (tk - i)
            out[k] = (u ** 3 * polys[i, 0] + 3 * u ** 2 * s * polys[i, 1]
                      + 3 * u * s ** 2 * polys[i, 2] + s ** 3 * polys[i, 3])
        assert np.abs(out - ref).max() < 1e-15


def test_case_load_rejects_missing_dofs(tmp_path):
    """没有 /design/dofs 的旧文件（或 mat_to_design 转出来的）要给出明确提示。"""
    import dataio
    monkeypatch_path = tmp_path / "design.h5"
    dataio.write_design(
        str(monkeypatch_path),
        in_P=np.zeros((3, 2)), in_V=np.zeros((3, 2)),
        ex_P=np.zeros((3, 2)), ex_V=np.zeros((3, 2)),
        port_id=np.zeros((1, 4)), port_is_analytic=np.zeros((1, 4)),
        port_coord=np.zeros((3, 4)), port_x_axis=np.zeros((3, 4)),
        incident_on=np.zeros((1, 4)), port_other_params=np.zeros((8, 4)),
    )
    with pytest.raises(ValueError, match="dofs"):
        design.Case.load(str(monkeypatch_path))
