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


def test_sensitivity_to_dofs_shape():
    """求解器给的 (n_ctrl+3, 2) 灵敏度应映射回 4*n_ctrl+4 长度。"""
    n = 6
    layout = design.Layout(n)
    ex_sens = np.arange((n + 3) * 2, dtype=float).reshape(n + 3, 2)
    in_sens = -ex_sens
    grad = design.sensitivity_to_dofs(ex_sens, in_sens)
    assert grad.shape == (layout.size,)


def test_case_load_roundtrip(tmp_path, monkeypatch):
    """design.h5 必须能独立还原出完整的 Case（不需要调用方手里有 DOF）。"""
    pytest.importorskip("Bezier")
    monkeypatch.chdir(tmp_path)

    geo = design.Geometry(w_wg=0.05, L_sim=0.7, L_wg_fixed=0.25)
    ph = design.physics_mod.Physics(freq=1.3, epsilon_wg=9.0, epsilon_amb=1.1,
                                    mu_wg=1.2, mu_amb=1.3)
    dofs = np.linspace(0.1, 0.9, design.Layout(5).size)
    case = design.make_case(dofs, physics=ph, geometry=geo, name="t", build=False)
    design.write_design(case)

    back = design.Case.load()
    assert np.allclose(back.dofs, dofs)
    assert back.physics == ph
    assert back.geometry == geo
    assert back.name == "t"


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
