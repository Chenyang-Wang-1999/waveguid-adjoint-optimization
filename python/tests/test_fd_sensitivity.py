"""伴随灵敏度的有限差分校核（可选、较慢）。

用有限差分校核伴随梯度：沿随机方向扰动设计变量，比较

    grad · delta                  （伴随给出，grad = -d|S1|^2/dDOF）
    f(dofs+delta) - f(dofs)       （有限差分），f = |S1|^2

默认不跑（要真的调 gmsh + 求解器，几十秒）。显式打开：

    WAVEGUIDE_FD_TEST=1 python3 -m pytest tests/test_fd_sensitivity.py -s

注意网格是**设计的函数**：微扰设计变量之后必须重划网格，否则曲线变了、
网格没变，物理量根本不会动（"网格前移"之前是靠 C++ 把节点吸附回新曲线
来模拟变形；现在这条路已经删掉了）。
"""

import os

import numpy as np
import pytest

import design
import optimizer
import paths
import solver

pytestmark = pytest.mark.skipif(
    not (os.path.exists(paths.SOLVER) and os.environ.get("WAVEGUIDE_FD_TEST")),
    reason="需要求解器，且要设 WAVEGUIDE_FD_TEST=1 才跑（较慢）",
)


def _transmission(case):
    """按当前设计变量重建网格，再跑一次正问题，返回 |S1|^2。"""
    case.ids = design.build_mesh(case)
    design.write_design(case)
    return float(np.abs(solver.run(case)["S_params"][1]) ** 2)


def _direction_scale(geo, n_ctrl):
    """各分量的"自然"量级：坐标取 w_wg/10，两个端点角取 0.05 rad。"""
    scale = np.full(design.Layout(n_ctrl).size, geo.w_wg / 10.0)
    scale[1] = 0.05
    scale[2 * n_ctrl + 3] = 0.05
    return scale


def test_gradient_vs_finite_difference(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    geo, layout = design.Geometry(), design.Layout(5)
    dofs0 = design.initial_dofs(geo, layout)
    case = design.make_case(dofs0, name="fdtest")

    grad, S = optimizer.dof_derivative(case)
    f0 = float(np.abs(S[1]) ** 2)
    scale = _direction_scale(geo, layout.n_ctrl)

    rng = np.random.default_rng(1)
    ratios = []
    for _ in range(2):
        delta = rng.standard_normal(len(dofs0)) * scale
        case.dofs = dofs0 + delta
        df = _transmission(case) - f0
        if abs(df) < 1e-12:
            continue
        dgrad = grad.dot(delta)
        ratios.append(dgrad / df)
        assert dgrad * df > 0, (
            "伴随梯度方向与有限差分不符：grad·δ=%g, Δf=%g" % (dgrad, df))

    print("\n[FD] grad·delta / Δf 实测比值：", ["%.3f" % r for r in ratios])
    assert ratios, "有限差分变化太小，无法判断方向"
