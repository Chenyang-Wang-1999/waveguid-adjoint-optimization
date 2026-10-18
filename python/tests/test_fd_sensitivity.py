"""伴随灵敏度的方向性校验（可选、较慢）。

用有限差分校核伴随梯度：沿随机方向扰动设计变量，比较

    grad · delta                  （伴随给出）
    f(dofs+delta) - f(dofs)       （有限差分），f = |S1|^2

默认不跑（要真的调 gmsh + 求解器，几十秒）。显式打开：

    WAVEGUIDE_FD_TEST=1 python3 -m pytest tests/test_fd_sensitivity.py -s

已知问题（**重构前就存在**，不是这次模块拆分引入的）：
``-DOF_sens`` 与 d|S1|^2/dDOF 之间差一个随方向变化的因子（实测 0.25 ~ 1.05），
所以本用例只断言**方向正确**（符号一致），并把比值打印出来供人工检查。
量级问题的可能来源：

* C++ 里 QOI 的右端项用的是 ``2*Re(conj(S)*dS/dA)``，与 |S|^2 的链式法则
  之间可能差一个常数 2；
* 边界节点由 ``mesh_snap_to_Bezier`` 吸附到曲线上，扰动设计变量后每个节点
  的参数 t 也会重新算，而网格灵敏度用的是"t 固定"的假设。

要查清需要单独做一组控制实验，不在本轮重构范围内。
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
    design.write_design(case)        # 只移动边界，不重划网格
    return float(np.abs(solver.run(case)["S_params"][1]) ** 2)


def test_gradient_direction(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    geo, layout = design.Geometry(), design.Layout(5)
    dofs0 = design.initial_dofs(geo, layout)
    case = design.make_case(dofs0, name="fdtest")

    grad, S = optimizer.dof_derivative(case)
    f0 = float(np.abs(S[1]) ** 2)

    rng = np.random.default_rng(1)
    ratios = []
    for _ in range(3):
        delta = rng.standard_normal(len(dofs0)) * 1e-4
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
