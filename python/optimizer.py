"""主优化器：目标函数与梯度的封装、梯度下降驱动。

注意：这一轮只做模块搬移，**优化算法与旧实现逐行等价**（同样 10 步、同样步长、
同样日志格式）。换成 L-BFGS-B 之类的改进是单独的一步。
"""

import os
from math import cos, pi

import numpy as np

import design
import physics as physics_mod
import solver

SENS_CSV = "sens_list.csv"
SPARAM_CSV = "S_params_list.csv"
DOF_CSV = "DOF_list.csv"


def last_dofs(path=DOF_CSV):
    """日志里最后一轮的设计变量；没有历史（文件缺失 / 为空 / 只有表头外的空行）时返回 None。

    空文件会让 ``np.loadtxt`` 返回 0 维、单行文件返回 1 维，所以这里把两种情况
    都归一化掉——"有没有历史"是一个只该在开始处判断一次的事。
    """
    try:
        if os.path.getsize(path) == 0:
            return None
    except OSError:
        return None
    arr = np.loadtxt(path, delimiter=",")
    arr = np.atleast_2d(arr)
    return arr[-1, :]


def dof_derivative(case):
    """一次伴随求解，返回 ``(-d|S1|^2/dDOF, S_params)``。

    符号约定沿用旧实现：``gradient_descent`` 里是 ``DOF += sens * alpha``，
    也就是对 |S1|^2 做梯度上升。
    """
    design.make_case(case.dofs, case.physics, case.geometry, case.name)
    result = solver.run(case)
    S = result["S_params"]
    sens = design.sensitivity_to_dofs(result["sens_ex"], result["sens_in"])
    return -sens.flatten(), S.flatten()


def _truncate_logs():
    for name in (SENS_CSV, SPARAM_CSV, DOF_CSV):
        open(name, "w").close()


def _append_logs(dofs, sens, S):
    with open(SENS_CSV, "a") as fp:
        np.savetxt(fp, sens.reshape((1, len(sens))), delimiter=',')
    with open(SPARAM_CSV, "a") as fp:
        np.savetxt(fp, S.reshape((1, len(S))), delimiter=',')
    with open(DOF_CSV, "a") as fp:
        np.savetxt(fp, dofs.reshape((1, len(dofs))), delimiter=',')


def gradient_descent(derivative, case, alpha, iterations=10):
    """固定步长的梯度上升（原样保留旧行为，含 11 次求解）。"""
    sens, S = derivative(case)
    for _ in range(iterations):
        case.dofs += sens * alpha
        case.dofs[case.dofs < 0.0] = 0.0

        sens, S = derivative(case)
        print("Transmission:", np.abs(S[1]) ** 2)
        _append_logs(case.dofs, sens, S)


def run(name="cross_waveguide", physics=None, geometry=None, n_ctrl=10,
        alpha0=1e-3, iterations=10):
    """跑一遍优化循环（等价于旧的 ``optimize_cross_waveguide2()``）。"""
    ph = physics or physics_mod.Physics()
    geo = geometry or design.Geometry()
    layout = design.Layout(n_ctrl)

    dofs = last_dofs()
    if dofs is not None:
        # 断点续跑：沿用最后一轮的设计变量，日志继续追加
        case = design.make_case(dofs, ph, geo, name, build=False)
    else:
        _truncate_logs()
        dofs = design.initial_dofs(geo, layout)
        case = design.make_case(dofs, ph, geo, name)

    alpha_base = geo.w_wg * np.ones_like(dofs)
    # 两个端点角（弧度）的量纲和坐标不同，单独放大步长
    alpha_base[1] *= (2 * pi / geo.w_wg) ** 2
    alpha_base[2 * n_ctrl + 3] *= (2 * pi / geo.w_wg) ** 2

    gradient_descent(dof_derivative, case, alpha_base * alpha0,
                     iterations=iterations)
    return case
