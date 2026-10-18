"""材料与频率参数，以及平板波导色散方程的求解。

原先散落在 ``original_mesh.py`` 的若干 float 参数，现在集中成一个 ``Physics``
对象；色散求解来自旧的 ``propagation_constant0.py``。
"""

from dataclasses import dataclass
from math import sqrt

import numpy as np
import scipy.optimize as opt


@dataclass(frozen=True)
class Physics:
    """一次仿真的物理参数（材料 + 频率）。"""

    freq: float = 1.55
    epsilon_wg: float = 3.4 ** 2
    epsilon_amb: float = 1.0
    mu_wg: float = 1.0
    mu_amb: float = 1.0

    @property
    def n_wg(self):
        return sqrt(self.epsilon_wg * self.mu_wg)

    @property
    def n_amb(self):
        return sqrt(self.epsilon_amb * self.mu_amb)

    @property
    def omega(self):
        return 2.0 * np.pi * self.freq


def _dispersion(beta, h, n0, n1, k0):
    """色散方程 g(beta) = kappa1*h/2 - arctan(p0/kappa1)。"""
    kappa1 = np.sqrt(k0 ** 2 * n1 ** 2 - beta ** 2)
    p0 = np.sqrt(beta ** 2 - k0 ** 2 * n0 ** 2)
    return kappa1 * h / 2 - np.arctan(p0 / kappa1)


def _ddispersion(beta, h, n0, n1, k0):
    kappa1 = np.sqrt(k0 ** 2 * n1 ** 2 - beta ** 2)
    p0 = np.sqrt(beta ** 2 - k0 ** 2 * n0 ** 2)
    return -h / 2 * beta / kappa1 - beta / (p0 * kappa1)


def beta_from_k0(n0, n1, h, k0_sweep):
    """对每个 k0 用牛顿迭代求传播常数 beta（沿用旧实现，含连续性初值）。"""
    betas = []
    beta_approx = (n0 + 8 * n1) * k0_sweep[0] / 9
    for k0 in k0_sweep:
        cost = lambda b: _dispersion(b, h, n0, n1, k0)
        dcost = lambda b: _ddispersion(b, h, n0, n1, k0)
        beta = opt.newton(cost, beta_approx, dcost)
        betas.append(beta)
        beta_approx = beta
        if n1 * k0 < beta or n0 * k0 > beta:
            print("Error! minmum or maximum frequency reached!")
    return betas


def beta(physics, width):
    """给定波导宽度和材料，返回基模传播常数 beta。"""
    return beta_from_k0(
        n0=physics.n_amb, n1=physics.n_wg, h=width,
        k0_sweep=[physics.freq * 2 * np.pi],
    )[0]
