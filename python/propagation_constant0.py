import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, pi

def cost_function_2(beta:float,
h:float, n0:float, n1:float, k0:float):
    # cost function f(beta) = kappa1*h - arctan(p0/kappa1)
    kappa1 = np.sqrt((k0**2) * (n1**2) - beta**2)
    p0 = np.sqrt(beta**2 - (k0**2)*(n0**2))
    return (kappa1*h/2 - np.arctan(p0/kappa1))

def diff_cost_function_2(beta:float,
h:float, n0:float, n1:float, k0:float):
    # diff of cost function
    kappa1 = np.sqrt((k0**2) * (n1**2) - beta**2)
    p0 = np.sqrt(beta**2 - (k0**2)*(n0**2))
    return (- h/2 * beta/kappa1-beta/(p0*kappa1))

def get_beta_from_k0(n0:float,n1:float,h:float, k0_sweep):
    beta_sweep = []
    beta_approx = (n0 + 8*n1) * k0_sweep[0] / 9
    for j in range(k0_sweep.__len__()):
        k0 = k0_sweep[j]
        cost_fun = lambda beta:cost_function_2(beta,h,n0,n1,k0)
        diff_cost_fun = lambda beta:diff_cost_function_2(beta,h,n0,n1,k0)
        beta = opt.newton(cost_fun,beta_approx,diff_cost_fun)
        beta_sweep.append(beta)
        beta_approx = beta
        if(n1*k0 < beta or n0*k0>beta):
            print("Error! minmum or maximum frequency reached!")
    return beta_sweep