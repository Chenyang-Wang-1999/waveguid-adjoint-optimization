from original_mesh import cross_waveguide, cross_waveguide_2
import numpy as np
import scipy.io as spio
from scipy import optimize
import os
from math import sin, cos, sqrt, pi

from paths import chdir_data, run_solver

# 求解器在 data/ 里读写 infile.pot / *.msh / *_result.mat
chdir_data()


def DOF_to_points(DOFs):
    N_ctrl = int((DOFs.__len__() - 2)/4+0.5)
    
    # convert DOF to points
    in_points_bd = np.array([[DOFs[0], 0]])
    in_points  = DOFs[1:2*N_ctrl+1] # all 2*N_ctrl points
    in_points = np.reshape(in_points, (N_ctrl, 2))
    in_points = np.vstack([in_points, in_points_bd])

    ex_points_bd = np.array([[0, DOFs[2*N_ctrl + 1]]])
    ex_points  = DOFs[2*N_ctrl + 2:] # all 2*N_ctrl points
    ex_points = np.reshape(ex_points, (N_ctrl, 2))
    ex_points = np.vstack([ex_points, ex_points_bd])

    return (in_points, ex_points)

def points_to_DOF(in_points, ex_points):
    N_ctrl = in_points.shape[0] - 1
    DOFs = np.array([])
    DOFs = np.hstack([DOFs, in_points[-1,0]])
    DOFs = np.hstack([DOFs, in_points[:-1,:].flatten()])
    DOFs = np.hstack([DOFs, ex_points[-1,1]])
    DOFs = np.hstack([DOFs, ex_points[:-1,:].flatten()])
    return DOFs

def DOF_to_points2(DOFs):
    N_ctrl = int((DOFs.__len__() - 4)/4+0.5)
    
    # convert DOF to points
    in_points_bd = np.array([[DOFs[0], 0]])
    in_end_theta = DOFs[1]
    in_points  = DOFs[2:2*N_ctrl+2] # all 2*N_ctrl points
    in_points = np.reshape(in_points, (N_ctrl, 2))
    in_points = np.vstack([in_points, in_points_bd])

    ex_points_bd = np.array([[0, DOFs[2*N_ctrl + 2]]])
    ex_end_theta = DOFs[2*N_ctrl + 3]
    ex_points  = DOFs[2*N_ctrl + 4:] # all 2*N_ctrl points
    ex_points = np.reshape(ex_points, (N_ctrl, 2))
    ex_points = np.vstack([ex_points, ex_points_bd])

    return (in_points, in_end_theta, ex_points, ex_end_theta)

def points_to_DOF2(in_points, in_end_theta, ex_points, ex_end_theta):
    N_ctrl = in_points.shape[0] - 1
    DOFs = np.array([])
    DOFs = np.hstack([DOFs, [in_points[-1,0], in_end_theta]])
    DOFs = np.hstack([DOFs, in_points[:-1,:].flatten()])
    DOFs = np.hstack([DOFs, [ex_points[-1,1], ex_end_theta]])
    DOFs = np.hstack([DOFs, ex_points[:-1,:].flatten()])
    return DOFs

def get_dof_derivative(DOFs, params:tuple):
    rerun, freq, w_wg, L_sim, L_fixed, epsilon= params
    if(rerun):
        cross_waveguide(DOFs, freq=freq, w_wg = w_wg, L_sim=L_sim, L_wg_fixed=L_fixed,
        epsilon_wg = epsilon)
        run_solver()
    data = spio.loadmat("cross_waveguide_result.mat")
    S_params = data["S_params_0"]
    ex_sens = data["ex_opt_sens"]
    in_sens = data["in_opt_sens"]
    DOF_sens = points_to_DOF(in_sens[1:,:],ex_sens[1:,:])
    return(-DOF_sens.flatten(), S_params.flatten())

def get_dof_derivative2(DOFs, params:tuple):
    rerun, freq, w_wg, L_sim, L_fixed, epsilon= params
    if(rerun):
        cross_waveguide_2(DOFs, freq=freq, w_wg = w_wg, L_sim=L_sim, L_wg_fixed=L_fixed,
        epsilon_wg = epsilon)
        run_solver()
    data = spio.loadmat("cross_waveguide_result.mat")
    S_params = data["S_params_0"]
    ex_sens = data["ex_opt_sens"]
    in_sens = data["in_opt_sens"]
    DOF_sens = points_to_DOF2(in_sens[1:-1,:],in_sens[-1,1],ex_sens[1:-1,:],ex_sens[-1,1])
    return(-DOF_sens.flatten(), S_params.flatten())
 

def get_S_params(DOFs, freq, w_wg, L_sim,L_fixed, epsilon):
    DOF_sens, S_params = get_dof_derivative(DOFs, (True, freq, w_wg, L_sim, L_fixed, epsilon))
    
    print("Transmission: ",(np.abs(S_params[1])**2))

    with open("sens_list.csv","a") as fp_sens:
        np.savetxt(fp_sens, DOF_sens.flatten(), delimiter=',')
    
    with open("S_params_list.csv","a") as fp_S_params:
        np.savetxt(fp_S_params, S_params.flatten(), delimiter=',')
    with open("DOF_list.csv","a") as fp_DOF:
        np.savetxt(fp_DOF, DOFs.flatten(), delimiter=',')
    return (-np.abs(S_params[1])**2)

def get_DOF_derivative(DOFs,freq, w_wg, L_sim, L_fixed, epsilon):
    DOF_sens, S_params = get_dof_derivative(DOFs, (True, freq, w_wg, L_sim, L_fixed, epsilon))
    return -DOF_sens

def gradient_descent(my_der_fun,DOF0, params, alpha):
    sens, S = my_der_fun(DOF0, params)
    DOF = DOF0
    
    while((sens**2).sum()/sens.__len__() > 1e-5):
        DOF += sens * alpha

        zero_points = DOF < 0.0
        DOF[zero_points] = 0

        sens, S = my_der_fun(DOF0, params)
        print("Transmission:", np.abs(S[1])**2)
        with open("sens_list.csv","a") as fp_sens:
            np.savetxt(fp_sens, sens.reshape((1, sens.__len__())), delimiter=',')
        
        with open("S_params_list.csv","a") as fp_S_params:
            np.savetxt(fp_S_params, S.reshape((1, S.__len__())), delimiter=',')
        with open("DOF_list.csv","a") as fp_DOF:
            np.savetxt(fp_DOF, DOF.reshape((1, DOF.__len__())), delimiter=',')

def optimize_cross_waveguide():
    freq = 1.55
    w_wg = 0.08
    L_sim = 1
    L_fixed = 0.3
    epsilon = 3.4**2

    N_ctrl = 10

    fp_sens = open("sens_list.csv","w")
    fp_sens.close()
    fp_S_params = open("S_params_list.csv","w")
    fp_S_params.close()
    fp_DOF = open("DOF_list.csv","w")
    fp_DOF.close()

    DOF0 = cross_waveguide(np.zeros((4*N_ctrl + 2)), freq=freq, w_wg = w_wg, L_sim=L_sim, L_wg_fixed=L_fixed,
            epsilon_wg = epsilon, use_original=True)
    
    lb = np.zeros_like(DOF0)
    ub = L_sim * np.ones_like(DOF0)
    bounds = optimize.Bounds(lb, ub)


    gradient_descent(DOF0, (True, freq, w_wg, L_sim, L_fixed, epsilon), w_wg/100)

    # Result = optimize.minimize(get_S_params, DOF0, (freq, w_wg, L_air, L_fixed, epsilon),
    #                 method = 'BFGS', jac=get_DOF_derivative, constraints=(bounds,))

def optimize_cross_waveguide2():
    freq = 1.55
    w_wg = 0.08
    L_opt = 0.2
    L_fixed = 0.3
    L_sim = (L_opt+L_fixed)/cos(pi/6)
    epsilon = 3.4**2

    N_ctrl = 10

    if(os.path.exists('DOF_list.csv')):
        DOFs = np.loadtxt("DOF_list.csv",delimiter=',')
        DOF0 = DOFs[-1,:]
    else:
        DOF0 = cross_waveguide_2(np.zeros((4*N_ctrl + 4)), freq=freq, w_wg = w_wg, L_sim=L_sim, L_wg_fixed=L_fixed,
            epsilon_wg = epsilon, use_original=True)
        fp_sens = open("sens_list.csv","w")
        fp_sens.close()
        fp_S_params = open("S_params_list.csv","w")
        fp_S_params.close()
        fp_DOF = open("DOF_list.csv","w")
        fp_DOF.close()


    
    lb = np.zeros_like(DOF0)
    ub = L_sim * np.ones_like(DOF0)
    bounds = optimize.Bounds(lb, ub)

    alpha_base = w_wg*np.ones_like(DOF0)
    alpha_base[1] *= (2*pi/w_wg)**2
    alpha_base[2*N_ctrl+3] *= (2*pi/w_wg)**2
    gradient_descent(get_dof_derivative2,DOF0, (True, freq, w_wg, L_sim, L_fixed, epsilon), alpha_base/500)

def test_sensitivity():
    # test waveguide generation
    N_ctrl = 5

    freq = 1.55
    w_wg = 0.08
    L_sim = 0.6
    L_fixed = 0.3
    epsilon = 3.4**2

    all_S_params = []
    all_sens = []
    all_DOFS = []
    os.system("rm *result.mat")
    DOFs = cross_waveguide_2(DOFs= np.zeros((4*N_ctrl + 4)), use_original=True, freq=freq,w_wg = w_wg,
        L_sim=L_sim, L_wg_fixed=L_fixed,epsilon_wg=epsilon)
    print(DOFs)
    all_DOFS.append(DOFs.flatten())
    
    params = (1, freq, w_wg, L_sim, L_fixed, epsilon)
    Sens, S_params = get_dof_derivative2(DOFs, params)

    all_S_params.append(S_params.flatten())
    all_sens.append(Sens.flatten())
    
    delta_DOF = np.random.rand(DOFs.shape[0]) * 0.01*w_wg
    delta_DOF[1] = np.random.rand()*2*pi*0.01
    delta_DOF[2*N_ctrl + 3] = np.random.rand()*2*pi*0.01
    DOFs += delta_DOF

    all_DOFS.append(DOFs.flatten())
    params = (1, freq, w_wg, L_sim, L_fixed, epsilon)
    Sens, S_params = get_dof_derivative2(DOFs, params)
    all_S_params.append(S_params.flatten())
    all_sens.append(Sens.flatten())

    data = {}
    data['all_DOFs'] = np.vstack(all_DOFS)
    data['all_S_params'] = np.vstack(all_S_params)
    data["all_sens"] = np.vstack(all_sens)

    spio.savemat('cross_waveguide_test_sense.mat',data)

    diff_FD = np.abs(data['all_S_params'][1,1])**2 - np.abs(data['all_S_params'][0,1])**2
    diff_sens = data['all_sens'].dot(data['all_DOFs'][1,:] - data['all_DOFs'][0,:])
    print(diff_sens)
    print(diff_FD)
    print(all_sens)

def test_DOF_conversion():
    N_ctrl = 5
    DOFs = np.random.rand(4*N_ctrl + 2)
    in_points, ex_points = DOF_to_points(DOFs)
    new_DOFs = points_to_DOF(in_points, ex_points)

    print(DOFs)
    print(new_DOFs-DOFs)



    
if __name__ == "__main__":
    optimize_cross_waveguide2()
