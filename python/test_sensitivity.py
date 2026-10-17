from original_mesh import cross_waveguide
import numpy as np
import scipy.io as spio

import os

from paths import chdir_data, run_solver

chdir_data()

def get_dof_derivative(DOFs, params:tuple):
    freq, w_wg, L_sim, L_fixed, epsilon, mu= params
    os.system("rm *result.mat")
    cross_waveguide(DOFs, freq=freq)
    run_solver()
    data = spio.loadmat("cross_waveguide_result.mat")
    S_params = data["S_params_0"]
    ex_sens = data["ex_opt_sens"]
    in_sens = data["in_opt_sens"]
    DOF_sens = points_to_DOF(in_sens[1:,:],ex_sens[1:,:])
    return(DOF_sens, S_params)
    

def test_sensitivity():
    # test waveguide generation
    N_ctrl = 5

    all_S_params = []
    all_sens = []
    all_DOFS = []
    os.system("rm *result.mat")
    DOFs = cross_waveguide(DOFs= np.zeros((4*N_ctrl + 2)), use_original=True, freq=0.1)

    all_DOFS.append(DOFs.flatten())
    
    Sens, S_params = get_dof_derivative(DOFs, (0.1,))

    all_S_params.append(S_params.flatten())
    all_sens.append(Sens.flatten())
    
    DOFs += np.random.rand(DOFs.shape[0]) * 0.01
    all_DOFS.append(DOFs.flatten())
    Sens, S_params = get_dof_derivative(DOFs, (0.1,))
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

if __name__ == "__main__":
    test_sensitivity()