import scipy.io as spio
import numpy as np
import os
import pygmsh,gmsh,pyvista
import Bezier
from math import sin,cos, sqrt, pi

def cross_waveguide_interp_2(
    DOFs:np.ndarray = np.array([0,0]),
    L_sim = 10,
    L_wg_fixed = 3,
    w_wg = 1,
    N_BPoints = 50,
):  
    # some vectors
    e1 = np.array([cos(pi/6), sin(pi/6)])
    e2 = np.array([-sin(pi/6), cos(pi/6)])
    # some points
    L_wg_cen = L_sim * cos(pi/6)
    wg_input_cen = L_wg_cen*e1
    wg_input_bd_ex = wg_input_cen + w_wg*e2 /2
    wg_input_bd_in = wg_input_cen - w_wg*e2 /2
    wg_fixed_bd_ex = wg_input_bd_ex - L_wg_fixed*e1
    wg_fixed_bd_in = wg_input_bd_in - L_wg_fixed*e1


    # DOFs: 4*N_ctrl + 2: 2*N_ctrl interior points, 2 boundary points
    # [x_in_bd, x_in_bk1, y_in_bk1, ..., y_ex_bd, x_ex_bk1, y_ex_bk1, ...]
    N_ctrl = int((DOFs.__len__() - 4)/4+0.5)
    
    # convert DOF to points
    in_points_bd = np.array([[DOFs[0], 0]])
    in_end_vec = np.array([[cos(DOFs[1]), sin(DOFs[1])]])
    in_points  = DOFs[2:2*N_ctrl+2] # all 2*N_ctrl points
    in_points = np.reshape(in_points, (N_ctrl, 2))
    in_points = np.vstack([wg_fixed_bd_in, in_points, in_points_bd])

    ex_points_bd = np.array([[0, DOFs[2*N_ctrl + 2]]])
    ex_end_vec = np.array([[cos(DOFs[2*N_ctrl + 3]),sin(DOFs[2*N_ctrl + 3])]])
    ex_points  = DOFs[2*N_ctrl + 4:] # all 2*N_ctrl points
    ex_points = np.reshape(ex_points, (N_ctrl, 2))
    ex_points = np.vstack([wg_fixed_bd_ex, ex_points, ex_points_bd])

    # print("1")
    # Number of Bezier points on each side
    ex_tangent = Bezier.Bezier_tangent_vector(ex_points, 0, 0, - e1, ex_end_vec.flatten(), True)
    in_tangent = Bezier.Bezier_tangent_vector(in_points, 0, 0, - e1, in_end_vec.flatten(), True)

    # print("2")

    t = np.linspace(0, N_ctrl + 1, N_BPoints)
    ex_line_pts = Bezier.Bezier_interp(ex_points, ex_tangent, t, 0)
    in_line_pts = Bezier.Bezier_interp(in_points, in_tangent, t, 0)

    np.savetxt("ex_line_pts.csv", ex_line_pts, delimiter=',')
    np.savetxt("in_line_pts.csv", in_line_pts, delimiter=',')


def cross_waveguide_interp(
    DOFs:np.ndarray = np.array([0,0]),
    L_sim = 10,
    L_wg_fixed = 3,
    w_wg = 1,
    N_BPoints = 20
):  
    # some vectors
    e1 = np.array([cos(pi/6), sin(pi/6)])
    e2 = np.array([-sin(pi/6), cos(pi/6)])
    # some points
    L_wg_cen = L_sim * cos(pi/6)
    wg_input_cen = L_wg_cen*e1
    wg_input_bd_ex = wg_input_cen + w_wg*e2 /2
    wg_input_bd_in = wg_input_cen - w_wg*e2 /2
    wg_fixed_bd_ex = wg_input_bd_ex - L_wg_fixed*e1
    wg_fixed_bd_in = wg_input_bd_in - L_wg_fixed*e1


    # DOFs: 4*N_ctrl + 2: 2*N_ctrl interior points, 2 boundary points
    # [x_in_bd, x_in_bk1, y_in_bk1, ..., y_ex_bd, x_ex_bk1, y_ex_bk1, ...]
    N_ctrl = int((DOFs.__len__() - 2)/4+0.5)
    
    # convert DOF to points
    in_points_bd = np.array([[DOFs[0], 0]])
    in_points  = DOFs[1:2*N_ctrl+1] # all 2*N_ctrl points
    in_points = np.reshape(in_points, (N_ctrl, 2))
    in_points = np.vstack([wg_fixed_bd_in, in_points, in_points_bd])

    ex_points_bd = np.array([[0, DOFs[2*N_ctrl + 1]]])
    ex_points  = DOFs[2*N_ctrl + 2:] # all 2*N_ctrl points
    ex_points = np.reshape(ex_points, (N_ctrl, 2))
    ex_points = np.vstack([wg_fixed_bd_ex, ex_points, ex_points_bd])

    # Number of Bezier points on each side
    ex_tangent = Bezier.Bezier_tangent_vector(ex_points, 0, 1, - e1, -e1, True)
    in_tangent = Bezier.Bezier_tangent_vector(in_points, 0, 1, - e1, -e1, True)

    t = np.linspace(0, N_ctrl + 1, N_BPoints)
    ex_line_pts = Bezier.Bezier_interp(ex_points, ex_tangent, t, 0)
    in_line_pts = Bezier.Bezier_interp(in_points, in_tangent, t, 0)

    np.savetxt("ex_line_pts.csv", ex_line_pts, delimiter=',')
    np.savetxt("in_line_pts.csv", in_line_pts, delimiter=',')


DOFs = np.loadtxt("DOFs.csv", delimiter=',')
cross_waveguide_interp_2(
    DOFs,
    L_wg_fixed = 0.3,
    
    L_sim = (0.5 + 0.3)/cos(pi/6) ,
    w_wg = 0.08,
    N_BPoints = 100
)
