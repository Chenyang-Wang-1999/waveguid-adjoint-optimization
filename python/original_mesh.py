import scipy.io as spio
import numpy as np
import os
import pygmsh,gmsh,pyvista
import Bezier
from math import sin,cos, sqrt, pi
from propagation_constant0 import get_beta_from_k0

from paths import chdir_data, run_solver

# 网格、.mat、infile.pot 都生成在 data/ 里
chdir_data()

def basic_config(pot_filename, varname_value_pair):
    pot_fp = open(pot_filename,"w")
    for pair in varname_value_pair:
        print("%s = %s"%(pair[0],pair[1]), file=pot_fp)

def cross_waveguide_2(
    DOFs:np.ndarray = np.array([0,0]),
    name = "cross_waveguide",
    L_sim = 10,
    L_wg_fixed = 3,
    w_wg = 1,
    N_BPoints = 50,
    epsilon_wg = 13,
    mu_wg = 1,
    epsilon_amb = 1,
    mu_amb = 1,
    freq = 1,
    use_original = False
):  
    print("Building geometry ...")
    #print(DOFs)
    mesh_size = 0.5*w_wg
    K = 10
    port_refine = 3
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
    
    if(use_original):
        x0 = w_wg/2/sin(pi/6)
        y0 = w_wg/2/cos(pi/6)
        ex_origin = np.linspace(wg_fixed_bd_ex,[0,y0], N_ctrl + 2)
        in_origin = np.linspace(wg_fixed_bd_in,[x0, 0], N_ctrl + 2)
        in_origin_dofs = in_origin[1:-1].flatten()
        ex_origin_dofs = ex_origin[1:-1].flatten()
        DOFs = np.hstack([[x0,7*pi/6],in_origin_dofs, [y0, 7*pi/6], ex_origin_dofs])
        return DOFs
    
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

    # mesh
    # if(os.path.exists(name+'.msh')):
    #     os.system("rm %s"%(name+'.msh'))
    os.system("rm %s*"%(name))
    # Initialize empty geometry using the build in kernel in GMSH
    # with pygmsh.geo.Geometry() as geom:
    geom = pygmsh.geo.Geometry()
    model = geom.__enter__()
    mx = np.array(
        [1, -1]
    )
    my = np.array(
        [-1,1]
    )
    mxy = -1

    origin = model.add_point((0,0,0), mesh_size=mesh_size)

    input_ex_pt = model.add_point(wg_input_bd_ex, mesh_size = mesh_size/port_refine)
    input_in_pt = model.add_point(wg_input_bd_in, mesh_size = mesh_size/port_refine)
    ex_opt_pts = [model.add_point(ex_line_pts[i,:], mesh_size = mesh_size) for i in range(N_BPoints)]
    in_opt_pts = [model.add_point(in_line_pts[i,:], mesh_size = mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mx = model.add_point(wg_input_bd_ex*mx, mesh_size = mesh_size/port_refine)
    input_in_pt_mx = model.add_point(wg_input_bd_in*mx, mesh_size = mesh_size/port_refine)
    ex_opt_pts_mx = [model.add_point(ex_line_pts[i,:]*mx, mesh_size = mesh_size) for i in range(N_BPoints)] 
    in_opt_pts_mx = [model.add_point(in_line_pts[i,:]*mx, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [in_opt_pts[-1]]

    input_ex_pt_my = model.add_point(wg_input_bd_ex*my, mesh_size = mesh_size/port_refine)
    input_in_pt_my = model.add_point(wg_input_bd_in*my, mesh_size = mesh_size/port_refine)
    ex_opt_pts_my = [model.add_point(ex_line_pts[i,:]*my, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [ex_opt_pts[-1]]
    in_opt_pts_my = [model.add_point(in_line_pts[i,:]*my, mesh_size = mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mxy = model.add_point(wg_input_bd_ex*mxy, mesh_size = mesh_size/port_refine)
    input_in_pt_mxy = model.add_point(wg_input_bd_in*mxy, mesh_size = mesh_size/port_refine)
    ex_opt_pts_mxy = [model.add_point(ex_line_pts[i,:]*mxy, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [ex_opt_pts_mx[-1]]
    in_opt_pts_mxy = [model.add_point(in_line_pts[i,:]*mxy, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [in_opt_pts_my[-1]]

    # lines
    curves = []
    ex_opt_m0 = [model.add_line(ex_opt_pts[i], ex_opt_pts[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_m0
    ex_fixed_m0 = [model.add_line(ex_opt_pts[0],input_ex_pt)]
    curves += ex_fixed_m0
    input_bd_m0 = [model.add_line(input_ex_pt, input_in_pt)]
    curves += input_bd_m0
    in_fixed_m0 = [model.add_line(input_in_pt, in_opt_pts[0])]
    curves += in_fixed_m0
    in_opt_m0 = [model.add_line(in_opt_pts[i], in_opt_pts[i+1]) for i in range(N_BPoints - 1)]
    curves += in_opt_m0
    curves += [model.add_line(in_opt_pts[-1], origin), model.add_line(origin, ex_opt_pts[-1])]
    m0_curves = curves
    wg1_loop = model.add_curve_loop(curves)
    # print(curves)
    wg1 = model.add_plane_surface(wg1_loop)

    

    curves = []
    ex_opt_mxy = [model.add_line(ex_opt_pts_mxy[i], ex_opt_pts_mxy[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_mxy
    ex_fixed_mxy = [model.add_line(ex_opt_pts_mxy[0],input_ex_pt_mxy)]
    curves += ex_fixed_mxy
    input_bd_mxy = [model.add_line(input_ex_pt_mxy, input_in_pt_mxy)]
    curves += input_bd_mxy
    in_fixed_mxy = [model.add_line(input_in_pt_mxy, in_opt_pts_mxy[0])]
    curves += in_fixed_mxy
    in_opt_mxy = [model.add_line(in_opt_pts_mxy[i], in_opt_pts_mxy[i+1]) for i in range(N_BPoints - 1)]
    # print(len(in_opt_pts_mxy))
    curves += in_opt_mxy
    curves += [model.add_line(in_opt_pts_mxy[-1], origin), model.add_line(origin, ex_opt_pts_mxy[-1])]
    # print(curves)
    mxy_curves = curves
    wg3_loop = model.add_curve_loop(curves)
    wg3 = model.add_plane_surface(wg3_loop)

    curves = []
    in_opt_mx = [model.add_line(in_opt_pts_mx[i], in_opt_pts_mx[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_mx
    in_fixed_mx = [model.add_line(in_opt_pts_mx[0],input_in_pt_mx)]
    curves += in_fixed_mx
    input_bd_mx = [model.add_line(input_in_pt_mx, input_ex_pt_mx)]
    curves += input_bd_mx
    ex_fixed_mx = [model.add_line(input_ex_pt_mx, ex_opt_pts_mx[0])]
    curves += ex_fixed_mx
    ex_opt_mx = [model.add_line(ex_opt_pts_mx[i], ex_opt_pts_mx[i+1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_mx
    curves += [ - mxy_curves[-1], - m0_curves[-2]]
    wg2_loop = model.add_curve_loop(curves)
    wg2 = model.add_plane_surface(wg2_loop)

    curves = []
    in_opt_my = [model.add_line(in_opt_pts_my[i], in_opt_pts_my[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_my
    in_fixed_my = [model.add_line(in_opt_pts_my[0],input_in_pt_my)]
    curves += in_fixed_my
    input_bd_my = [model.add_line(input_in_pt_my, input_ex_pt_my)]
    curves += input_bd_my
    ex_fixed_my = [model.add_line(input_ex_pt_my, ex_opt_pts_my[0])]
    curves += ex_fixed_my
    ex_opt_my = [model.add_line(ex_opt_pts_my[i], ex_opt_pts_my[i+1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_my
    curves += [- m0_curves[-1], -mxy_curves[-2]]
    # print(curves)
    wg4_loop = model.add_curve_loop(curves)
    wg4 = model.add_plane_surface(wg4_loop)

    sim_region_pts = [model.add_point((L_sim*cos(2*pi*n/6),L_sim*sin(2*pi*n/6)),mesh_size = K*mesh_size) for n in range(6)]

    curves = []
    port_m0_air1 = [model.add_line(sim_region_pts[0],input_in_pt)]
    curves += port_m0_air1
    curves += in_fixed_m0
    curves += in_opt_m0
    curves += in_opt_mx
    curves += in_fixed_mx
    port_mx_air1 = [model.add_line(input_in_pt_mx, sim_region_pts[0])]
    curves += port_mx_air1
    
    air1 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[1],sim_region_pts[2])]
    port_my_air2 = [model.add_line(sim_region_pts[2],input_ex_pt_my)]
    curves += port_my_air2
    curves += ex_fixed_my
    curves += ex_opt_my
    curves += ex_opt_m0
    curves += ex_fixed_m0
    port_m0_air2 = [model.add_line(input_ex_pt, sim_region_pts[1])]
    curves += port_m0_air2
    
    air2 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    port_mxy_air3 = [model.add_line(sim_region_pts[3],input_in_pt_mxy)]
    curves += port_mxy_air3
    curves += in_fixed_mxy
    curves += in_opt_mxy
    curves += in_opt_my
    curves += in_fixed_my
    port_my_air3 = [model.add_line(input_in_pt_my, sim_region_pts[3])]
    curves += port_my_air3
    air3 = model.add_plane_surface(model.add_curve_loop(curves))
    
    curves = []
    curves += [model.add_line(sim_region_pts[4],sim_region_pts[5])]
    port_mx_air4 = [model.add_line(sim_region_pts[5],input_ex_pt_mx)]
    curves += port_mx_air4
    curves += ex_fixed_mx
    curves += ex_opt_mx
    curves += ex_opt_mxy
    curves += ex_fixed_mxy
    port_mxy_air4 = [model.add_line(input_ex_pt_mxy, sim_region_pts[4])]
    curves += port_mxy_air4
    air4 = model.add_plane_surface(model.add_curve_loop(curves))

    model.synchronize()
    model.add_physical([wg1,wg2,wg3,wg4], "waveguide")
    model.add_physical([air1, air2, air3, air4],"air")
    model.add_physical(input_bd_m0 + port_m0_air1 + port_m0_air2,"input")
    model.add_physical(input_bd_mxy + port_mxy_air3 + port_mxy_air4, "output1")
    model.add_physical(input_bd_mx + port_mx_air1 + port_mx_air4, "output21")
    model.add_physical(input_bd_my + port_my_air2 + port_my_air3, "output22")
    model.add_physical(ex_opt_m0 + ex_opt_mx + ex_opt_my + ex_opt_mxy, "ex_opt")
    model.add_physical(in_opt_m0 + in_opt_mx + in_opt_my + in_opt_mxy, "in_opt")
    # gmsh.model.mesh.setOrder(2)
    mesh = model.generate_mesh(dim=2,order=2)
    gmsh.write(name+'.msh')
    gmsh.clear()
    geom.__exit__()
        # # model.add_physical(ex_opt_lines, "ex opt lines")
    
    # mesh.write(name+".msh",file_format = 'gmsh22', binary=False)

    # Write matlab file
    port_id = np.array([3,4,5,6], dtype=int)
    port_analytic = np.array([1,1,1,1], dtype=int)
    port_coords = np.array(
        [
            wg_input_cen,
            wg_input_cen*mxy, 
            wg_input_cen*mx,
            wg_input_cen*my
        ]
    )
    port_coords = np.vstack([port_coords.transpose(), np.zeros((1,4))])
    
    e3 = np.array([cos(pi/3), sin(pi/3)])
    port_x_axis = np.array(
        [
            e2,
            -e2,
            e3,
            -e3
        ]
    )
    port_x_axis = np.vstack([port_x_axis.transpose(), np.zeros((1,4))])

    incident_on = np.array([1,0,0,0], dtype=int)

    # get beta
    # print(sqrt(epsilon_amb*mu_amb))
    # print(sqrt(epsilon_wg*mu_wg))
    # print(w_wg)
    # print([freq*2*pi])
    beta = get_beta_from_k0(n0=sqrt(epsilon_amb*mu_amb), n1 = sqrt(epsilon_wg*mu_wg),
            h = w_wg, k0_sweep=[freq*2*pi])[0]
    print("beta:%f"%(beta))
    port_other_params = np.array(
        [
            [1, 1, 1, 1],   # amplitude
            [w_wg, w_wg, w_wg, w_wg], # width
            [0, 0, 0, 0],
            [beta, beta, beta, beta],
            [epsilon_wg, epsilon_wg, epsilon_wg, epsilon_wg],
            [epsilon_amb, epsilon_amb, epsilon_amb, epsilon_amb],
            [mu_wg, mu_wg, mu_wg, mu_wg],
            [mu_amb, mu_amb, mu_amb, mu_amb]
        ]
    )
    data = {}
    data["port_id"] = port_id
    data["port_x_axis"] = port_x_axis
    data["port_coord"] = port_coords
    data["incident_on"] = incident_on
    data["port_other_params"] = port_other_params
    data["port_is_analytic"] = port_analytic
    data["in_Bezier_P"] = in_points
    data["in_Bezier_V"] = in_tangent
    data["ex_Bezier_P"] = ex_points
    data["ex_Bezier_V"] = ex_tangent
    spio.savemat(name + '.mat', data)

    # varname_value pair
    varname_value_pair = [
        ["mat_file", name+'.mat'],
        ["result_mat_name", name+"_result.mat"],
        ["mesh_file", name+'.msh'],
        ["out_mesh_file", 'out_' + name +'.msh'],
        ["port_id_varname", "port_id"],
        ["port_analytic_varname", "port_is_analytic"],
        ["port_coords_varname", "port_coord"],
        ["port_other_params_varname", "port_other_params"],
        ["port_x_axis_varname", "port_x_axis"],
        ["port_incidence_on_varname","incident_on"],
        ["start_point_mode","0"],
        ["end_point_mode","0"],
        ["result_mat",name+'_result.mat'],
        ["epsilon_wg","%f"%(epsilon_wg)],
        ["epsilon_amb","%f"%(epsilon_amb)],
        ["mu_wg","%f"%(mu_wg)],
        ["mu_amb","%f"%(mu_amb)],
        ['freq',"%f"%(freq)]
    ]
    basic_config("infile.pot", varname_value_pair)

def cross_waveguide(
    DOFs:np.ndarray = np.array([0,0]),
    name = "cross_waveguide",
    L_sim = 10,
    L_wg_fixed = 3,
    w_wg = 1,
    N_BPoints = 50,
    epsilon_wg = 13,
    mu_wg = 1,
    epsilon_amb = 1,
    mu_amb = 1,
    freq = 1,
    use_original = False
):  
    print("Building geometry ...")
    #print(DOFs)
    mesh_size = 0.2*w_wg
    K = 10
    port_refine = 3
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
    
    if(use_original):
        x0 = w_wg/2/sin(pi/6)
        y0 = w_wg/2/cos(pi/6)
        ex_origin = np.linspace(wg_fixed_bd_ex,[0,y0], N_ctrl + 2)
        in_origin = np.linspace(wg_fixed_bd_in,[x0, 0], N_ctrl + 2)
        in_origin_dofs = in_origin[1:-1].flatten()
        ex_origin_dofs = ex_origin[1:-1].flatten()
        DOFs = np.hstack([[x0],in_origin_dofs, [y0], ex_origin_dofs])
        return DOFs
    
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

    # mesh
    # if(os.path.exists(name+'.msh')):
    #     os.system("rm %s"%(name+'.msh'))
    os.system("rm %s*"%(name))
    # Initialize empty geometry using the build in kernel in GMSH
    # with pygmsh.geo.Geometry() as geom:
    geom = pygmsh.geo.Geometry()
    model = geom.__enter__()
    mx = np.array(
        [1, -1]
    )
    my = np.array(
        [-1,1]
    )
    mxy = -1

    origin = model.add_point((0,0,0), mesh_size=mesh_size)

    input_ex_pt = model.add_point(wg_input_bd_ex, mesh_size = mesh_size/port_refine)
    input_in_pt = model.add_point(wg_input_bd_in, mesh_size = mesh_size/port_refine)
    ex_opt_pts = [model.add_point(ex_line_pts[i,:], mesh_size = mesh_size) for i in range(N_BPoints)]
    in_opt_pts = [model.add_point(in_line_pts[i,:], mesh_size = mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mx = model.add_point(wg_input_bd_ex*mx, mesh_size = mesh_size/port_refine)
    input_in_pt_mx = model.add_point(wg_input_bd_in*mx, mesh_size = mesh_size/port_refine)
    ex_opt_pts_mx = [model.add_point(ex_line_pts[i,:]*mx, mesh_size = mesh_size) for i in range(N_BPoints)] 
    in_opt_pts_mx = [model.add_point(in_line_pts[i,:]*mx, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [in_opt_pts[-1]]

    input_ex_pt_my = model.add_point(wg_input_bd_ex*my, mesh_size = mesh_size/port_refine)
    input_in_pt_my = model.add_point(wg_input_bd_in*my, mesh_size = mesh_size/port_refine)
    ex_opt_pts_my = [model.add_point(ex_line_pts[i,:]*my, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [ex_opt_pts[-1]]
    in_opt_pts_my = [model.add_point(in_line_pts[i,:]*my, mesh_size = mesh_size) for i in range(N_BPoints)]

    input_ex_pt_mxy = model.add_point(wg_input_bd_ex*mxy, mesh_size = mesh_size/port_refine)
    input_in_pt_mxy = model.add_point(wg_input_bd_in*mxy, mesh_size = mesh_size/port_refine)
    ex_opt_pts_mxy = [model.add_point(ex_line_pts[i,:]*mxy, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [ex_opt_pts_mx[-1]]
    in_opt_pts_mxy = [model.add_point(in_line_pts[i,:]*mxy, mesh_size = mesh_size) for i in range(N_BPoints-1)] + [in_opt_pts_my[-1]]

    # lines
    curves = []
    ex_opt_m0 = [model.add_line(ex_opt_pts[i], ex_opt_pts[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_m0
    ex_fixed_m0 = [model.add_line(ex_opt_pts[0],input_ex_pt)]
    curves += ex_fixed_m0
    input_bd_m0 = [model.add_line(input_ex_pt, input_in_pt)]
    curves += input_bd_m0
    in_fixed_m0 = [model.add_line(input_in_pt, in_opt_pts[0])]
    curves += in_fixed_m0
    in_opt_m0 = [model.add_line(in_opt_pts[i], in_opt_pts[i+1]) for i in range(N_BPoints - 1)]
    curves += in_opt_m0
    curves += [model.add_line(in_opt_pts[-1], origin), model.add_line(origin, ex_opt_pts[-1])]
    m0_curves = curves
    wg1_loop = model.add_curve_loop(curves)
    # print(curves)
    wg1 = model.add_plane_surface(wg1_loop)

    

    curves = []
    ex_opt_mxy = [model.add_line(ex_opt_pts_mxy[i], ex_opt_pts_mxy[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += ex_opt_mxy
    ex_fixed_mxy = [model.add_line(ex_opt_pts_mxy[0],input_ex_pt_mxy)]
    curves += ex_fixed_mxy
    input_bd_mxy = [model.add_line(input_ex_pt_mxy, input_in_pt_mxy)]
    curves += input_bd_mxy
    in_fixed_mxy = [model.add_line(input_in_pt_mxy, in_opt_pts_mxy[0])]
    curves += in_fixed_mxy
    in_opt_mxy = [model.add_line(in_opt_pts_mxy[i], in_opt_pts_mxy[i+1]) for i in range(N_BPoints - 1)]
    # print(len(in_opt_pts_mxy))
    curves += in_opt_mxy
    curves += [model.add_line(in_opt_pts_mxy[-1], origin), model.add_line(origin, ex_opt_pts_mxy[-1])]
    # print(curves)
    mxy_curves = curves
    wg3_loop = model.add_curve_loop(curves)
    wg3 = model.add_plane_surface(wg3_loop)

    curves = []
    in_opt_mx = [model.add_line(in_opt_pts_mx[i], in_opt_pts_mx[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_mx
    in_fixed_mx = [model.add_line(in_opt_pts_mx[0],input_in_pt_mx)]
    curves += in_fixed_mx
    input_bd_mx = [model.add_line(input_in_pt_mx, input_ex_pt_mx)]
    curves += input_bd_mx
    ex_fixed_mx = [model.add_line(input_ex_pt_mx, ex_opt_pts_mx[0])]
    curves += ex_fixed_mx
    ex_opt_mx = [model.add_line(ex_opt_pts_mx[i], ex_opt_pts_mx[i+1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_mx
    curves += [ - mxy_curves[-1], - m0_curves[-2]]
    wg2_loop = model.add_curve_loop(curves)
    wg2 = model.add_plane_surface(wg2_loop)

    curves = []
    in_opt_my = [model.add_line(in_opt_pts_my[i], in_opt_pts_my[i-1]) for i in range(N_BPoints - 1, 0, -1)]
    curves += in_opt_my
    in_fixed_my = [model.add_line(in_opt_pts_my[0],input_in_pt_my)]
    curves += in_fixed_my
    input_bd_my = [model.add_line(input_in_pt_my, input_ex_pt_my)]
    curves += input_bd_my
    ex_fixed_my = [model.add_line(input_ex_pt_my, ex_opt_pts_my[0])]
    curves += ex_fixed_my
    ex_opt_my = [model.add_line(ex_opt_pts_my[i], ex_opt_pts_my[i+1]) for i in range(N_BPoints - 1)]
    curves += ex_opt_my
    curves += [- m0_curves[-1], -mxy_curves[-2]]
    # print(curves)
    wg4_loop = model.add_curve_loop(curves)
    wg4 = model.add_plane_surface(wg4_loop)

    sim_region_pts = [model.add_point((L_sim*cos(2*pi*n/6),L_sim*sin(2*pi*n/6)),mesh_size = K*mesh_size) for n in range(6)]

    curves = []
    port_m0_air1 = [model.add_line(sim_region_pts[0],input_in_pt)]
    curves += port_m0_air1
    curves += in_fixed_m0
    curves += in_opt_m0
    curves += in_opt_mx
    curves += in_fixed_mx
    port_mx_air1 = [model.add_line(input_in_pt_mx, sim_region_pts[0])]
    curves += port_mx_air1
    
    air1 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    curves += [model.add_line(sim_region_pts[1],sim_region_pts[2])]
    port_my_air2 = [model.add_line(sim_region_pts[2],input_ex_pt_my)]
    curves += port_my_air2
    curves += ex_fixed_my
    curves += ex_opt_my
    curves += ex_opt_m0
    curves += ex_fixed_m0
    port_m0_air2 = [model.add_line(input_ex_pt, sim_region_pts[1])]
    curves += port_m0_air2
    
    air2 = model.add_plane_surface(model.add_curve_loop(curves))

    curves = []
    port_mxy_air3 = [model.add_line(sim_region_pts[3],input_in_pt_mxy)]
    curves += port_mxy_air3
    curves += in_fixed_mxy
    curves += in_opt_mxy
    curves += in_opt_my
    curves += in_fixed_my
    port_my_air3 = [model.add_line(input_in_pt_my, sim_region_pts[3])]
    curves += port_my_air3
    air3 = model.add_plane_surface(model.add_curve_loop(curves))
    
    curves = []
    curves += [model.add_line(sim_region_pts[4],sim_region_pts[5])]
    port_mx_air4 = [model.add_line(sim_region_pts[5],input_ex_pt_mx)]
    curves += port_mx_air4
    curves += ex_fixed_mx
    curves += ex_opt_mx
    curves += ex_opt_mxy
    curves += ex_fixed_mxy
    port_mxy_air4 = [model.add_line(input_ex_pt_mxy, sim_region_pts[4])]
    curves += port_mxy_air4
    air4 = model.add_plane_surface(model.add_curve_loop(curves))

    model.synchronize()
    model.add_physical([wg1,wg2,wg3,wg4], "waveguide")
    model.add_physical([air1, air2, air3, air4],"air")
    model.add_physical(input_bd_m0 + port_m0_air1 + port_m0_air2,"input")
    model.add_physical(input_bd_mxy + port_mxy_air3 + port_mxy_air4, "output1")
    model.add_physical(input_bd_mx + port_mx_air1 + port_mx_air4, "output21")
    model.add_physical(input_bd_my + port_my_air2 + port_my_air3, "output22")
    model.add_physical(ex_opt_m0 + ex_opt_mx + ex_opt_my + ex_opt_mxy, "ex_opt")
    model.add_physical(in_opt_m0 + in_opt_mx + in_opt_my + in_opt_mxy, "in_opt")
    # gmsh.model.mesh.setOrder(2)
    mesh = model.generate_mesh(dim=2,order=2)
    gmsh.write(name+'.msh')
    gmsh.clear()
    geom.__exit__()
        # # model.add_physical(ex_opt_lines, "ex opt lines")
    
    # mesh.write(name+".msh",file_format = 'gmsh22', binary=False)

    # Write matlab file
    port_id = np.array([3,4,5,6], dtype=int)
    port_analytic = np.array([1,1,1,1], dtype=int)
    port_coords = np.array(
        [
            wg_input_cen,
            wg_input_cen*mxy, 
            wg_input_cen*mx,
            wg_input_cen*my
        ]
    )
    port_coords = np.vstack([port_coords.transpose(), np.zeros((1,4))])
    
    e3 = np.array([cos(pi/3), sin(pi/3)])
    port_x_axis = np.array(
        [
            e2,
            -e2,
            e3,
            -e3
        ]
    )
    port_x_axis = np.vstack([port_x_axis.transpose(), np.zeros((1,4))])

    incident_on = np.array([1,0,0,0], dtype=int)

    # get beta
    # print(sqrt(epsilon_amb*mu_amb))
    # print(sqrt(epsilon_wg*mu_wg))
    # print(w_wg)
    # print([freq*2*pi])
    beta = get_beta_from_k0(n0=sqrt(epsilon_amb*mu_amb), n1 = sqrt(epsilon_wg*mu_wg),
            h = w_wg, k0_sweep=[freq*2*pi])[0]
    print("beta:%f"%(beta))
    port_other_params = np.array(
        [
            [1, 1, 1, 1],   # amplitude
            [w_wg, w_wg, w_wg, w_wg], # width
            [0, 0, 0, 0],
            [beta, beta, beta, beta],
            [epsilon_wg, epsilon_wg, epsilon_wg, epsilon_wg],
            [epsilon_amb, epsilon_amb, epsilon_amb, epsilon_amb],
            [mu_wg, mu_wg, mu_wg, mu_wg],
            [mu_amb, mu_amb, mu_amb, mu_amb]
        ]
    )
    data = {}
    data["port_id"] = port_id
    data["port_x_axis"] = port_x_axis
    data["port_coord"] = port_coords
    data["incident_on"] = incident_on
    data["port_other_params"] = port_other_params
    data["port_is_analytic"] = port_analytic
    data["in_Bezier_P"] = in_points
    data["in_Bezier_V"] = in_tangent
    data["ex_Bezier_P"] = ex_points
    data["ex_Bezier_V"] = ex_tangent
    spio.savemat(name + '.mat', data)

    # varname_value pair
    varname_value_pair = [
        ["mat_file", name+'.mat'],
        ["result_mat_name", name+"_result.mat"],
        ["mesh_file", name+'.msh'],
        ["out_mesh_file", 'out_' + name +'.msh'],
        ["port_id_varname", "port_id"],
        ["port_analytic_varname", "port_is_analytic"],
        ["port_coords_varname", "port_coord"],
        ["port_other_params_varname", "port_other_params"],
        ["port_x_axis_varname", "port_x_axis"],
        ["port_incidence_on_varname","incident_on"],
        ["result_mat",name+'_result.mat'],
        ["epsilon_wg","%f"%(epsilon_wg)],
        ["epsilon_amb","%f"%(epsilon_amb)],
        ["mu_wg","%f"%(mu_wg)],
        ["mu_amb","%f"%(mu_amb)],
        ['freq',"%f"%(freq)]
    ]
    basic_config("infile.pot", varname_value_pair)

def test_waveguide_gen():
    # test waveguide generation
    N_ctrl = 5
    DOFs = cross_waveguide(DOFs= np.zeros((4*N_ctrl + 2)), use_original=True)
    DOFs += np.random.rand(DOFs.shape[0]) * 0
    cross_waveguide(DOFs, freq=0.1)
    
def refine_mesh():
    test_waveguide_gen()
    pot_filename = 'infile.pot'
    varname_value_pair = [
        ['mesh_file','cross_waveguide.msh'],
        ['epsilon_wg','13'],
        ['epsilon_air','1'],
        ['freq', '0.5']
    ]
    basic_config(pot_filename, varname_value_pair)

    os.system("./mesh_refinement.exe -ksp_monitor")

def post_process():
    data = spio.loadmat("cross_waveguide_result.mat")
    Ez = np.zeros_like(data['Ez_0'])
    Power_out = np.zeros_like(data['power_out_0'])
    for j in range(int(data['N_process'][0]+0.5)):
        Ez += data['Ez_%d'%(j)]
        Power_out += data['power_out_%d'%(j)]
    new_data = {}
    new_data['Ez'] = Ez
    new_data['power_out'] = Power_out
    spio.savemat("cross_waveguide_post.mat",new_data)

def visualization():
    # mesh = pyvista.read("out_cross_waveguide.msh")
    # data = spio.loadmat("cross_waveguide_post.mat")
    mesh = pyvista.read("remote_mesh.msh")
    data = spio.loadmat("remote_data.mat")
    mesh['E_norm'] = np.abs(data['Ez'])
    mesh['Ez'] = data['Ez'].real
    mesh.plot(scalars = 'E_norm')
    mesh.plot(scalars = 'Ez')
    print(np.abs(data['power_out'])**2)

def test_optimize():
    test_waveguide_gen()
    run_solver("")
    post_process()

def test():
    my_DOF = [1, 1,2, 2,2, 3,2, 4,2, 4, 2,1, 2,2, 3,2, 5,2]
    my_points, my_expoints = cross_waveguide(my_DOF)
    print(my_points)
    print(my_expoints)

if __name__ == "__main__":
    # test_optimize()
    visualization()
