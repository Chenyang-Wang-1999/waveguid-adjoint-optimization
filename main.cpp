// libmesh basic
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_base.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/elem.h>
#include <libmesh/boundary_info.h>

// libmesh linalg
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/enum_solver_package.h>

// libmesh fem
#include <libmesh/fe.h>
#include <libmesh/parameters.h>
#include <libmesh/fem_context.h>
#include <libmesh/dof_map.h>

// quadrature
#include <libmesh/quadrature_gauss.h>

// libmesh io

#include <libmesh/getpot.h>
#include <MATio>

// C libraries
#include <math.h>
#include <iostream>
#include <Eigen/Core>
#include <mpi.h>

// Debug mode
#define NO_MESH_REFINEMENT
// #define DEBUG
#ifdef DEBUG
    #define WEXTRAP
#endif // DEBUG

// #define DEBUG_UNEXPECTED_MPI_ABORT


// My functions
#include "EM_mesh_flagging.h"
#include "2D_waveguide_model_params.h"
#include "Bezier.h"




using namespace libMesh;

// Function declaration
void mesh_snap_to_Bezier(Eigen::MatrixXd &Bezier_p, Eigen::MatrixXd &Bezier_V, MeshBase & mesh, subdomain_id_type opt_bd_id);
void integrate_along_side(
                      EquationSystems & es,
                      std::vector<subdomain_id_type> & selected_domain,
                      Eigen::VectorXcd & integration_value,
                      Number (*integrand)( Parameters&,
                        FEMContext&),
                      unsigned int parent_domain_id=0
                    );
Number N_poav_integrand( Parameters& params,
                FEMContext& context);
Number mode_decompose(Parameters & params, FEMContext & context );
void get_data_vec(System &sys, Eigen::VectorXcd *sol, int varid);
void assemble_waveguide(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name));
void assemble_QOI_rhs(System & system, const Parameters & params);

void assemble_sensitivity_K(System & system, Parameters & params,
        SparseMatrix<Number> & dKdp, unsigned int point_id, unsigned int dim_id,
        boundary_id_type boundary_id, const Eigen::MatrixXd & Bezier_p, const Eigen::MatrixXd & Bezier_v);

int main(int argc, char ** argv)
{
    // Initialize
    LibMeshInit init(argc, argv);

    GetPot infile("infile.pot");

    WaveguideParams curr_params("infile.pot");

    // create mesh
    Mesh mesh(init.comm());

    // load init mesh
    // std::string mesh_file = infile("mesh_file", std::string("mesh.msh"));
    mesh.read(curr_params.mesh_file);
    BoundaryInfo mesh_bd= mesh.get_boundary_info();
    // refine
    MeshRefinement my_refinement(mesh);

    // set flagging
    double epsilon_wg = infile("epsilon_wg", 13.0);
    double epsilon_air = infile("epsilon_amb",1.0);
    double mu_wg = infile("mu_wg", 1.0);
    double mu_amb = infile("mu_amb", 1.0);
    double freq = infile("freq",1.0);
    unsigned int save_field = infile("save_field",0);
    std::string result_mat_name = infile("result_mat_name",std::string());
    subdomain_id_type id_wg = mesh.get_id_by_name("waveguide");
    subdomain_id_type id_input = mesh_bd.get_id_by_name("input");

    // subdomain_id_type id_output1 = mesh_bd.get_id_by_name("output1");
    // subdomain_id_type id_output21 = mesh_bd.get_id_by_name("output21");
    subdomain_id_type id_ex_opt = mesh_bd.get_id_by_name("ex_opt");
    subdomain_id_type id_in_opt = mesh_bd.get_id_by_name("in_opt");

    Eigen::MatrixXd in_Bezier_P, in_Bezier_V, ex_Bezier_P, ex_Bezier_V;
    matio::MatioFile matfp(curr_params.mat_file);
    matfp.read_mat("in_Bezier_P", in_Bezier_P);
    matfp.read_mat("in_Bezier_V", in_Bezier_V);

    matfp.read_mat("ex_Bezier_P", ex_Bezier_P);
    matfp.read_mat("ex_Bezier_V", ex_Bezier_V);
    matfp.close();

    std::vector<subdomain_id_type> model_domain_id = {id_wg};
    std::vector<double> model_index = {1.0/(8*freq/8)> sqrt(epsilon_wg) ? 1.0/(8*freq/8): sqrt(epsilon_wg)};

    std::vector<subdomain_id_type> port_id;

    for(unsigned int j = 0 ; j < curr_params.port_id.size(); j ++)
    {
        port_id.push_back(curr_params.port_id(j));
    }

    EMMeshFlagging my_model(&mesh,model_domain_id, model_index, sqrt(epsilon_air),
            port_id, 1.0/20);
    my_model.set_freq(freq);

#ifndef NO_MESH_REFINEMENT
    unsigned int refinement_times = 0;
    bool whether_continue;
    do
    {
        // refine
        my_refinement.flag_elements_by(my_model);
        whether_continue = my_refinement.refine_elements();
        refinement_times ++;
        std::cout << "Curr refinement: " << refinement_times <<'\n';
    }while(whether_continue);
#endif // NO_MESH_REFINEMENT

    mesh.add_node_datum<double>("Bezier_t");


    mesh_snap_to_Bezier(ex_Bezier_P, ex_Bezier_V, mesh, id_ex_opt);

    mesh_snap_to_Bezier(in_Bezier_P, in_Bezier_V, mesh, id_in_opt);

    // output
    mesh.write(curr_params.out_mesh_file);

    // MPI_Barrier(init.comm().get());
    MPI_Barrier(MPI_COMM_WORLD);

    // eq systems
        // eq systems
    EquationSystems eq_systems(mesh);
    eq_systems.add_system<LinearImplicitSystem>("waveguide");
    eq_systems.get_system("waveguide").add_variable("Ez", SECOND);
    eq_systems.get_system("waveguide").attach_assemble_function(assemble_waveguide);
    eq_systems.parameters.set<double>("epsilon_wg") = epsilon_wg;
    eq_systems.parameters.set<double>("epsilon_amb") = epsilon_air;
    eq_systems.parameters.set<double>("mu_wg") = mu_wg;
    eq_systems.parameters.set<double>("mu_amb") = mu_amb;
    eq_systems.parameters.set<double>("freq") = freq;
    eq_systems.parameters.set<WaveguideParams*>("waveguide_params") = &curr_params;
    eq_systems.parameters.set<subdomain_id_type>("wg_id") = id_wg;
    unsigned short Ez_id = eq_systems.get_system("waveguide").variable_number("Ez");
    std::cout << "Id of Ez: "<<Ez_id <<'\n';
    FEMContext context(eq_systems.get_system("waveguide"),{Ez_id});

    eq_systems.init();
    eq_systems.print_info();

    // solve
    eq_systems.get_system("waveguide").solve();
    // assemble_waveguide(eq_systems,"waveguide");


    // MPI_Barrier(init.comm().get());
    MPI_Barrier(MPI_COMM_WORLD);


    // Get S params and Solve adjoint problem
    // Get S params
    Eigen::VectorXcd S_params(4);
    std::vector<subdomain_id_type> selected_domains = {3,4,5,6};
    integrate_along_side(eq_systems, selected_domains, S_params, mode_decompose);

    // Solve adjoint problem
    unsigned int selected_port = 1; // output 1
    eq_systems.parameters.set<Complex>("curr_S") = S_params(selected_port);
    eq_systems.parameters.set<subdomain_id_type>("QOI_port") = selected_port;
    assemble_QOI_rhs(eq_systems.get_system("waveguide"), eq_systems.parameters);
    eq_systems.get_system("waveguide").assemble_before_solve = false;
    eq_systems.get_system("waveguide").adjoint_solve();

    // Get all sensitivity
    Eigen::MatrixXd ex_opt_sens(ex_Bezier_P.rows(), ex_Bezier_P.cols());
    Eigen::MatrixXd in_opt_sens(in_Bezier_P.rows(), in_Bezier_P.cols());


    // std::unique_ptr< NumericVector< Number > > adjoint_copy = eq_systems.get_system("waveguide").get_adjoint_solution(0).clone();
    System & curr_system = eq_systems.get_system("waveguide");
    curr_system.add_vector("middle_vec");
    curr_system.add_matrix("dKdp");
    NumericVector<Number> & middle_vec = curr_system.get_vector("middle_vec");
    for(unsigned int dim_id = 0; dim_id < ex_Bezier_P.cols(); dim_id++)
    {
        for(unsigned int point_id = 0; point_id < ex_Bezier_P.rows(); point_id++)
        {
            // SparseMatrix<Number> dKdp(init.comm());
            // NumericVector<Number> middle_vec(init.comm());

            curr_system.get_vector("middle_vec").zero();
            curr_system.get_matrix("dKdp").zero();
            assemble_sensitivity_K(eq_systems.get_system("waveguide"), eq_systems.parameters, curr_system.get_matrix("dKdp"), point_id, dim_id,
                    id_ex_opt, ex_Bezier_P, ex_Bezier_V);
            curr_system.get_matrix("dKdp").vector_mult(curr_system.get_vector("middle_vec"), eq_systems.get_system("waveguide").get_vector("Ez"));

            ex_opt_sens(point_id, dim_id) = std::real(eq_systems.get_system("waveguide").get_adjoint_solution(0).dot(middle_vec));
        }
    }

    for(unsigned int dim_id = 0; dim_id < in_Bezier_P.cols(); dim_id++)
    {
        for(unsigned int point_id = 0; point_id < in_Bezier_P.rows(); point_id++)
        {
            curr_system.get_vector("middle_vec").zero();
            curr_system.get_matrix("dKdp").zero();
            assemble_sensitivity_K(eq_systems.get_system("waveguide"), eq_systems.parameters, curr_system.get_matrix("dKdp"), point_id, dim_id,
                    id_in_opt, in_Bezier_P, in_Bezier_V);
            curr_system.get_matrix("dKdp").vector_mult(curr_system.get_vector("middle_vec"), eq_systems.get_system("waveguide").get_vector("Ez"));

            in_opt_sens(point_id, dim_id) = std::real(eq_systems.get_system("waveguide").get_adjoint_solution(0).dot(middle_vec));
        }
    }

    /*
     *
    std::vector<Complex> sol_vec;
    eq_systems.build_solution_vector(sol_vec);
    std::cout<<sol_vec.size()<<'\n';
    //*
    Eigen::Map<Eigen::VectorXcd> sol_VectorXcd(sol_vec.data(),sol_vec.size());
    //*/
    processor_id_type curr_id = mesh.processor_id();

    static Eigen::VectorXcd sol(eq_systems.get_system("waveguide").n_dofs());
    std::stringstream varname_stream;
    varname_stream << "Ez_" << curr_id;
    std::string varname = varname_stream.str();
    if(save_field)
    {
        // sol.setZero();
        get_data_vec(eq_systems.get_system("waveguide"), &sol, Ez_id);

    }




    std::stringstream integrate_stream;
    integrate_stream << "S_params_" << curr_id;
    std::string integrate_name = integrate_stream.str();
    // Output
    for(unsigned int i = 0; i < mesh.n_processors(); i++)
    {
        // save to mat file
        // std::cout << "i: "<<i << "current_id: " << curr_id << '\n';
        if(i == curr_id)
        {
            matio::MatioFile matfp(result_mat_name);
            if(save_field)
            {
                matfp.write_mat(varname.data(), sol,true);
            }
            matfp.write_mat(integrate_name.data(), S_params, true);
            matfp.close();
        }
        // MPI_Barrier(init.comm().get());
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(curr_id == 0)
    {
        matio::MatioFile matfp(result_mat_name);
        Eigen::VectorXi N_processor(1);
        N_processor(0) = mesh.n_processors();
        matfp.write_mat("N_process", N_processor, true);
        matfp.write_mat("ex_opt_sens", ex_opt_sens, true);
        matfp.write_mat("in_opt_sens", in_opt_sens, true);
        matfp.close();
    }
    // MPI_Finalize();
    return 0;
}

void mesh_snap_to_Bezier(Eigen::MatrixXd &Bezier_p, Eigen::MatrixXd &Bezier_V, MeshBase & mesh, subdomain_id_type opt_bd_id)
{
    BoundaryInfo & mesh_boundary = mesh.get_boundary_info();
    for(Node* node: mesh.local_node_ptr_range())
    {
        if(mesh_boundary.has_boundary_id(node, opt_bd_id))
        {
            Eigen::VectorXd node_pt(2);
            bool x_posi = ((*node)(0) > 0), y_posi = ((*node)(1) > 0);
            // convert to absolute value
            node_pt(0) = x_posi?(*node)(0):(-(*node)(0));
            node_pt(1) = y_posi?(*node)(1):(-(*node)(1));

            double t = Nearest_t(Bezier_p, Bezier_V, node_pt, 1e-5);

            // append t
            node->set_extra_datum(0, t);
            Eigen::VectorXd t_vec(1);
            t_vec(0) = t;
            Eigen::MatrixXd curr_pt;

            curr_pt = Bezier_interp(Bezier_p, Bezier_V, t_vec, 0);

            // reset node
            (*node)(0) = x_posi?curr_pt(0,0):-curr_pt(0,0);
            (*node)(1) = y_posi?curr_pt(0,1):-curr_pt(0,1);
        }
    }
}



void integrate_along_side(
                      EquationSystems & es,
                      std::vector<subdomain_id_type> & selected_domain,
                      Eigen::VectorXcd & integration_value,
                      Number (*integrand)( Parameters&,
                        FEMContext&),
                      unsigned int parent_domain_id
                    )
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  System & sys = es.get_system("waveguide");

  FEMContext context(sys);
  const MeshBase & mesh = sys.get_mesh();

  // Get a constant reference to the mesh object.
  //const MeshBase & mesh = es.get_mesh();
  //const BoundaryInfo & mesh_boundary = mesh.get_boundary_info();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // initialize the solution vector
  integration_value.resize(selected_domain.size());
  integration_value.setZero();

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      context.pre_fe_reinit(sys,elem);
      context.elem_fe_reinit();
      bool use_curr_elem = true;
      if(parent_domain_id)
      {
        use_curr_elem = (context.get_elem().subdomain_id() == parent_domain_id);
      }

      if(use_curr_elem)
      {

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (auto side : elem->side_index_range())
        {
         // get boundary id
            context.side = side;
            context.side_fe_reinit();
            for(unsigned int integral_id = 0; integral_id < selected_domain.size(); integral_id ++)
            {
                if(context.has_side_boundary_id(selected_domain[integral_id]))
                {
                    // has current id
                    integration_value(integral_id) += integrand(es.parameters,context);
                    //std::cout << "Current id: " << selected_domain[integral_id];
                    // std::cout << integration_value(integral_id) <<'\n';
                }
            }
        }
      }
  // All done!
  }
}

Number N_poav_integrand( Parameters& params,
                FEMContext& context)
{
    subdomain_id_type wg_id = params.get<subdomain_id_type>("wg_id");
    double epsilon, mu;
    if(context.get_elem().subdomain_id() == wg_id)
    {
        epsilon = params.get<double>("epsilon_wg");
        mu = params.get<double>("mu_wg");
    }
    else
    {
        epsilon = params.get<double>("epsilon_amb");
        mu = params.get<double>("mu_amb");
    }
    double freq = params.get<double> ("freq");
    double omega = 2*M_PI*freq;

    // get boundary fe
    FEBase* fe_face;
    fe_face = context.get_side_fe(0); // var0: Ez

    // The XYZ locations (in physical space) of the
    // quadrature points on the face.  This is where
    // we will interpolate the boundary value function.
    const std::vector<Point> & q_point = fe_face->get_xyz();
    const std::vector<Real> & JxW = fe_face -> get_JxW();
    const std::vector<Point> & normal = fe_face -> get_normals();

    // Compute the shape function values on the element
    // face.

    Number result = 0.0;

    Point averged_normal_vec(0,0,0);
    for(unsigned int i = 0; i < normal.size(); i++)
    {
        averged_normal_vec += normal[i];
    }
    averged_normal_vec /= normal.size();

    for(unsigned int qp = 0; qp < JxW.size(); qp ++)
    {
        Number Ez_qp = context.side_value(0,qp);
        Gradient Ez_grad = context.side_gradient(0, qp);

        Number Hx_qp = - (Complex(1j)/(omega*mu)) * Ez_grad(1);
        Number Hy_qp = (Complex(1j)/(omega*mu)) * Ez_grad(0);

        Real Poavx_qp = - (std::conj(Ez_qp)*Hy_qp).real()/2;
        Real Poavy_qp = (std::conj(Ez_qp)*Hx_qp).real()/2;

        Real nx = averged_normal_vec(0);
        Real ny = averged_normal_vec(1);
        Real normal_norm = sqrt(pow(nx,2) + pow(ny,2));
        Real nPoav = Poavx_qp*nx/normal_norm + Poavy_qp*ny/normal_norm;
        result += nPoav*JxW[qp];
        /* debug
        std::cout << "Ez: " << Ez_qp <<'\n';
        std::cout << "Hx: " << Hx_qp <<'\n';
        std::cout << "Hy: " << Hy_qp << '\n';
        std::cout << "Poavx: " << Poavx_qp <<'\n';
        std::cout << "Poavy: " << Poavy_qp << '\n';
        std::cout << "current_normal: "<< nx <<", "<< ny <<'\n';
        std::cout << "norm of normal: "<< normal_norm <<'\n';
        std::cout << "nPoav: " << nPoav <<'\n';
        std::cout <<"result: "<<result <<'\n';
        //*/
    }
    return result;
}

Number mode_decompose(Parameters & params, FEMContext & context )
{
    // decompose the field to the corresponding port mode
    subdomain_id_type wg_id = params.get<subdomain_id_type>("wg_id");
    double epsilon, mu;
    if(context.get_elem().subdomain_id() == wg_id)
    {
        epsilon = params.get<double>("epsilon_wg");
        mu = params.get<double>("mu_wg");
    }
    else
    {
        epsilon = params.get<double>("epsilon_amb");
        mu = params.get<double>("mu_amb");
    }
    double freq = params.get<double> ("freq");
    double omega = 2*M_PI*freq;

    WaveguideParams * wg_params = params.get<WaveguideParams*>("waveguide_params");

    // get boundary fe
    FEBase* fe_face;
    fe_face = context.get_side_fe(0); // var0: Ez

    // The XYZ locations (in physical space) of the
    // quadrature points on the face.  This is where
    // we will interpolate the boundary value function.
    const std::vector<Point> & q_point = fe_face->get_xyz();
    const std::vector<Real> & JxW = fe_face -> get_JxW();
    const std::vector<Point> & normal = fe_face -> get_normals();

    // Compute the shape function values on the element
    // face.

    Number result = 0.0;

    for(unsigned int curr_port = 0; curr_port < wg_params ->n_port; curr_port ++)
    {
        if(context.has_side_boundary_id(wg_params->port_id(curr_port)))
        {
            for(unsigned int qp = 0; qp < JxW.size(); qp ++)
            {
                Number Ez_qp = context.side_value(0,qp);
                Gradient Ez_grad = context.side_gradient(0, qp);

                Number Hx_qp = - (Complex(1j)/(omega*mu)) * Ez_grad(1);
                Number Hy_qp = (Complex(1j)/(omega*mu)) * Ez_grad(0);

                Eigen::Vector3d pt_vec = {q_point[qp](0), q_point[qp](1), 0};
                Number mode_Dz_qp, mode_Bx_qp, mode_By_qp;
                wg_params ->analytic_E_field(pt_vec, curr_port, omega, mode_Dz_qp, true, WaveguideModes::Dz);
                wg_params ->analytic_E_field(pt_vec, curr_port, omega, mode_Bx_qp, true, WaveguideModes::Bx);
                wg_params ->analytic_E_field(pt_vec, curr_port, omega, mode_By_qp, true, WaveguideModes::By);

                result += (std::conj(mode_Dz_qp)*Ez_qp + std::conj(mode_Bx_qp)*Hx_qp + std::conj(mode_By_qp)*Hy_qp)*JxW[qp]/4;
                /* debug
                std::cout << "Ez: " << Ez_qp <<'\n';
                std::cout << "Hx: " << Hx_qp <<'\n';
                std::cout << "Hy: " << Hy_qp << '\n';
                std::cout << "Poavx: " << Poavx_qp <<'\n';
                std::cout << "Poavy: " << Poavy_qp << '\n';
                std::cout << "current_normal: "<< nx <<", "<< ny <<'\n';
                std::cout << "norm of normal: "<< normal_norm <<'\n';
                std::cout << "nPoav: " << nPoav <<'\n';
                std::cout <<"result: "<<result <<'\n';
                //*/
            }
        }
    }
    return result;
}

void get_data_vec(System &sys, Eigen::VectorXcd *sol, int varid)
{
    // get data vector by FEM context
    const MeshBase & mesh = sys.get_mesh();
    const DofMap & dof_map = sys.get_dof_map();
    const BoundaryInfo & mesh_boundary = mesh.get_boundary_info();
    NumericVector<Number> & sys_sol = *sys.solution;
    sol->setZero();

    // localize
    int i = 0;
    for(const auto node: mesh.local_node_ptr_range())
    {
        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices(node, dof_indices, varid);
        std::vector<Number> node_sol(dof_indices.size());
        sys_sol.get(dof_indices, node_sol);

        // average
        (*sol)(node->id()) = 0;
        for(i = 0; i<dof_indices.size(); i++)
        {
            (*sol)(node->id()) += node_sol[i];

        }
        (*sol)(node->id())/=dof_indices.size();
    }
}

void assemble_waveguide(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "waveguide");

  // Get system parameters
  WaveguideParams * waveguide_params = es.parameters.get<WaveguideParams *>("waveguide_params");
  double epsilon_wg = es.parameters.get<double>("epsilon_wg");
  double mu_wg = es.parameters.get<double>("mu_wg");
  double epsilon_amb = es.parameters.get<double>("epsilon_amb");
  double mu_amb = es.parameters.get<double>("mu_amb");
  double freq = es.parameters.get<double>("freq");
  double omega = 2*M_PI*freq;
  subdomain_id_type wg_id = es.parameters.get<subdomain_id_type>("wg_id");

  // waveguide propagation constants
  // std::cout << waveguide_params->n_port;
  std::vector<double> beta(waveguide_params->n_port);
  std::complex<double> curr_value;
  int port_index;

  #ifndef DEBUG_UNEXPECTED_MPI_ABORT
  for(port_index = 0; port_index < waveguide_params->n_port; port_index ++)
  {
    beta[port_index] = waveguide_params->port_other_params(3, port_index);
  }

  // FEM parameters
  Real FEM_beta_wg = - pow(omega, 2) * epsilon_wg * mu_wg;
  Real FEM_beta_amb = - pow(omega,2) * epsilon_amb * mu_amb;
  Complex FEM_gamma;
  Complex FEM_q;


  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  const BoundaryInfo & mesh_boundary = mesh.get_boundary_info();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("waveguide");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // The global system matrix
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element ranges are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // range will iterate from the first to the last element on
  // the local processor.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the
  // active_local_element_ptr_range.


  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Cache the number of degrees of freedom on this element, for
      // use as a loop bound later.  We use cast_int to explicitly
      // convert from size() (which may be 64-bit) to unsigned int
      // (which may be 32-bit but which is definitely enough to count
      // *local* degrees of freedom.
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // With one variable, we should have the same number of degrees
      // of freedom as shape functions.
      libmesh_assert_equal_to (n_dofs, phi.size());

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.
      Ke.resize (n_dofs, n_dofs);

      Fe.resize (n_dofs);

      // \int J(x) dW ->\sum JxW(qp)

      // Now loop over the quadrature points.  This handles
      // the numeric integration.

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test functions (i) against
          // the trial functions (j).
          for (unsigned int i=0; i != n_dofs; i++)
          {
          /*
            std::cout << "current_dof_indices: " << dof_indices[i] <<'\n';
            std::cout << "node_indices: "<< elem->node_id(i) << '\n';
            //*/
            for (unsigned int j=0; j != n_dofs; j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                if(elem->subdomain_id() == wg_id)
                {
                    Ke(i,j) += JxW[qp]*FEM_beta_wg*(phi[i][qp]*phi[j][qp]);
                }
                else
                {
                    Ke(i,j) += JxW[qp]*FEM_beta_amb*(phi[i][qp]*phi[j][qp]);
                }
                //std::cout << dphi[i][qp];
              }
          }
        }


        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.

        for (auto side : elem->side_index_range())
        {
         // get boundary id

          if (elem->neighbor_ptr(side) == nullptr)
            {
              // std::cout << qface.n_points()<<'\n';
              // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = fe_face->get_JxW();



              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point> & qface_point = fe_face->get_xyz();


              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);

              // Some shape functions will be 0 on the face, but for
              // ease of indexing and generality of code we loop over
              // them anyway
              libmesh_assert_equal_to (n_dofs, phi_face.size());

              // which boundary?
              bool is_port = false;
              is_port = false;
              for(port_index = 0; port_index < waveguide_params->n_port; port_index ++)
              {
	       if(mesh_boundary.has_boundary_id(elem, side, waveguide_params->port_id(port_index)))
                {
            #ifdef DEBUG
                std::cout << "current id:" << waveguide_params->port_id(port_index) <<'\n';
            #endif // DEBUG
                    is_port = true;
                    // is current port
                    // FEM parameters
                    FEM_gamma = - 1j*beta[port_index];

                    // Assemble K, the same
                    for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // Matrix contribution of the L2 projection.

                      for (unsigned int i=0; i != n_dofs; i++)
                      {
                        for (unsigned int j=0; j != n_dofs; j++)
                        {
                          Ke(i,j) += JxW_face[qp]*FEM_gamma*phi_face[i][qp]*phi_face[j][qp];
                        }
                      // Right-hand-side contribution of the L2
                      }
                    }
                    // Assemble F
                    // whether the incidence is on
                  if(waveguide_params ->port_incidence_on(port_index))
                  {
              #ifdef DEBUG
                    std::cout << waveguide_params -> port_id(port_index) << '\n';
              #endif // DEBUG
                    for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                        Eigen::Vector3d pt = {qface_point[qp](0), qface_point[qp](1), 0.0};
                        std::complex<double> E0;
                        waveguide_params ->analytic_E_field(pt, port_index, omega, E0, true, WaveguideModes::Ez);
                        // FEM_q = (-(2.0*1j)*beta[port_index])*E0;
                    #ifdef DEBUG
                        std::cout << E0<<'\n';
                    #endif // DEBUG
                        FEM_q = FEM_gamma*2*E0;

                        // std::cout << E0 <<','<<FEM_q<<'\n';
                        for(unsigned int i=0; i != n_dofs; i++)
                        {
                            Fe(i) += JxW_face[qp]*FEM_q *phi_face[i][qp];
                        }
                    }
                  }
                  break;
                }
              }

              if(! is_port)
              {
                // not port
                // Loop over the face quadrature points for integration.
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // The location on the boundary of the current


                        /* PEC
                      // The penalty value.  \frac{1}{\epsilon}
                      // in the discussion above.
                      const Real penalty = 1.e20;

                      // Matrix contribution of the L2 projection.
                      for (unsigned int i=0; i != n_dofs; i++)
                        for (unsigned int j=0; j != n_dofs; j++)
                        {
                          Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                        }
                      // Right-hand-side contribution of the L2
                      //*/

                      // Scattering BC
                      FEM_gamma = - Complex(1j)*sqrt(epsilon_amb*mu_amb)*omega;

                        // Assemble K, the same
                        for (unsigned int qp=0; qp<qface.n_points(); qp++)
                        {
                          // Matrix contribution of the L2 projection.

                          for (unsigned int i=0; i != n_dofs; i++)
                          {
                            for (unsigned int j=0; j != n_dofs; j++)
                            {
                              Ke(i,j) += JxW_face[qp]*FEM_gamma*phi_face[i][qp]*phi_face[j][qp];
                            }
                          // Right-hand-side contribution of the L2
                          }
                        }

                    }
              }


              }
            }


      //std::cout << "4\n";
      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The  SparseMatrix::add_matrix()
      // and  NumericVector::add_vector() members do this for us.
      matrix.add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
      //std::cout << "5\n";

    }
    #endif // DEBUG_UNEXPECTED_MPI_ABORT
  // All done!
}

void assemble_QOI_rhs(System & system, const Parameters & params)
{
    NumericVector<Number> * adjoint_rhs;
    if(system.n_qois())
    {
        // qoi has been defined
        adjoint_rhs = &system.get_adjoint_rhs(0);
    }
    else
    {
        adjoint_rhs = &system.add_adjoint_rhs(0);
    }

    // get parameters
    WaveguideParams* wg_params = params.get<WaveguideParams*>("waveguide_params");
    subdomain_id_type which_port = params.get<subdomain_id_type>("QOI_port");
    double freq = params.get<double>("freq");
    Number S = params.get<Number>("curr_S");
    double omega = freq*2*M_PI;
    subdomain_id_type port_id = wg_params ->port_id(which_port);

    // get mesh
    MeshBase & mesh = system.get_mesh();
    BoundaryInfo & mesh_boundary = mesh.get_boundary_info();

    const unsigned int dim = mesh.mesh_dimension();

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH);
    fe_face->attach_quadrature_rule (&qface);

    DenseVector<Number> Fe;
    std::vector<dof_id_type> dof_indices;


    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices (elem, dof_indices);
            const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());
        Fe.resize (n_dofs);

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.

        for (auto side : elem->side_index_range())
        {
         // get boundary id

          if (elem->neighbor_ptr(side) == nullptr)
          {
              if(mesh_boundary.has_boundary_id(elem, side, port_id))
              {
                fe_face->reinit(elem, side);
                // Integration: \int 1/4 (D_mode^*\cdot E + H_mode^* \cdot B)
                // Bx = - i /omega \partial_y E^z
                // By = i / omega  \partial_x E^z
                const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
                const std::vector<std::vector<Real>> & dphidx_face = fe_face->get_dphidx();
                const std::vector<std::vector<Real>> & dphidy_face = fe_face->get_dphidy();
                const std::vector<Real> & JxW_face = fe_face->get_JxW();
                const std::vector<Point> & qface_point = fe_face->get_xyz();

                // Integrate
                for(unsigned int i = 0; i < n_dofs; i++)
                {
                    for(unsigned int qp = 0; qp < JxW_face.size(); qp++)
                    {
                    //IEEE TRANSACTIONS ON MAGNETICS, VOL. 46, NO. 8, AUGUST 2010
                        Eigen::Vector3d pt_vec = {qface_point[qp](0), qface_point[qp](1), qface_point[qp](2)};
                        Complex mode_Dz, mode_Hx, mode_Hy;
                        wg_params->analytic_E_field(pt_vec, which_port, omega, mode_Dz, true, WaveguideModes::Dz);
                        wg_params->analytic_E_field(pt_vec, which_port, omega, mode_Hx, true, WaveguideModes::Hx);
                        wg_params->analytic_E_field(pt_vec, which_port, omega, mode_Hy, true, WaveguideModes::Hy);

                        Number dBxdA = - Complex(1j)/omega * dphidy_face[i][qp];
                        Number dBydA = Complex(1j)/omega * dphidx_face[i][qp];

                        Number dSdA = (std::conj(mode_Dz)*phi_face[i][qp] +
                                    std::conj(mode_Hx)*dBxdA +
                                    std::conj(mode_Hy)*dBydA
                                    )/4;
                        Fe(i) += 2*std::real(std::conj(S)*dSdA) + 2*Complex(1j)*std::imag(std::conj(S)*dSdA);
                    }
                }

              }
          }
        }


      //std::cout << "4\n";
      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_vector(Fe, dof_indices);

      adjoint_rhs->add_vector    (Fe, dof_indices);
      //std::cout << "5\n";

    }
}

void assemble_sensitivity_K(System & system, Parameters & params,
        SparseMatrix<Number> & dKdp, unsigned int point_id, unsigned int dim_id,
        boundary_id_type boundary_id, const Eigen::MatrixXd & Bezier_p, const Eigen::MatrixXd & Bezier_v)
{
    Eigen::MatrixXd Sens_v = Bezier_V_sensitivity(point_id, dim_id, Bezier_p, 0, 1, Bezier_v.row(0), Bezier_v.row(0), true);
    Eigen::MatrixXd Sens_p(Bezier_p.rows(), Bezier_p.cols());
    Sens_p.setZero();
    Sens_p(point_id, dim_id) = 1;

    // loop around the geometry, and integrate around the boundary
    // get parameters
    subdomain_id_type which_boundary = boundary_id;
    subdomain_id_type wg_id = params.get<subdomain_id_type>("wg_id");
    double freq = params.get<double>("freq");
    double epsilon_wg = params.get<double>("epsilon_wg");
    double epsilon_amb = params.get<double>("epsilon_amb");
    double mu_wg = params.get<double>("mu_wg");
    double mu_amb = params.get<double>("mu_amb");
    double omega = freq*2*M_PI;

    // FEM_beta = - pow(omega, 2) * epsilon * mu
    double Delta_FEM_beta = - pow(omega, 2) * epsilon_wg * mu_wg + pow(omega, 2) * epsilon_amb * mu_amb;

    // get mesh
    MeshBase & mesh = system.get_mesh();
    BoundaryInfo & mesh_boundary = mesh.get_boundary_info();

    const unsigned int dim = mesh.mesh_dimension();

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH);
    fe_face->attach_quadrature_rule (&qface);

    DenseMatrix<Number> Ke;
    std::vector<dof_id_type> dof_indices;

    for(Elem* elem: mesh.active_local_element_ptr_range())
    {
        if(elem->subdomain_id() == wg_id)
        {
            for(auto side:elem->side_index_range())
            {
                if(mesh_boundary.has_boundary_id(elem, side, which_boundary))
                {
                    dof_map.dof_indices (elem, dof_indices);
                    const unsigned int n_dofs =
                        cast_int<unsigned int>(dof_indices.size());
                    Ke.resize (n_dofs, n_dofs);

                    // reinitialize and get fe data
                    fe_face->reinit(elem, side);
                    const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
                    const std::vector<std::vector<Real>> & dphidx_face = fe_face->get_dphidx();
                    const std::vector<std::vector<Real>> & dphidy_face = fe_face->get_dphidy();
                    const std::vector<Point> & normal_face = fe_face->get_normals();
                    const std::vector<Real> & JxW_face = fe_face->get_JxW();
                    const std::vector<Point> & qface_point = fe_face->get_xyz();

                    // Bezier interpolation
                    Eigen::VectorXd t_list(n_dofs);
                    for(unsigned int i = 0; i < n_dofs; i++)
                    {
                        t_list(i) = elem->node_ptr(i) ->get_extra_datum<double>(0);
                    }

                    Eigen::MatrixXd P_sens = Bezier_interp(Sens_p, Sens_v, t_list, 0);
                    Eigen::MatrixXd V_tan = Bezier_interp(Bezier_p, Bezier_v, t_list, 1);

                    // normalize
                    for(unsigned int i = 0; i<n_dofs; i++)
                    {
                        double norm_tangent = V_tan.row(i).norm();
                        V_tan.row(i) /= norm_tangent;
                    }

                    // assemble matrix
                    /* For FEM
                            Ke(i,j) = JxW[qp]*(dphi[i][qp]*dphi[j][qp]) + JxW[qp]*FEM_beta*(phi[i][qp]*phi[j][qp]);
                            FEM_beta = - pow(omega, 2) * epsilon * mu
                       For sensitivity
                            dKe(i,j) = \Delta FEM_beta*(phi[i][qp]*phi[j][qp]) dS
                            where \Delta FEM_beta = {FEM_beta_wg - FEM_beta_amb  ,normalxtangent > 0 (+ indicates enlarging the waveguide)
                                                    {- ...
                    */

                    for(unsigned int i = 0; i<n_dofs; i++)
                    {
                        for(unsigned int qp = 0; qp < JxW_face.size(); qp++)
                        {
                            bool enlarge = (normal_face[i](0)*V_tan(i,1) - normal_face[i](1)*V_tan(i, 0)>0); // normal(i) x  V_tan(i)
                            double volume_elem = P_sens(i, 0) * V_tan(i, 1) - P_sens(i, 1)*V_tan(i, 0);     // P_sens(i) x V_tan(i)
                            double curr_FEM_beta = enlarge? Delta_FEM_beta: - Delta_FEM_beta;
                            for(unsigned int j = 0; j < n_dofs; j++)
                            {
                                Ke(i,j) += curr_FEM_beta * (phi_face[i][qp]*phi_face[j][qp]) * volume_elem * JxW_face[qp];

                            }
                        }
                    }// for(i:dofs)

                    dKdp.add_matrix(Ke, dof_indices);

                }
            }
        }
    }



}
