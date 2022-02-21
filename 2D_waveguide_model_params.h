#include <string>
#include <Eigen/Core>
#include <complex.h>

#ifndef Two_D_WAVEGUIDE_MODEL_PARAMS_H_INCLUDED
#define Two_D_WAVEGUIDE_MODEL_PARAMS_H_INCLUDED

class WaveguideParams
{
public:
    WaveguideParams(std::string getpot_filename);
    std::string mat_file; // .mat file
    std::string mesh_file; //.msh file
    std::string out_mesh_file; //output mesh file
    int n_port; // number of ports
    Eigen::RowVectorXi port_id; // domain id of the ports
    Eigen::Matrix<double, 3, Eigen::Dynamic> port_coords; //coordinates of the reference point of the port
    Eigen::Matrix<double, 3, Eigen::Dynamic> port_x_axis; // unit vector of x axis of the port
    Eigen::RowVectorXi port_analytic; // whether use analytic expressions
    Eigen::RowVectorXi port_incidence_on; // whether the port incidence is on
    Eigen::MatrixXd port_other_params; //other parameters of the port, (amplitude, width, mode_index, beta, epsilon_wg, epsilon_amb, mu_wg, mu_amb)
    void analytic_E_field(Eigen::Vector3d pt, int curr_port, double omega, std::complex<double>& value, bool normalize=true);
};

#endif // Two_D_WAVEGUIDE_MODEL_PARAMS_H_INCLUDED
