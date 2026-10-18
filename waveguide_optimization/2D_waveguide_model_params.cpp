#include "2D_waveguide_model_params.h"
#include "h5io.h"
#include <Eigen/Core>
#include <libmesh/getpot.h>
#include <string.h>
#include <math.h>
#include <iostream>

// #define DEBUG

WaveguideParams::WaveguideParams(std::string getpot_filename)
{
    // read getpot
    GetPot input_file(getpot_filename);
    design_file = input_file("design_file",std::string("design.h5"));
    result_file = input_file("result_file",std::string("result.h5"));
    mesh_file = input_file("mesh_file",std::string());
    out_mesh_file = input_file("out_mesh_file",std::string());
    // n_port = input_file("n_port",0);
    std::string port_id_varname = input_file("port_id_varname",std::string("port_id"));
    std::string port_analytic_varname = input_file("port_analytic_varname",std::string("is_analytic"));
    std::string port_coords_varname = input_file("port_coords_varname",std::string("coords"));
    std::string port_other_params_varname = input_file("port_other_params_varname",std::string("other_params"));
    std::string port_x_axis_varname = input_file("port_x_axis_varname", std::string("x_axis"));
    std::string port_incidence_on_varname = input_file("port_incidence_on_varname", std::string("incidence_on"));



    // 端口表来自 design.h5
    {
        h5io::File design(design_file, false);
        const std::string base = "ports/";
        port_analytic = h5io::read_matrix_int(design.id(), base + port_analytic_varname).row(0);
        port_id = h5io::read_matrix_int(design.id(), base + port_id_varname).row(0);
        port_coords = h5io::read_matrix(design.id(), base + port_coords_varname);
        port_other_params = h5io::read_matrix(design.id(), base + port_other_params_varname);
        port_x_axis = h5io::read_matrix(design.id(), base + port_x_axis_varname);
        port_incidence_on = h5io::read_matrix_int(design.id(), base + port_incidence_on_varname).row(0);
    }

    n_port = port_id.size();

    #ifdef DEBUG
        std::cout << port_id_varname.data() <<'\n';
        std::cout << port_coords << '\n';
        std::cout << "port_id:" << port_id << '\n';
        std::cout << "port_incidence_on:" << port_incidence_on << '\n';
        std::cout << "port_analytic:" << port_analytic << '\n';
    #endif // DEBUG
}//WaveguideParams::WaveguideParams


void WaveguideParams::analytic_E_field(Eigen::Vector3d pt, int curr_port, double omega, std::complex<double>& value, bool normalize,
    WaveguideModes::FieldComp which_field)
{
    // get parameters
    // (amplitude, width, mode_index, beta, epsilon_wg, epsilon_amb, mu_wg, mu_amb)
    double w = port_other_params(1, curr_port); // width
    double beta = port_other_params(3, curr_port); // beta
    double epsilon_wg = port_other_params(4, curr_port); // epsilon_wg
    double epsilon_amb = port_other_params(5, curr_port); // epsilon_amb
    double mu_wg = port_other_params(6, curr_port); // mu_wg
    double mu_amb = port_other_params(7, curr_port); // mu_amb
    unsigned int mode_index = (unsigned int)(port_other_params(2, curr_port) + 0.5); // mode index
    // other parameters
    double kappa = sqrt(beta*beta - omega*omega*epsilon_amb*mu_amb);
    double k = sqrt(omega*omega*epsilon_wg*mu_wg - beta*beta);
    double mode_energy;
    double mode_flux;
    double omega_sq = omega*omega;
    double beta_sq = beta*beta;

    // Get y coordinate
    Eigen::Vector3d cen = port_coords.col(curr_port);
    Eigen::Vector3d rel_vec = pt - cen;
    Eigen::Vector3d y_axis = port_x_axis.col(curr_port);
    Eigen::Vector3d x_axis = {y_axis(1), - y_axis(0), 0};
    double x = rel_vec.dot(port_x_axis.col(curr_port));


    if(which_field == WaveguideModes::Ex || which_field == WaveguideModes::Ey || which_field == WaveguideModes::Hz||
        which_field == WaveguideModes::Dx || which_field == WaveguideModes::Dy || which_field == WaveguideModes::Bz)
    {
        value = 0.0;
    }
    else
    {

        std::complex<double> Ez, Hx_loc, Hy_loc;
        // even mode or odd mode
        // fundamental mode is 0
        if(mode_index % 2)
        {
            // odd mode
            mode_energy = pow(sin(k*w/2), 2)*(epsilon_amb + beta_sq/(omega_sq*mu_amb))*exp(- kappa*w)/(4*kappa)
                        + (epsilon_wg + beta_sq/(omega_sq*mu_wg))*(w - sin(k*w)/k)/8;
            mode_flux = beta*pow(sin(k*w/2), 2)*exp(-kappa*w)/(2*kappa*omega*mu_amb)
                        + beta*(w - sin(k*w)/k)/(4*omega*mu_wg);


            double C_factor;

            if(normalize)
            {
                C_factor = 1.0/sqrt(mode_energy);
            }
            else
            {
                double amplitude = port_other_params(0, curr_port);
                C_factor = amplitude / sqrt(mode_flux);
            }

            // calculate the field
            if(x < - w/2)
            {
                Ez = - C_factor * sin(k*w/2)*exp(kappa * x);
                Hx_loc = std::complex<double>(1j)/(omega*mu_amb)* kappa * C_factor * sin(k*w/2)*exp(kappa * x);
                Hy_loc = - beta/(omega * mu_amb) * Ez;
            }
            else if(x > w/2)
            {
                Ez = C_factor * sin(k*w/2) * exp(- kappa*x);
                Hx_loc = std::complex<double>(1j)/(omega*mu_amb) * kappa * C_factor * sin(k*w/2) * exp(- kappa*x);
                Hy_loc = - beta/(omega * mu_amb) * Ez;
            }
            else
            {
                Ez = C_factor * sin(k*x);
                Hx_loc = - std::complex<double>(1j)/(omega*mu_wg) * C_factor * k * cos(k*x);
                Hy_loc = - beta/(omega * mu_wg) * Ez;
            }

        }
        else
        {
            // even mode
            mode_energy = pow(cos(k*w/2), 2)*(epsilon_amb + beta_sq/(omega_sq*mu_amb))*exp(-kappa*w)/(4*kappa)
                            + (epsilon_wg + beta_sq/(omega_sq * mu_wg))*(w + sin(k*w)/k)/8;
            mode_flux = beta*pow(cos(k*w/2),2) * exp(- kappa * w)/ (2 * kappa * omega * mu_amb)
                            + beta * (w + sin(k*w)/k) / (4*omega*mu_wg);
            /* debug
            std::cout << "x = " << x <<'\n';
            std::cout << "mode_energy:"<<mode_energy<<'\n';
            std::cout << "mode_flux:"<<mode_flux<<'\n';
            std::cout << "beta:"<<beta<<'\n';
            std::cout << "omega:" << omega<<'\n';
            std::cout << omega*omega*epsilon_wg*mu_wg - beta*beta << '\n';
            std::cout << beta*beta - omega*omega*epsilon_amb*mu_amb << '\n';
        //*/
        // end debug
            double C_factor;
            if(normalize)
            {
                C_factor = 1/sqrt(mode_energy);
            }
            else
            {
                double amplitude = port_other_params(0, curr_port);
                C_factor = amplitude / sqrt(mode_flux);
            }

            // calculate the field
            if(x < - w/2)
            {
                Ez = C_factor * cos(k*w/2)*exp(kappa * x);
                Hx_loc = - std::complex<double>(1j)/(omega*mu_amb)* kappa * C_factor * cos(k*w/2)*exp(kappa * x);
                Hy_loc = - beta/(omega * mu_amb) * Ez;
            }
            else if(x > w/2)
            {
                Ez = C_factor * cos(k*w/2) * exp(- kappa*x);
                Hx_loc = std::complex<double>(1j)/(omega*mu_amb)* kappa * C_factor * cos(k*w/2)*exp( - kappa * x);
                Hy_loc = - beta/(omega * mu_amb) * Ez;
            }
            else
            {
                Ez = C_factor * cos(k*x);
                Hx_loc = std::complex<double>(1j)/(omega*mu_wg) * C_factor * k* sin(k*x);
                Hy_loc = - beta/(omega*mu_wg) * Ez;
            }
        }

        double epsilon_loc, mu_loc;
        if(x < -w/2 || x > w/2)
        {
            epsilon_loc = epsilon_amb;
            mu_loc = mu_amb;
        }
        else
        {
            epsilon_loc = epsilon_wg;
            mu_loc = mu_wg;
        }

        value = 1.0;
        switch(which_field)
        {
        case WaveguideModes::Dz:
            value *= epsilon_loc;
        case WaveguideModes::Ez:
            value *= Ez;
            break;

        case WaveguideModes::Bx:
            value *= mu_loc;
        case WaveguideModes::Hx:
            value *= Hx_loc * x_axis(0) + Hy_loc * y_axis(0);
            break;

        case WaveguideModes::By:
            value *= mu_loc;
        case WaveguideModes::Hy:
            value *= Hx_loc * x_axis(1) + Hy_loc * y_axis(1);
            break;
        }
    }
}
