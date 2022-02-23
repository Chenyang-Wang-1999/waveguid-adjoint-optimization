#include "Bezier.h"
#include <vector>
#include "TriDiagonal.h"
#include <Eigen/Core>
#include <iostream>

/*
#define DEBUG

//*/
// #define WEXTRAP

Eigen::MatrixXd Bezier_tangent_vector(const Eigen::MatrixXd &P, unsigned int start_point_mode, unsigned int end_point_mode,
        Eigen::VectorXd vec_start, Eigen::VectorXd vec_end, bool autoscale)
{
    unsigned int N_points = P.rows();   // number of points
    unsigned int dim = P.cols();        // dimension of the points

    // set indices
    std::vector<double> a(N_points - 1), b(N_points), c(N_points - 1);
    std::vector<Eigen::VectorXd> rhs(N_points), sol(N_points);

    // initialize
    rhs[0] = P.row(0);
    sol[0] = P.row(0);

    // assemble the bulk
    for(unsigned int i=1; i<N_points - 1; i++)
    {
        rhs[i] = 3*P.row(i + 1) - 3*P.row(i - 1);
        a[i-1] = 1.0;
        b[i] = 4.0;
        c[i] = 1.0;
    }
    sol[N_points - 1] = P.row(0);

    // assemble the boundary
    // start point
    if(start_point_mode == 0)   // given_vector
    {
        b[0] = 1;
        c[0] = 0;
        if(autoscale)
        {
            double vec_norm = vec_start.norm();
            double tan_vec_norm = (P.row(1) - P.row(0)).norm();
            rhs[0] = vec_start*tan_vec_norm/vec_norm;
        }
        else
        {
            rhs[0] = vec_start;
        }
    }
    else if(start_point_mode == 1) // zero second diff
    {
        b[0] = 2;
        c[0] = 1;
        rhs[0] = 3*(P.row(1) - P.row(0));
    }
    else
    {
        std::cout << "Not supported";
    }

    // end points
    if(end_point_mode == 0) //given_vector
    {
        a[N_points - 2] = 0;
        b[N_points - 1] = 1;
        if(autoscale)
        {
            double vec_norm = vec_end.norm();
            double tan_vec_norm = (P.row(N_points - 1) - P.row(N_points - 2)).norm();
            rhs[N_points - 1] = vec_end*tan_vec_norm / vec_norm;
        }
        else
        {
            rhs[N_points - 1] = vec_end;
        }
    }
    else if(end_point_mode == 1)    // zero second diff
    {
        a[N_points - 2] = 1;
        b[N_points - 1] = 2;
        rhs[N_points - 1] = 3*(P.row(N_points - 1) - P.row(N_points - 2));
    }

    // std::cout<<"previous pass\n";
    // solve
    Solve_tridiagonal(a, b, c, rhs, sol);
    // std::cout << "tridiag pass\n";
    // assemble Eigen::Vector
    // Eigen::Map<Eigen::MatrixXd> V(sol.data()->data(), N_points, dim);

    //*
    Eigen::MatrixXd V(N_points,dim);
    for(unsigned int i=0; i< N_points; i++)
    {
        V.row(i) = sol[i];
    }
    //*/
    return V;
}

//*
Eigen::MatrixXd Bezier_interp(const Eigen::MatrixXd &P, const Eigen::MatrixXd &V, Eigen::VectorXd &t,  unsigned int diff_order)
{
    unsigned int dim = P.cols();
    unsigned int N_interp_points = t.size();
    unsigned int N_ctrl_points = P.rows();
    // initialize interp points
    Eigen::MatrixXd B(N_interp_points, dim);

    for(unsigned int i = 0; i< N_interp_points; i++)
    {
        if(t(i) < 0)
        {
            #ifdef WEXTRAP
                std::cout << "Warning: Extrapolating for t = " << t(i) << ", snap to 0 \n";
            #endif // WEXTRAP
            t(i) = 0.0;
        }
        else if (t(i) > N_ctrl_points - 1)
        {
            #ifdef WEXTRAP
                std::cout << "Warning: Extrapolating for t = " << t(i) << ", snap to the boundary \n";
            #endif // WEXTRAP
            t(i) = N_ctrl_points - 1;
        }
        int id = (int)t(i);    // id ~ id + 1
        double loc_t = t(i) - id;    // get local t
        double loc_inv_t = 1 - loc_t;
        //std::cout << "id: "<<id <<", local t: "<<loc_t << ", local inv:" << loc_inv_t<<'\n';
        if(-1e-6 < loc_t && loc_t< 1e-6)
        {
            // std::cout << "Is control point\n";
            id = (int)(t(i) + 0.5);
            if(diff_order == 0)
            {
                B.row(i) = P.row(id);
            }
            else if(diff_order == 1)
            {
                B.row(i) = V.row(id);
            }
            else if(diff_order == 2)
            {
                if(id == P.rows() - 1)// right boundary
                {
                    B.row(i) = 6*P.row(id - 1) - 6*P.row(id) + 2* V.row(id-1) + 4*V.row(id);
                }
                else
                {
                    B.row(i) = 6*P.row(id + 1) - 6*P.row(id) - 4*V.row(id) - 2*V.row(id+1);
                }
            }

        }
        else
        {
            Eigen::VectorXd Cp(dim);
            Cp = P.row(id) + V.row(id)/3;
            Eigen::VectorXd Cn(dim);
            Cn = P.row(id + 1) - V.row(id+1)/3;
            //std::cout << P.row(id)<<'\n';
            //std::cout << Cp.transpose() <<'\n';
            //std::cout << Cn.transpose()<<'\n';
            //std::cout << P.row(id+1)<<'\n';
            if(diff_order == 0)
            {
                B.row(i) =  P.row(id) * loc_inv_t*loc_inv_t*loc_inv_t +
                        3 * Cp.transpose() * loc_inv_t * loc_inv_t * loc_t +
                        3 * Cn.transpose() * loc_inv_t * loc_t * loc_t +
                        P.row(id + 1) * loc_t * loc_t * loc_t;

            }
            else if(diff_order == 1)
            {
                B.row(i) = V.row(id) * loc_inv_t * loc_inv_t
                            - 2*(3*P.row(id) - 3*P.row(id+1) + V.row(id) + V.row(id+1))*loc_t*loc_inv_t
                            +V.row(id+1)*loc_t * loc_t;
            }
            else if (diff_order == 2)
            {
                B.row(i) = 2*(6*Cp.transpose() - 6*Cn.transpose() + V.row(id) + V.row(id+1))*loc_t - 6*Cp.transpose() + 6*Cn.transpose() - 2*V.row(id);
            }
        }

    }
    return B;
}
//*/

Eigen::MatrixXd Bezier_V_sensitivity(unsigned int id_P, unsigned int id_dim,
        const Eigen::MatrixXd &P, unsigned int start_point_mode, unsigned int end_point_mode,
        Eigen::VectorXd vec_start, Eigen::VectorXd vec_end, bool autoscale)
{
    unsigned int N_points = P.rows();
    unsigned int dim = P.cols();

    // Point sensitivity
    Eigen::MatrixXd P_sens(N_points, dim);
    P_sens.setZero();
    P_sens(id_P, id_dim) = 1;

    // start and end vectors
    if(start_point_mode == 0)
    {
        if(id_P == 0 || id_P == 1)
        {
            if(autoscale)
            {
                double vec_norm = vec_start.norm();
                double tan_vec_norm = (P.row(1) - P.row(0)).norm();
                vec_start *= ( P(id_P,id_dim) - P((id_P+1)%2, id_dim) )/vec_norm/tan_vec_norm;
            }
        }
        else
        {
            vec_start.setZero();
        }
    }

    if(end_point_mode == 0)
    {
        if(id_P == N_points - 2 || id_P == N_points - 1)
        {
            if(autoscale)
            {
                double vec_norm = vec_end.norm();
                double tan_vec_norm = (P.row(N_points - 1) - P.row(N_points - 2)).norm();
                vec_end *= ( P(id_P,id_dim) - P((id_P+3-N_points)%2 + N_points - 2, id_dim) )/vec_norm/tan_vec_norm;
            }
        }
        else
        {
            vec_end.setZero();
        }
    }

    return Bezier_tangent_vector(P_sens, start_point_mode, end_point_mode, vec_start, vec_end, false);
}

double Nearest_t(const Eigen::MatrixXd &P_ctrl, const Eigen::MatrixXd &V, Eigen::VectorXd pt, double AbsTol)
{
    Eigen::VectorXd t_interp = Eigen::VectorXd::LinSpaced(100, 0, P_ctrl.rows() - 1);
    Eigen::MatrixXd P = Bezier_interp(P_ctrl, V, t_interp, 0);

    unsigned int N_ctrl = P.rows();
    unsigned int dim = P.cols();

    double nearest_t;

    // 1. Get the approximate value by line segments
    Eigen::VectorXd distance(N_ctrl);

    for(unsigned int curr_row = 0; curr_row < N_ctrl; curr_row ++)
    {
        distance(curr_row) = (P.row(curr_row) - pt.transpose()).norm();
    }

    Eigen::Index min_id;
    double min_dist = distance.minCoeff(& min_id);

    nearest_t = t_interp(min_id);

    // 2. Newton to get
    Eigen::VectorXd T_test(1);
    T_test(0) = nearest_t;
    Eigen::MatrixXd P_t;
    P_t = Bezier_interp(P_ctrl,V,T_test,0);
    double curr_error = (pt.transpose() - P_t.row(0)).norm();
    for(int j = 0; j<1000; j++)
    {
        Eigen::MatrixXd D1_t, D2_t;
        D1_t = Bezier_interp(P_ctrl, V, T_test, 1);
        D2_t = Bezier_interp(P_ctrl, V, T_test, 2);
        double First_diff, Second_diff;
        First_diff = (D1_t.row(0)).dot(P_t.row(0).transpose() - pt);
        Second_diff = (D2_t.row(0)).dot(P_t.row(0).transpose() - pt) + D1_t.row(0).dot(D1_t.row(0));
        T_test(0) += - First_diff/Second_diff;
        P_t = Bezier_interp(P_ctrl,V,T_test,0);

        // error estimation
        curr_error = First_diff>0?First_diff:-First_diff;
#ifdef DEBUG
        std::cout << "Current error:"<< curr_error <<'\n';
#endif // DEBUG
        if(curr_error < AbsTol)
        {
            break;
        }
    }

    return T_test(0);
}
