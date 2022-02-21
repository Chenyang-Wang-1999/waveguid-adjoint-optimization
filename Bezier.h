#ifndef BEZIER_INCLUDED
#define BEZIER_INCLUDED

#ifdef EXPORT_PYLIB
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py=pybind11;
#endif // EXPORT_PYLIB

#include <Eigen/Core>

using namespace std;

Eigen::MatrixXd Bezier_tangent_vector(const Eigen::MatrixXd &P, unsigned int start_point_mode, unsigned int end_point_mode,
        Eigen::VectorXd vec_start=Eigen::Vector3d({0,0,1}), Eigen::VectorXd vec_end=Eigen::Vector3d({0,0,1}), bool autoscale=true);

Eigen::MatrixXd Bezier_interp(const Eigen::MatrixXd &P, const Eigen::MatrixXd &V, Eigen::VectorXd &t,  unsigned int diff_order);

double Nearest_t(const Eigen::MatrixXd &P_ctrl, const Eigen::MatrixXd &V, Eigen::VectorXd pt, double AbsTol=1e-2);

#ifdef EXPORT_PYLIB
PYBIND11_MODULE(Bezier, m) {

    m.def("Bezier_tangent_vector", &Bezier_tangent_vector);
    m.def("Bezier_interp", &Bezier_interp);
    m.def("Nearest_t", &Nearest_t);
}
#endif // EXPORT_PYLIB

#endif // BEZIER_INCLUDED
