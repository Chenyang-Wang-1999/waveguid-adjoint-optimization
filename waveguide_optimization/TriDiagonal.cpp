#include "TriDiagonal.h"


Eigen::VectorXd Solve_tridiagonal(TridiagonalMatrix &M, Eigen::VectorXd &r)
{

    // Solve tridiagonal equation MX=b
    Eigen::VectorXd X(M.Size);
    Eigen::VectorXd beta(M.Size), rho(M.Size); // middle variant
    int i;
    double l;
    // Error check
    if(r.size() != M.Size)
    {
        std::cout << "Error! Matrix size dose not match the RHS! \n";
        exit(0);
    }
    beta(0) = M.b(0);
    rho(0) = r(0);
    // upper trianglize
    for(i=1; i < M.Size ; i++)
    {
        l = M.a(i-1)/beta(i-1);
        beta(i) = M.b(i) - l*M.c(i-1);
        rho(i) = r(i) - l*rho(i-1);
    }

    // solve equation
    X(M.Size - 1) = rho(M.Size - 1)/beta(M.Size - 1);
    for(i = M.Size-2; i >=0 ; i--)
    {
        X(i) = (rho(i)-M.c(i)*X(i+1))/beta(i);
    }
    return X;
}// Solve_tridiagonal

