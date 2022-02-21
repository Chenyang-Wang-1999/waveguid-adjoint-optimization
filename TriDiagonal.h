#include <vector>
#include <iostream>
#ifndef TRIDIAGONAL_H_INCLUDED
#define TRIDIAGONAL_H_INCLUDED

template<typename T>
void Solve_tridiagonal(const std::vector<double> & a, const std::vector<double> &b, const std::vector<double> & c,
                        const std::vector<T> & r, std::vector<T> &sol)
{
    // size of the system
    unsigned int sys_size = r.size();
    if(a.size() != sys_size - 1 || b.size()!=sys_size || c.size()!=sys_size-1)
    {
        std::cout << "Matrix size does not match!\n";
    }

    // Middle variables
    std::vector<double> beta(sys_size);
    std::vector<T> rho(sys_size);
    double l;

    beta[0] = b[0];
    rho[0] = r[0];

    // upper trianglize
    for(unsigned int i=1; i < sys_size ; i++)
    {
        // std::cout << i <<'\n';
        l = a[i-1]/beta[i-1];
        beta[i] = b[i] - l*c[i-1];
        rho[i] = r[i] - l*rho[i-1];
    }
    // std::cout << "pass 1 \n";
    // solve equation
    sol[sys_size - 1] = rho[sys_size - 1]/beta[sys_size-1];
    for(int i = sys_size-2; i >=0 ; i--)
    {
        sol[i] = (rho[i]-c[i]*sol[i+1])/beta[i];
    }
    // std::cout << "pass 2 \n";
}// Solve_tridiagonal

#endif // TRIDIAGONAL_H_INCLUDED

