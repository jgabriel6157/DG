#include "NewtonSolver.hxx"
#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <cmath>
#include <iostream>
#include <functional>

NewtonSolver::NewtonSolver(const Mesh& mesh) 
                         : integrator(mesh) {}

Vector NewtonSolver::solve(Vector alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, double tolerance, int maxIteration, std::function<double(int,double)> basisFunction, int quadratureOrder)
{
    Vector F(3);
    Matrix J(3,3);
    Vector G(3);
    double norm = 1;
    int count = 0;
    F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder);
    while (norm > tolerance)
    {
        count+=1;
        J = createJ(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder);
        G = J.CalculateInverse()*F;
        for (int i=0; i<3; i++)
        {
            alpha[i]-=G[i];
        }
        F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder);
        norm = F.CalculateNorm(1);
        if (count > maxIteration)
        {
            norm = 0;
        }
    }

    return alpha;
}

Vector NewtonSolver::createF(Vector alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder)
{
    Vector F(3);
    Vector tilde(1);
    
    for (int m=0; m<3; m++)
    {
        for (int i=0; i<quadratureOrder; i++)
        {

            double integral = integrator.integrate(alpha, basisFunction, m, roots[i]);

            switch (m)
            {
            case 0:
                tilde = rho;
                break;
            case 1:
                tilde = u;
                break;
            case 2:
                tilde = rt;
                break;
            default:
                std::cout << "******************issue with F******************" << "\n";
                break;
            }

            double moment = SpecialFunctions::computeMoment(tilde, basisFunction, 1, roots[i]);

            F[m] += weights[i]*basisFunction(0,roots[i])*(nu*integral-moment)*dx/2.0;
        }
    }

    return F;
}

Matrix NewtonSolver::createJ(Vector alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder)
{
    Matrix J(3,3);

    for (int m=0; m<3; m++)
    {
        for (int n=0; n<3; n++)
        {
            for (int i=0; i<quadratureOrder; i++)
            {
                double integral = integrator.integrate(alpha, basisFunction, m+n, roots[i]);

                if (n==0)
                {
                    J(m,n) += weights[i]*basisFunction(0,roots[i])*basisFunction(0,roots[i])*(nu*integral)*dx/2.0;
                }
                else
                {
                    J(m,n) -= weights[i]*basisFunction(0,roots[i])*basisFunction(0,roots[i])*(nu*integral)*dx/2.0;
                }
            }
        }
    }

    return J;
}