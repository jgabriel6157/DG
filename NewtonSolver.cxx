#include "NewtonSolver.hxx"
#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <cmath>
#include <iostream>
#include <functional>

NewtonSolver::NewtonSolver(const Mesh& mesh) 
                         : integrator(mesh) {}

Matrix NewtonSolver::solve(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, double tolerance, int maxIteration, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax, bool test)
{
    Vector F(3*lMax);
    Vector FNorm(3*lMax);
    Matrix J(3*lMax,3*lMax);
    Vector G(3*lMax);
    double norm = 1;
    int count = 0;
    F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
    if (test==true)
    {
        // F.Print();
        // std::cout << rho[0] << "\n";
        // std::cout << u[0] << "\n";
        // std::cout << rt[0] << "\n";
        // alpha.Print();
        // std::cout << F.CalculateNorm(1)/rho[0] << "\n";
    }
    // F.Print();
    // alpha.Print();
    while (norm > tolerance)
    {
        count+=1;
        // std::cout << count << "\n";
        J = createJ(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        G = J.CalculateInverse()*F;
        for (int i=0; i<3; i++)
        {
            for (int l=0; l<lMax; l++)
            {
                alpha(i,l)-=G[i+l*3];
            }
        }
        F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        FNorm = F;
        for (int l=0; l<lMax; l++)
        {
            for (int m=0; m<3; m++)
            {
                double divisor;
                switch (m)
                {
                case 0:
                    divisor = rho[0];
                    break;
                case 1:
                    divisor = fabs(u[0]);
                    break;
                case 2:
                    divisor = rt[0];
                    break;
                default:
                    std::cout << "******************issue with F******************" << "\n";
                    break;
                }
                divisor = rt[0];
                if (divisor > 1)
                {
                    FNorm[m+l*3] = F[m+l*3]/divisor;
                }
            }
        }
        // norm = F.CalculateNorm(1)/rho[0];
        norm = FNorm.CalculateNorm(1);
        if (test==true)
        {
            // alpha.Print();
            std::cout << norm << "\n";
        }
        if (norm != norm)
        {
            std::cout << "norm is NaN" << "\n";
        }
        if (count > maxIteration)
        {
            std::cout << "alpha did not converge after " << maxIteration << " iterations. Norm is "<< norm <<"\n"; 
            norm = 0;
        }
    }

    // std::cout << count << "\n";
    // FNorm.Print();

    // if (test==true)
    // {
    //     alpha.Print();
    // }

    return alpha;
}

Vector NewtonSolver::createF(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
{
    Vector F(3*lMax);
    Vector tilde(lMax);
    
    for (int l=0; l<lMax; l++)
    {
        for (int m=0; m<3; m++)
        {
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
            for (int i=0; i<quadratureOrder; i++)
            {

                double integral = integrator.integrate(alpha, basisFunction, m, roots[i], lMax);

                double moment = SpecialFunctions::computeMoment(tilde, basisFunction, lMax, roots[i]);

                if (m==1)
                {
                    // std::cout << i << "\n";
                    // std::cout << integral << "\n";
                    // std::cout << moment << "\n";
                }

                F[m+l*3] += weights[i]*basisFunction(l,roots[i])*nu*(integral-moment)*dx/2.0;
            }
        }
    }

    return F;
}

Matrix NewtonSolver::createJ(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
{
    Matrix J(3*lMax,3*lMax);

    for (int l=0; l<lMax; l++)
    {
        for (int m=0; m<3; m++)
        {
            for (int p=0; p<lMax; p++)
            {
                for (int n=0; n<3; n++)
                {
                    for (int i=0; i<quadratureOrder; i++)
                    {
                        double integral = integrator.integrate(alpha, basisFunction, m+n, roots[i], lMax);

                        if (n==0)
                        {
                            J(m+l*3,n+p*3) += weights[i]*basisFunction(p,roots[i])*basisFunction(l,roots[i])*(nu*integral)*dx/2.0;
                        }
                        else
                        {
                            J(m+l*3,n+p*3) -= weights[i]*basisFunction(p,roots[i])*basisFunction(l,roots[i])*(nu*integral)*dx/2.0;
                        }
                    }
                }
            }
        }
    }

    return J;
}