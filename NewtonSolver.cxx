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
    Vector F_new(3*lMax);
    Matrix J(3*lMax,3*lMax);
    Vector G(3*lMax);
    double norm = 1;
    int count = 0;
    F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
    if (test==true)
    {
        alpha.Print();
        // rho.Print();
        // u.Print();
        // rt.Print();
        // F.Print();
        std::cout << F.CalculateNorm(1) << "\n";
    }
    while (norm > tolerance)
    {
        count+=1;
        if (test==true)
        {
            std::cout << count << "\n";
        }
        J = createJ(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        if (test == true)
        {
            G = J.Transpose()*F;

            for (int i=0; i<3; i++)
            {
                for (int l=0; l<lMax; l++)
                {
                    alpha(i,l)-=0.001*G[i+l*3];
                }
            }
        }
        else
        {
            G = J.CalculateInverse()*F;


            double alpha_step = 1.0;
            double beta = 0.5;
            double c = 1e-4;

            while (true)
            {
                Matrix alpha_new = alpha;

                for (int i=0; i<3; i++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        alpha_new(i,l)-=alpha_step*G[i+l*3];
                    }
                }
                F_new = createF(alpha_new, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);

                // std::cout << F_new.CalculateNorm(1) << "\n";
                // std::cout << c * alpha_step * (G*F) << "\n";
                if (test==true)
                {
                    std::cout << F_new.CalculateNorm(1) << "\n";
                }

                if (F_new.CalculateNorm(1) <= F.CalculateNorm(1) + c * alpha_step * (G*F))
                {
                    alpha = alpha_new;
                    break;
                }
                else
                {
                    alpha_step*=beta;
                }

                if (alpha_step < 1e-10)
                {
                    std::cout << "Line search has failed" << "\n";
                    break;
                }
            }
        }

        F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        norm = F.CalculateNorm(1);
        if (test==true)
        {
            alpha.Print();
            std::cout << norm << "\n";
        }
        if (norm != norm)
        {
            std::cout << "norm is NaN" << "\n";
        }
        if (count > maxIteration)
        {
            norm = 0;
            std::cout << "alpha did not converge after " << maxIteration << " iterations\n"; 
        }
    }

    // std::cout << count << "\n";
    // F.Print();

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