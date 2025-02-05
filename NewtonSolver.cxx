#include "NewtonSolver.hxx"
#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <cmath>
#include <iostream>
#include <functional>
#include <chrono>

NewtonSolver::NewtonSolver(const Mesh& mesh) 
                         : integrator(mesh) {}

// Matrix NewtonSolver::solve(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, double tolerance, int maxIteration, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax, bool test)
// {
//     Vector F(3*lMax);
//     Matrix J(3*lMax,3*lMax);
//     Vector G(3*lMax);
//     double norm = 1;
//     int count = 0;
//     F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
//     while (norm > tolerance)
//     {
//         count+=1;
//         // std::cout << count << "\n";
//         J = createJ(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
//         G = J.CalculateInverse()*F;
//         for (int i=0; i<3; i++)
//         {
//             for (int l=0; l<lMax; l++)
//             {
//                 alpha(i,l)-=G[i+l*3];
//             }
//         }
//         F = createF(alpha, nu, rho, u, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
//         norm = F.CalculateNorm(1);
//         if (test==true)
//         {
//             // alpha.Print();
//             std::cout << norm << "\n";
//         }
//         if (norm != norm)
//         {
//             std::cout << "norm is NaN" << "\n";
//         }
//         if (count > maxIteration)
//         {
//             std::cout << "alpha did not converge after " << maxIteration << " iterations. Norm is "<< norm <<"\n";
//             norm = 0;
//         }
//     }

//     return alpha;
// }

Matrix NewtonSolver::solve(Matrix alpha, double nu, Vector rho, Vector ux, Vector uy, Vector uz, Vector rt, double dx, Vector roots, Vector weights, double tolerance, int maxIteration, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax, bool test)
{
    Vector F(5*lMax);
    Matrix J(5*lMax,5*lMax);
    Vector G(5*lMax);
    double norm = 1;
    int count = 0;
    F = createF(alpha, nu, rho, ux, uy, uz, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
    if (test)
    {
        F.Print();
        alpha.Print();
    }
    while (norm > tolerance)
    {
        count+=1;
        // std::cout << count << "\n";
        J = createJ(alpha, nu, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        G = J.CalculateInverse()*F;
        for (int i=0; i<5; i++)
        {
            for (int l=0; l<lMax; l++)
            {
                alpha(i,l)-=G[i+l*5];
            }
        }
        F = createF(alpha, nu, rho, ux, uy, uz, rt, dx, roots, weights, basisFunction, quadratureOrder, lMax);
        norm = F.CalculateNorm(1);
        if (test)
        {
            F.Print();
            alpha.Print();
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

    return alpha;
}

// Vector NewtonSolver::createF(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
// {
//     Vector F(3*lMax);
//     Vector MomentF(3*lMax);
//     Vector tilde(lMax);
    
//     for (int l=0; l<lMax; l++)
//     {
//         for (int m=0; m<3; m++)
//         {
//             switch (m)
//             {
//             case 0:
//                 tilde = rho;
//                 break;
//             case 1:
//                 tilde = u;
//                 break;
//             case 2:
//                 tilde = rt;
//                 break;
//             default:
//                 std::cout << "******************issue with F******************" << "\n";
//                 break;
//             }
//             for (int i=0; i<quadratureOrder; i++)
//             {
//                 double integral = integrator.integrate(alpha, basisFunction, m, roots[i], lMax);

//                 double moment = SpecialFunctions::computeMoment(tilde, basisFunction, lMax, roots[i]);

//                 F[m+l*3] += weights[i]*basisFunction(l,roots[i])*nu*(integral-moment)*dx/2.0;
//             }
//         }
//     }

//     return F;
// }

// Matrix NewtonSolver::createJ(Matrix alpha, double nu, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
// {
//     Matrix J(3*lMax,3*lMax);

//     for (int l=0; l<lMax; l++)
//     {
//         for (int m=0; m<3; m++)
//         {
//             for (int p=0; p<lMax; p++)
//             {
//                 for (int n=0; n<3; n++)
//                 {
//                     for (int i=0; i<quadratureOrder; i++)
//                     {
//                         double integral = integrator.integrate(alpha, basisFunction, m+n, roots[i], lMax);

//                         if (n==0)
//                         {
//                             J(m+l*3,n+p*3) += weights[i]*basisFunction(p,roots[i])*basisFunction(l,roots[i])*(nu*integral)*dx/2.0;
//                         }
//                         else
//                         {
//                             J(m+l*3,n+p*3) -= weights[i]*basisFunction(p,roots[i])*basisFunction(l,roots[i])*(nu*integral)*dx/2.0;
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return J;
// }

Vector NewtonSolver::createF(Matrix alpha, double nu, Vector rho, Vector ux, Vector uy, Vector uz, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
{
    Vector F(5*lMax);
    Vector tilde(lMax);

    for (int i=0; i<quadratureOrder; i++)
    {
        Vector integrals = integrator.integrateAlphaMoments(alpha,basisFunction,roots[i],lMax);

        for (int m=0; m<5; m++)
        {
            switch (m)
            {
            case 0:
                tilde = rho;
                break;
            case 1:
                tilde = ux;
                break;
            case 2:
                tilde = uy;
                break;
            case 3:
                tilde = uz;
                break;
            case 4:
                tilde = rt;
                break;
            default:
                std::cout << "******************issue with F******************" << "\n";
                break;
            }
            
            double moment = SpecialFunctions::computeMoment(tilde, basisFunction, lMax, roots[i]);

            double integral = integrals[m];

            for (int l=0; l<lMax; l++)
            {
                F[m+l*5] += weights[i]*basisFunction(l,roots[i])*nu*(integral-moment)*dx/2.0;
            }
        }
    }

    return F;
}

Matrix NewtonSolver::createJ(Matrix alpha, double nu, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax)
{
    Matrix J(5*lMax,5*lMax);
    Matrix Vindex(5,5);

    Vindex(0,0) = 0; Vindex(0,1) = 1; Vindex(0,2) = 2; Vindex(0,3) = 3; Vindex(0,4) = 4;
    Vindex(1,0) = 1; Vindex(1,1) = 5; Vindex(1,2) = 6; Vindex(1,3) = 7; Vindex(1,4) = 8;
    Vindex(2,0) = 2; Vindex(2,1) = 6; Vindex(2,2) = 9; Vindex(2,3) = 10; Vindex(2,4) = 11;
    Vindex(3,0) = 3; Vindex(3,1) = 7; Vindex(3,2) = 10; Vindex(3,3) = 12; Vindex(3,4) = 13;
    Vindex(4,0) = 4; Vindex(4,1) = 8; Vindex(4,2) = 11; Vindex(4,3) = 13; Vindex(4,4) = 14;

    double jacob = nu*dx/2.0;
    for (int i=0; i<quadratureOrder; i++)
    {
        Vector integrals = integrator.getVelocityIntegrals(alpha, basisFunction, roots[i], lMax);
        // double prefactor = jacob*weights[i];
        for (int l=0; l<lMax; l++)
        {
            // prefactor*=basisFunction(l,roots[i]);
            for (int p=0; p<lMax; p++)
            {
                // prefactor*=basisFunction(p,roots[i]);
                for (int m=0; m<5; m++)
                {
                    for (int n=0; n<5; n++)
                    {
                        double integral = integrals[static_cast<int>(Vindex(m,n))];
                        if (m != 4)
                        {
                            // J(m+l*5,n+p*5) += prefactor*integral;
                            J(m+l*5,n+p*5) += jacob*weights[i]*basisFunction(l,roots[i])*basisFunction(p,roots[i])*integral;
                        }
                        else
                        {
                            // J(m+l*5,n+p*5) -= prefactor*integral;
                            J(m+l*5,n+p*5) -= jacob*weights[i]*basisFunction(l,roots[i])*basisFunction(p,roots[i])*integral;
                        }
                    }
                }
            }
        }
    }

    return J;
}