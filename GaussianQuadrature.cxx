#include <iostream>
#include <cmath>
#include <cassert>
#include "GaussianQuadrature.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <functional>

//Perform Gaussian quadrature integration for two specified functions of order n_func1, n_func2 to a certain quadrature order
double GaussianQuadrature::integrate(std::function<double(int, double)> func1, int n_func1, std::function<double(int, double)> func2, int n_func2, int quadratureOrder,
                        Vector roots, Vector weights)
{
    double integral = 0;

    for (int i=0; i<quadratureOrder; i++)
    {
        integral += weights[i]*func1(n_func1,roots[i])*func2(n_func2,roots[i]);
    }

    return integral;
}

double GaussianQuadrature::integrate(std::function<double(int, double)> func1, int n_func1, Matrix alpha, double vx, int lMax, int quadratureOrder, Vector roots, Vector weights)
{
    double integral = 0;

    Vector alphaSum(lMax);

    for (int l=0; l<lMax; l++)
    {
        alphaSum[l] = alpha(0,l)-vx*alpha(1,l)-pow(vx,2)*alpha(2,l);
    }

    for (int i=0; i<quadratureOrder; i++)
    {
        double exponent = 0;
        for (int l=0; l<lMax; l++)
        {
            exponent += func1(l,roots[i])*alphaSum[l];
        }
        integral += weights[i]*func1(n_func1,roots[i])*exp(exponent);
    }

    return integral;
}

//Calculate weights for Gaussian quadrature
Vector GaussianQuadrature::calculateWeights(int quadratureOrder, Vector roots)
{
    Vector weights(quadratureOrder);

    double root;
    double w;
    for (int i=0; i<quadratureOrder; i++)
    {
        
        root = roots[i];
        w = 2.0/((1.0-pow(root,2.0))*pow(SpecialFunctions::legendreDerivative(quadratureOrder,root),2.0));
        
        weights[i] = w;
    }

    return weights;
}