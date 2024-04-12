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
    double y = 0;

    for (int i=0; i<quadratureOrder; i++)
    {
        y += weights[i]*func1(n_func1,roots[i])*func2(n_func2,roots[i]);
    }

    return y;
}

//Perform Gaussian quadrature integration for three specified functions of order n_func1, n_func2, n_func3 to a certain quadrature order
double GaussianQuadrature::integrate(std::function<double(int, double)> func1, int n_func1, std::function<double(int, double)> func2, int n_func2, 
                        std::function<double(int, double)> func3, int n_func3, int quadratureOrder, Vector roots, Vector weights)
{
    double y = 0;

    for (int i=0; i<quadratureOrder; i++)
    {
        y += weights[i]*func1(n_func1,roots[i])*func2(n_func2,roots[i])*func3(n_func3,roots[i]);
    }

    return y;
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
        
        // weights.push_back(w);
        weights[i] = w;
    }

    return weights;
}