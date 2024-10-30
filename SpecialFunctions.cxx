#include <iostream>
#include <cmath>
#include <cassert>
#include <functional>
#include "SpecialFunctions.hxx"
#include "Vector.hxx"
#include "Matrix.hxx"

//Function to calculate Legendre polynomial of order n at point x
double SpecialFunctions::legendre(int n, double x)
{
    assert (n>=0);
    switch (n)
    {
        case 0:
            return 1;
            break;
        case 1:
            return x;
            break;
        default:
            return ((2.0*n-1.0)*x*legendre(n-1,x)-(n-1)*legendre(n-2,x))/n;
            break;
    }
}

//Function to calculate orthonormal Legendre polynomial of order n at point x
double SpecialFunctions::legendreOrthonormal(int n, double x)
{
    assert (n>=0);
    return sqrt((2.0*n+1.0)/2.0)*legendre(n,x);
}

//Function to calculate derivate of the Legendre polynomial of order n at point x
double SpecialFunctions::legendreDerivative(int n, double x)
{
    assert (n>=0);
    switch (n)
    {
        case 0:
            return 0;
            break;
        case 1:
            return 1;
            break;
        default:
            return ((2.0*n-1.0)*(x*legendreDerivative(n-1,x)+legendre(n-1,x))-(n-1)*legendreDerivative(n-2,x))/n;
    }
}

//Function to calculate derivative of orthonormal Legendre polynomial of order n at point x
double SpecialFunctions::legendreOrthonormalDerivative(int n, double x)
{
    assert(n>=0);
    return sqrt((2.0*n+1.0)/2.0)*legendreDerivative(n,x);
}

//Function to calculate quadratic basis function of 'order' n at point x
double SpecialFunctions::quadratic(int n, double x)
{
    assert(n>=0);
    assert(n<3);
    switch (n)
    {
    case 0:
        return -x*(1.0-x)/2.0;
        break;
    case 1:
        return (1.0-x)*(1.0+x);
        break;
    case 2:
        return x*(1+x)/2;
        break;
    default:
        return 0;
        break;
    }
}

//Function to calculate derivative of the quadratic basis function of 'order' n at point x
double SpecialFunctions::quadraticDerivative(int n, double x)
{
    assert(n>=0);
    assert(n<3);
    switch (n)
    {
    case 0:
        return x-1.0/2.0;
        break;
    case 1:
        return -2.0*x;
        break;
    case 2:
        return x+1.0/2.0;
        break;
    default:
        return 0;
        break;
    }
}

//Function to calculate linear basis function of 'order' n at point x
double SpecialFunctions::linear(int n, double x)
{
    assert(n>=0);
    assert(n<2);
    switch (n)
    {
    case 0:
        return (1.0-x)/2.0;
        break;
    case 1:
        return (1.0+x)/2.0;
        break;
    default:
        return 0;
        break;
    }
}

//Function to calculate derivative of linear basis function of 'order' n at point x
double SpecialFunctions::linearDerivative(int n, double x)
{
    assert(n>=0);
    assert(n<2);
    switch (n)
    {
    case 0:
        return -1.0/2.0;
        break;
    case 1:
        return 1.0/2.0;
        break;
    default:
        return 0;
        break;
    }
}

//Square pulse (top hat) function centered at pi with length 2
double SpecialFunctions::topHat(double x)
{
    // if ((x<M_PI-1.0)||(x>M_PI+1.0))
    if ((x<-1.0)||(x>1.0))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

//Gaussian pulse centered at pi
double SpecialFunctions::gaussianPulse(double x)
{
    // return exp(-1.0*pow(x-1.0*M_PI,2.0));
    return exp(-50.0*pow(x-0.5,2.0));
}

//Gaussian pulse centered at +-0.5 pi
double SpecialFunctions::twinGaussianPulse(double x)
{
    // return exp(-5.0*pow(x-0.5*M_PI,2.0))+exp(-5.0*pow(x+0.5*M_PI,2.0));
    // return 1e18*(exp(-(1.0e-8)*pow(x-10000.0,2.0))+exp(-(1.0e-8)*pow(x+10000.0,2.0)));
    return 1e18*(exp(-(1.0e-8)*pow(x-00000,2.0)));
    // return 1e18*(exp(-(1.0)*pow(x-1.0,2.0))+exp(-(1.0)*pow(x+1.0,2.0)));
}

//Constant functions that is = 1 for all x
double SpecialFunctions::constantFunction(double x)
{
    return 1;
}

//Newton Raphson method to find root of Legendre polynomial of root n with initial guess x0
double SpecialFunctions::newtonRaphson(int n, double x0)
{
    double x_root = x0;

    while (fabs(legendre(n,x_root))>1.0e-10)
    {
        x_root -= legendre(n,x_root)/legendreDerivative(n,x_root);
    }

    return x_root;
}

//Find roots of Legendre polynomial of order n using Newton Raphson method 
Vector SpecialFunctions::legendreRoots(int n)
{
    Vector roots(n);

    for (int i=0; i<n; i++)
    {
        double x0 = (1.0-1.0/(8.0*pow(n,2.0))+1.0/(8.0*pow(n,3.0)))*cos(M_PI*(4.0*(i+1.0)-1.0)/(4.0*n+2.0)); //intial guess
        double x = newtonRaphson(n,x0);
        roots[i] = x;
    }

    return roots;
}

//Return sign of value (0 returns 1)
double SpecialFunctions::sign(double x)
{
    if (x>0)
    {
        return 1.0;
    }
    else if (x<0)
    {
        return -1.0;
    }
    else
    {
        return 1;
    }
}

//Return smaller value between a,b
double SpecialFunctions::min(double a, double b)
{
    if (a<b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

//Return minmod of a,b,c
double SpecialFunctions::minmod(double a, double b, double c)
{
    double s = (sign(a)+sign(b)+sign(c))/3.0;
    if (fabs(s)==1)
    {
        return s*fabs(min(min(a,b),c));
    }
    else
    {
        return 0;
    }
}

double SpecialFunctions::computeMoment(Vector moment, std::function<double(int,double)> basisFunction, int lMax, double x)
{
    double momentValue = 0;

    for (int l=0; l<lMax; l++)
    {
        momentValue += moment[l]*basisFunction(l,x);
    }

    return momentValue;
}

double SpecialFunctions::computeMaxwellian(double rho, double u, double rt, double vx)
{
    // double vt2 = rt*9.58134e7;
    double vt2 = rt;
    return rho*exp(-pow(vx-u,2)/(2.0*vt2))/pow(2.0*M_PI*vt2,0.5);
}

double SpecialFunctions::getF(Matrix uPre, int lMax, std::function<double(int,double)> basisFunction, int j, double x)
{
    double f = 0;

    for (int l=0; l<lMax; l++)
    {
        f += uPre(l,j)*basisFunction(l,x); //Should probably be 2.0*(x-xj)/dx
    }

    return f;
}