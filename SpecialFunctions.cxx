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
    return exp(-5.0*pow(x-0.1*M_PI,2.0))+exp(-5.0*pow(x+0.6*M_PI,2.0));
}

//Constant functions that is = 1 for all x
double SpecialFunctions::constantFunction(double x)
{
    return 1;
}

//Initial condition for ionization and cx test cases
double SpecialFunctions::inelasticICx(double x)
{
    double n0 = 5.0;
    double Lx = 40.0;
    double value;
    if (x>20)
    {
        value = pow(cosh(-(Lx/2.0-(x-20.0)-2.0)/2.0),-2)+1e-6;
    }
    else
    {
        value = pow(cosh((Lx/2.0+(x-20.0)-2.0)/2.0),-2)+1e-6;
    }
    value *= n0;

    return value;
}

double SpecialFunctions::inelasticICvx(double x)
{
    double rho = 1.0;
    double u = 0.0;
    double rt = 2.0;
    double vt2 = rt;
    return rho*exp(-pow(x-u,2)/(2.0*vt2))/sqrt(2.0*M_PI*vt2);
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
    double vt2 = rt;
    return rho*exp(-pow(vx-u,2)/(2.0*vt2))/pow(2.0*M_PI*vt2,0.5);
}

double SpecialFunctions::computeMaxwellian3(double rho, double ux, double uy, double uz, double rt, double vx, double vy, double vz)
{
    double vt2 = rt;
    return rho*exp(-(pow(vx-ux,2)+pow(vy-uy,2)+pow(vz-uz,2))/(2.0*vt2))/pow(2.0*M_PI*vt2,1.5);
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

double SpecialFunctions::computeSigmav(double T, double E)
{
    Matrix a(9,9);
    a(0,0) = -1.829079581680e+01; a(0,1) = 1.640252721210e-01;   a(0,2) = 3.364564509137e-02;
    a(1,0) = 2.169137615703e-01;  a(1,1) = -1.106722014459e-01;  a(1,2) = -1.382158680424e-03;
    a(2,0) = 4.307131243894e-02;  a(2,1) = 8.948693624917e-03;   a(2,2) = -1.209480567154e-02;
    a(3,0) = -5.754895093075e-04; a(3,1) = 6.062141761233e-03;   a(3,2) = 1.075907881928e-03;
    a(4,0) = -1.552077120204e-03; a(4,1) = -1.210431587568e-03;  a(4,2) = 8.297212635856e-04;
    a(5,0) = -1.876800283030e-04; a(5,1) = -4.052878751584e-05;  a(5,2) = -1.907025662962e-04;
    a(6,0) = 1.125490270962e-04;  a(6,1) = 2.875900435985e-05;   a(6,2) = 1.338839628570e-05;
    a(7,0) = -1.238982763007e-05; a(7,1) = -2.616998139678e-06;  a(7,2) = -1.171762874107e-07;
    a(8,0) = 4.163596197181e-07;  a(8,1) = 7.558092849125e-08;   a(8,2) = -1.328404104165e-08;

    a(0,3) = 9.530225559189e-03;  a(0,4) = -8.519413589968e-04; a(0,5) = -1.247583860943e-03;
    a(1,3) = 7.348786286628e-03;  a(1,4) = -6.343059502294e-04; a(1,5) = -1.919569450380e-04;
    a(2,3) = -3.675019470470e-04; a(2,4) = 1.039643390686e-03;  a(2,5) = -1.553840717902e-04;
    a(3,3) = -8.119301728339e-04; a(3,4) = 8.911036876068e-06;  a(3,5) = 3.175388949811e-05;
    a(4,3) = 1.361661816974e-04;  a(4,4) = -1.008928628425e-04; a(4,5) = 1.080693990468e-05;
    a(5,3) = 1.141663041636e-05;  a(5,4) = 1.775681984457e-05;  a(5,5) = -3.149286923815e-06;
    a(6,3) = -4.340802793033e-06; a(6,4) = -7.003521917385e-07; a(6,5) = 2.318308730487e-07;
    a(7,3) = 3.517971869029e-07;  a(7,4) = -4.928692832866e-08; a(7,5) = 1.756388998863e-10;
    a(8,3) = -9.170850253981e-09; a(8,4) = 3.208853883734e-09;  a(8,5) = -3.952740758950e-10;

    a(0,6) = 3.014307545716e-04;  a(0,7) = -2.499323170044e-05; a(0,8) = 6.932627237765e-07;
    a(1,6) = 4.075019351738e-05;  a(1,7) = -2.850044983009e-06; a(1,8) = 6.966822400446e-08;
    a(2,6) = 2.670827249272e-06;  a(2,7) = 7.695300597935e-07;  a(2,8) = -3.783302281524e-08;
    a(3,6) = -4.515123641755e-06; a(3,7) = 2.187439283954e-07;  a(3,8) = -2.911233951880e-09;
    a(4,6) = 5.106059413591e-07;  a(4,7) = -1.299275586093e-07; a(4,8) = 5.117133050290e-09;
    a(5,6) = 3.105491554749e-08;  a(5,7) = 2.274394089017e-08;  a(5,8) = -1.130988250912e-09;
    a(6,6) = -6.030983538280e-09; a(6,7) = -1.755944926274e-09; a(6,8) = 1.005189187279e-10;
    a(7,6) = -1.446756795654e-10; a(7,7) = 7.143183138281e-11;  a(7,8) = -3.989884105603e-12;
    a(8,6) = 2.739558475782e-11;  a(8,7) = -1.693040208927e-12; a(8,8) = 6.388219930167e-14;

    double E_min = 0.1;
    double E_max = 2e4;
    double T_min = 0.1;
    double T_max = 2e4;

    E = std::max(E,E_min);
    E = std::min(E,E_max);
    T = std::max(T,T_min);
    T = std::min(T,T_max);

    double logT = log(T);
    double logE = log(E);

    double result = 0;
    for (int i=0; i<9; i++)
    {
        for (int j=0; j<9; j++)
        {
            result += a(j,i)*pow(logE,i)*pow(logT,j);
        }
    }

    return exp(result)*1e-6;
}