#ifndef SPECIALFUNCTIONSHEADERDEF
#define SPECIALFUNCTIONSHEADERDEF

#include "Vector.hxx"
#include "Matrix.hxx"
#include <functional>

class SpecialFunctions
{
private:
    static double newtonRaphson(int n, double x0);
public:
    SpecialFunctions() {};
    static double legendre(int n, double x);
    static double legendreOrthonormal(int n, double x);
    static double legendreDerivative(int n, double x);
    static double legendreOrthonormalDerivative(int n, double x);
    static double quadratic(int n, double x);
    static double quadraticDerivative(int n, double x);
    static double linear(int n, double x);
    static double linearDerivative(int n, double x);

    static Vector legendreRoots(int n);

    static double topHat(double x);
    static double gaussianPulse(double x);
    static double twinGaussianPulse(double x);
    static double constantFunction(double x);
    static double inelasticICx(double x);
    static double inelasticICvx(double x);

    static double sign(double x);
    static double min(double a, double b);
    static double minmod(double a, double b, double c);

    //compute the value of your moment at spatial point x
    static double computeMoment(Vector moment, std::function<double(int,double)> basisFunction, int lMax, double x);

    //compute the value of f_eq (Maxwellian) from the moments and vx
    static double computeMaxwellian(double rho, double u, double rt, double vx);

    //compute the reconstructed f(x,t)
    static double getF(Matrix uPre, int lMax, std::function<double(int,double)> basisFunction, int j, double x);
};



#endif