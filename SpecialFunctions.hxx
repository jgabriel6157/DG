#ifndef SPECIALFUNCTIONSHEADERDEF
#define SPECIALFUNCTIONSHEADERDEF

#include "Vector.hxx"

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
    static double constantFunction(double x);

    static double sign(double x);
    static double min(double a, double b);
    static double minmod(double a, double b, double c);
};



#endif