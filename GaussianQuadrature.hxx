#ifndef GAUSSIANQUADRATUREHEADERDEF
#define GAUSSIANQUADRATUREHEADERDEF

#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <functional>

class GaussianQuadrature
{
public:
    static double integrate(std::function<double(int, double)> func1, int n_func1, std::function<double(int, double)> func2, int n_func2, int quadratureOrder,
                            Vector roots, Vector weights);
    static double integrate(std::function<double(int, double)> func1, int n_func1, Matrix alpha, double vx, int lMax, int quadratureOrder, Vector roots, Vector weights);
    static Vector calculateWeights(int quadratureOrder, Vector roots);
};






#endif