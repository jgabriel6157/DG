#ifndef NEWTONSOLVERHEADERDEF
#define NEWTONSOLVERHEADERDEF

#include <functional>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "Mesh.hxx"
#include "NewtonCotes.hxx"
#include "SpecialFunctions.hxx"

class NewtonSolver
{
    public:
    NewtonSolver(const Mesh& mesh);

    Matrix solve(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, double tolerance, int maxIteration, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax);

private:
    NewtonCotes integrator;

    Vector createF(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax);

    Matrix createJ(Matrix alpha, double nu, Vector rho, Vector u, Vector rt, double dx, Vector roots, Vector weights, std::function<double(int,double)> basisFunction, int quadratureOrder, int lMax);
};










#endif