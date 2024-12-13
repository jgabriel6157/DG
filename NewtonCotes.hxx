#ifndef NEWTONCOTES_H
#define NEWTONCOTES_H

#include <functional>
#include "Mesh.hxx"
#include "Vector.hxx"
#include "Matrix.hxx"

class NewtonCotes 
{
public:
    NewtonCotes(const Mesh& mesh);

    Vector integrate(Matrix M, int lMax, int power);

    Vector integrate3f(Matrix M, int lMax);
    Vector integrate3vxf(Matrix M, int lMax);
    Vector integrate3vyf(Matrix M, int lMax);
    Vector integrate3vzf(Matrix M, int lMax);
    Vector integrate3v2f(Matrix M, int lMax);
    Matrix integrateMoments(Matrix M, int lMax);

    double integrate(Matrix f, int lMax, std::function<double(int,double)> basisFunction, double x);

    double integrate(Matrix alpha, std::function<double(int,double)> basisFunction, int power, double x, int lMax);

private:
    const Mesh& mesh;

    double testMaxwellian(Matrix alpha, std::function<double(int,double)> basisFunction, int power, double x, int k, int lMax);

    Vector computeWeights(int nv);


};

#endif // NEWTONCOTES_H