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

    double integrate(Vector alpha, std::function<double(int,double)> basisFunction, int power, double x);

private:
    const Mesh& mesh;

    double testMaxwellian(Vector alpha, std::function<double(int,double)> basisFunction, int power, double x, int k);


};

#endif // NEWTONCOTES_H