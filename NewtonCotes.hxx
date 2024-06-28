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

private:
    const Mesh& mesh;


};

#endif // NEWTONCOTES_H