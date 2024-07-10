#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include <cmath>
#include <iostream>

NewtonCotes::NewtonCotes(const Mesh& mesh) : mesh(mesh) {}

Vector NewtonCotes::integrate(Matrix M, int lMax, int power) 
{
    int nvx = mesh.getNVX();
    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        integral[l] = M(l,0)*pow(mesh.getVelocity(0),power) + M(l,nvx-1)*pow(mesh.getVelocity(nvx-1),power);

        for (int i = 1; i < nvx - 1; i += 2) 
        {
            double vx = mesh.getVelocity(i);
            integral[l] += 4*M(l,i)*pow(vx,power);
        }
        for (int i = 2; i < nvx - 1; i += 2) 
        {
            double vx = mesh.getVelocity(i);
            integral[l] += 2 * M(l,i)*pow(vx,power);
        }

        integral[l] *= mesh.getDVX() / 3.0;
    }
    return integral;
}