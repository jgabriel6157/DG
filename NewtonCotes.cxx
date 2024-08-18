#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <cmath>
#include <iostream>

NewtonCotes::NewtonCotes(const Mesh& mesh) : mesh(mesh) {}

Vector NewtonCotes::integrate(Matrix M, int lMax, int power) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        integral[l] = M(l,0)*pow(mesh.getVelocity(0),power) + M(l,nvx-1)*pow(mesh.getVelocity(nvx-1),power);

        for (int k = 1; k < nvx - 1; k += 2) 
        {
            double vx = mesh.getVelocity(k);
            integral[l] += 4.0*M(l,k)*pow(vx,power);
        }
        for (int k = 2; k < nvx - 1; k += 2) 
        {
            double vx = mesh.getVelocity(k);
            integral[l] += 2.0*M(l,k)*pow(vx,power);
        }

        integral[l] *= dvx / 3.0;
    }
    return integral;
}

double NewtonCotes::integrate(Matrix f, int lMax, std::function<double(int,double)> basisFunction, double x)
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double integral;

    integral = fabs(SpecialFunctions::getF(f,lMax,basisFunction,0,x))*log(fabs(SpecialFunctions::getF(f,lMax,basisFunction,0,x)));
    integral += fabs(SpecialFunctions::getF(f,lMax,basisFunction,nvx-1,x))*log(fabs(SpecialFunctions::getF(f,lMax,basisFunction,nvx-1,x)));

    for (int k=1; k<nvx-1; k+=2)
    {
        double val = fabs(SpecialFunctions::getF(f,lMax,basisFunction,k,x));
        integral += 4.0*val*log(val);
    }
    for (int k=2; k<nvx-1; k+=2)
    {
        double val = fabs(SpecialFunctions::getF(f,lMax,basisFunction,k,x));
        integral += 2.0*val*log(val);
    }

    integral *= -dvx/3.0;
    
    return integral;
}

double NewtonCotes::integrate(Matrix alpha, std::function<double(int,double)> basisFunction, int power, double x, int lMax)
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double integral;

    integral = testMaxwellian(alpha,basisFunction,power,x,0,lMax)+testMaxwellian(alpha,basisFunction,power,x,nvx-1,lMax);

    for (int k=1; k<nvx-1; k+=2)
    {
        integral += 4.0*testMaxwellian(alpha,basisFunction,power,x,k,lMax);
    }
    for (int k=2; k<nvx-1; k+=2)
    {
        integral += 2.0*testMaxwellian(alpha,basisFunction,power,x,k,lMax);
    }

    integral *= dvx/3.0;

    return integral;

}

double NewtonCotes::testMaxwellian(Matrix alpha, std::function<double(int,double)> basisFunction, int power, double x, int k, int lMax)
{
    double vx = mesh.getVelocity(k);
    double exponent = 0;

    for (int l=0; l<lMax; l++)
    {
        exponent += basisFunction(l,x)*(alpha(0,l)-vx*alpha(1,l)-pow(vx,2)*alpha(2,l));
    }

    return pow(vx,power)*exp(exponent);
}