#include "NewtonCotes.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include <cmath>
#include <iostream>

NewtonCotes::NewtonCotes(const Mesh& mesh) : mesh(mesh) {}

//only integrates for vx
Vector NewtonCotes::integrate(Matrix M, int lMax, int power) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        integral[l] = M(l,0)*pow(mesh.getVelocityX(0),power) + M(l,nvx-1)*pow(mesh.getVelocityX(nvx-1),power);

        for (int k = 1; k < nvx - 1; k += 2) 
        {
            double vx = mesh.getVelocityX(k);
            integral[l] += 4.0*M(l,k)*pow(vx,power);
        }
        for (int k = 2; k < nvx - 1; k += 2) 
        {
            double vx = mesh.getVelocityX(k);
            integral[l] += 2.0*M(l,k)*pow(vx,power);
        }

        integral[l] *= dvx / 3.0;
    }
    return integral;
}

Vector NewtonCotes::integrate3f(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    integral[l] += weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy);
                }
            }
        }
        integral[l] *= dvx*dvy*dvz/27.0;
    }
    return integral;
}

Vector NewtonCotes::integrate3vxf(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            double vx = mesh.getVelocityX(kx);
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    integral[l] += weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy)*vx;
                }
            }
        }
        integral[l] *= dvx*dvy*dvz/27.0;
    }
    return integral;
}

Vector NewtonCotes::integrate3vyf(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                double vy = mesh.getVelocityY(ky);
                for (int kz=0; kz<nvz; kz++)
                {
                    integral[l] += weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy)*vy;
                }
            }
        }
        integral[l] *= dvx*dvy*dvz/27.0;
    }
    return integral;
}

Vector NewtonCotes::integrate3vzf(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    double vz = mesh.getVelocityZ(kz);
                    integral[l] += weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy)*vz;
                }
            }
        }
        integral[l] *= dvx*dvy*dvz/27.0;
    }
    return integral;
}

Vector NewtonCotes::integrate3v2f(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            double vx2 = pow(mesh.getVelocityX(kx),2);
            for (int ky=0; ky<nvy; ky++)
            {
                double vy2 = pow(mesh.getVelocityY(ky),2);
                for (int kz=0; kz<nvz; kz++)
                {
                    double vz2 = pow(mesh.getVelocityZ(kz),2);
                    integral[l] += weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy)*(vx2+vy2+vz2);
                }
            }
        }
        integral[l] *= dvx*dvy*dvz/27.0;
    }
    return integral;
}

Matrix NewtonCotes::integrateMoments(Matrix M, int lMax) 
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Matrix integral(lMax, 5);
    for (int l=0; l<lMax; l++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            double vx = mesh.getVelocityX(kx);
            for (int ky=0; ky<nvy; ky++)
            {
                double vy = mesh.getVelocityY(ky);
                for (int kz=0; kz<nvz; kz++)
                {
                    double vz = mesh.getVelocityZ(kz);
                    double val = weightsX[kx]*weightsY[ky]*weightsZ[kz]*M(l,kz+ky*nvz+kx*nvz*nvy);
                    integral(l,0) += val;
                    integral(l,1) += val*vx;
                    integral(l,2) += val*vy;
                    integral(l,3) += val*vz;
                    integral(l,4) += val*(vx*vx+vy*vy+vz*vz);
                }
            }
        }
    }
    integral = integral * (dvx*dvy*dvz/27.0);
    return integral;
}

double NewtonCotes::integrate(Matrix f, int lMax, std::function<double(int,double)> basisFunction, double x)
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);
    double integral;

    for (int kx=0; kx<nvx; kx++)
    {
        for (int ky=0; ky<nvy; ky++)
        {
            for (int kz=0; kz<nvz; kz++)
            {
                double val = fabs(SpecialFunctions::getF(f,lMax,basisFunction,kz+ky*nvz+kx*nvz*nvx,x));
                integral += weightsX[kx]*weightsY[ky]*weightsZ[kz]*val*log(val);
            }
        }
    }
    integral *= -dvx*dvy*dvz/27.0;
    
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
    double vx = mesh.getVelocityX(k); //NOTE: needs to be updated to reflect 3V
    double exponent = 0;

    for (int l=0; l<lMax; l++)
    {
        exponent += basisFunction(l,x)*(alpha(0,l)-vx*alpha(1,l)-pow(vx,2)*alpha(2,l));
    }

    return pow(vx,power)*exp(exponent);
}

Vector NewtonCotes::computeWeights(int nv)
{
    Vector weights(nv); 
    weights[0] = 1.0;
    weights[nv-1] = 1.0;
    for (int k=1; k<nv-1; k+=2)
    {
        weights[k] = 4.0;
    }
    for (int k=2; k<nv-1; k+=2)
    {
        weights[k] = 2.0;
    }

    return weights;
}