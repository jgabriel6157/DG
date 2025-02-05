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

Vector NewtonCotes::integrateAlphaMoments(Matrix alpha, std::function<double(int,double)> basisFunction, double x, int lMax)
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

    Vector integral(5);
    for (int kx=0; kx<nvx; kx++)
    {
        double vx = mesh.getVelocityX(kx);
        for (int ky=0; ky<nvy; ky++)
        {
            double vy = mesh.getVelocityY(ky);
            for (int kz=0; kz<nvz; kz++)
            {
                double vz = mesh.getVelocityZ(kz);
                double feq = testMaxwellian(alpha, basisFunction, x, vx, vy, vz, lMax);
                double val = weightsX[kx]*weightsY[ky]*weightsZ[kz]*feq;
                integral[0] += val;
                integral[1] += val*vx;
                integral[2] += val*vy;
                integral[3] += val*vz;
                integral[4] += val*(vx*vx+vy*vy+vz*vz);
            }
        }
    }
    integral = integral * (dvx*dvy*dvz/27.0);
    return integral;
}

Vector NewtonCotes::getVelocityIntegrals(Matrix alpha, std::function<double(int,double)> basisFunction, double x, int lMax)
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

    Vector integral(15);
    for (int kx=0; kx<nvx; kx++)
    {
        double vx = mesh.getVelocityX(kx);
        double vx2 = vx*vx;
        for (int ky=0; ky<nvy; ky++)
        {
            double vy = mesh.getVelocityY(ky);
            double vy2 = vy*vy;
            double vxvy = vx*vy;
            for (int kz=0; kz<nvz; kz++)
            {
                double vz = mesh.getVelocityZ(kz);
                double vz2 = vz*vz;
                double v2 = vx2+vy2+vz2;
                double feq = testMaxwellian(alpha, basisFunction, x, vx, vy, vz, lMax);
                double val = weightsX[kx]*weightsY[ky]*weightsZ[kz]*feq;
                integral[0] += val;
                integral[1] += val*vx;
                integral[2] += val*vy;
                integral[3] += val*vz;
                integral[4] -= val*v2; //val*(vx*vx+vy*vy+vz*vz);
                integral[5] += val*vx2; //val*vx*vx;
                integral[6] += val*vxvy; //val*vx*vy;
                integral[7] += val*vx*vz;
                integral[8] -= val*vx*v2; //val*vx*(vx*vx+vy*vy+vz*vz);
                integral[9] += val*vy2; //val*vy*vy;
                integral[10] += val*vy*vz;
                integral[11] -= val*vy*v2; //val*vy*(vx*vx+vy*vy+vz*vz);
                integral[12] += val*vz2; //val*vz*vz;
                integral[13] -= val*vz*v2; //val*vz*(vx*vx+vy*vy+vz*vz);
                integral[14] += val*v2*v2; //val*(vx*vx+vy*vy+vz*vz)*(vx*vx+vy*vy+vz*vz);
            }
        }
    }
    integral = integral * (dvx*dvy*dvz/27.0);
    return integral;
}

Vector NewtonCotes::integrate3fnCXavg(Matrix M, int lMax, double Ti, double ui)
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    double J = dvx*dvy*dvz/27.0;
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    
    for (int kx=0; kx<nvx; kx++)
    {
        double vx = mesh.getVelocityX(kx);
        double vx_ui2 = (vx-ui)*(vx-ui);
        for (int ky=0; ky<nvy; ky++)
        {
            double vy2 = mesh.getVelocityY(ky)*mesh.getVelocityY(ky);
            double weightXY = weightsX[kx]*weightsY[ky];
            for (int kz=0; kz<nvz; kz++)
            {
                double vz2 = mesh.getVelocityZ(kz)*mesh.getVelocityZ(kz);
                double weightXYZ = weightXY*weightsZ[kz];
                double E = 0.5*(vx_ui2+vy2+vz2);
                double sigmavg = weightXYZ*SpecialFunctions::computeSigmav(Ti,E)*(1e18)/(9822.766369779);
                for (int l=0; l<lMax; l++)
                {
                    integral[l] += M(l,kz+ky*nvz+kx*nvz*nvy)*sigmavg;
                }
            }
        }
    }
    integral = integral*J;
    return integral;
}

Vector NewtonCotes::integrate3fnCX(Matrix M, int lMax, double vx, double vy, double vz)
{
    double dvx = mesh.getDVX();
    int nvx = mesh.getNVX();
    double dvy = mesh.getDVY();
    int nvy = mesh.getNVY();
    double dvz = mesh.getDVZ();
    int nvz = mesh.getNVZ();
    double J = dvx*dvy*dvz/27.0;
    Vector weightsX = computeWeights(nvx);
    Vector weightsY = computeWeights(nvy);
    Vector weightsZ = computeWeights(nvz);

    Vector integral(lMax);
    for (int kx=0; kx<nvx; kx++)
    {
        double vxRel2 = (vx-mesh.getVelocityX(kx))*(vx-mesh.getVelocityX(kx));
        for (int ky=0; ky<nvy; ky++)
        {
            double vyRel2 = (vy-mesh.getVelocityY(ky))*(vy-mesh.getVelocityY(ky));
            double weightXY = weightsX[kx]*weightsY[ky];
            for (int kz=0; kz<nvz; kz++)
            {
                double vzRel2 = (vz-mesh.getVelocityZ(kz))*(vz-mesh.getVelocityZ(kz));
                double weightXYZ = weightXY*weightsZ[kz];
                double relVelocity = sqrt(vxRel2+vyRel2+vzRel2);
                double sigma = SpecialFunctions::computeSigma(0.5*relVelocity*relVelocity)*(1e18);
                double sigmaVelocity = weightXYZ*relVelocity*sigma;
                for (int l=0; l<lMax; l++)
                {
                    integral[l] += M(l,kz+ky*nvz+kx*nvz*nvy)*sigmaVelocity;
                }
            }
        }
    }
    integral = integral*J;
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

// double NewtonCotes::integrate(Matrix alpha, std::function<double(int,double)> basisFunction, int power, double x, int lMax)
// {
//     double dvx = mesh.getDVX();
//     int nvx = mesh.getNVX();
//     double integral;

//     integral = testMaxwellian(alpha,basisFunction,power,x,0,lMax)+testMaxwellian(alpha,basisFunction,power,x,nvx-1,lMax);

//     for (int k=1; k<nvx-1; k+=2)
//     {
//         integral += 4.0*testMaxwellian(alpha,basisFunction,power,x,k,lMax);
//     }
//     for (int k=2; k<nvx-1; k+=2)
//     {
//         integral += 2.0*testMaxwellian(alpha,basisFunction,power,x,k,lMax);
//     }

//     integral *= dvx/3.0;

//     return integral;

// }

double NewtonCotes::testMaxwellian(Matrix alpha, std::function<double(int,double)> basisFunction, double x, double vx, double vy, double vz, int lMax)
{
    double exponent = 0;

    for (int l=0; l<lMax; l++)
    {
        exponent += basisFunction(l,x)*(alpha(0,l)+vx*alpha(1,l)+vy*alpha(2,l)+vz*alpha(3,l)-(vx*vx+vy*vy+vz*vz)*alpha(4,l));
    }

    return exp(exponent);
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