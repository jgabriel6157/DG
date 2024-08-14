#include <functional>
#include <cmath>
#include <iostream>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include "GaussianQuadrature.hxx"
#include "Solver.hxx"
#include "Mesh.hxx"
#include "NewtonCotes.hxx"

//constructor
Solver::Solver(const Mesh& mesh, double dt, double a, int lMax, double alpha) 
    : mesh(mesh), integrator(mesh),dt(dt), a(a), lMax(lMax), alpha(alpha),
      M_invS(lMax,lMax), M_invT(lMax,lMax*lMax), M_invF1Minus(lMax,lMax), M_invF0Minus(lMax,lMax), M_invF1Plus(lMax,lMax), M_invF0Plus(lMax,lMax),
      uPre(lMax,(mesh.getNX()+2)*mesh.getNVX()), uIntermediate(lMax,(mesh.getNX()+2)*mesh.getNVX()), uPost(lMax,(mesh.getNX()+2)*mesh.getNVX()) {}

//deconstructor
Solver::~Solver() {}

void Solver::createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    Matrix M(lMax,lMax);
    Matrix S(lMax,lMax);
    Matrix T(lMax,lMax*lMax);
    Matrix F1Minus(lMax,lMax);
    Matrix F0Minus(lMax,lMax);
    Matrix F1Plus(lMax,lMax);
    Matrix F0Plus(lMax,lMax);   
    Matrix M_inv(lMax,lMax);

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M(i,j) = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,quadratureOrder,roots,weights)/2;
            S(i,j) = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,quadratureOrder,roots,weights);
            for (int k=0; k<lMax; k++)
            {
                T(i,j+k*lMax) = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,basisFunction,k,quadratureOrder,roots,weights)/2;
            }
            F1Minus(i,j) = (basisFunction(i,1))*(basisFunction(j,1));
            F0Minus(i,j) = (basisFunction(i,-1))*(basisFunction(j,1));
            F1Plus(i,j) = (basisFunction(i,1))*(basisFunction(j,-1));
            F0Plus(i,j) = (basisFunction(i,-1))*(basisFunction(j,-1));
            if (fabs(M(i,j)) < 1e-10)
            {
                M(i,j) = 0;
            }
            if (fabs(S(i,j)) < 1e-10)
            {
                S(i,j) = 0;
            }
        }
    }
    
    M_inv = M.CalculateInverse();

    M_invS = M_inv*S;
    M_invT = M_inv*T;
    M_invF1Minus = M_inv*F1Minus;
    M_invF0Minus = M_inv*F0Minus;
    M_invF1Plus = M_inv*F1Plus;
    M_invF0Plus = M_inv*F0Plus;
}

//initialize using the Least Squares method
void Solver::initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunctionX, std::function<double(double)> inputFunctionVX)
{
    const auto& cells = mesh.getCells();

    double nx = mesh.getNX();
    double nvx = mesh.getNVX();
    for (int j=0; j<nx; j++)
    {
        double dx = cells[j].dx;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;
        
        for (int k=0; k<nvx; k++)
        {
            double vx = mesh.getVelocity(k);
            Vector uInitialize(lMax);
            
            double x;
            Vector y(10);
            Matrix bigX(10,lMax);

            for (int i=0; i<10; i++)
            {
                x = leftVertex+i*dx/9.0;
                y[i] = inputFunctionX(x)*inputFunctionVX(vx);
                for (int l=0; l<lMax; l++)
                {
                    bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
                }
            }

            uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;

            for (int l=0; l<lMax; l++)
            {
                uPre(l,k+j*nvx) = uInitialize[l];
            }
        }
    }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor, std::function<double(int,double)> basisFunction)
{
    double nu = 1.0;
    double A = 2.91e-14;
    double P = 0;
    double X = 0.232;
    double K = 0.39;
    double U = 13.6/20.0;
    double sigma_iz = A*(1+P*sqrt(U))*pow(U,K)*exp(-U)/(X+U);
    
    double ne = 5.0e18;
    double Crec = 4.98e18;
    double cs = sqrt(20.0*9.58134e7);

    double nu_cx = 2.2e-14;
    double ni = ne;
    double ui = 0;
    double Ti = 20;

    const auto& cells = mesh.getCells();

    double nx = mesh.getNX();
    double nvx = mesh.getNVX();
    for (int j=0; j<nx; j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double dx = cells[j].dx;

        Matrix fj(lMax,nvx);
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<lMax; l++)
            {
                fj(l,k) = uPre(l,k+j*nvx);
            }
        }

        Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
        Vector u = integrator.integrate(fj, lMax, 1); //u tilde
        Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde

        Vector rho_i(lMax);
        rho_i[0] = ni;

        for (int k=0; k<nvx; k++)
        {
            // std::cout << j << "\n";
            // std::cout << k << "\n";
            double vx = mesh.getVelocity(k);
            double fluxFactorMinus = (1.0+SpecialFunctions::sign(vx))/2.0;
            double fluxFactorPlus = (1.0-SpecialFunctions::sign(vx))/2.0;

            // calculate Ghost cells
            Vector fL = fitMaxwellian(basisFunction, Crec, cs, 2.0, vx, j);
            Vector fR = fitMaxwellian(basisFunction, Crec, -cs, 2.0, vx, j);

            // fit ion distribution function
            Vector fi = fitMaxwellian(basisFunction, ni, 0, 20, vx, j);

            Matrix M_invS1(lMax,lMax*lMax);
            Matrix M_invS2(lMax,lMax*lMax);

            Vector f_tilde(lMax);
            for (int l=0; l<lMax; l++)
            {
                uBefore(l,k+nx*nvx) = fL[l];
                uBefore(l,k+(nx+1)*nvx) = fR[l];
                f_tilde[l] = uBefore(l,k+j*nvx);
            }
            Vector fCX = fitCX(basisFunction, ni, ui, Ti, rho, f_tilde, k, j);
            // std::cout << "\n" << k << ":\n";
            for (int l=0; l<lMax; l++)
            {
                uAfter(l,k+j*nvx)=0;
                for (int i=0; i<lMax; i++)
                {
                    uAfter(l,k+j*nvx)+=M_invS(l,i)*uBefore(i,k+j*nvx);
                    uAfter(l,k+j*nvx)-=fluxFactorMinus*M_invF1Minus(l,i)*uBefore(i,k+j*nvx);
                    uAfter(l,k+j*nvx)+=fluxFactorMinus*M_invF0Minus(l,i)*uBefore(i,k+leftNeighborIndex*nvx);
                    uAfter(l,k+j*nvx)-=fluxFactorPlus*M_invF1Plus(l,i)*uBefore(i,k+rightNeighborIndex*nvx);
                    uAfter(l,k+j*nvx)+=fluxFactorPlus*M_invF0Plus(l,i)*uBefore(i,k+j*nvx);
                }
                uAfter(l,k+j*nvx)*=vx;
                uAfter(l,k+j*nvx)/=dx;
                // uAfter(l,k+j*nvx)-=ne*uBefore(l,k+j*nvx)*sigma_iz;
                // for (int i=0; i<lMax; i++)
                // {
                //     for (int m=0; m<lMax; m++)
                //     {
                //         M_invS1(l,i)+=M_invT(l,i+m*lMax)*rho_i[m];
                //         M_invS2(l,i)+=M_invT(l,i+m*lMax)*rho[m];
                //     }
                //     uAfter(l,k+j*nvx)-=nu_cx*M_invS1(l,i)*uBefore(i,k+j*nvx);
                //     uAfter(l,k+j*nvx)+=nu_cx*M_invS2(l,i)*fi[i];
                // }
                uAfter(l,k+j*nvx)-=nu_cx*fCX[l];
                // uAfter(l,k+j*nvx) -= ni*nu_cx*uBefore(l,k+j*nvx);
                // uAfter(l,k+j*nvx) += nu_cx*fCX[l];
                uAfter(l,k+j*nvx)*=dt;
                uAfter(l,k+j*nvx)+=uBefore(l,k+j*nvx);
                
                uAfter(l,k+j*nvx)*=timesFactor;
                uAfter(l,k+j*nvx)+=plusFactor*uPre(l,k+j*nvx); //Note that uPre != uBefore
            }
            // std::cout << fCX[0] << "\n";
            // std::cout << uAfter(0,k+j*nvx) << "\n";
        }
    }
}

void Solver::advance(std::function<double(int,double)> basisFunction)
{
    //First stage of solver
    advanceStage(uPre, uPost, 0.0, 1.0, basisFunction);
    // Second stage of solver
    advanceStage(uPost, uIntermediate, 3.0/4.0, 1.0/4.0, basisFunction);
    //Third stage of solver
    advanceStage(uIntermediate, uPost, 1.0/3.0, 2.0/3.0, basisFunction);
    
    uPre = uPost;
}

void Solver::slopeLimiter()
{
    const auto& cells = mesh.getCells();
    double u1Lim;

    for (int j=0; j<mesh.getNX(); j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        u1Lim = SpecialFunctions::minmod(uPre(1,j),uPre(0,rightNeighborIndex)-uPre(0,j),uPre(0,j)-uPre(0,leftNeighborIndex));
        if (u1Lim-uPre(1,j)>1e-10)
        {
            uPre(1,j) = u1Lim;
            for (int l=2; l<lMax; l++)
            {
                uPre(l,j) = 0;
            }
        }
    }
}

const double Solver::getSolution(int l, int j)
{
    return uPre(l,j);
}

const double Solver::getError(int tMax, std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction)
{
    const auto& cells = mesh.getCells();
    double error = 0;
    double solutionSum = 0;
    for (int j=0; j<mesh.getNX(); j++)
    {
        double y[10] = {0}, sol[10], x[10];
        double dx = cells[j].dx;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;
        for (int i=0; i<10; i++)
        {
            x[i] = leftVertex+i*dx/9.0;
            for (int l=0; l<lMax; l++)
            {
                y[i] += uPre(l,j)*basisFunction(l,(2.0/dx)*(x[i]-xj));
            }
            // sol[i] = inputFunction(x[i]);
            sol[i] = inputFunction(x[i]-2.0*M_PI*tMax*dt);
            // sol[i] = sin(x[i]-2.0*M_PI*tMax*dt);
            // sol[i] = exp(-1.0*pow(x[i]-1.0*M_PI-2.0*M_PI*(tMax*dt),2.0));
            // if ((x[i]<M_PI-1.0)||(x[i]>M_PI+1.0))
            // {
            //     sol[i] = 0;
            // }
            // else
            // {
            //     sol[i] = 1;
            // }
            error+=pow(y[i]-sol[i],2.0);
            solutionSum+=pow(sol[i],2.0);
        }
    }
    return sqrt(error/solutionSum);
}

double Solver::getF(std::function<double(int,double)> basisFunction, int lMax, int j, int k, double x)
{
    double f = 0;
    int nvx = mesh.getNVX();

    for (int l=0; l<lMax; l++)
    {
        f += uPre(l,k+j*nvx)*basisFunction(l,x);
    }

    return f;
}

Vector Solver::getMoments(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    Vector moments(3);
    double mass = 0;
    double momentum = 0;
    double energy = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    double nvx = mesh.getNVX();
    for (int j=0; j<mesh.getNX(); j++)
    {
        Matrix fj(lMax,nvx);
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<lMax; l++)
            {
                fj(l,k) = uPre(l,k+j*nvx);
            }
        }

        Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
        Vector u = integrator.integrate(fj, lMax, 1); //u tilde
        Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde
        for (int i=0; i<quadratureOrder; i++)
        {
            mass += weights[i]*computeMoment(rho, basisFunction, lMax, roots[i]);
            momentum += weights[i]*computeMoment(u, basisFunction, lMax, roots[i]);
            energy += weights[i]*computeMoment(rt, basisFunction, lMax, roots[i]);
        }
    }
    moments[0] = mass;
    moments[1] = momentum;
    moments[2] = energy;
    return moments;
}


double Solver::computeMoment(Vector moment, std::function<double(int,double)> basisFunction, int lMax, double x)
{
    double momentValue = 0;

    for (int l=0; l<lMax; l++)
    {
        momentValue += moment[l]*basisFunction(l,x);
    }

    return momentValue;
}

//initialize using the Least Squares method
Vector Solver::fitMaxwellian(std::function<double(int,double)> basisFunction, Vector rho, Vector u, Vector rt, double vx, int j)
{
    const auto& cells = mesh.getCells();

    double dx = cells[j].dx;
    double leftVertex = cells[j].vertices[0];
    double xj = leftVertex+dx/2.0;
    
    Vector uInitialize(lMax);
    
    double x;
    Vector y(10);
    Matrix bigX(10,lMax);

    for (int i=0; i<10; i++)
    {
        x = leftVertex+i*dx/9.0;
        double density = computeMoment(rho, basisFunction,lMax,x); //Might need to be 2.0*(x-xj)/dx instead of x
        double meanVelocity = computeMoment(u, basisFunction,lMax,x)/density; //Might need to be 2.0*(x-xj)/dx instead of x
        double temperature = (computeMoment(rt, basisFunction,lMax,x)-density*pow(meanVelocity,2))/density; //Might need to be 2.0*(x-xj)/dx instead of x
        y[i] = SpecialFunctions::computeMaxwellian(density,meanVelocity,temperature,vx);
        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}

Vector Solver::fitMaxwellian(std::function<double(int,double)> basisFunction, double density, double meanVelocity, double temperature, double vx, int j)
{
    const auto& cells = mesh.getCells();

    double dx = cells[j].dx;
    double leftVertex = cells[j].vertices[0];
    double xj = leftVertex+dx/2.0;
    
    Vector uInitialize(lMax);
    
    double x;
    Vector y(10);
    Matrix bigX(10,lMax);

    for (int i=0; i<10; i++)
    {
        x = leftVertex+i*dx/9.0;
        y[i] = SpecialFunctions::computeMaxwellian(density,meanVelocity,temperature,vx);
        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}

Vector Solver::getDensity(int j)
{
    int nvx = mesh.getNVX();
    Matrix fj(lMax,nvx);
    for (int k=0; k<nvx; k++)
    {
        for (int l=0; l<lMax; l++)
        {
            fj(l,k) = uPre(l,k+j*nvx);
        }
    }

    return integrator.integrate(fj, lMax, 0); //rho tilde
}

Vector Solver::fitCX(std::function<double(int,double)> basisFunction, double density_i, double meanVelocity_i, double temperature_i, Vector rho_n, Vector f_tilde, int k, int j)
{
    const auto& cells = mesh.getCells();

    double dx = cells[j].dx;
    double leftVertex = cells[j].vertices[0];
    double xj = leftVertex+dx/2.0;

    double vx = mesh.getVelocity(k);
    
    Vector uInitialize(lMax);

    int res = 10;
    double x;
    Vector y(res);
    Matrix bigX(res,lMax);

    // std::cout << vx << ":\n\n";

    for (int i=0; i<res; i++)
    {
        x = leftVertex+i*dx/(res-1.0);
        double density_n = computeMoment(rho_n, basisFunction, lMax, 2.0*(x-xj)/dx);
        double f_n = computeMoment(f_tilde, basisFunction, lMax, 2.0*(x-xj)/dx);
        y[i] = density_i * f_n - density_n * SpecialFunctions::computeMaxwellian(density_i,meanVelocity_i,temperature_i,vx);

        // y[i] = density_n*SpecialFunctions::computeMaxwellian(density_i,meanVelocity_i,temperature_i,vx);

        // std::cout << x << "       " << y[i] << "\n";
        // std::cout << x << "\n";
        // std::cout << density_i << "\n";
        // std::cout << f_n << "\n";
        // std::cout << density_n << "\n";
        // std::cout << SpecialFunctions::computeMaxwellian(density_i,meanVelocity_i,temperature_i,vx) << "\n\n";
        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}