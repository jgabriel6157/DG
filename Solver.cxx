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
      M_invS(lMax,lMax), M_invF1Minus(lMax,lMax), M_invF0Minus(lMax,lMax), M_invF1Plus(lMax,lMax), M_invF0Plus(lMax,lMax),
      uPre(lMax,mesh.getNX()*mesh.getNVX()), uIntermediate(lMax,mesh.getNX()*mesh.getNVX()), uPost(lMax,mesh.getNX()*mesh.getNVX()) {}

//deconstructor
Solver::~Solver() {}

void Solver::createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    Matrix M(lMax,lMax);
    Matrix S(lMax,lMax);
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

        for (int k=0; k<nvx; k++)
        {
            double vx = mesh.getVelocity(k);
            double fluxFactorMinus = (1.0+SpecialFunctions::sign(vx))/2.0;
            double fluxFactorPlus = (1.0-SpecialFunctions::sign(vx))/2.0;

            //calculate feq
            Vector feq = fitMaxwellian(basisFunction, rho, u, rt, vx, j);
            for (int l=0; l<lMax; l++)
            {
                // std::cout << feq[l] << "\n";
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
                uAfter(l,k+j*nvx)+=nu*(feq[l]-uBefore(l,k+j*nvx));
                uAfter(l,k+j*nvx)*=dt;
                uAfter(l,k+j*nvx)+=uBefore(l,k+j*nvx);
                
                uAfter(l,k+j*nvx)*=timesFactor;
                uAfter(l,k+j*nvx)+=plusFactor*uPre(l,k+j*nvx); //Note that uPre != uBefore
            }
        }
    }
}

void Solver::advance(std::function<double(int,double)> basisFunction)
{
    //First stage of solver
    advanceStage(uPre, uPost, 0.0, 1.0, basisFunction);
    //Second stage of solver
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
        for (int i=0; i<100; i++)
        {
            x[i] = leftVertex+i*dx/99.0;
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

double Solver::getF(Matrix& uPre, std::function<double(int,double)> basisFunction, int lMax, int j, double x)
{
    double f = 0;

    for (int l=0; l<lMax; l++)
    {
        f += uPre(l,j)*basisFunction(l,x);
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

double Solver::computeMaxwellian(double rho, double u, double rt, double vx)
{
    return rho*exp(-pow(vx-u,2)/(2.0*rt))/pow(2.0*M_PI*rt,0.5);
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
        double density = computeMoment(rho, basisFunction,lMax,x);
        double meanVelocity = computeMoment(u, basisFunction,lMax,x)/density;
        double temperature = (computeMoment(rt, basisFunction,lMax,x)-density*pow(meanVelocity,2))/density;
        y[i] = computeMaxwellian(density,meanVelocity,temperature,vx);
        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}