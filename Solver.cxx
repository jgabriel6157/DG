#include <functional>
#include <cmath>
#include <iostream>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include "GaussianQuadrature.hxx"
#include "Solver.hxx"
#include "Mesh.hxx"

//constructor
Solver::Solver(const Mesh& mesh, double dt, double a, int lMax, double alpha) 
    : mesh(mesh), dt(dt), a(a), lMax(lMax), alpha(alpha),
      M_invS(lMax,lMax), M_invF1(lMax,lMax), M_invF2(lMax,lMax), M_invF3(lMax,lMax), M_invF4(lMax,lMax),
      uPre(lMax,mesh.getNumCells()), uIntermediate(lMax,mesh.getNumCells()), uPost(lMax,mesh.getNumCells()) {}

//deconstructor
Solver::~Solver() {}

void Solver::createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    Matrix M(lMax,lMax);
    Matrix S(lMax,lMax);
    Matrix F1(lMax,lMax);
    Matrix F2(lMax,lMax);
    Matrix F3(lMax,lMax);
    Matrix F4(lMax,lMax);   
    Matrix M_inv(lMax,lMax);
    double fluxFactorPlus = (1.0+SpecialFunctions::sign(a)*(1.0-alpha))/2.0;
    double fluxFactorMinus = (1.0-SpecialFunctions::sign(a)*(1.0-alpha))/2.0;

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M(i,j) = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,quadratureOrder,roots,weights)/2;
            S(i,j) = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,quadratureOrder,roots,weights);
            F1(i,j) = fluxFactorPlus*(basisFunction(i,1))*(basisFunction(j,1));
            F2(i,j) = fluxFactorPlus*(basisFunction(i,-1))*(basisFunction(j,1));
            F3(i,j) = fluxFactorMinus*(basisFunction(i,1))*(basisFunction(j,-1));
            F4(i,j) = fluxFactorMinus*(basisFunction(i,-1))*(basisFunction(j,-1));
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
    M_invF1 = M_inv*F1;
    M_invF2 = M_inv*F2;
    M_invF3 = M_inv*F3;
    M_invF4 = M_inv*F4;
}

//initialize using the Least Squares method
void Solver::initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction)
{
    const auto& cells = mesh.getCells();

    for (int j=0; j<mesh.getNumCells(); j++)
    {
        double dx = cells[j].cellLength;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;
        Vector uInitialize(lMax);
        
        double x;
        Vector y(10);
        Matrix bigX(10,lMax);

        for (int i=0; i<10; i++)
        {
            x = leftVertex+i*dx/9.0;
            y[i] = inputFunction(x);
            for (int l=0; l<lMax; l++)
            {
                bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
            }
        }

        uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;

        for (int l=0; l<lMax; l++)
        {
            uPre(l,j) = uInitialize[l];
        }
    }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor)
{
    const auto& cells = mesh.getCells();

    for (int j=0; j<mesh.getNumCells(); j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double dx = cells[j].cellLength;
        for (int l=0; l<lMax; l++)
        {
            uAfter(l,j)=0;
            for (int i=0; i<lMax; i++)
            {
                uAfter(l,j)+=M_invS(l,i)*uBefore(i,j);
                uAfter(l,j)-=M_invF1(l,i)*uBefore(i,j);
                uAfter(l,j)+=M_invF2(l,i)*uBefore(i,leftNeighborIndex);
                uAfter(l,j)-=M_invF3(l,i)*uBefore(i,rightNeighborIndex);
                uAfter(l,j)+=M_invF4(l,i)*uBefore(i,j);
            }
            uAfter(l,j)*=a;
            uAfter(l,j)*=dt;
            uAfter(l,j)/=dx;
            uAfter(l,j)+=uBefore(l,j);
            
            uAfter(l,j)*=timesFactor;
            uAfter(l,j)+=plusFactor*uPre(l,j); //Note that uPre != uBefore
        }
    }
}

void Solver::advance()
{
    //First stage of solver
    advanceStage(uPre, uPost, 0.0, 1.0);
    //Second stage of solver
    advanceStage(uPost, uIntermediate, 3.0/4.0, 1.0/4.0);
    //Third stage of solver
    advanceStage(uIntermediate, uPost, 1.0/3.0, 2.0/3.0);
    
    uPre = uPost;
}

void Solver::slopeLimiter()
{
    const auto& cells = mesh.getCells();
    double u1Lim;

    for (int j=0; j<mesh.getNumCells(); j++)
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
    for (int j=0; j<mesh.getNumCells(); j++)
    {
        double y[10] = {0}, sol[10], x[10];
        double dx = cells[j].cellLength;
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

double Solver::getF(Matrix& uPre, std::function<double(int,double)> basisFunction, int lMax, int j, double x)
{
    double f = 0;

    for (int l=0; l<lMax; l++)
    {
        f += uPre(l,j)*basisFunction(l,x);
    }

    return f;
}

double Solver::getMass(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    double mass = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    for (int j=0; j<mesh.getNumCells(); j++)
    {
        for (int i=0; i<quadratureOrder; i++)
        {
            mass += weights[i]*getF(uPre, basisFunction, lMax, j, roots[i]);
        }
    }
    return mass;
}


