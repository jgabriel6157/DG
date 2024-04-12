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
    : mesh(mesh), dt(dt), a(a), lMax(lMax), alpha(alpha), M_invT(lMax,lMax*lMax),
      M_invS(lMax,lMax), M_invF1Minus(lMax,lMax), M_invF0Minus(lMax,lMax), M_invF1Plus(lMax,lMax), M_invF0Plus(lMax,lMax),
      uPre(lMax,mesh.getNVX()), uIntermediate(lMax,mesh.getNVX()), uPost(lMax,mesh.getNVX()), vxWeights(lMax,mesh.getNVX()),
      gPre(lMax,mesh.getNVX()), gPost(lMax,mesh.getNVX()) {}

//deconstructor
Solver::~Solver() {}

void Solver::createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    Matrix M(lMax,lMax);
    Matrix T(lMax,lMax*lMax);
    Matrix S(lMax,lMax);
    Matrix F1Minus(lMax,lMax);
    Matrix F0Minus(lMax,lMax);
    Matrix F1Plus(lMax,lMax);
    Matrix F0Plus(lMax,lMax);   
    Matrix M_inv(lMax,lMax);
    double fluxFactorPlus = (1.0+SpecialFunctions::sign(a)*(1.0-alpha))/2.0;
    double fluxFactorMinus = (1.0-SpecialFunctions::sign(a)*(1.0-alpha))/2.0;

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M(i,j) = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,quadratureOrder,roots,weights)/2;
            S(i,j) = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,quadratureOrder,roots,weights);
            for (int k=0; k<lMax; k++)
            {
                T(i,j+k*lMax) = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,basisFunction,k,quadratureOrder,roots,weights);
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

    // std::cout << M_invF0Minus(0,0) << "\n";
    // std::cout << M_invF1Minus(0,0) << "\n";
    // std::cout << M_invF0Plus(0,0) << "\n";
    // std::cout << M_invF1Plus(0,0) << "\n";
    // std::cout << M_invS(0,0) << "\n";
    // std::cout << M_invT(0,0) << "\n";
}

//initialize using the Least Squares method
void Solver::initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction)
{
    const auto& cells = mesh.getCells();

    for (int j=0; j<mesh.getNVX(); j++)
    {
        double dx = cells[j].dvx;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;
        Vector uInitialize(lMax);
        Vector uInitializeW(lMax);
        
        double x;
        Vector y(10);
        Vector yW(10);
        Matrix bigX(10,lMax);

        for (int i=0; i<10; i++)
        {
            x = leftVertex+i*dx/9.0;
            y[i] = inputFunction(x);
            yW[i] = x;
            for (int l=0; l<lMax; l++)
            {
                bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
            }
        }

        uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
        uInitializeW = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*yW;

        for (int l=0; l<lMax; l++)
        {
            uPre(l,j) = uInitialize[l];
            vxWeights(l,j) = uInitializeW[l];
        }
    }
}

void Solver::advanceVelocity(Matrix& uBefore, Matrix& uAfter, Matrix& gBefore, Matrix& gAfter, int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    const auto& cells = mesh.getCells();
    const int nvx = mesh.getNVX();
    double mass = getMass(quadratureOrder, basisFunction);
    double momentum = getMomentum(quadratureOrder, basisFunction);
    double energy = getEnergy(quadratureOrder, basisFunction);

    double fMax = getF(uPre, basisFunction, lMax, mesh.getNVX()-1,1);
    double fMin = getF(uPre, basisFunction, lMax, 0,-1);
    double temperature = getTemperature(mass, momentum, energy, fMax, fMin);
    double velocity = getVelocity(mass, momentum, temperature, fMax, fMin);

    for (int j=0; j<nvx; j++)
    // for (int j=25; j<26; j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double dvx = cells[j].dvx;
        
        for (int l=0; l<lMax; l++)
        {
            gAfter(l,j) = 0;
            for (int i=0; i<lMax; i++)
            {
                gAfter(l,j)-=M_invS(l,i)*uBefore(i,j);
                // std::cout << gAfter(l,j) << "\n";
                if (j==0)
                {
                    gAfter(l,j)-=M_invF0Plus(l,i)*uBefore(i,j);
                }
                else
                {
                    gAfter(l,j)-=M_invF0Minus(l,i)*uBefore(i,leftNeighborIndex);
                }
                // std::cout << gAfter(l,j) << "\n";
                gAfter(l,j)+=M_invF1Minus(l,i)*uBefore(i,j);
                // std::cout << gAfter(l,j) << "\n";
            }
            gAfter(l,j)/=dvx;
            gAfter(l,j)*=sqrt(temperature);
            // std::cout << gAfter(l,j) << "\n";
        }
    }
    alpha = fmax(fabs(cells[0].vertices[0]-velocity),fabs(cells[nvx-1].vertices[1]-velocity));

    for (int j=0; j<nvx; j++)
    // for (int j=25; j<26; j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double leftVertex = cells[j].vertices[0];
        double rightVertex = cells[j].vertices[1];
        double dvx = cells[j].dvx;
        if (j==25)
        {
            // std::cout << alpha << "\n";
            // std::cout << temperature << "\n";
            // std::cout << gAfter(0,26) << "\n";
            // std::cout << gAfter(0,25) << "\n";
        }

        for(int l=0; l<lMax; l++)
        {
            uAfter(l,j) = 0;
            for (int i=0; i<lMax; i++)
            {
                if (j==nvx-1)
                {
                    // uAfter(l,j) += M_invF1Minus(l,i)*(0.5*(alpha)*uBefore(i,j));
                    // uAfter(l,j) += 0;
                }
                else
                {
                    uAfter(l,j) += M_invF1Plus(l,i)*(0.5*(rightVertex-velocity+alpha)*uBefore(i,rightNeighborIndex)+sqrt(temperature)*gAfter(i,rightNeighborIndex));
                    uAfter(l,j) += M_invF1Minus(l,i)*(0.5*(rightVertex-velocity-alpha)*uBefore(i,j));
                }
                // std::cout << uAfter(l,j) << "\n";
                // uAfter(l,j) += M_invF1Minus(l,i)*(0.5*(rightVertex-velocity-alpha)*uBefore(i,j));
                // std::cout << uAfter(l,j) << "\n";
                // uAfter(l,j) -= M_invF0Plus(l,i)*(0.5*(leftVertex-velocity+alpha)*uBefore(i,j)+sqrt(temperature)*gAfter(i,j));
                // std::cout << uAfter(l,j) << "\n";
                if (j==0)
                {
                    // uAfter(l,j) -= M_invF0Plus(l,i)*(0.5*(leftVertex-velocity-alpha)*uBefore(i,j));
                }
                else
                {
                    uAfter(l,j) -= M_invF0Plus(l,i)*(0.5*(leftVertex-velocity+alpha)*uBefore(i,j)+sqrt(temperature)*gAfter(i,j));
                    uAfter(l,j) -= M_invF0Minus(l,i)*(0.5*(leftVertex-velocity-alpha)*uBefore(i,leftNeighborIndex));
                }
                // std::cout << uAfter(l,j) << "\n";
                Matrix M_invV(lMax,lMax);
                for (int k=0; k<lMax; k++)
                {
                    M_invV(l,i)+=M_invT(l,k+i*lMax)*vxWeights(k,j);
                }
                uAfter(l,j) -= M_invV(l,i)*uBefore(i,j);
                // std::cout << uAfter(l,j) << "\n";
                uAfter(l,j) += M_invS(l,i)*(velocity*uBefore(i,j)-sqrt(temperature)*gAfter(i,j));
                // std::cout << uAfter(l,j) << "\n";
            }
            uAfter(l,j)*=dt;
            uAfter(l,j)/=dvx;
            // std::cout << uAfter(l,j) << "\n";
            uAfter(l,j)+=uBefore(l,j);
            // std::cout << uAfter(l,j) << "\n";
        }   
    }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor)
{
    const auto& cells = mesh.getCells();

    for (int j=0; j<mesh.getNVX(); j++)
    {
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double dx = cells[j].dvx;
        for (int l=0; l<lMax; l++)
        {
            uAfter(l,j)=0;
            for (int i=0; i<lMax; i++)
            {
                uAfter(l,j)+=M_invS(l,i)*uBefore(i,j);
                uAfter(l,j)-=M_invF1Minus(l,i)*uBefore(i,j);
                uAfter(l,j)+=M_invF0Minus(l,i)*uBefore(i,leftNeighborIndex);
                uAfter(l,j)-=M_invF1Plus(l,i)*uBefore(i,rightNeighborIndex);
                uAfter(l,j)+=M_invF0Plus(l,i)*uBefore(i,j);
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

void Solver::advance(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    // //First stage of solver
    // advanceStage(uPre, uPost, 0.0, 1.0);
    // //Second stage of solver
    // advanceStage(uPost, uIntermediate, 3.0/4.0, 1.0/4.0);
    // //Third stage of solver
    // advanceStage(uIntermediate, uPost, 1.0/3.0, 2.0/3.0);

    advanceVelocity(uPre,uPost,gPre,gPost,quadratureOrder,basisFunction);
    
    uPre = uPost;
}

void Solver::slopeLimiter()
{
    const auto& cells = mesh.getCells();
    double u1Lim;

    for (int j=0; j<mesh.getNVX(); j++)
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
    for (int j=0; j<mesh.getNVX(); j++)
    {
        double y[10] = {0}, sol[10], x[10];
        double dx = cells[j].dvx;
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
    const auto& cells = mesh.getCells();

    for (int j=0; j<mesh.getNVX(); j++)
    {
        Cell cell = cells[j];
        double dvx = cell.dvx;

        for (int i=0; i<quadratureOrder; i++)
        {
            mass += weights[i]*getF(uPre, basisFunction, lMax, j, roots[i])*dvx/2;
        }
    }
    return mass;
}

double Solver::getMomentum(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    double momentum = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    const auto& cells = mesh.getCells();
    double nvx = mesh.getNVX();

        for (int k=0; k<nvx; k++)
        {
            Cell cell = cells[k];
            double dvx = cell.dvx;
            double vxStart = cell.vertices[0];
            double vxj = vxStart+dvx/2.0;

            for (int l=0; l<quadratureOrder; l++)
            {
                double vx = vxj+((dvx/2.0)*roots[l]);
                momentum += vx*weights[l]*getF(uPre, basisFunction, lMax, k, roots[l])*dvx/2;
            }
        }
    return momentum;

}

double Solver::getEnergy(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    double energy = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    const auto& cells = mesh.getCells();
    double nvx = mesh.getNVX();

    for (int k=0; k<nvx; k++)
    {
        Cell cell = cells[k];
        double dvx = cell.dvx;
        double vxStart = cell.vertices[0];
        double vxj = vxStart+dvx/2.0;

        for (int l=0; l<quadratureOrder; l++)
        {
            double vx = vxj+((dvx/2.0)*roots[l]);
            energy += (1.0/2.0)*pow(vx,2)*weights[l]*getF(uPre, basisFunction, lMax, k, roots[l])*dvx/2;
        }
    }
    return energy;

}

double Solver::getEntropy(int quadratureOrder, std::function<double(int,double)> basisFunction)
{
    double entropy = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    for (int k=0; k<mesh.getNVX(); k++)
    {
        for (int l=0; l<quadratureOrder; l++)
        {
            entropy -= weights[l]*fabs(getF(uPre, basisFunction, lMax, k, roots[l]))*std::log(fabs(getF(uPre, basisFunction, lMax, k, roots[l])));
        }
    }
    return entropy;
}

double Solver::getTemperature(double mass, double momentum, double energy, double fMax, double fMin)
{
    const auto& cells = mesh.getCells();
    double nvx = mesh.getNVX();
    double vMax = cells[nvx-1].vertices[1];
    double vMin = cells[0].vertices[0];

    double A = fMax-fMin;
    double B = vMax*fMax-vMin*fMin;

    return (mass*2*energy-pow(momentum,2))/(pow(mass,2)-(mass*B)+(momentum*A));
}

double Solver::getVelocity(double mass, double momentum, double temperature, double fMax, double fMin)
{
    // const auto& cells = mesh.getCells();
    // double nvx = mesh.getNVX();
    // double vMax = cells[nvx-1].vertices[1];
    // double vMin = cells[0].vertices[0];

    double A = fMax-fMin;
    // double B = vMax*fMax-vMin*fMin;

    // return (momentum+A*(mass*energy-pow(momentum,2))/(pow(mass,2)-(mass*B)+(momentum*A)))/mass;
    return (momentum+temperature*A)/mass;
}
