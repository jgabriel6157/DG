#include <functional>
#include <cmath>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include "Solver.hxx"

//constructor
Solver::Solver(double dx, double dt, double a, int jMax, int lMax, Matrix& M_invS, Matrix& M_invF1, Matrix& M_invF2, Matrix& M_invF3, Matrix& M_invF4)
    : dx(dx), dt(dt), a(a), jMax(jMax), lMax(lMax), M_invS(M_invS), M_invF1(M_invF1), M_invF2(M_invF2), M_invF3(M_invF3), M_invF4(M_invF4),
      uPre(lMax,jMax), uIntermediate(lMax,jMax), uPost(lMax,jMax) {}

//deconstructor
Solver::~Solver() {}

//initialize using the Least Squares method
void Solver::initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction)
{
    for (int j=0; j<jMax; j++)
    {
        double xj = j*dx+dx/2.0;
        Vector uInitialize(lMax);
        
        double x;
        Vector y(10);
        Matrix bigX(10,lMax);
        Matrix bigXT(lMax,10);
        Matrix bigXprod(lMax,lMax);
        Matrix bigXprodInv(lMax,lMax);
        Matrix hugeX(lMax,10);

        for (int i=0; i<10; i++)
        {
            x = xj-dx/2.0+i*dx/9.0;
            y[i] = inputFunction(x);
            for (int l=0; l<lMax; l++)
            {
                bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
                bigXT(l,i) = basisFunction(l,2.0*(x-xj)/dx);
            }
        }

        bigXprod = bigXT*bigX;

        bigXprodInv = bigXprod.CalculateInverse();

        hugeX = bigXprodInv*bigXT;

        uInitialize = hugeX*y;

        for (int l=0; l<lMax; l++)
        {
            uPre(l,j) = uInitialize[l];
        }
    }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor)
{
    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            uAfter(l,j)=0;
            for (int i=0; i<lMax; i++)
            {
                uAfter(l,j)+=M_invS(l,i)*uBefore(i,j);
                uAfter(l,j)-=M_invF1(l,i)*uBefore(i,j);
                if (j==0)
                {
                    uAfter(l,j)+=M_invF2(l,i)*uBefore(i,jMax-1);
                }
                else
                {
                    uAfter(l,j)+=M_invF2(l,i)*uBefore(i,j-1);
                }
                if (j==jMax-1)
                {
                    uAfter(l,j)-=M_invF3(l,i)*uBefore(i,0);
                }
                else
                {
                    uAfter(l,j)-=M_invF3(l,i)*uBefore(i,j+1);
                }
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
    double u1Lim;

    for (int j=0; j<jMax; j++)
    {
        if (j==0)
        {
            u1Lim = SpecialFunctions::minmod(uPre(1,j),uPre(0,j+1)-uPre(0,j),uPre(0,j)-uPre(0,jMax-1));
            if (u1Lim-uPre(1,j)>1e-10)
            {
                uPre(1,j) = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPre(l,j) = 0;
                }
            }
        }
        else if (j==jMax-1)
        {
            u1Lim = SpecialFunctions::minmod(uPre(1,j),uPre(0,0)-uPre(0,j),uPre(0,j)-uPre(0,j-1));
            if (u1Lim-uPre(1,j)>1e-10)
            {
                uPre(1,j) = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPre(l,j) = 0;
                }
            }
        }
        else
        {
            u1Lim = SpecialFunctions::minmod(uPre(1,j),uPre(0,j+1)-uPre(0,j),uPre(0,j)-uPre(0,j-1));
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
}

const double Solver::getSolution(int l, int j)
{
    return uPre(l,j);
}

const double Solver::getError(int tMax, std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction)
{
    double error = 0;
    double solutionSum = 0;
    for (int j=0; j<jMax; j++)
    {
        double y[10] = {0}, sol[10], x[10];
        for (int i=0; i<10; i++)
        {
            x[i] = j*dx+i*dx/9.0;
            for (int l=0; l<lMax; l++)
            {
                y[i] += uPre(l,j)*basisFunction(l,(2.0/dx)*(x[i]-(j*dx+dx/2.0)));
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






