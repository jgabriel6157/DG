#include <functional>
#include <cmath>
#include <iostream>
#include <cassert>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include "GaussianQuadrature.hxx"
#include "Solver.hxx"
#include "Mesh.hxx"
#include "NewtonCotes.hxx"
#include "NewtonSolver.hxx"

//constructor
Solver::Solver(const Mesh& mesh, double dt, double a, int lMax) 
    : mesh(mesh), integrator(mesh), newtonSolver(mesh), dt(dt), a(a), lMax(lMax), alphaDomain(3*mesh.getNX(), lMax),
      M_invDiag(lMax), M_invS(lMax,lMax), M_invT(lMax,lMax*lMax), M_invF1Minus(lMax,lMax), M_invF0Minus(lMax,lMax), M_invF1Plus(lMax,lMax), M_invF0Plus(lMax,lMax),
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

    for (int l=0; l<lMax; l++)
    {
        M_invDiag[l] = M_inv(l,l);
    }

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
                // x = leftVertex+i*dx/10.0;
                // if (x<0.5)
                // {
                //    y[i] = SpecialFunctions::computeMaxwellian(1.0,0.0,1.0,vx);
                // }
                // else
                // {
                //    y[i] = SpecialFunctions::computeMaxwellian(0.125,0.0,0.8,vx);
                // }
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
    // for (int k=0; k<nvx; k++)
    // {
    //     for (int l=0; l<lMax; l++)
    //     {
    //         uPre(l,k+nx*nvx) = uPre(l,k+0*nvx); //Left BC
    //         uPre(l,k+(nx+1)*nvx) = uPre(l,k+(nx-1)*nvx); //Right BC
    //     }
    // }
}

void Solver::initializeAlpha(std::function<double(int,double)> basisFunction)
{
    const auto& cells = mesh.getCells();

    double nx = mesh.getNX();
    double nvx = mesh.getNVX();
    for (int j=0; j<nx; j++)
    {
        double dx = cells[j].dx;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;

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

        Vector uInitialize(3*lMax);
        Vector y(10*nvx);
        Matrix bigX(10*nvx,3*lMax);
        
        for (int k=0; k<nvx; k++)
        {
            // std::cout << j << "\n";
            // std::cout << k << "\n";
            double vx = mesh.getVelocity(k);
            double x;

            for (int i=0; i<10; i++)
            {
                x = leftVertex+i*dx/9.0;
                double density = computeMoment(rho, basisFunction,lMax,2.0*(x-xj)/dx);
                double meanVelocity = computeMoment(u, basisFunction,lMax,2.0*(x-xj)/dx)/density;
                double temperature = (computeMoment(rt, basisFunction,lMax,2.0*(x-xj)/dx)-density*pow(meanVelocity,2))/density;
                // temperature /= (9.58134e7);
                // if (temperature < 0)
                // {
                //     temperature = 0.1;
                // }
                double arg = SpecialFunctions::computeMaxwellian(density,meanVelocity,temperature,vx);
                // if (arg < 1e-10)
                // {
                //     arg = 1e-10;
                // }
                y[i+k*10] = log(arg);
                // if (j==31 && k==0)
                // {
                //     std::cout << "x = " << x << "\n";
                //     std::cout << density << "\n";
                //     std::cout << meanVelocity << "\n";
                //     std::cout << temperature << "\n";
                //     std::cout << arg << "\n";
                //     std::cout << y[i+k*10] << "\n";
                // }
                for (int m=0; m<3; m++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        if (m==0)
                        {
                            bigX(i+k*10,m+l*3) = basisFunction(l,2.0*(x-xj)/dx);
                        }
                        else
                        {
                            bigX(i+k*10,m+l*3) = -basisFunction(l,2.0*(x-xj)/dx)*pow(vx,m);
                        }
                    }
                }
            }
        }

        uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;

        // if (j==31)
        // {
        //     rho.Print();
        //     u.Print();
        //     rt.Print();
        //     // ((bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()).Print();
        //     y.Print();
        // }

        for (int m=0; m<3; m++)
        {
            for (int l=0; l<lMax; l++)
            {
                alphaDomain(m+j*3,l) = uInitialize[m+l*3];
            }
        }
    }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor, std::function<double(int,double)> basisFunction, int quadratureOrder)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    double nu = 100.0;

    double A = 2.91e-14;
    double P = 0;
    double X = 0.232;
    double K = 0.39;
    double U = 13.6/20.0;
    double sigma_iz = A*(1+P*sqrt(U))*pow(U,K)*exp(-U)/(X+U);

    sigma_iz *= 1e18;
    sigma_iz /= 9822.766369779;
    
    double ne = 5.0;
    double Crec = 4.98;
    double cs = sqrt(20.0);

    double nu_cx = (2.2e-14)*(1e18)/(9822.766369779);
    double ni = ne;
    double ui = cs;
    double Ti = 20;

    const auto& cells = mesh.getCells();

    // initializeAlpha(basisFunction);

    double nx = mesh.getNX();
    double nvx = mesh.getNVX();
    for (int j=0; j<nx; j++)
    {
        // std::cout << j << "\n";
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

        Matrix feq_tot(lMax,nvx);

        Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
        Vector u = integrator.integrate(fj, lMax, 1); //u tilde
        Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde

        Vector rho_i(lMax);
        rho_i[0] = ni;

        // if (j<32)
        // {
        //     ui = -cs;
        // }
        // else
        // {
        //     ui = cs;
        // }

        Matrix alpha(3,lMax);
        for (int m=0; m<3; m++)
        {
            for (int l=0; l<lMax; l++)
            {
                alpha(m,l) = alphaDomain(m+j*3,l);
            }
        }
        bool test = false;
        if (j==0)
        {
            test = false;
        }
        // std::cout << "Calculate Alphas" << "\n";
        // alpha = newtonSolver.solve(alpha, nu, rho, u, rt, dx, roots, weights, pow(10,-13), 100, basisFunction, quadratureOrder, lMax, test);
        // std::cout << "Alphas calculated" << "\n";
        for (int m=0; m<3; m++)
        {
            for (int l=0; l<lMax; l++)
            {
                alphaDomain(m+j*3,l) = alpha(m,l);
            }
        }

        for (int k=0; k<nvx; k++)
        {
            // std::cout << k << "\n";
            double vx = mesh.getVelocity(k);
            double fluxFactorMinus = (1.0+SpecialFunctions::sign(vx))/2.0;
            double fluxFactorPlus = (1.0-SpecialFunctions::sign(vx))/2.0;

            // calculate Ghost cells
            Vector fL = fitMaxwellian(basisFunction, Crec, cs, 10.0, vx, j);
            Vector fR = fitMaxwellian(basisFunction, Crec, -cs, 10.0, vx, j);

            // fit ion distribution function
            // Vector fi = fitMaxwellian(basisFunction, ni, ui, Ti, vx, j);

            // Matrix M_invS1(lMax,lMax*lMax);
            // Matrix M_invS2(lMax,lMax*lMax);

            Vector f_tilde(lMax);
            for (int l=0; l<lMax; l++)
            {
                uBefore(l,k+nx*nvx) = fL[l];
                uBefore(l,k+(nx+1)*nvx) = fR[l];
                // uBefore(l,k+nx*nvx) = uBefore(l,k+0*nvx); //Left BC
                // uBefore(l,k+(nx+1)*nvx) = uBefore(l,k+(nx-1)*nvx); //Right BC
                f_tilde[l] = uBefore(l,k+j*nvx);
            }
            Vector fCX = fitCX(basisFunction, ni, ui, Ti, rho, f_tilde, k, j);
            // std::cout << "\n" << k << ":\n";
            //calculate feq
            // Vector feq = fitMaxwellian(basisFunction, alpha, vx, j);
            for (int l=0; l<lMax; l++)
            {
                // feq_tot(l,k) = (nu/2.0)*GaussianQuadrature::integrate(basisFunction,l,alpha,vx,lMax,quadratureOrder,roots,weights);
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
                // uAfter(l,k+j*nvx)-=ne*uBefore(l,k+j*nvx)*sigma_iz; //This line for ionization
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
                uAfter(l,k+j*nvx)-=nu_cx*fCX[l]; //This line for CX
                // uAfter(l,k+j*nvx) -= ni*nu_cx*uBefore(l,k+j*nvx);
                // uAfter(l,k+j*nvx) += nu_cx*fCX[l];
                // uAfter(l,k+j*nvx)+=nu*(feq[l]-uBefore(l,k+j*nvx));
                // uAfter(l,k+j*nvx)+=nu*M_invDiag[l]*GaussianQuadrature::integrate(basisFunction,l,alpha,vx,lMax,quadratureOrder,roots,weights)/2.0;
                // uAfter(l,k+j*nvx)-=nu*uBefore(l,k+j*nvx);
                uAfter(l,k+j*nvx)*=dt;
                uAfter(l,k+j*nvx)+=uBefore(l,k+j*nvx);
                
                uAfter(l,k+j*nvx)*=timesFactor;
                uAfter(l,k+j*nvx)+=plusFactor*uPre(l,k+j*nvx); //Note that uPre != uBefore
            }
            // std::cout << fCX[0] << "\n";
            // std::cout << uAfter(0,k+j*nvx) << "\n";
        }
        // Vector rho_eq = integrator.integrate(feq_tot,lMax,0);
        // Vector u_eq = integrator.integrate(feq_tot,lMax,1);
        // Vector rt_eq = integrator.integrate(feq_tot,lMax,2);
        // for (int l=0; l<lMax; l++)
        // {
            // std::cout << "l = " << l << "\n";
            // std::cout << rho_eq[l] << "\n";
            // std::cout << pow(M_invDiag[l],-1)*rho[l] << "\n";
        // }
        // std::cout << (rho_eq[0]-rho[0]) << "\n";
    }
}

void Solver::advance(std::function<double(int,double)> basisFunction, int quadratureOrder)
{
    //First stage of solver
    // std::cout << "Solver stage 1" << "\n";
    advanceStage(uPre, uPost, 0.0, 1.0, basisFunction, quadratureOrder);
    //Second stage of solver
    // std::cout << "Solver stage 2" << "\n";
    advanceStage(uPost, uIntermediate, 3.0/4.0, 1.0/4.0, basisFunction, quadratureOrder);
    //Third stage of solver
    // std::cout << "Solver stage 3" << "\n";
    advanceStage(uIntermediate, uPost, 1.0/3.0, 2.0/3.0, basisFunction, quadratureOrder);

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
    // assert(uPre(l,j)==uPre(l,j)); //Might be better ways to ensure value is not NaN
    // if (uPre(l,j)!=uPre(l,j))
    // {
    //     std::cout << "NaN at l=" << l << "; k+j*nvx=" << j << "\n";
    // }
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
    Vector moments(4);
    double mass = 0;
    double momentum = 0;
    double energy = 0;
    double entropy = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    const auto& cells = mesh.getCells();

    double nvx = mesh.getNVX();
    for (int j=0; j<mesh.getNX(); j++)
    {
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
        for (int i=0; i<quadratureOrder; i++)
        {
            mass += weights[i]*SpecialFunctions::computeMoment(rho, basisFunction, lMax, roots[i])*dx/2.0;
            momentum += weights[i]*SpecialFunctions::computeMoment(u, basisFunction, lMax, roots[i])*dx/2.0;
            energy += weights[i]*SpecialFunctions::computeMoment(rt, basisFunction, lMax, roots[i])*dx/2.0;
            entropy += weights[i]*integrator.integrate(fj, lMax, basisFunction, roots[i])*dx/2.0;
        }
    }
    moments[0] = mass; //return rho
    moments[1] = momentum/mass; //return u
    moments[2] = energy; //return E
    moments[3] = entropy; //return S
    return moments;
}

//compute the value of your moment at normalized point x
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
Vector Solver::fitMaxwellian(std::function<double(int,double)> basisFunction, Matrix alpha, double vx, int j)
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
        double exponent = 0;
        for (int l=0; l<lMax; l++)
        {
            exponent += basisFunction(l,2.0*(x-xj)/dx)*(alpha(0,l)-vx*alpha(1,l)-pow(vx,2)*alpha(2,l));
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
        y[i] = exp(exponent);
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}

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
        double ui = meanVelocity_i*SpecialFunctions::sign(x-cells.back().vertices[1]/2.0);
        y[i] = density_i * f_n - density_n * SpecialFunctions::computeMaxwellian(density_i,ui,temperature_i,vx);

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

Vector Solver::getMoment(int j, int power)
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

    return integrator.integrate(fj, lMax, power); //rho tilde
}
