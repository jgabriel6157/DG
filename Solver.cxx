#include <functional>
#include <cmath>
#include <iostream>
#include <cassert>
#include <chrono>

#include "Matrix.hxx"
#include "Vector.hxx"
#include "SpecialFunctions.hxx"
#include "GaussianQuadrature.hxx"
#include "Solver.hxx"
#include "Mesh.hxx"
#include "NewtonCotes.hxx"
#include "NewtonSolver.hxx"
#include "FunctionMapper.hxx"

//constructor
Solver::Solver(const Mesh& mesh, double dt, int lMax, std::function<double(int,double)> basisFunction, int quadratureOrder,
               bool ionization, bool cx, bool bgk, int bc) 
    : mesh(mesh), integrator(mesh), newtonSolver(mesh), dt(dt), lMax(lMax), basisFunction(basisFunction), quadratureOrder(quadratureOrder), 
      ionization(ionization), cx(cx), bgk(bgk), bc(bc), alphaDomain(3*mesh.getNX(), lMax),
      M_invDiag(lMax), M_invS(lMax,lMax), M_invF1Minus(lMax,lMax), M_invF0Minus(lMax,lMax), M_invF1Plus(lMax,lMax), M_invF0Plus(lMax,lMax),
      uPre(lMax,(mesh.getNX()+2)*mesh.getNVX()*mesh.getNVY()*mesh.getNVZ()), uIntermediate(lMax,(mesh.getNX()+2)*mesh.getNVX()*mesh.getNVY()*mesh.getNVZ()), 
      uPost(lMax,(mesh.getNX()+2)*mesh.getNVX()*mesh.getNVY()*mesh.getNVZ()) {}

//deconstructor
Solver::~Solver() {}

void Solver::createMatrices()
{
    auto basisFunctionDerivative = FunctionMapper::getDerivative(FunctionMapper::getFunctionName(basisFunction));
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

    for (int l=0; l<lMax; l++)
    {
        M_invDiag[l] = M_inv(l,l);
    }

    M_invS = M_inv*S;
    M_invF1Minus = M_inv*F1Minus;
    M_invF0Minus = M_inv*F0Minus;
    M_invF1Plus = M_inv*F1Plus;
    M_invF0Plus = M_inv*F0Plus;
}

//initialize using the Least Squares method
void Solver::initialize(std::function<double(double, double, double, double)> inputFunction)
{
    const auto& cells = mesh.getCells();

    double nx = mesh.getNX();
    double nvx = mesh.getNVX();
    double nvy = mesh.getNVY();
    double nvz = mesh.getNVZ();
    for (int j=0; j<nx; j++)
    {
        double dx = cells[j].dx;
        double leftVertex = cells[j].vertices[0];
        double xj = leftVertex+dx/2.0;
        
        for (int kx=0; kx<nvx; kx++)
        {
            double vx = mesh.getVelocityX(kx);
            for (int ky = 0; ky<nvy; ky++)
            {
                double vy = mesh.getVelocityY(ky);
                for (int kz = 0; kz<nvz; kz++)
                {
                    double vz = mesh.getVelocityZ(kz);
                    
                    Vector uInitialize(lMax);
                    
                    double x;
                    Vector y(10);
                    Matrix bigX(10,lMax);

                    for (int i=0; i<10; i++)
                    {
                        x = leftVertex+i*dx/9.0;
                        // x = leftVertex+(i+1)*dx/11.0;
                        y[i] = inputFunction(x,vx,vy,vz);
                        for (int l=0; l<lMax; l++)
                        {
                            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
                        }
                    }

                    uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;

                    for (int l=0; l<lMax; l++)
                    {
                        uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx) = uInitialize[l];
                    }
                }
            }
        }
    }
}

void Solver::initializeAlpha()
{
    // const auto& cells = mesh.getCells();

    // double nx = mesh.getNX();
    // double nvx = mesh.getNVX();
    // for (int j=0; j<nx; j++)
    // {
    //     double dx = cells[j].dx;
    //     double leftVertex = cells[j].vertices[0];
    //     double xj = leftVertex+dx/2.0;

    //     Matrix fj(lMax,nvx);
    //     for (int k=0; k<nvx; k++)
    //     {
    //         for (int l=0; l<lMax; l++)
    //         {
    //             fj(l,k) = uPre(l,k+j*nvx);
    //         }
    //     }

    //     Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
    //     Vector u = integrator.integrate(fj, lMax, 1); //u tilde
    //     Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde

    //     Vector uInitialize(3*lMax);
    //     Vector y(10*nvx);
    //     Matrix bigX(10*nvx,3*lMax);
        
    //     for (int k=0; k<nvx; k++)
    //     {
    //         double vx = mesh.getVelocity(k);
    //         double x;

    //         for (int i=0; i<10; i++)
    //         {
    //             x = leftVertex+(i+1)*dx/11.0; //Ignores end points which are more likely to be NaN in certain edge cases
    //             // x = leftVertex+(i+2)*dx/13.0; //Ignores end points which are more likely to be NaN in certain edge cases
    //             double density = SpecialFunctions::computeMoment(rho, basisFunction,lMax,2.0*(x-xj)/dx);
    //             double meanVelocity = SpecialFunctions::computeMoment(u, basisFunction,lMax,2.0*(x-xj)/dx)/density;
    //             double temperature = (SpecialFunctions::computeMoment(rt, basisFunction,lMax,2.0*(x-xj)/dx)-density*pow(meanVelocity,2))/density;

    //             double arg = SpecialFunctions::computeMaxwellian(density,meanVelocity,temperature,vx);
    //             // std::cout << "i = " << i << "\n";
    //             // std::cout << density << "\n";
    //             // std::cout << meanVelocity << "\n";
    //             // std::cout << temperature << "\n";
    //             // std::cout << arg << "\n";

    //             y[i+k*10] = log(arg);

    //             for (int m=0; m<3; m++)
    //             {
    //                 for (int l=0; l<lMax; l++)
    //                 {
    //                     if (m==0)
    //                     {
    //                         bigX(i+k*10,m+l*3) = basisFunction(l,2.0*(x-xj)/dx);
    //                     }
    //                     else
    //                     {
    //                         bigX(i+k*10,m+l*3) = -basisFunction(l,2.0*(x-xj)/dx)*pow(vx,m);
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;

    //     for (int m=0; m<3; m++)
    //     {
    //         for (int l=0; l<lMax; l++)
    //         {
    //             alphaDomain(m+j*3,l) = uInitialize[m+l*3];
    //             if (uInitialize[m+l*3]!=uInitialize[m+l*3])
    //             {
    //                 std::cout << j << "\n";
    //             }
    //             assert(uInitialize[m+l*3]==uInitialize[m+l*3]);
    //         }
    //     }
    // }
}

void Solver::advanceStage(Matrix& uBefore, Matrix& uAfter, double plusFactor, double timesFactor)
{
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    double nu = 1000.0;

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

    int nx = mesh.getNX();
    int nvx = mesh.getNVX();
    int nvy = mesh.getNVY();
    int nvz = mesh.getNVZ();
    #pragma omp parallel for schedule(dynamic)
    for (int j=0; j<nx; j++)
    {
        // std::cout << j << "\n";
        int leftNeighborIndex = cells[j].neighbors[0];
        int rightNeighborIndex = cells[j].neighbors[1];
        double dx = cells[j].dx;

        Matrix fj(lMax,nvx*nvy*nvz);
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        fj(l,kz+ky*nvz+kx*nvz*nvy) = uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                    }
                }
            }
        }

        // Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
        // Vector u = integrator.integrate(fj, lMax, 1); //u tilde
        // Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde
        Vector rho(lMax);
        Vector ux(lMax);
        Vector uy(lMax);
        Vector uz(lMax);
        Vector rt(lMax);
        Matrix moments = integrator.integrateMoments(fj, lMax);
        for (int l=0; l<lMax; l++)
        {
            rho[l] = moments(l,0);
            ux[l] = moments(l,1);
            uy[l] = moments(l,2);
            uz[l] = moments(l,3);
            rt[l] = moments(l,4);
        }
        // Vector rho = integrator.integrate3f(fj, lMax); //rho tilde
        // Vector ux = integrator.integrate3vxf(fj, lMax); //ux tilde
        // Vector uy = integrator.integrate3vyf(fj, lMax); //uy tilde
        // Vector uz = integrator.integrate3vzf(fj, lMax); //uz tilde
        // Vector rt = integrator.integrate3v2f(fj, lMax); //rt tilde

        Vector rho_i(lMax);
        rho_i[0] = ni;

        Matrix alpha(3,lMax);
        if (bgk)
        {
            // for (int m=0; m<3; m++)
            // {
            //     for (int l=0; l<lMax; l++)
            //     {
            //         alpha(m,l) = alphaDomain(m+j*3,l);
            //     }
            // }
            // bool test = false;
            // if (j==0)
            // {
            //     test = false;
            // }
            // // std::cout << "Calculate Alphas" << "\n";
            // alpha = newtonSolver.solve(alpha, nu, rho, u, rt, dx, roots, weights, pow(10,-13), 100, basisFunction, quadratureOrder, lMax, test);
            // // std::cout << "Alphas calculated" << "\n";
            // for (int m=0; m<3; m++)
            // {
            //     for (int l=0; l<lMax; l++)
            //     {
            //         alphaDomain(m+j*3,l) = alpha(m,l);
            //     }
            // }
        }

        for (int kx=0; kx<nvx; kx++)
        {
            double vx = mesh.getVelocityX(kx);
            double fluxFactorMinus = (1.0+SpecialFunctions::sign(vx))/2.0;
            double fluxFactorPlus = (1.0-SpecialFunctions::sign(vx))/2.0;

            for(int ky=0; ky<nvy; ky++)
            {
                double vy = mesh.getVelocityY(ky);
                for(int kz=0; kz<nvz; kz++)
                {
                    double vz = mesh.getVelocityZ(kz);

                    // calculate Ghost cells
                    // Vector fL = fitMaxwellian(Crec, cs, 2.0, vx, j);
                    // Vector fR = fitMaxwellian(Crec, -cs, 2.0, vx, j);
                    Vector fL = fitMaxwellian3(Crec, cs, 2.0, vx, vy, vz, j);
                    Vector fR = fitMaxwellian3(Crec, -cs, 2.0, vx, vy, vz, j);

                    Vector f_tilde(lMax);
                    for (int l=0; l<lMax; l++)
                    {
                        if (bc==1)
                        {
                            uBefore(l,kz+ky*nvz+kx*nvz*nvy+nx*nvz*nvy*nvx) = fL[l];
                            uBefore(l,kz+ky*nvz+kx*nvz*nvy+(nx+1)*nvz*nvy*nvx) = fR[l];
                        }
                        else if (bc==2)
                        {
                            uBefore(l,kz+ky*nvz+kx*nvz*nvy+nx*nvz*nvy*nvx) = uBefore(l,kz+ky*nvz+kx*nvz*nvy+0*nvz*nvy*nvx); //Left BC
                            uBefore(l,kz+ky*nvz+kx*nvz*nvy+(nx+1)*nvz*nvy*nvx) = uBefore(l,kz+ky*nvz+kx*nvz*nvy+(nx-1)*nvz*nvy*nvx); //Right BC
                        }
                        f_tilde[l] = uBefore(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                    }

                    Vector fCX(lMax);
                    if (cx)
                    {
                        fCX = fitCX(ni, ui, Ti, rho, f_tilde, kx, ky, kz, j);
                    }

                    for (int l=0; l<lMax; l++)
                    {
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)=0;
                        for (int i=0; i<lMax; i++)
                        {
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)+=M_invS(l,i)*uBefore(i,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)-=fluxFactorMinus*M_invF1Minus(l,i)*uBefore(i,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)+=fluxFactorMinus*M_invF0Minus(l,i)*uBefore(i,kz+ky*nvz+kx*nvz*nvy+leftNeighborIndex*nvz*nvy*nvx);
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)-=fluxFactorPlus*M_invF1Plus(l,i)*uBefore(i,kz+ky*nvz+kx*nvz*nvy+rightNeighborIndex*nvz*nvy*nvx);
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)+=fluxFactorPlus*M_invF0Plus(l,i)*uBefore(i,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                        }
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)*=vx;
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)/=dx;
                        if (ionization)
                        {
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)-=ne*uBefore(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)*sigma_iz; //This line for ionization
                        }
                        if (cx)
                        {
                            uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)-=nu_cx*fCX[l]; //This line for CX
                        }
                        if (bgk)
                        {
                            // uAfter(l,k+j*nvx)+=nu*M_invDiag[l]*GaussianQuadrature::integrate(basisFunction,l,alpha,vx,lMax,quadratureOrder,roots,weights)/2.0; //BGK
                            // uAfter(l,k+j*nvx)-=nu*uBefore(l,k+j*nvx); //BGK
                        }
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)*=dt;
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)+=uBefore(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                        
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)*=timesFactor;
                        uAfter(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx)+=plusFactor*uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx); //Note that uPre != uBefore
                    }
                }
            }
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

void Solver::slopeLimiter() // Needs updating
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
    assert(uPre(l,j)==uPre(l,j)); //Might be better ways to ensure value is not NaN
    if (uPre(l,j)!=uPre(l,j))
    {
        std::cout << "NaN at l=" << l << "; k+j*nvx=" << j << "\n";
    }
    return uPre(l,j);
}

Vector Solver::getMoments()
{
    Vector moments(6);
    double mass = 0;
    double momentumX = 0;
    double momentumY = 0;
    double momentumZ = 0;
    double energy = 0;
    double entropy = 0;
    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);
    const auto& cells = mesh.getCells();

    double nvx = mesh.getNVX();
    double nvy = mesh.getNVY();
    double nvz = mesh.getNVZ();
    for (int j=0; j<mesh.getNX(); j++)
    {
        double dx = cells[j].dx;
        Matrix fj(lMax,nvz*nvy*nvx);
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        fj(l,kz+ky*nvz+kx*nvz*nvy) = uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                    }
                }
            }
        }

        // Vector rho = integrator.integrate(fj, lMax, 0); //rho tilde
        // Vector u = integrator.integrate(fj, lMax, 1); //u tilde
        // Vector rt = integrator.integrate(fj, lMax, 2); //rt tilde
        Vector rho(lMax);
        Vector ux(lMax);
        Vector uy(lMax);
        Vector uz(lMax);
        Vector rt(lMax);
        Matrix moments = integrator.integrateMoments(fj, lMax);
        for (int l=0; l<lMax; l++)
        {
            rho[l] = moments(l,0);
            ux[l] = moments(l,1);
            uy[l] = moments(l,2);
            uz[l] = moments(l,3);
            rt[l] = moments(l,4);
        }
        for (int i=0; i<quadratureOrder; i++)
        {
            mass += weights[i]*SpecialFunctions::computeMoment(rho, basisFunction, lMax, roots[i])*dx/2.0;
            momentumX += weights[i]*SpecialFunctions::computeMoment(ux, basisFunction, lMax, roots[i])*dx/2.0;
            momentumY += weights[i]*SpecialFunctions::computeMoment(uy, basisFunction, lMax, roots[i])*dx/2.0;
            momentumZ += weights[i]*SpecialFunctions::computeMoment(uz, basisFunction, lMax, roots[i])*dx/2.0;
            energy += weights[i]*SpecialFunctions::computeMoment(rt, basisFunction, lMax, roots[i])*dx/2.0;
            entropy += weights[i]*integrator.integrate(fj, lMax, basisFunction, roots[i])*dx/2.0;
        }
    }
    moments[0] = mass; //return rho
    moments[1] = momentumX/mass; //return ux
    moments[2] = momentumY/mass; //return uy
    moments[3] = momentumZ/mass; //return uz
    moments[4] = energy; //return E
    moments[5] = entropy; //return S
    return moments;
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

Vector Solver::getRho(int j)
{
    int nvx = mesh.getNVX();
    int nvy = mesh.getNVY();
    int nvz = mesh.getNVZ();
    Matrix fj(lMax,nvz*nvy*nvx);
    for (int kx=0; kx<nvx; kx++)
    {
        for (int ky=0; ky<nvy; ky++)
        {
            for (int kz=0; kz<nvz; kz++)
            {
                for (int l=0; l<lMax; l++)
                {
                    fj(l,kz+ky*nvz+kx*nvz*nvy) = uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                }
            }
        }
    }

    return integrator.integrate3f(fj,lMax);
}

Vector Solver::getU(int j, int coord)
{
    assert(coord>=0);
    assert(coord<3);
    int nvx = mesh.getNVX();
    int nvy = mesh.getNVY();
    int nvz = mesh.getNVZ();
    Matrix fj(lMax,nvz*nvy*nvx);
    for (int kx=0; kx<nvx; kx++)
    {
        for (int ky=0; ky<nvy; ky++)
        {
            for (int kz=0; kz<nvz; kz++)
            {
                for (int l=0; l<lMax; l++)
                {
                    fj(l,kz+ky*nvz+kx*nvz*nvy) = uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                }
            }
        }
    }

    if (coord==0)
    {
        return integrator.integrate3vxf(fj,lMax);
    }
    else if (coord==1)
    {
        return integrator.integrate3vyf(fj,lMax);
    }
    else
    {
        return integrator.integrate3vzf(fj,lMax);
    }
}

Vector Solver::getE(int j)
{
    int nvx = mesh.getNVX();
    int nvy = mesh.getNVY();
    int nvz = mesh.getNVZ();
    Matrix fj(lMax,nvz*nvy*nvx);
    for (int kx=0; kx<nvx; kx++)
    {
        for (int ky=0; ky<nvy; ky++)
        {
            for (int kz=0; kz<nvz; kz++)
            {
                for (int l=0; l<lMax; l++)
                {
                    fj(l,kz+ky*nvz+kx*nvz*nvy) = uPre(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx);
                }
            }
        }
    }

    return integrator.integrate3v2f(fj,lMax);
}

//initialize using the Least Squares method
Vector Solver::fitMaxwellian(Matrix alpha, double vx, int j)
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

Vector Solver::fitMaxwellian(Vector rho, Vector u, Vector rt, double vx, int j)
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
        double density = SpecialFunctions::computeMoment(rho, basisFunction,lMax,2.0*(x-xj)/dx);
        double meanVelocity = SpecialFunctions::computeMoment(u, basisFunction,lMax,2.0*(x-xj)/dx)/density;
        double temperature = (SpecialFunctions::computeMoment(rt, basisFunction,lMax,2.0*(x-xj)/dx)-density*pow(meanVelocity,2))/density;
        y[i] = SpecialFunctions::computeMaxwellian(density,meanVelocity,temperature,vx);
        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;
}

Vector Solver::fitMaxwellian(double density, double meanVelocity, double temperature, double vx, int j)
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

Vector Solver::fitMaxwellian3(double density, double meanVelocity, double temperature, double vx, double vy, double vz, int j)
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
        y[i] = SpecialFunctions::computeMaxwellian3(density,meanVelocity,temperature,vx,vy,vz);
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

Vector Solver::fitCX(double density_i, double meanVelocity_i, double temperature_i, Vector rho_n, Vector f_tilde, int kx, int ky, int kz, int j)
{
    const auto& cells = mesh.getCells();

    double dx = cells[j].dx;
    double leftVertex = cells[j].vertices[0];
    double xj = leftVertex+dx/2.0;

    double vx = mesh.getVelocityX(kx);
    double vy = mesh.getVelocityY(ky);
    double vz = mesh.getVelocityZ(kz);
    
    Vector uInitialize(lMax);

    int res = 10;
    double x;
    Vector y(res);
    Matrix bigX(res,lMax);

    for (int i=0; i<res; i++)
    {
        x = leftVertex+i*dx/(res-1.0);
        double density_n = SpecialFunctions::computeMoment(rho_n, basisFunction, lMax, 2.0*(x-xj)/dx);
        double f_n = SpecialFunctions::computeMoment(f_tilde, basisFunction, lMax, 2.0*(x-xj)/dx);
        double ui = meanVelocity_i*SpecialFunctions::sign(x-cells.back().vertices[1]/2.0);
        y[i] = density_i * f_n - density_n * SpecialFunctions::computeMaxwellian3(density_i,ui,temperature_i,vx,vy,vz);

        for (int l=0; l<lMax; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    return uInitialize = (bigX.Transpose()*bigX).CalculateInverse()*bigX.Transpose()*y;   
}
