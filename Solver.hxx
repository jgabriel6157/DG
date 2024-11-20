#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF

#include <functional>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "Mesh.hxx"
#include "NewtonCotes.hxx"
#include "NewtonSolver.hxx"

class Solver
{
private:
    Mesh mesh;
    NewtonCotes integrator;
    NewtonSolver newtonSolver;
    double dt;
    int lMax;
    Matrix alphaDomain;
    Vector M_invDiag;
    Matrix M_invS;
    Matrix M_invF1Minus;
    Matrix M_invF0Minus;
    Matrix M_invF1Plus;
    Matrix M_invF0Plus;
    Matrix uPre;
    Matrix uIntermediate;
    Matrix uPost;

    void advanceStage(Matrix& uPre, Matrix& uPost, double plusFactor, double timesFactor, std::function<double(int,double)> basisFunction, int quadratureOrder);

public:
    //constructor 
    Solver(const Mesh& mesh, double dt, int lMax);
    
    //deconstructor
    ~Solver();

    //create the mass, stiffness and flux matricies as used by the solver
    void createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder);

    //initialize using the Least Squares method
    void initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunctionX, std::function<double(double)> inputFunctionVX);

    void initializeAlpha(std::function<double(int,double)> basisFunction);

    //advance time step using 3rd order SSP RK
    void advance(std::function<double(int,double)> basisFunction, int quadratureOrder);

    //use minmod slope limiter
    void slopeLimiter();

    //get solution at (l,j)
    const double getSolution(int l, int j);

    //compute the mass, momentum and energy from f
    Vector getMoments(int quadratureOrder, std::function<double(int,double)> basisFunction);

    //compute the density from f
    Vector getMoment(int j, int power);

    //fit Maxwellian
    Vector fitMaxwellian(std::function<double(int,double)> basisFunction, Vector rho, Vector u, Vector rt, double vx, int j);

    Vector fitMaxwellian(std::function<double(int,double)> basisFunction, double density, double meanVelocity, double temperature, double vx, int j);

    //compute the density from f
    Vector getDensity(int j);

    Vector fitCX(std::function<double(int,double)> basisFunction, double density, double meanVelocity, double temperature, Vector rho_n, Vector f_tilde, int k, int j);
    //fit Maxwellian
    Vector fitMaxwellian(std::function<double(int,double)> basisFunction, Matrix alpha, double vx, int j);

};



#endif