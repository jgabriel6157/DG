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
    std::function<double(int,double)> basisFunction;
    int quadratureOrder;
    bool ionization;
    bool cx;
    bool bgk;
    int bc;
    Matrix alphaDomain;
    Vector M_invDiag;
    Matrix M_invS;
    Matrix M_invT;
    Matrix M_invF1Minus;
    Matrix M_invF0Minus;
    Matrix M_invF1Plus;
    Matrix M_invF0Plus;
    Matrix uPre;
    Matrix uIntermediate;
    Matrix uPost;
    Matrix fSource;
    Matrix fi;

    void advanceStage(Matrix& uPre, Matrix& uPost, double plusFactor, double timesFactor);

public:
    //constructor 
    Solver(const Mesh& mesh, double dt, int lMax, std::function<double(int,double)> basisFunction, int quadratureOrder,
           bool ionization, bool cx, bool bgk, int bc);
    
    //deconstructor
    ~Solver();

    //create the mass, stiffness and flux matricies as used by the solver
    void createMatrices();

    //initialize using the Least Squares method
    void initialize(std::function<double(double, double, double, double)> inputFunction);

    void initializeSource();

    void initializeIons();

    void initializeAlpha();

    //advance time step using 3rd order SSP RK
    void advance();

    //use minmod slope limiter
    void slopeLimiter();

    //get solution at (l,j)
    const double getSolution(int l, int j);

    //compute the mass, momentum and energy from f
    Vector getMoments();

    //compute the density from f
    Vector getMoment(int j, int power);

    Vector getRho(int j);

    Vector getU(int j, int coord);

    Vector getE(int j);

    //fit Maxwellian
    Vector fitMaxwellian(Vector rho, Vector u, Vector rt, double vx, int j);

    Vector fitMaxwellian(double density, double meanVelocity, double temperature, double vx, int j);

    Vector fitMaxwellian3(double density, double meanVelocityX, double meanVelocityY, double meanVelocityZ, double temperature, double vx, double vy, double vz);

    //compute the density from f
    Vector getDensity(int j);

    Vector fitCX(double density_i, double meanVelocity_i, double temperature_i, Vector rho_n, Vector f_tilde, int kx, int ky, int kz, int j);

    Vector fitCX(double density_i, double meanVelocity_i, double temperature_i, Vector f_tilde, Vector fnCXavg, double vx, double vy, double vz, int j);
    
    //fit Maxwellian
    Vector fitMaxwellian(Matrix alpha, double vx, int j);

};



#endif