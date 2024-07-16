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
    double a;
    int lMax;
    Vector alpha;
    Matrix M_invS;
    Matrix M_invF1Minus;
    Matrix M_invF0Minus;
    Matrix M_invF1Plus;
    Matrix M_invF0Plus;
    Matrix uPre;
    Matrix uIntermediate;
    Matrix uPost;

    void advanceStage(Matrix& uPre, Matrix& uPost, double plusFactor, double timesFactor, std::function<double(int,double)> basisFunction, int quadratureOrder);

    //compute the reconstructed f(x,t)
    double getF(Matrix& uPre, std::function<double(int,double)> basisFunction, int lMax, int j, double x);

public:
    //constructor 
    Solver(const Mesh& mesh, double dt, double a, int lMax, Vector alpha);
    
    //deconstructor
    ~Solver();

    //create the mass, stiffness and flux matricies as used by the solver
    void createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder);

    //initialize using the Least Squares method
    void initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunctionX, std::function<double(double)> inputFunctionVX);

    //advance time step using 3rd order SSP RK
    void advance(std::function<double(int,double)> basisFunction, int quadratureOrder);

    //use minmod slope limiter
    void slopeLimiter();

    //get solution at (l,j)
    const double getSolution(int l, int j);

    //get error of solution
    const double getError(int tMax, std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction);

    //compute the mass, momentum and energy from f
    Vector getMoments(int quadratureOrder, std::function<double(int,double)> basisFunction);

    //fit Maxwellian
    Vector fitMaxwellian(std::function<double(int,double)> basisFunction, Vector alpha, double vx, int j);

};



#endif