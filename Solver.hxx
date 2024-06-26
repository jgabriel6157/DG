#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF

#include <functional>
#include "Matrix.hxx"
#include "Vector.hxx"
#include "Mesh.hxx"

class Solver
{
private:
    Mesh mesh;
    // double dx;
    double dt;
    double a;
    // int jMax;
    int lMax;
    double alpha;
    Matrix M_invS;
    Matrix M_invF1;
    Matrix M_invF2;
    Matrix M_invF3;
    Matrix M_invF4;
    Matrix uPre;
    Matrix uIntermediate;
    Matrix uPost;

    void advanceStage(Matrix& uPre, Matrix& uPost, double plusFactor, double timesFactor);

    //compute the reconstructed f(x,t)
    double getF(Matrix& uPre, std::function<double(int,double)> basisFunction, int lMax, int j, double x);

public:
    //constructor 
    // Solver(double dx, double dt, double a, int jMax, int lMax, double alpha);
    Solver(const Mesh& mesh, double dt, double a, int lMax, double alpha);
    
    //deconstructor
    ~Solver();

    //create the mass, stiffness and flux matricies as used by the solver
    void createMatrices(std::function<double(int,double)> basisFunction, std::function<double(int,double)> basisFunctionDerivative, int quadratureOrder);

    //initialize using the Least Squares method
    void initialize(std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction);

    //advance time step using 3rd order SSP RK
    void advance();

    //use minmod slope limiter
    void slopeLimiter();

    //get solution at (l,j)
    const double getSolution(int l, int j);

    //get error of solution
    const double getError(int tMax, std::function<double(int,double)> basisFunction, std::function<double(double)> inputFunction);

    //compute the total mass from f
    double getMass(int quadratureOrder, std::function<double(int,double)> basisFunction);

};



#endif