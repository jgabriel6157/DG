#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF
#include <functional>
#include "Matrix.hxx"
#include "Vector.hxx"

class Solver
{
private:
    double dx;
    double dt;
    double a;
    int jMax;
    int lMax;
    Matrix M_invS;
    Matrix M_invF1;
    Matrix M_invF2;
    Matrix M_invF3;
    Matrix M_invF4;
    Matrix uPre;
    Matrix uIntermediate;
    Matrix uPost;

    void advanceStage(Matrix& uPre, Matrix& uPost, double plusFactor, double timesFactor);
public:
    //constructor 
    Solver(double dx, double dt, double a, int jMax, int lMax, Matrix& M_invS, Matrix& M_invF1, Matrix& M_invF2, Matrix& M_invF3, Matrix& M_invF4);
    
    //deconstructor
    ~Solver();

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

};



#endif