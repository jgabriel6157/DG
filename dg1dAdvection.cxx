#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <chrono>
#include <functional>
#include <map>
#include <vector>

class SpecialFunctions
{
public:
    //Constructor
    SpecialFunctions() {}

    //Function to calculate Legendre polynomial of order n at point x
    static double legendre(int n, double x)
    {
        assert (n>=0);
        switch (n)
        {
            case 0:
                return 1;
                break;
            case 1:
                return x;
                break;
            default:
                return ((2.0*n-1.0)*x*legendre(n-1,x)-(n-1)*legendre(n-2,x))/n;
                break;
        }
    }

    //Function to calculate orthonormal Legendre polynomial of order n at point x
    static double legendreOrthonormal(int n, double x)
    {
        assert (n>=0);
        return sqrt((2.0*n+1.0)/2.0)*legendre(n,x);
    }

    //Function to calculate derivate of the Legendre polynomial of order n at point x
    static double legendreDerivative(int n, double x)
    {
        assert (n>=0);
        switch (n)
        {
            case 0:
                return 0;
                break;
            case 1:
                return 1;
                break;
            default:
                return ((2.0*n-1.0)*(x*legendreDerivative(n-1,x)+legendre(n-1,x))-(n-1)*legendreDerivative(n-2,x))/n;
        }
    }

    //Function to calculate derivative of orthonormal Legendre polynomial of order n at point x
    static double legendreOrthonormalDerivative(int n, double x)
    {
        assert(n>=0);
        return sqrt((2.0*n+1.0)/2.0)*legendreDerivative(n,x);
    }

    //Function to calculate quadratic basis function of 'order' n at point x
    static double quadratic(int n, double x)
    {
        assert(n>=0);
        assert(n<3);
        switch (n)
        {
        case 0:
            return -x*(1.0-x)/2.0;
            break;
        case 1:
            return (1.0-x)*(1.0+x);
            break;
        case 2:
            return x*(1+x)/2;
            break;
        default:
            return 0;
            break;
        }
    }

    //Function to calculate derivative of the quadratic basis function of 'order' n at point x
    static double quadraticDerivative(int n, double x)
    {
        assert(n>=0);
        assert(n<3);
        switch (n)
        {
        case 0:
            return x-1.0/2.0;
            break;
        case 1:
            return -2.0*x;
            break;
        case 2:
            return x+1.0/2.0;
            break;
        default:
            return 0;
            break;
        }
    }

    //Function to calculate linear basis function of 'order' n at point x
    static double linear(int n, double x)
    {
        assert(n>=0);
        assert(n<2);
        switch (n)
        {
        case 0:
            return (1.0-x)/2.0;
            break;
        case 1:
            return (1.0+x)/2.0;
            break;
        default:
            return 0;
            break;
        }
    }

    //Function to calculate derivative of linear basis function of 'order' n at point x
    static double linearDerivative(int n, double x)
    {
        assert(n>=0);
        assert(n<2);
        switch (n)
        {
        case 0:
            return -1.0/2.0;
            break;
        case 1:
            return 1.0/2.0;
            break;
        default:
            return 0;
            break;
        }
    }

    //Square pulse (top hat) function centered at pi with length 2
    static double topHat(double x)
    {
        if ((x<M_PI-1.0)||(x>M_PI+1.0))
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }

    //Gaussian pulse centered at pi
    static double gaussianPulse(double x)
    {
        return exp(-1.0*pow(x-1.0*M_PI,2.0));
    }

    //Find roots of Legendre polynomial of order n using Newton Raphson method 
    static std::vector<double> legendreRoots(int n)
    {
        std::vector<double> roots;

        for (int i=0; i<n; i++)
        {
            double x0 = (1.0-1.0/(8.0*pow(n,2.0))+1.0/(8.0*pow(n,3.0)))*cos(M_PI*(4.0*(i+1.0)-1.0)/(4.0*n+2.0)); //intial guess
            double x = newtonRaphson(n,x0);
            roots.push_back(x);
        }

        return roots;
    }

private:
    //Newton Raphson method to find root of Legendre polynomial of root n with initial guess x0
    static double newtonRaphson(int n, double x0)
    {
        double x_root = x0;

        while (fabs(legendre(n,x_root))>1.0e-10)
        {
            x_root -= legendre(n,x_root)/legendreDerivative(n,x_root);
        }

        return x_root;
    }

};

class GaussianQuadrature
{
public:
    //Perform Gaussian quadrature integration for two specified functions of order n_func1, n_func2 to a certain quadrature order
    static double integrate(std::function<double(int, double)> func1, int n_func1, std::function<double(int, double)> func2, int n_func2, int quadratureOrder,
                            std::vector<double> roots, std::vector<double> weights)
    {
        double y = 0;

        for (int i=0; i<quadratureOrder; i++)
        {
            y += weights[i]*func1(n_func1,roots[i])*func2(n_func2,roots[i]);
        }

        return y;
    }

    //Calculate weights for Gaussian quadrature
    static std::vector<double>calculateWeights(int quadratureOrder, std::vector<double> roots)
    {
        std::vector<double> weights;

        double root;
        double w;
        for (int i=0; i<quadratureOrder; i++)
        {
            
            root = roots[i];
            w = 2.0/((1.0-pow(root,2.0))*pow(SpecialFunctions::legendreDerivative(quadratureOrder,root),2.0));
            
            weights.push_back(w);
        }

        return weights;
    }
};

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

double sign(double x);
double minmod(double a, double b, double c);
double min(double a, double b);

void LeastSquares(double* uInitialize, double xj, double dx, std::function<double(int,double)> basisFunction, int order, std::function<double(double)> inputFunction);
void MatrixMultiply(double** M1, double** M2, double** M, int rowM1, int middleSize, int columnM2);
void MatrixMultiply(double** M, double* A1, double* A, int rowM1, int middleSize);
void MatrixInverse(double** M, double** M_inverse, int order);

double getError(double** u, int jMax, std::function<double(int,double)> basisFunction, int lMax, double dx, int tMax, double dt, std::function<double(double)> inputFunction);

void firstOrderEulerPlusTimes(double** uPre, double** uPost, double** M_invS, double** M_invF1, double** M_invF2, double** M_invF3, double** M_invF4,
                              double** uPlus, double plusFactor, double timesFactor, double dx, double dt, double a, int jMax, int lMax);
void slopeLimiter(double** uPost, int jMax, int lMax);

int main(int argc, char* argv[])
{
    std::string argName[11] = {"a","jMax","lMax","tMax","quadratureOrder","length","dt","basis","test","alpha","input"};
    std::string argString[11];
    readFile("input.txt",argName,argString,11);

    double a = assignDouble(argString[0]);
    int jMax = assignInt(argString[1]);
    int lMax = assignInt(argString[2]);
    int tMax = assignInt(argString[3]);
    int quadratureOrder = assignInt(argString[4]);
    double length = assignDouble(argString[5]);
    double dt = assignDouble(argString[6]);
    std::string basis = argString[7];
    bool test = assignBool(argString[8]);
    double alpha = assignDouble(argString[9]);
    std::string input = argString[10];
    
    double dx = length/jMax;
    double fluxFactorPlus = (1.0+sign(a)*(1.0-alpha))/2.0;
    double fluxFactorMinus = (1.0-sign(a)*(1.0-alpha))/2.0;

    double** M = new double* [lMax];
    double** M_inv = new double* [lMax];
    double** S = new double* [lMax];
    double** F1 = new double* [lMax];
    double** F2 = new double* [lMax];
    double** F3 = new double* [lMax];
    double** F4 = new double* [lMax];
    double** M_invS = new double* [lMax];
    double** M_invF1 = new double* [lMax];
    double** M_invF2 = new double* [lMax];
    double** M_invF3 = new double* [lMax];
    double** M_invF4 = new double* [lMax];
    double** uPre = new double* [lMax];
    double** uPost = new double* [lMax];
    double** uIntermediate = new double* [lMax];
    for (int l=0; l<lMax; l++)
    {
        M[l] = new double [lMax];
        M_inv[l] = new double [lMax];
        S[l] = new double [lMax];
        F1[l] = new double [lMax];
        F2[l] = new double [lMax];
        F3[l] = new double [lMax];
        F4[l] = new double [lMax];
        M_invS[l] = new double [lMax];
        M_invF1[l] = new double [lMax];
        M_invF2[l] = new double [lMax];
        M_invF3[l] = new double [lMax];
        M_invF4[l] = new double [lMax];
        uPre[l] = new double [jMax];
        uPost[l] = new double [jMax];
        uIntermediate[l] = new double [jMax];
    }
    double* uInitialize = new double [lMax];
    double xj;
    
    std::ofstream write_output("Output.csv");
    assert(write_output.is_open());

    //Define map for basis string to basis function
    std::map<std::string, std::function<double(int, double)>> basisFunctionMap = 
    {
        {"legendre", SpecialFunctions::legendre},
        {"legendreDerivative", SpecialFunctions::legendreDerivative},
        {"legendreOrthonormal", SpecialFunctions::legendreOrthonormal},
        {"legendreOrthonormalDerivative", SpecialFunctions::legendreOrthonormalDerivative},
        {"quadratic", SpecialFunctions::quadratic},
        {"quadraticDerivative", SpecialFunctions::quadraticDerivative},
        {"linear", SpecialFunctions::linear},
        {"linearDerivative", SpecialFunctions::linearDerivative}
    };
    auto basisFunctionIterator = basisFunctionMap.find(basis);
    auto basisFunctionIteratorDerivative = basisFunctionMap.find(basis+"Derivative");
    if (basisFunctionIterator == basisFunctionMap.end())
    {
        std::cerr << "Unknown function: " << basis << "\n";
        return 1;
    }
    std::function<double(int, double)> basisFunction = basisFunctionIterator->second; //Grab basis function
    std::function<double(int, double)> basisFunctionDerivative = basisFunctionIteratorDerivative->second; //Grab basis function derivative

    //Define map for input string to input function
    std::map<std::string, std::function<double(double)>> inputFunctionMap = 
    {
        {"sin", sin},
        {"topHat", SpecialFunctions::topHat},
        {"pulse", SpecialFunctions::gaussianPulse}
    };
    auto inputFunctionIterator = inputFunctionMap.find(input);
    if (inputFunctionIterator == inputFunctionMap.end())
    {
        std::cerr << "Unknown function: " << input << "\n";
        return 1;
    }
    std::function<double(double)> inputFunction = inputFunctionIterator->second; //Grab input function

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> roots = SpecialFunctions::legendreRoots(quadratureOrder);
    std::vector<double> weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M[i][j] = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,quadratureOrder,roots,weights)/2;
            S[i][j] = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,quadratureOrder,roots,weights);
            F1[i][j] = fluxFactorPlus*(basisFunction(i,1))*(basisFunction(j,1));
            F2[i][j] = fluxFactorPlus*(basisFunction(i,-1))*(basisFunction(j,1));
            F3[i][j] = fluxFactorMinus*(basisFunction(i,1))*(basisFunction(j,-1));
            F4[i][j] = fluxFactorMinus*(basisFunction(i,-1))*(basisFunction(j,-1));
            if (fabs(M[i][j]) < 1e-10)
            {
                M[i][j] = 0;
            }
            if (fabs(S[i][j]) < 1e-10)
            {
                S[i][j] = 0;
            }
        }
    }
    
    MatrixInverse(M,M_inv,lMax);
    for (int l=0; l<lMax; l++)
    {
        delete[] M[l];
    }
    delete[] M;

    MatrixMultiply(M_inv,S,M_invS,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F1,M_invF1,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F2,M_invF2,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F3,M_invF3,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F4,M_invF4,lMax,lMax,lMax);

    for (int l=0; l<lMax; l++)
    {
        delete[] M_inv[l];
        delete[] S[l];
        delete[] F1[l];
        delete[] F2[l];
        delete[] F3[l];
        delete[] F4[l];
    }
    delete[] M_inv;
    delete[] S;
    delete[] F1;
    delete[] F2;
    delete[] F3;
    delete[] F4;

    for (int j=0; j<jMax; j++)
    {
        xj = j*dx+dx/2.0;
        LeastSquares(uInitialize, xj, dx, basisFunction, lMax, inputFunction);
        for (int l=0; l<lMax; l++)
        {
            uPre[l][j] = uInitialize[l];
        }
    }
    delete[] uInitialize;

    for (int t=0; t<tMax; t++)
    {
        firstOrderEulerPlusTimes(uPre,uPost,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,0.0,1.0,dx,dt,a,jMax,lMax);

        firstOrderEulerPlusTimes(uPost,uIntermediate,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,3.0/4.0,1.0/4.0,dx,dt,a,jMax,lMax);

        firstOrderEulerPlusTimes(uIntermediate,uPost,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,1.0/3.0,2.0/3.0,dx,dt,a,jMax,lMax);

        // slopeLimiter(uPost,jMax,lMax);

        for (int j=0; j<jMax; j++)
        {
            for (int l=0; l<lMax; l++)
            {
                uPre[l][j] = uPost[l][j];
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(stop-start);
    std::cout << duration.count() << " ms" << "\n";

    if (test)
    {
        std::cout << getError(uPre, jMax, basisFunction, lMax, dx, tMax, dt, inputFunction) << "\n";
    }

    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            write_output << uPre[l][j] << "\n";
        }
    }

    write_output.close();

    for (int l=0; l<lMax; l++)
    {
        delete[] M_invS[l];
        delete[] M_invF1[l];
        delete[] M_invF2[l];
        delete[] M_invF3[l];
        delete[] M_invF4[l];
        delete[] uPre[l];
        delete[] uPost[l];
    }
    delete[] M_invS;
    delete[] M_invF1;
    delete[] M_invF2;
    delete[] M_invF3;
    delete[] M_invF4;
    delete[] uPre;
    delete[] uPost;

    return 0;
}

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables)
{
    std::ifstream inputFile(filename);
    assert(inputFile.is_open());

    std::string word,value,dum1,valueCheck;

    while (inputFile >> word)
    {
        for (int i=0;i<numberOfVariables;i++)
        {
            if (word==argName[i])
            {
                inputFile >> word >> value;
                argString[i] = value;
            }
        }
    } 

    inputFile.close();
}

int assignInt(std::string varString)
{
    int number = 1;
    double const pi = M_PI;
    std::string valueCheck = "*";
    std::string value = varString;
    std::string dum1;
    while (valueCheck.find("*") != std::string::npos)
    {
        valueCheck = value;
        dum1 = value.substr(0,value.find("*"));
        if (dum1 == "pi")
        {
            number*=pi;
        }
        else
        {
            number*=stoi(dum1);
        }
        value = value.substr(value.find("*")+1);
    }
    return number;
}

double assignDouble(std::string varString)
{
    double number = 1.0;
    double const pi = M_PI;
    std::string valueCheck = "*";
    std::string value = varString;
    std::string dum1;
    while (valueCheck.find("*") != std::string::npos)
    {
        valueCheck = value;
        dum1 = value.substr(0,value.find("*"));
        if (dum1 == "pi")
        {
            number*=pi;
        }
        else
        {
            number*=stod(dum1);
        }
        value = value.substr(value.find("*")+1);
    }
    return number;
}

bool assignBool(std::string varString)
{
    if (varString=="true")
    {
        return true;
    }
    else if (varString=="false")
    {
        return false;
    }
    else
    {
        std::cout << "Invalid input for bool, returning false" << "\n";
        return false;
    }
}

void MatrixMultiply(double** M1, double** M2, double** M, int rowM1, int middleSize, int columnM2)
{
    for (int i=0; i<rowM1; i++)
    {
        for (int j=0; j<columnM2; j++)
        {
            M[i][j] = 0;        
            for (int k=0; k<middleSize; k++)
            {
                M[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

void MatrixMultiply(double** M, double* A1, double* A, int rowM1, int middleSize)
{
    for (int i=0; i<rowM1; i++)
    {
        A[i] = 0;
        for (int j=0; j<middleSize; j++)
        {
            A[i] += M[i][j]*A1[j];
        }
    }
}

void MatrixInverse(double** M, double** M_inverse, int order)
{
    //Take inverse of Matrix M through Gaussian Jordan elimination
    double ratio;
    double M_inverse_big[order][order*2];
    
    // Augmenting Identity Matrix of Order n 
    for(int i=0; i<order; i++)
    {
        for(int j=0; j<order; j++)
        {
            M_inverse_big[i][j] = M[i][j];
            if (i==j)
            {
                M_inverse_big[i][j+order] = 1;
            }
            else
            {
                M_inverse_big[i][j+order] = 0;
            }
        }
    }
    // Applying Gauss Jordan Elimination 
    for(int i=0; i<order; i++)
    {
        assert(M_inverse_big[i][i] != 0.0);
        for(int j=0; j<order; j++)
        {
            if(i!=j)
            {
                ratio = M_inverse_big[j][i]/M_inverse_big[i][i];
                for(int k=0; k<2*order; k++)
                {
                    M_inverse_big[j][k] = M_inverse_big[j][k] - ratio*M_inverse_big[i][k];
                }
            }
        }
    }
    // Row Operation to Make Principal Diagonal to 1 
    for(int i=0; i<order;i++)
    {
        for(int j=order; j<2*order; j++)
        {
            M_inverse[i][j-order] = M_inverse_big[i][j]/M_inverse_big[i][i];
        }
    }
}

double sign(double x)
{
    //Return sign of value (0 = 0)
    if (x>0)
    {
        return 1.0;
    }
    else if (x<0)
    {
        return -1.0;
    }
    else
    {
        return 0;
    }
}

double minmod(double a, double b, double c)
{
    double s = (sign(a)+sign(b)+sign(c))/3.0;
    if (fabs(s)==1)
    {
        return s*fabs(min(min(a,b),c));
    }
    else
    {
        return 0;
    }
}

double min(double a, double b)
{
    if (a<b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

void LeastSquares(double* uInitialize, double xj, double dx, std::function<double(int,double)> basisFunction, int order, std::function<double(double)> inputFunction)
{
    //Compute initial condition using Least Squares method
    double x;
    double* y = new double [10];
    double** bigX = new double* [10];
    double** bigXT = new double* [order];
    double** bigXprod = new double* [order];
    double** bigXprodInv = new double* [order];
    double** hugeX = new double* [order];
    for (int i=0; i<10; i++)
    {
        bigX[i] = new double [order];
    }
    for (int l=0; l<order; l++)
    {
        bigXT[l] = new double [10];
        bigXprod[l] = new double [order];
        bigXprodInv[l] = new double [order];
        hugeX[l] = new double [10];
    }

    for (int i=0; i<10; i++)
    {
        x = xj-dx/2.0+i*dx/9.0;
        y[i] = inputFunction(x);
        for (int l=0; l<order; l++)
        {
            bigX[i][l] = basisFunction(l,2.0*(x-xj)/dx);
            bigXT[l][i] = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    MatrixMultiply(bigXT,bigX,bigXprod,order,10,order);
    for (int i=0; i<10; i++)
    {
        delete[] bigX[i];
    }
    delete[] bigX;

    MatrixInverse(bigXprod,bigXprodInv,order);
    for (int l=0; l<order; l++)
    {
        delete[] bigXprod[l];
    }
    delete[] bigXprod;

    MatrixMultiply(bigXprodInv,bigXT,hugeX,order,order,10);
    for (int l=0; l<order; l++)
    {
        delete[] bigXprodInv[l];
        delete[] bigXT[l];
    }
    delete[] bigXprodInv;
    delete[] bigXT;

    MatrixMultiply(hugeX,y,uInitialize,order,10);
    for (int l=0; l<order; l++)
    {
        delete[] hugeX[l];
    }
    delete[] hugeX;
    delete[] y;
}

double getError(double** u, int jMax, std::function<double(int,double)> basisFunction, int lMax, double dx, int tMax, double dt, std::function<double(double)> inputFunction)
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
                y[i]+=u[l][j]*basisFunction(l,(2.0/dx)*(x[i]-(j*dx+dx/2.0)));
            }
            sol[i] = inputFunction(x[i]);
            // sol[i] = inputFunction(x[i]-2.0*M_PI*tMax*dt);
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

void firstOrderEulerPlusTimes(double** uPre, double** uPost, double** M_invS, double** M_invF1, double** M_invF2, double** M_invF3, double** M_invF4,
                              double** uPlus, double plusFactor, double timesFactor, double dx, double dt, double a, int jMax, int lMax)
{
    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            uPost[l][j]=0;
            for (int i=0; i<lMax; i++)
            {
                uPost[l][j]+=M_invS[l][i]*uPre[i][j];
                uPost[l][j]-=M_invF1[l][i]*uPre[i][j];
                if (j==0)
                {
                    uPost[l][j]+=M_invF2[l][i]*uPre[i][jMax-1];
                }
                else
                {
                    uPost[l][j]+=M_invF2[l][i]*uPre[i][j-1];
                }
                if (j==jMax-1)
                {
                    uPost[l][j]-=M_invF3[l][i]*uPre[i][0];
                }
                else
                {
                    uPost[l][j]-=M_invF3[l][i]*uPre[i][j+1];
                }
                uPost[l][j]+=M_invF4[l][i]*uPre[i][j];
            }
            uPost[l][j]*=a;
            uPost[l][j]*=dt;
            uPost[l][j]/=dx;
            uPost[l][j]+=uPre[l][j];
            
            uPost[l][j]*=timesFactor;
            uPost[l][j]+=plusFactor*uPlus[l][j];
        }
    }
}

void slopeLimiter(double** uPost, int jMax, int lMax)
{
    double u1Lim;

    for (int j=0; j<jMax; j++)
    {
        if (j==0)
        {
            u1Lim = minmod(uPost[1][j],uPost[0][j+1]-uPost[0][j],uPost[0][j]-uPost[0][jMax-1]);
            if (u1Lim-uPost[1][j]>1e-10)
            {
                uPost[1][j] = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost[l][j] = 0;
                }
            }
        }
        else if (j==jMax-1)
        {
            u1Lim = minmod(uPost[1][j],uPost[0][0]-uPost[0][j],uPost[0][j]-uPost[0][j-1]);
            if (u1Lim-uPost[1][j]>1e-10)
            {
                uPost[1][j] = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost[l][j] = 0;
                }
            }
        }
        else
        {
            u1Lim = minmod(uPost[1][j],uPost[0][j+1]-uPost[0][j],uPost[0][j]-uPost[0][j-1]);
            if (u1Lim-uPost[1][j]>1e-10)
            {
                uPost[1][j] = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost[l][j] = 0;
                }
            }
        }
    }
}
