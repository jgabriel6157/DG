#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <chrono>
#include <functional>
#include <map>

#include "Vector.hxx"
#include "Matrix.hxx"
#include "SpecialFunctions.hxx"
#include "GaussianQuadrature.hxx"

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

double sign(double x);
double minmod(double a, double b, double c);
double min(double a, double b);

void LeastSquares(Vector& uInitialize, double xj, double dx, std::function<double(int,double)> basisFunction, int order, std::function<double(double)> inputFunction);

double getError(Matrix u, int jMax, std::function<double(int,double)> basisFunction, int lMax, double dx, int tMax, double dt, std::function<double(double)> inputFunction);

void firstOrderEulerPlusTimes(Matrix uPre, Matrix& uPost, Matrix M_invS, Matrix M_invF1, Matrix M_invF2, Matrix M_invF3, Matrix M_invF4,
                              Matrix uPlus, double plusFactor, double timesFactor, double dx, double dt, double a, int jMax, int lMax);
void slopeLimiter(Matrix& uPost, int jMax, int lMax);

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

    Matrix M(lMax,lMax);
    Matrix M_inv(lMax,lMax);
    Matrix S(lMax,lMax);
    Matrix F1(lMax,lMax);
    Matrix F2(lMax,lMax);
    Matrix F3(lMax,lMax);
    Matrix F4(lMax,lMax);
    Matrix M_invS(lMax,lMax);
    Matrix M_invF1(lMax,lMax);
    Matrix M_invF2(lMax,lMax);
    Matrix M_invF3(lMax,lMax);
    Matrix M_invF4(lMax,lMax);
    Matrix uPre(lMax,jMax);
    Matrix uPost(lMax,jMax);
    Matrix uIntermediate(lMax,jMax);
    Vector uInitialize(lMax);
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

    Vector roots = SpecialFunctions::legendreRoots(quadratureOrder);
    Vector weights = GaussianQuadrature::calculateWeights(quadratureOrder, roots);

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M(i,j) = GaussianQuadrature::integrate(basisFunction,i,basisFunction,j,quadratureOrder,roots,weights)/2;
            S(i,j) = GaussianQuadrature::integrate(basisFunctionDerivative,i,basisFunction,j,quadratureOrder,roots,weights);
            F1(i,j) = fluxFactorPlus*(basisFunction(i,1))*(basisFunction(j,1));
            F2(i,j) = fluxFactorPlus*(basisFunction(i,-1))*(basisFunction(j,1));
            F3(i,j) = fluxFactorMinus*(basisFunction(i,1))*(basisFunction(j,-1));
            F4(i,j) = fluxFactorMinus*(basisFunction(i,-1))*(basisFunction(j,-1));
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
    M_invF1 = M_inv*F1;
    M_invF2 = M_inv*F2;
    M_invF3 = M_inv*F3;
    M_invF4 = M_inv*F4;

    for (int j=0; j<jMax; j++)
    {
        xj = j*dx+dx/2.0;
        LeastSquares(uInitialize, xj, dx, basisFunction, lMax, inputFunction);
        for (int l=0; l<lMax; l++)
        {
            uPre(l,j) = uInitialize[l];
        }
    }

    for (int t=0; t<tMax; t++)
    {
        firstOrderEulerPlusTimes(uPre,uPost,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,0.0,1.0,dx,dt,a,jMax,lMax);

        firstOrderEulerPlusTimes(uPost,uIntermediate,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,3.0/4.0,1.0/4.0,dx,dt,a,jMax,lMax);

        firstOrderEulerPlusTimes(uIntermediate,uPost,M_invS,M_invF1,M_invF2,M_invF3,M_invF4,uPre,1.0/3.0,2.0/3.0,dx,dt,a,jMax,lMax);

        // slopeLimiter(uPost,jMax,lMax);

        uPre = uPost;
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
            write_output << uPre(l,j) << "\n";
        }
    }

    write_output.close();

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

void LeastSquares(Vector& uInitialize, double xj, double dx, std::function<double(int,double)> basisFunction, int order, std::function<double(double)> inputFunction)
{
    //Compute initial condition using Least Squares method
    double x;
    Vector y(10);
    Matrix bigX(10,order);
    Matrix bigXT(order,10);
    Matrix bigXprod(order,order);
    Matrix bigXprodInv(order,order);
    Matrix hugeX(order,10);

    for (int i=0; i<10; i++)
    {
        x = xj-dx/2.0+i*dx/9.0;
        y[i] = inputFunction(x);
        for (int l=0; l<order; l++)
        {
            bigX(i,l) = basisFunction(l,2.0*(x-xj)/dx);
            bigXT(l,i) = basisFunction(l,2.0*(x-xj)/dx);
        }
    }

    bigXprod = bigXT*bigX;

    bigXprodInv = bigXprod.CalculateInverse();

    hugeX = bigXprodInv*bigXT;

    uInitialize = hugeX*y;
}

double getError(Matrix u, int jMax, std::function<double(int,double)> basisFunction, int lMax, double dx, int tMax, double dt, std::function<double(double)> inputFunction)
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
                y[i] += u(l,j)*basisFunction(l,(2.0/dx)*(x[i]-(j*dx+dx/2.0)));
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

void firstOrderEulerPlusTimes(Matrix uPre, Matrix& uPost, Matrix M_invS, Matrix M_invF1, Matrix M_invF2, Matrix M_invF3, Matrix M_invF4,
                              Matrix uPlus, double plusFactor, double timesFactor, double dx, double dt, double a, int jMax, int lMax)
{
    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            uPost(l,j)=0;
            for (int i=0; i<lMax; i++)
            {
                uPost(l,j)+=M_invS(l,i)*uPre(i,j);
                uPost(l,j)-=M_invF1(l,i)*uPre(i,j);
                if (j==0)
                {
                    uPost(l,j)+=M_invF2(l,i)*uPre(i,jMax-1);
                }
                else
                {
                    uPost(l,j)+=M_invF2(l,i)*uPre(i,j-1);
                }
                if (j==jMax-1)
                {
                    uPost(l,j)-=M_invF3(l,i)*uPre(i,0);
                }
                else
                {
                    uPost(l,j)-=M_invF3(l,i)*uPre(i,j+1);
                }
                uPost(l,j)+=M_invF4(l,i)*uPre(i,j);
            }
            uPost(l,j)*=a;
            uPost(l,j)*=dt;
            uPost(l,j)/=dx;
            uPost(l,j)+=uPre(l,j);
            
            uPost(l,j)*=timesFactor;
            uPost(l,j)+=plusFactor*uPlus(l,j);
        }
    }
}

void slopeLimiter(Matrix& uPost, int jMax, int lMax)
{
    double u1Lim;

    for (int j=0; j<jMax; j++)
    {
        if (j==0)
        {
            u1Lim = minmod(uPost(1,j),uPost(0,j+1)-uPost(0,j),uPost(0,j)-uPost(0,jMax-1));
            if (u1Lim-uPost(1,j)>1e-10)
            {
                uPost(1,j) = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost(l,j) = 0;
                }
            }
        }
        else if (j==jMax-1)
        {
            u1Lim = minmod(uPost(1,j),uPost(0,0)-uPost(0,j),uPost(0,j)-uPost(0,j-1));
            if (u1Lim-uPost(1,j)>1e-10)
            {
                uPost(1,j) = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost(l,j) = 0;
                }
            }
        }
        else
        {
            u1Lim = minmod(uPost(1,j),uPost(0,j+1)-uPost(0,j),uPost(0,j)-uPost(0,j-1));
            if (u1Lim-uPost(1,j)>1e-10)
            {
                uPost(1,j) = u1Lim;
                for (int l=2; l<lMax; l++)
                {
                    uPost(l,j) = 0;
                }
            }
        }
    }
}
