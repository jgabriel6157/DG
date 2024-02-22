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
#include "Solver.hxx"

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

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

    bool slopeLimit = true;
    
    double dx = length/jMax;
    double fluxFactorPlus = (1.0+SpecialFunctions::sign(a)*(1.0-alpha))/2.0;
    double fluxFactorMinus = (1.0-SpecialFunctions::sign(a)*(1.0-alpha))/2.0;

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

    Solver solver(dx, dt, a, jMax, lMax, M_invS, M_invF1, M_invF2, M_invF3, M_invF4);

    solver.initialize(basisFunction, inputFunction);

    for (int t=0; t<tMax; t++)
    {
        solver.advance();

        if (slopeLimit)
        {
            solver.slopeLimiter();
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(stop-start);
    std::cout << duration.count() << " ms" << "\n";

    if (test)
    {
        std::cout << solver.getError(tMax, basisFunction, inputFunction) << "\n";
    }

    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            write_output << solver.getSolution(l,j) << "\n";
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
