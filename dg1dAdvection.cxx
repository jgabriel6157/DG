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
// #include "FunctionMapper.hxx"

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

int main(int argc, char* argv[])
{
    std::string argName[12] = {"a","jMax","lMax","tMax","quadratureOrder","length","dt","basis","test","alpha","input","slopeLimiter"};
    std::string argString[12];
    readFile("input.txt",argName,argString,12);

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
    bool slopeLimit = assignBool(argString[11]);

    // auto basisFunction = FunctionMapper::getFunction<FunctionMapper::FunctionType1>(basis);
    // auto basisFunctionDerivative = FunctionMapper::getFunction<FunctionMapper::FunctionType1>(basis+"Derivative");
    // auto inputFunction = FunctionMapper::getFunction<FunctionMapper::FunctionType2>(input);
    
    double dx = length/jMax;
    
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
        {"sin", [](double x){return sin(x);}},
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

    Solver solver(dx,dt,a,jMax,lMax,alpha);

    solver.createMatrices(basisFunction, basisFunctionDerivative, quadratureOrder);

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
