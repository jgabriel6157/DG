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
#include "FunctionMapper.hxx"
#include "Mesh.hxx"

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

int main(int argc, char* argv[])
{
    std::string argName[15] = {"a","jMax","lMax","tMax","quadratureOrder","length","dt","basis","test","alpha","input","slopeLimiter","nout","nvx","maxVX"};
    std::string argString[15];
    readFile("input.txt",argName,argString,15);

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
    int nout = assignInt(argString[12]);
    int nvx = assignInt(argString[13]);
    double domainMaxVX = assignDouble(argString[14]);

    FunctionMapper::initializeMap();

    auto basisFunction = FunctionMapper::getFunction<std::function<double(int,double)>>(basis);
    auto basisFunctionDerivative = FunctionMapper::getFunction<FunctionMapper::FunctionType1>(basis+"Derivative");
    auto inputFunction = FunctionMapper::getFunction<FunctionMapper::FunctionType2>(input);
    
    lMax+=1;
    Mesh mesh(jMax, nvx, length, domainMaxVX);
    int outputTimeStep = tMax/nout;

    std::ofstream write_output("Output.csv");
    assert(write_output.is_open());

    auto start = std::chrono::high_resolution_clock::now();
    Solver solver(mesh, dt, a, lMax, alpha);
    
    solver.createMatrices(basisFunction, basisFunctionDerivative, quadratureOrder);
    solver.initialize(basisFunction, inputFunction, SpecialFunctions::constantFunction);
    for (int t=0; t<tMax; t++)
    {
        solver.advance();

        // if (slopeLimit)
        // {
        //     solver.slopeLimiter();
        // }

        // if (t%outputTimeStep==0)
        // {
        //     for (int j=0; j<jMax; j++)
        //     {
        //         for (int l=0; l<lMax; l++)
        //         {
        //             write_output << solver.getSolution(l,j) << "\n";
        //         }
        //     }
        // }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(stop-start);
    std::cout << duration.count() << " ms" << "\n";

    // if (test)
    // {
    //     std::cout << solver.getError(tMax, basisFunction, inputFunction) << "\n";
    // }
    for (int k=0; k<nvx; k++)
    {
        for (int j=0; j<jMax; j++)
        {
            for (int lvx=0; lvx<lMax; lvx++)
            {
                for (int lx=0; lx<lMax; lx++)
                {
                    write_output << solver.getSolution(lx+lvx*lMax,j+k*jMax) << "\n";
                }
            }
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
