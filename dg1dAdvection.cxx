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

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);

int main(int argc, char* argv[])
{
    std::string argName[12] = {"jMax","lMax","tMax","quadratureOrder","length","dt","basis","input","slopeLimiter","nout","nvx","maxVX"};
    std::string argString[12];
    readFile("input.txt",argName,argString,12);

    int jMax = assignInt(argString[0]);
    int lMax = assignInt(argString[1]);
    int tMax = assignInt(argString[2]);
    int quadratureOrder = assignInt(argString[3]);
    double length = assignDouble(argString[4]);
    double dt = assignDouble(argString[5]);
    std::string basis = argString[6];
    std::string input = argString[7];
    bool slopeLimit = assignBool(argString[8]);
    int nout = assignInt(argString[9]);
    int nvx = assignInt(argString[10]);
    double domainMaxVX = assignDouble(argString[11]);

    FunctionMapper::initializeMap();

    auto basisFunction = FunctionMapper::getFunction<std::function<double(int,double)>>(basis);
    auto basisFunctionDerivative = FunctionMapper::getFunction<FunctionMapper::FunctionType1>(basis+"Derivative");
    auto inputFunction = FunctionMapper::getFunction<FunctionMapper::FunctionType2>(input);
    
    lMax+=1;
    Mesh mesh(jMax, nvx, length, domainMaxVX);
    if (nout>tMax)
    {
        std::cout << "Invalid nout (nout > tMax)! nout set equal to tMax" << "\n";
        nout = tMax;
    }
    int outputTimeStep = tMax/nout;
    
    std::ofstream write_output("Output.csv");
    assert(write_output.is_open());

    std::ofstream write_density("Density.csv");
    assert(write_density.is_open());

    std::ofstream write_velocity("Velocity.csv");
    assert(write_velocity.is_open());

    std::ofstream write_temperature("Temperature.csv");
    assert(write_temperature.is_open());

    std::ofstream write_moments("Moments.csv");
    assert(write_moments.is_open());

    auto start = std::chrono::high_resolution_clock::now();

    Solver solver(mesh, dt, lMax, basisFunction, quadratureOrder);

    solver.createMatrices();

    solver.initialize(SpecialFunctions::inelasticICx, SpecialFunctions::inelasticICvx);

    // solver.initialize(SpecialFunctions::gaussianPulse, inputFunction);

    std::cout << "initialization complete" << "\n";
    solver.initializeAlpha();
    std::cout << "alpha initialization complete" << "\n";
    for (int j=0; j<jMax; j++)
    {
        Vector rho = solver.getMoment(j,0);
        Vector u = solver.getMoment(j,1);
        Vector rt = solver.getMoment(j,2);
        for (int l=0; l<lMax; l++)
        {
            write_density << rho[l] << "\n";
            write_velocity << u[l] << "\n";
            write_temperature << rt[l] << "\n";
        }
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<lMax; l++)
            {
                write_output << solver.getSolution(l,k+j*nvx) << "\n";
            }
        }
    }

    Vector moments = solver.getMoments();
    double M0 = moments[0];
    double U0 = moments[1];
    double E0 = moments[2];
    double S0 = moments[3];
    // std::cout << M0 << "\n";
    // std::cout << U0 << "\n";
    // std::cout << E0 << "\n";
    // std::cout << S0 << "\n";
    std::cout << "start" << "\n";
    for (int t=0; t<tMax; t++)
    {
        solver.advance();

        if ((t+1)%outputTimeStep==0)
        {
            std::cout << "t = " << t << "\n";
            Vector moments = solver.getMoments();
            write_moments << (moments[0]-M0)/M0 << "\n";
            write_moments << (moments[1]-U0)/U0 << "\n";
            write_moments << (moments[2]-E0)/E0 << "\n";
            write_moments << (moments[3]-S0)/fabs(S0) << "\n";
            // std::cout << (moments[0]-M0)/M0 << "\n";
            // std::cout << (moments[1]-U0)/U0 << "\n";
            // std::cout << (moments[2]-E0)/E0 << "\n";
            // std::cout << (moments[3]-S0)/fabs(S0) << "\n";
            for (int j=0; j<jMax; j++)
            {
                Vector rho = solver.getMoment(j,0);
                Vector u = solver.getMoment(j,1);
                Vector rt = solver.getMoment(j,2);
                for (int l=0; l<lMax; l++)
                {
                    write_density << rho[l] << "\n";
                    write_velocity << u[l] << "\n";
                    write_temperature << rt[l] << "\n";
                }
                for (int k=0; k<nvx; k++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        write_output << solver.getSolution(l,k+j*nvx) << "\n";
                    }
                }
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(stop-start);
    std::cout << duration.count() << " s" << "\n";

    write_output.close();
    write_density.close();
    write_velocity.close();
    write_temperature.close();
    write_moments.close();

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
