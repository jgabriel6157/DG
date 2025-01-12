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
#include "Parser.hxx"

#include "muParser.h"

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
void readOutput(std::string filename, double* values);
int assignInt(std::string varString);
double assignDouble(std::string varString);
bool assignBool(std::string varString);
int assignBC(std::string varString);

int main(int argc, char* argv[])
{
    std::string argName[21] = {"jMax","lMax","tMax","quadratureOrder","length","dt","basis","input","slopeLimiter","nout","nvx","maxVX","nvy","maxVY","nvz","maxVZ",
                               "ionization","cx","bgk","bc","resume"};
    std::string argString[21];
    readFile("input.txt",argName,argString,21);

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
    int nvy = assignInt(argString[12]);
    double domainMaxVY = assignDouble(argString[13]);
    int nvz = assignInt(argString[14]);
    double domainMaxVZ = assignDouble(argString[15]);
    bool ionization = assignBool(argString[16]);
    bool cx = assignBool(argString[17]);
    bool bgk = assignBool(argString[18]);
    int bc = assignBC(argString[19]);
    bool resume = assignBool(argString[20]);

    FunctionMapper::initializeMap();
    auto basisFunction = FunctionMapper::getFunction<std::function<double(int,double)>>(basis);

    Parser functionParser;
    functionParser.setExpression(input);
    auto inputFunction = functionParser.getFunction();
    
    lMax+=1;
    Mesh mesh(jMax, nvx, nvy, nvz, length, domainMaxVX, domainMaxVY, domainMaxVZ,bc);
    if (nout>tMax)
    {
        std::cout << "Invalid nout (nout > tMax)! nout set equal to tMax" << "\n";
        nout = tMax;
    }
    int outputTimeStep = tMax/nout;

    double* flattenedOutput = new double [jMax*nvx*nvy*nvz*lMax];  
    if (resume)
    {
        readOutput("lastOutput.csv",flattenedOutput);
    }


    std::ofstream write_output;
    std::ofstream write_density;
    std::ofstream write_velocity_x;
    std::ofstream write_velocity_y;
    std::ofstream write_velocity_z;
    std::ofstream write_temperature;
    std::ofstream write_moments;
    if (resume)
    {
        write_output.open("Output.csv", std::ios::app);
        write_density.open("Density.csv", std::ios::app);
        write_velocity_x.open("VelocityX.csv", std::ios::app);
        write_velocity_y.open("VelocityY.csv", std::ios::app);
        write_velocity_z.open("VelocityZ.csv", std::ios::app);
        write_temperature.open("Temperature.csv", std::ios::app);
        write_moments.open("Moments.csv", std::ios::app);
    }
    else
    {
        write_output.open("Output.csv");
        write_density.open("Density.csv");
        write_velocity_x.open("VelocityX.csv");
        write_velocity_y.open("VelocityY.csv");
        write_velocity_z.open("VelocityZ.csv");
        write_temperature.open("Temperature.csv");
        write_moments.open("Moments.csv");
    }

    assert(write_output.is_open());
    assert(write_density.is_open());
    assert(write_velocity_x.is_open());
    assert(write_velocity_y.is_open());
    assert(write_velocity_z.is_open());
    assert(write_temperature.is_open());
    assert(write_moments.is_open());

    auto start = std::chrono::high_resolution_clock::now();

    Solver solver(mesh, dt, lMax, basisFunction, quadratureOrder, ionization, cx, bgk, bc);

    solver.createMatrices();
    if (resume)
    {
        solver.resume(inputFunction, flattenedOutput);
        delete[] flattenedOutput;
    }
    else
    {
        solver.initialize(inputFunction);
    }
    std::cout << "initialization complete" << "\n";

    solver.initializeSource();
    std::cout << "Source initialization complete" << "\n";

    solver.initializeIons();
    std::cout << "Ion initialization complete" << "\n";
    
    if (bgk)
    {
        solver.initializeAlpha();
        std::cout << "alpha initialization complete" << "\n";
    }

    if (!resume)
    {
        for (int j=0; j<jMax; j++)
        {
            Vector rho = solver.getRho(j);
            Vector ux = solver.getU(j,0);
            Vector uy = solver.getU(j,1);
            Vector uz = solver.getU(j,2);
            Vector rt = solver.getE(j);
            for (int l=0; l<lMax; l++)
            {
                write_density << rho[l] << "\n";
                write_velocity_x << ux[l] << "\n";
                write_velocity_y << uy[l] << "\n";
                write_velocity_z << uz[l] << "\n";
                write_temperature << rt[l] << "\n";
            }
            // for (int kx=0; kx<nvx; kx++)
            // {
            //     for (int ky=0; ky<nvy; ky++)
            //     {
            //         for (int kz=0; kz<nvz; kz++)
            //         {
            //             for (int l=0; l<lMax; l++)
            //             {
            //                 write_output << solver.getSolution(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx) << "\n";
            //             }
            //         }
            //     }
            // }
        }
    }
    // Vector moments = solver.getMoments();
    // double M0 = moments[0];
    // double UX0 = moments[1];
    // double UY0 = moments[2];
    // double UZ0 = moments[3];
    // double E0 = moments[4];
    // double S0 = moments[5];
    std::cout << "start" << "\n";
    auto startLoop = std::chrono::high_resolution_clock::now();
    double lastTime = 0;
    for (int t=0; t<tMax; t++)
    {
        solver.advance();

        if ((t+1)%outputTimeStep==0)
        {
            // Vector moments = solver.getMoments();
            // write_moments << (moments[0]-M0)/M0 << "\n";
            // write_moments << (moments[1]-UX0)/UX0 << "\n";
            // write_moments << (moments[2]-UY0)/UY0 << "\n";
            // write_moments << (moments[3]-UZ0)/UZ0 << "\n";
            // write_moments << (moments[4]-E0)/E0 << "\n";
            // write_moments << (moments[5]-S0)/fabs(S0) << "\n";
            // if (bc==0)
            // {
            //     std::cout << (moments[0]-M0)/M0 << "\n";
            //     std::cout << (moments[1]-UX0)/UX0 << "\n";
            //     std::cout << (moments[2]-UY0)/UY0 << "\n";
            //     std::cout << (moments[3]-UZ0)/UZ0 << "\n";
            //     std::cout << (moments[4]-E0)/E0 << "\n";
            //     std::cout << (moments[5]-S0)/fabs(S0) << "\n";
            // }
            for (int j=0; j<jMax; j++)
            {
                Vector rho = solver.getRho(j);
                Vector ux = solver.getU(j,0);
                Vector uy = solver.getU(j,1);
                Vector uz = solver.getU(j,2);
                Vector rt = solver.getE(j);
                for (int l=0; l<lMax; l++)
                {
                    write_density << rho[l] << "\n";
                    write_velocity_x << ux[l] << "\n";
                    write_velocity_y << uy[l] << "\n";
                    write_velocity_z << uz[l] << "\n";
                    write_temperature << rt[l] << "\n";
                }
                // for (int kx=0; kx<nvx; kx++)
                // {
                //     for (int ky=0; ky<nvy; ky++)
                //     {
                //         for (int kz=0; kz<nvz; kz++)
                //         {
                //             for (int l=0; l<lMax; l++)
                //             {
                //                 write_output << solver.getSolution(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx) << "\n";
                //             }
                //         }
                //     }
                // }
            }
            auto stopLoop = std::chrono::high_resolution_clock::now();
            auto durationLoop = std::chrono::duration<double>(stopLoop-startLoop);
            double timePerIter = (durationLoop.count()-lastTime)/outputTimeStep;
            lastTime = durationLoop.count();
            std::cout << "t = " << t << "\n";
            std::cout << "ETA: " << timePerIter*(tMax-t) << " s\n"; 
        }
    }

    std::ofstream write_lastOutput;
    write_lastOutput.open("lastOutput.csv");
    assert(write_lastOutput.is_open());
    for (int j=0; j<jMax; j++)
    {
        for (int kx=0; kx<nvx; kx++)
        {
            for (int ky=0; ky<nvy; ky++)
            {
                for (int kz=0; kz<nvz; kz++)
                {
                    for (int l=0; l<lMax; l++)
                    {
                        write_lastOutput << solver.getSolution(l,kz+ky*nvz+kx*nvz*nvy+j*nvz*nvy*nvx) << "\n";
                    }
                }
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(stop-start);
    std::cout << duration.count() << " s" << "\n";

    write_lastOutput.close();
    write_output.close();
    write_density.close();
    write_velocity_x.close();
    write_velocity_y.close();
    write_velocity_z.close();
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

void readOutput(std::string filename, double* values)
{
    std::ifstream inputFile(filename);
    assert(inputFile.is_open());
    double value;
    int i=0;
    while (inputFile >> value)
    {
        values[i] = value;
        i++;
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

int assignBC(std::string varString)
{
    int bc = 0;
    if (varString=="periodic")
    {
        bc = 0;
    }
    else if (varString=="source")
    {
        bc = 1;
    }
    else if (varString=="copy")
    {
        bc = 2;
    }
    else
    {
        std::cout << "Invalid Boundary Condition, defaulting to periodic boundary conditions" << "\n";
    }
    return bc;
}