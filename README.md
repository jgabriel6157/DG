# DG Code
## Overview
DG code is the temporary name for this code which solves the 1D advection equation du/dt+adu/dx=0 for arbitrary order and grid size along with a python script to plot results.

## Features
- Solves 1X0V advection equation to arbitrary order
- Uniform, structured mesh
- Periodic boundary conditions
- Simple minmod slope limiter
- Lax-Friedrichs (Rusanov) numerical flux
- Input parameters specified via input.txt
- Results can be animated using plot_results.py

## Requirements
- C++ compiler (supporting C++17 or higher)
- Python

## Usage
1. Compile the code using the provided Makefile:
```
make
```
2. Run the solver executable:
```
./dg1d
```
The time it took the simulation to run and error (if test = true in input.txt) are output
3. After running the solver, animate the results using the python script:
```
python3 plot_results.py
```
This will generate an animation of the solution.

## Input File Format
The input file (input.txt) must contain the following parameters:
- Advection speed (a)
- Number of cells (jMax)
- Order of solution (lMax)
- Number of time steps (tMax)
- Order of Gaussian Quadrature (quadratureOrder)
- Length of domain (length)
- Time step (dt)
- Basis function (basis)
    - Only options are legendre, quadratic (lMax = 2), linear (lMax = 1)
- If the error should be output (test)
- What value alpha to use for Lax-Friedrichs numerical flux
    - alpha = 0 corresponds to upwind flux
    - alpha = 1 corresponds to centered flux
- Initial function (input)
    - Only options are sin (a sine wave), pulse (a Gaussian pulse), topHat (a discontinuous top hat function)
- If slope limiter should be used (slopeLimiter)
- Number of outputed time steps (nout)

## Code Structure
- dg1dAdvection.cxx: Entry point of the program, reads input parameters and orchestrates the simulation.
- Solver.cxx: Defines the solver class responsible for solving the advection equation and other relevant tasks
- Mesh.cxx: Defines the mesh representing the computational grid
- Matrix/Vector.cxx: Defines the Matrix/Vector types and operator definitions
- SpecialFunctions.cxx: Defines all non-native functions
- GaussianQuadrature.cxx: Performs Gaussian quadrature numerical integration
- FunctionMapper.cxx: Maps string of function name to function

## Author
-[Jack Gabriel](https://github.com/jgabriel6157)


