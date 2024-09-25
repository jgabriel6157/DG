import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import fnmatch
import numpy as np
import os
import random
import time

def getFunction(basis,n,x):
    if basis=='legendre':
        return LegendreP(n,x)
    elif basis=='legendreOrthonormal':
        return LegendrePorthonormal(n,x)
    elif basis=='quadratic':
        return Quadratic(n,x)
    elif basis=='linear':
        return linear(n,x)
    else:
        return 0

def LegendreP(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    else:
        return ((2.0*n-1.0)*x*LegendreP(n-1,x)-(n-1)*LegendreP(n-2,x))/n
    
def LegendrePorthonormal(n,x):
    return np.sqrt((2.0*n+1.0)/2.0)*LegendreP(n,x)

def Quadratic(n,x):
    if n==0:
        return -x*(1.0-x)/2.0
    elif n==1:
        return (1.0-x)*(1.0+x)
    else:
        return x*(1+x)/2
    
def linear(n,x):
    if n==0:
        return (1.0-x)/2.0
    else:
        return (1+x)/2

def assignFloat(varString):
    number = 1.0
    valueCheck = "*"
    value = varString
    while("*" in valueCheck):
        valueCheck = value
        if "*" in value:
            dum1 = value[0:value.index("*")]
            value = value[value.index("*")+1:]
        else:
            dum1 = value
        if dum1 == "pi":
            number*=np.pi
        else:
            number*=float(dum1)
        
        
    return number

fileName = 'Output.csv'
# fileNameSol = 'OutputSol.csv'
inputFile = open('input.txt','r')

while True:
    inputParam = inputFile.readline()

    if inputParam[0:4]=='jMax':        
        jMax = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:6]=='length':
        value = inputParam[inputParam.index('=')+2:-1]
        length = assignFloat(value)
    elif inputParam[0:4]=='lMax':
        lMax = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:5]=='basis':
        basis = inputParam[inputParam.index('=')+2:-1]
    elif inputParam[0:4]=='tMax':
        tMax = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:2]=='dt':
        dt = float(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:4]=='nout':
        nout = int(inputParam[inputParam.index('=')+2:-1])
    if not inputParam:
        break
nout+=1
lMax+=1
inputFile.close()

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
# valuesSol = pd.read_csv(fileNameSol,header=None)
# valuesSol = valuesSol[0].to_numpy()
# m = jMax*lMax*nvx*(nout-1)
m=288
dx = length/jMax
# dvx = 1.0/nvx
u = np.zeros((lMax,jMax))
# uSol = np.zeros((lMax,jMax,nvx))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for j in range(jMax):
    for lx in range(lMax):
        u[lx,j] = values[m]
        # uSol[lx,j,k] = valuesSol[m]
        m = m+1
print(m)
# vx = 0
# for j in range(jMax):
#     xj = j*dx+dx/2
#     y = np.zeros(10)
#     x = np.zeros(10)
#     for i in range(10):
#         x[i] = j*dx+i*dx/9.0
#         for l in range(lMax):
#             y[i] += u[l][j][vx]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
#     plt.plot(x,y,color='red')
sumNum = 0
sumDem = 0

# Number of quadrature points
nQuad = 10
# Precompute quadrature points and weights on the reference interval [-1, 1]
quadPoints, quadWeights = np.polynomial.legendre.leggauss(nQuad)

for vx in [2*np.pi]:
    sumNum = 0
    sumDem = 0
    for j in range(jMax):
        xj = j * dx + dx / 2
        y_offset = vx
        y = np.zeros(nQuad)
        sol = np.zeros(nQuad)
        
        # Map quadrature points from [-1, 1] to [j*dx, (j+1)*dx]
        x = 0.5 * dx * (quadPoints + 1) + j * dx
        
        for i in range(nQuad):
            for l in range(lMax):
                # Apply basis function and sum contributions
                y[i] += u[l][j] * getFunction(basis, l, (2.0 / dx) * (x[i] - xj))
            
            # Define the solution to compare against
            # sol[i] = max(
            #     np.exp(-(x[i] - np.pi - tMax * dt * y_offset) ** 2),
            #     np.exp(-(x[i] + np.pi - tMax * dt * y_offset) ** 2),
            #     np.exp(-(x[i] - 3 * np.pi - tMax * dt * y_offset) ** 2)
            # )
            sol[i] = np.sin(x[i] - tMax * dt * y_offset)
        
        # Use Gaussian quadrature weights in the error calculation
        for i in range(nQuad):
            # L2 error requires squaring the difference for sumNum
            sumNum += (y[i] - sol[i]) ** 2 * quadWeights[i] * (0.5 * dx)  # Account for the dx scaling in the transformation
            # Also square the solution for sumDem
            sumDem += sol[i] ** 2 * quadWeights[i] * (0.5 * dx)
        
        # Plotting the results (optional)
        ax.plot(x, [y_offset] * len(x), y-sol, color='red')
    print(np.sqrt(sumNum / sumDem))

# Print the L2 error by taking the square root of the ratio
print(np.sqrt(sumNum / sumDem))

# for vx in range(nvx):
#     for j in range(jMax):
#         xj = j*dx+dx/2
#         y = np.zeros(10)
#         x = np.zeros(10)
#         sol = np.zeros(10)
#         y_offset = -domainMaxVX + vx*dvx
#         for i in range(10):
#             x[i] = j*dx+i*dx/9.0
#             for l in range(lMax):
#                 y[i] += u[l][j][vx]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
#                 # sol[i] += uSol[l][j][vx]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
#             sol[i] = max(np.exp(-(x[i]-np.pi-tMax*dt*y_offset)**2),np.exp(-(x[i]+np.pi-tMax*dt*y_offset)**2),np.exp(-(x[i]-3*np.pi-tMax*dt*y_offset)**2))
#             sumNum+=abs(y[i]-sol[i])
#             sumDem+=sol[i]
#             # sol[i] = np.exp(-(x[i]-np.pi-tMax*dt*y_offset)**2)
#         ax.plot(x,[y_offset]*len(x),y,color='red')
#         # ax.plot(x,[y_offset]*len(x),sol,color='black')
#         # ax.plot(x,[y_offset]*len(x),y-sol,color='red')
# print(sumNum/sumDem)
# for j in range(jMax):
#     xj = j*dx+dx/2
#     for k in range(nvx):
#         y = np.zeros((10,10))
#         x = np.zeros((10,10))
#         vx = np.zeros((10,10))
#         vx_center = -domainMaxVX+k*dvx+dvx/2
#         # vx_center = domainMaxVX+k*dvx+dvx/2
#         for i in range(10):
#             for n in range(10):
#                 x[i,n] = j*dx+i*dx/9.0
#                 vx[i,n] = -domainMaxVX+k*dvx+n*dvx/9.0
#                 # vx[i,n] = domainMaxVX+k*dvx+n*dvx/9.0
#                 for lx in range(lMax):
#                     for lvx in range(lMax):
#                         y[i,n] += u[lx,lvx,j,k]*getFunction(basis,lx,(2/dx)*(x[i,n]-xj))*getFunction(basis,lvx,(2/dvx)*(vx[i,n]-vx_center))
#         ax.plot_wireframe(x,vx,y)
plt.tight_layout()
plt.show()