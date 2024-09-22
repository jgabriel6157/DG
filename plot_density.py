import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import fnmatch
import numpy as np
from scipy.integrate import quad
from scipy.integrate import simpson
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

fileName = 'Density.csv'
fileNameSol = 'Density1.csv'
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
    elif inputParam[0:3]=='nvx':
        nvx = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:5]=='maxVX':
        domainMaxVX = assignFloat(inputParam[inputParam.index('=')+2:-1])
    if not inputParam:
        break
nout+=1
lMax+=1
inputFile.close()

A = 0.291e-7
P = 0
U = 13.6 / 20
X = 0.232
K = 0.39
n0 = 5*10**18  
Tw = 2  
csL = np.sqrt(20)  
csR = -np.sqrt(20)  
Crec = 4.98*10**18  
Lz = length  
    
def sigma(A, P, U, X, K):
    return A * (1 + P * np.sqrt(U)) * U**K * np.exp(-U) / (X + U) * 1e-6

def integrand_1(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csL) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * sigma(A, P, U, X, K) * np.sqrt(1.672e-27 / 1.602e-19) * (z) / v)

def integrand_2(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csR) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * sigma(A, P, U, X, K) * np.sqrt(1.672e-27 / 1.602e-19) * (z - Lz) / v)

def f(z):
    integral_1, _ = quad(integrand_1, 0, np.inf, args=(z,))
    integral_2, _ = quad(integrand_2, -np.inf, 0, args=(z,))
    return integral_1 + integral_2

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
valuesSol = pd.read_csv(fileNameSol,header=None)
valuesSol = valuesSol[0].to_numpy()
m = 33600
dx = length/jMax
dvx = 2*domainMaxVX/(nvx-1)
# dvx = 1.0/nvx
u = np.zeros((lMax,jMax))
uSol = np.zeros((lMax,jMax))

fig = plt.figure()
ax = fig.gca()
ax.set_yscale('log')

for j in range(jMax):
    for lx in range(lMax):
        u[lx,j] = values[m]
        # uSol[lx,j] = valuesSol[m]
        m = m+1

l2Norm = 0
solution = 0

res = 9
for j in range(jMax):
    xj = j*dx+dx/2
    y = np.zeros(res)
    x = np.zeros(res)
    sol = np.zeros(res)
    error = np.zeros(res)
    y2 = np.zeros(res)
    for i in range(res):
        x[i] = j*dx+i*dx/(res-1)
        sol[i] = f(x[i])
        for l in range(lMax):
            y[i] += u[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
            # sol[i] += uSol[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
        error[i] = (y[i]-sol[i])**2
        y2[i] = y[i]**2
    l2Norm+=simpson(error,x=x)
    solution+=simpson(y2,x=x)
    plt.plot(x,y,color='red')
    plt.plot(x,sol,'k:')

print(np.sqrt(l2Norm/solution))

plt.show()
