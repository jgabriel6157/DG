import pandas as pd
import matplotlib.pyplot as plt
import fnmatch
import numpy as np
import os
import random

def getFunction(basis,n,x):
    if basis=='Legendre':
        return LegendreP(n,x)
    elif basis=='LegendreOrthonormal':
        return LegendrePorthonormal(n,x)
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

    if not inputParam:
        break

inputFile.close()

values = pd.read_csv('/home/jack/Documents/DG1D/Output.csv',header=None)
values = values[0].to_numpy()
k = 0
dx = length/jMax
u = np.zeros((lMax,jMax))
for j in range(jMax):
    for l in range(lMax):
        u[l][j] = values[k]
        k=k+1


for j in range(jMax):
    y = np.zeros(10)
    sol = np.zeros(10)
    x = np.zeros(10)
    for i in range(10):
        x[i] = j*dx+i*dx/9.0
        for l in range(lMax):
            y[i] = y[i] + u[l][j]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
        sol[i] = np.sin(x[i])
    plt.plot(x,y,color='red')
    # plt.xlim((0,1))
    # plt.plot(x,sol,color='k',linestyle='--')

plt.show()