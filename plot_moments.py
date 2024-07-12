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

fileNameDensity = 'OutputDensity.csv'
fileNameVelocity = 'OutputVelocity.csv'
fileNameTemperature = 'OutputTemperature.csv'
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
    elif inputParam[0:3]=='nvx':
        nvx = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:5]=='maxVX':
        domainMaxVX = assignFloat(inputParam[inputParam.index('=')+2:-1])
    if not inputParam:
        break
nout+=1
lMax+=1
inputFile.close()

valuesDensity = pd.read_csv(fileNameDensity,header=None)
valuesDensity = valuesDensity[0].to_numpy()
valuesVelocity = pd.read_csv(fileNameVelocity,header=None)
valuesVelocity = valuesVelocity[0].to_numpy()
valuesTemperature = pd.read_csv(fileNameTemperature,header=None)
valuesTemperature = valuesTemperature[0].to_numpy()


m = 0
dx = length/jMax

uDensity = np.zeros((lMax,jMax))
uVelocity = np.zeros((lMax,jMax))
uTemperature = np.zeros((lMax,jMax))

fig, axs = plt.subplots(2,2)

for j in range(jMax):
    for l in range(lMax):
        uDensity[l,j] = valuesDensity[m]
        uVelocity[l,j] = valuesVelocity[m]
        uTemperature[l,j] = valuesTemperature[m]
        m = m+1

for j in range(jMax):
    xj = j*dx+dx/2
    density = np.zeros(10)
    velocity = np.zeros(10)
    temperature = np.zeros(10)
    x = np.zeros(10)
    for i in range(10):
        x[i] = j*dx+i*dx/9.0
        for l in range(lMax):
            density[i] += uDensity[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
            velocity[i] += uVelocity[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
            temperature[i] += uTemperature[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
    axs[0,0].plot(x,density,color='red')
    axs[0,0].set_title('Density')
    axs[0,1].plot(x,velocity,color='red')
    axs[0,1].set_title('Velocity')
    axs[1,0].plot(x,temperature,color='red')
    axs[1,0].set_title('Temperature')

plt.show()