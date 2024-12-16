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

fileNameDensity = 'Density.csv'
fileNameVelocity = 'VelocityX.csv'
fileNameTemperature = 'Temperature.csv'
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

fig,ax = plt.subplots()
ax.set_yscale('log')
lines = [ax.plot([], [], lw=2,color='red')[0] for _ in range(jMax)]
plt.xlim(0,length)
plt.ylim(1e8,1e19)
# plt.ylim(-0.25,1.25)
# plt.ylim(0,40)

valuesDensity = pd.read_csv(fileNameDensity,header=None)
valuesDensity = valuesDensity[0].to_numpy()
valuesVelocity = pd.read_csv(fileNameVelocity,header=None)
valuesVelocity = valuesVelocity[0].to_numpy()
valuesTemperature = pd.read_csv(fileNameTemperature,header=None)
valuesTemperature = valuesTemperature[0].to_numpy()
k = 0
dx = length/jMax
rho = np.zeros((lMax,jMax,nout))
rhou = np.zeros((lMax,jMax,nout))
rt = np.zeros((lMax,jMax,nout))
for t in range(nout):
    for j in range(jMax):
        for l in range(lMax):
            rho[l][j][t] = valuesDensity[k]
            rhou[l][j][t] = valuesVelocity[k]
            rt[l][j][t] = valuesTemperature[k]
            k=k+1

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def generate_data(t,j):
    y = np.zeros(10)
    sol = np.zeros(10)
    x = np.zeros(10)
    for i in range(10):
        x[i] = j*dx+i*dx/9.0
        density = 0
        velocity = 0
        temperature = 0
        for l in range(lMax):
            density += rho[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocity += rhou[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            temperature += rt[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
        velocity/=density
        temperature = (temperature-density*velocity**2)/(3*density)
        y[i] = density*1e18
    return x,y

x = np.zeros((jMax,10))
y = np.zeros((jMax,10))
def animate(t):
    for j in range(jMax):
        x[j],y[j] = generate_data(t,j)
    for j, line in enumerate(lines):
        line.set_data(x[j],y[j])
    ax.set_title(f"Timestep: {t}")
    return lines

ani = FuncAnimation(fig, animate, frames=nout, init_func=init, repeat=False, interval = 100)

plt.show()