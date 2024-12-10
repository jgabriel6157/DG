import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import pandas as pd
from scipy.interpolate import interp1d

# Define constants
# gamma = 5.0 / 3  # Ratio of specific heats
gamma = 3.0
time = 0.1  # Simulation time (seconds)
x0 = 0.5  # Initial discontinuity position

# Left state (region 1)
P1 = 1.0  # Pressure
rho1 = 1.0  # Density

# Right state (region 5)
P5 = 0.1  # Pressure
rho5 = 0.125  # Density

cs1 = sqrt(gamma * P1 / rho1)  # Speed of sound in region 1
cs5 = sqrt(gamma * P5 / rho5)  # Speed of sound in region 5

# Calculate post-shock states iteratively
G = (gamma - 1) / (gamma + 1)
beta = (gamma - 1) / (2 * gamma)

P3 = 0.5 * (P5 + P1)
while True:
    u4 = (P3 - P5) * sqrt((1 - G) / (rho5 * (P3 + G * P5)))
    u3 = u4
    P3_new = (P1 ** beta - u3 * sqrt(G**2 * rho1 / ((1 - G**2) * P1 ** (1 / gamma)))) ** (1 / beta)
    if abs(P3_new - P3) < 1e-10:
        break
    P3 = 0.5 * (P3 + P3_new)

P4 = P3
rho3 = rho1 * (P3 / P1) ** (1 / gamma)
rho4 = rho5 * (P4 + G * P5) / (P5 + G * P4)

# Positions of discontinuities
pos12 = x0 - cs1 * time
pos34 = x0 + u3 * time
v_shock = u3 * ((rho4 / rho5) / ((rho4 / rho5) - 1))
pos45 = x0 + v_shock * time

cs2 = cs1 - ((gamma - 1) / 2) * u3
pos23 = x0 + (u3 - cs2) * time

# Solution on a grid
nx = 400
xSol = np.linspace(0, 1, nx)
rhoSol = np.zeros_like(xSol)
velocity = np.zeros_like(xSol)
pressure = np.zeros_like(xSol)

for i, xi in enumerate(xSol):
    if xi < pos12:
        rhoSol[i] = rho1
        velocity[i] = 0
        pressure[i] = P1
    elif xi < pos23:
        c = G * ((x0 - xi) / time) + (1 - G) * cs1
        rhoSol[i] = rho1 * (c / cs1) ** (2 / (gamma - 1))
        velocity[i] = (1 - c / cs1) * 2 / (gamma - 1) * cs1
        pressure[i] = P1 * (c / cs1) ** (2 * gamma / (gamma - 1))
    elif xi < pos34:
        rhoSol[i] = rho3
        velocity[i] = u3
        pressure[i] = P3
    elif xi < pos45:
        rhoSol[i] = rho4
        velocity[i] = u4
        pressure[i] = P4
    else:
        rhoSol[i] = rho5
        velocity[i] = 0
        pressure[i] = P5

temperature = pressure / rhoSol

rho_interp = interp1d(xSol,rhoSol,kind='cubic')
vel_interp = interp1d(xSol,velocity,kind='cubic')
pressure_interp = interp1d(xSol,pressure,kind='cubic')
temp_interp = interp1d(xSol,temperature,kind='cubic')

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
fileNameVelocity = 'Velocity.csv'
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

x = np.zeros(10*jMax)
densitySim = np.zeros(10*jMax)
velocitySim = np.zeros(10*jMax)
temperatureSim = np.zeros(10*jMax)
pressureSim = np.zeros(10*jMax)
t = 100
for j in range(jMax):
    for i in range(10):
        x[i+j*10] = j*dx+i*dx/9.0
        densityFoo = 0
        velocityFoo = 0
        temperatureFoo = 0
        for l in range(lMax):
            densityFoo += rho[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            velocityFoo += rhou[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            temperatureFoo += rt[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
        velocityFoo/=densityFoo
        temperatureFoo = (temperatureFoo-densityFoo*velocityFoo**2)/densityFoo
        densitySim[i+j*10] = densityFoo
        velocitySim[i+j*10] = velocityFoo
        temperatureSim[i+j*10] = temperatureFoo
        pressureSim[i+j*10] = densityFoo*temperatureFoo

# Plotting
plt.figure(figsize=(12, 8))

plt.subplot(4, 1, 1)
plt.plot(xSol, rhoSol, label="Solution", color="blue", linestyle="-")
plt.plot(x,densitySim, label="Simulation", color="blue", linestyle="--")
plt.ylabel("Density")
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(xSol, velocity, label="Solution", color="orange", linestyle="-")
plt.plot(x,velocitySim, label="Simulation", color="orange", linestyle="--")
plt.ylabel("Velocity")
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(xSol, temperature, label="Solution", color="green", linestyle="-")
plt.plot(x,temperatureSim, label="Simulation", color="green", linestyle="--")
plt.ylabel("Temperature")
plt.xlabel("Position")
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(xSol, temperature*rhoSol, label="Solution", color="red", linestyle="-")
plt.plot(x,temperatureSim*densitySim, label="Simulation", color="red", linestyle="--")
plt.ylabel("Pressure")
plt.xlabel("Position")
plt.legend()

plt.tight_layout()
plt.show()