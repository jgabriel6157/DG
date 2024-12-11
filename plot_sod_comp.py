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

fileNameDensity = 'Density1e1_2.csv'
fileNameVelocity = 'Velocity1e1_2.csv'
fileNameTemperature = 'Temperature1e1_2.csv'
fileNameDensity2 = 'Density1e2_2.csv'
fileNameVelocity2 = 'Velocity1e2_2.csv'
fileNameTemperature2 = 'Temperature1e2_2.csv'
fileNameDensity3 = 'Density1e3_2.csv'
fileNameVelocity3 = 'Velocity1e3_2.csv'
fileNameTemperature3 = 'Temperature1e3_2.csv'
fileNameDensity4 = 'Density1e3_2.csv'
fileNameVelocity4 = 'Velocity1e3_2.csv'
fileNameTemperature4 = 'Temperature1e3_2.csv'
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
valuesDensity2 = pd.read_csv(fileNameDensity2,header=None)
valuesDensity2 = valuesDensity2[0].to_numpy()
valuesVelocity2 = pd.read_csv(fileNameVelocity2,header=None)
valuesVelocity2 = valuesVelocity2[0].to_numpy()
valuesTemperature2 = pd.read_csv(fileNameTemperature2,header=None)
valuesTemperature2 = valuesTemperature2[0].to_numpy()
valuesDensity3 = pd.read_csv(fileNameDensity3,header=None)
valuesDensity3 = valuesDensity3[0].to_numpy()
valuesVelocity3 = pd.read_csv(fileNameVelocity3,header=None)
valuesVelocity3 = valuesVelocity3[0].to_numpy()
valuesTemperature3 = pd.read_csv(fileNameTemperature3,header=None)
valuesTemperature3 = valuesTemperature3[0].to_numpy()
valuesDensity4 = pd.read_csv(fileNameDensity4,header=None)
valuesDensity4 = valuesDensity4[0].to_numpy()
valuesVelocity4 = pd.read_csv(fileNameVelocity4,header=None)
valuesVelocity4 = valuesVelocity4[0].to_numpy()
valuesTemperature4 = pd.read_csv(fileNameTemperature4,header=None)
valuesTemperature4 = valuesTemperature4[0].to_numpy()
k = 0
dx = length/jMax
rho = np.zeros((lMax,jMax,nout))
rhou = np.zeros((lMax,jMax,nout))
rt = np.zeros((lMax,jMax,nout))
rho2 = np.zeros((lMax,jMax,nout))
rhou2 = np.zeros((lMax,jMax,nout))
rt2 = np.zeros((lMax,jMax,nout))
rho3 = np.zeros((lMax,jMax,nout))
rhou3 = np.zeros((lMax,jMax,nout))
rt3 = np.zeros((lMax,jMax,nout))
rho4 = np.zeros((lMax,jMax,nout))
rhou4 = np.zeros((lMax,jMax,nout))
rt4 = np.zeros((lMax,jMax,nout))
for t in range(nout):
    for j in range(jMax):
        for l in range(lMax):
            rho[l][j][t] = valuesDensity[k]
            rhou[l][j][t] = valuesVelocity[k]
            rt[l][j][t] = valuesTemperature[k]
            rho2[l][j][t] = valuesDensity2[k]
            rhou2[l][j][t] = valuesVelocity2[k]
            rt2[l][j][t] = valuesTemperature2[k]
            rho3[l][j][t] = valuesDensity3[k]
            rhou3[l][j][t] = valuesVelocity3[k]
            rt3[l][j][t] = valuesTemperature3[k]
            rho4[l][j][t] = valuesDensity4[k]
            rhou4[l][j][t] = valuesVelocity4[k]
            rt4[l][j][t] = valuesTemperature4[k]
            k=k+1

x = np.zeros(10*jMax)
densitySim = np.zeros(10*jMax)
velocitySim = np.zeros(10*jMax)
temperatureSim = np.zeros(10*jMax)
pressureSim = np.zeros(10*jMax)
densitySim2 = np.zeros(10*jMax)
velocitySim2 = np.zeros(10*jMax)
temperatureSim2 = np.zeros(10*jMax)
pressureSim2 = np.zeros(10*jMax)
densitySim3 = np.zeros(10*jMax)
velocitySim3 = np.zeros(10*jMax)
temperatureSim3 = np.zeros(10*jMax)
pressureSim3 = np.zeros(10*jMax)
densitySim4 = np.zeros(10*jMax)
velocitySim4 = np.zeros(10*jMax)
temperatureSim4 = np.zeros(10*jMax)
pressureSim4 = np.zeros(10*jMax)
t = 100
for j in range(jMax):
    for i in range(10):
        x[i+j*10] = j*dx+i*dx/9.0
        densityFoo = 0
        velocityFoo = 0
        temperatureFoo = 0
        densityFoo2 = 0
        velocityFoo2 = 0
        temperatureFoo2 = 0
        densityFoo3 = 0
        velocityFoo3 = 0
        temperatureFoo3 = 0
        densityFoo4 = 0
        velocityFoo4 = 0
        temperatureFoo4 = 0
        for l in range(lMax):
            densityFoo += rho[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            velocityFoo += rhou[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            temperatureFoo += rt[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            densityFoo2 += rho2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            velocityFoo2 += rhou2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            temperatureFoo2 += rt2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            densityFoo3 += rho3[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            velocityFoo3 += rhou3[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            temperatureFoo3 += rt3[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            densityFoo4 += rho4[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            velocityFoo4 += rhou4[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
            temperatureFoo4 += rt4[l][j][t]*getFunction(basis,l,(2/dx)*(x[i+j*10]-(j*dx+dx/2)))
        velocityFoo/=densityFoo
        temperatureFoo = (temperatureFoo-densityFoo*velocityFoo**2)/densityFoo
        velocityFoo2/=densityFoo2
        temperatureFoo2 = (temperatureFoo2-densityFoo2*velocityFoo2**2)/densityFoo2
        velocityFoo3/=densityFoo3
        temperatureFoo3 = (temperatureFoo3-densityFoo3*velocityFoo3**2)/densityFoo3
        velocityFoo4/=densityFoo4
        temperatureFoo4 = (temperatureFoo4-densityFoo4*velocityFoo4**2)/densityFoo4
        densitySim[i+j*10] = densityFoo
        velocitySim[i+j*10] = velocityFoo
        temperatureSim[i+j*10] = temperatureFoo
        pressureSim[i+j*10] = densityFoo*temperatureFoo
        densitySim2[i+j*10] = densityFoo2
        velocitySim2[i+j*10] = velocityFoo2
        temperatureSim2[i+j*10] = temperatureFoo2
        pressureSim2[i+j*10] = densityFoo2*temperatureFoo2
        densitySim3[i+j*10] = densityFoo3
        velocitySim3[i+j*10] = velocityFoo3
        temperatureSim3[i+j*10] = temperatureFoo3
        pressureSim3[i+j*10] = densityFoo3*temperatureFoo3
        densitySim4[i+j*10] = densityFoo4
        velocitySim4[i+j*10] = velocityFoo4
        temperatureSim4[i+j*10] = temperatureFoo4
        pressureSim4[i+j*10] = densityFoo4*temperatureFoo4
# Plotting
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(x,densitySim, label="Kn = 1e-1", color="blue", linestyle="-")
plt.plot(x,densitySim2, label="Kn = 1e-2", color="orange", linestyle="-")
plt.plot(x,densitySim3, label="Kn = 1e-3", color="green", linestyle="-")
plt.plot(xSol, rhoSol, label="Solution", color="black", linestyle="--")
plt.ylabel("Density")
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(x,velocitySim, label="Kn = 1e-1", color="blue", linestyle="-")
plt.plot(x,velocitySim2, label="Kn = 1e-2", color="orange", linestyle="-")
plt.plot(x,velocitySim3, label="Kn = 1e-3", color="green", linestyle="-")
plt.plot(xSol, velocity, label="Solution", color="black", linestyle="--")
plt.ylabel("Velocity")
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(x,temperatureSim, label="Kn = 1e-1", color="blue", linestyle="-")
plt.plot(x,temperatureSim2, label="Kn = 1e-2", color="orange", linestyle="-")
plt.plot(x,temperatureSim3, label="Kn = 1e-3", color="green", linestyle="-")
plt.plot(xSol, temperature, label="Solution", color="black", linestyle="--")
plt.ylabel("Temperature")
plt.xlabel("Position")
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(x,pressureSim, label="Kn = 1e-1", color="blue", linestyle="-")
plt.plot(x,pressureSim2, label="Kn = 1e-2", color="orange", linestyle="-")
plt.plot(x,pressureSim3, label="Kn = 1e-3", color="green", linestyle="-")
plt.plot(xSol, temperature*rhoSol, label="Solution", color="black", linestyle="--")
plt.ylabel("Pressure")
plt.xlabel("Position")
plt.legend()

plt.tight_layout()
plt.show()