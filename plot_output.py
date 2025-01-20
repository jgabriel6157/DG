import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

fig = plt.figure()
ax = fig.gca()

fileNameDensity = 'DensityCXapprox.csv'
fileNameVelocityX = 'VelocityXCXapprox.csv'
fileNameVelocityY = 'VelocityYCXapprox.csv'
fileNameVelocityZ = 'VelocityZCXapprox.csv'
fileNameTemperature = 'TemperatureCXapprox.csv'
fileNameDensity2 = 'DensityGkeyll.csv'
fileNameVelocityX2 = 'VelocityXGkeyll.csv'
fileNameVelocityY2 = 'VelocityYGkeyll.csv'
fileNameVelocityZ2 = 'VelocityZGkeyll.csv'
fileNameTemperature2 = 'TemperatureGkeyll.csv'
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
valuesVelocityX = pd.read_csv(fileNameVelocityX,header=None)
valuesVelocityX = valuesVelocityX[0].to_numpy()
valuesVelocityY = pd.read_csv(fileNameVelocityY,header=None)
valuesVelocityY = valuesVelocityY[0].to_numpy()
valuesVelocityZ = pd.read_csv(fileNameVelocityZ,header=None)
valuesVelocityZ = valuesVelocityZ[0].to_numpy()
valuesTemperature = pd.read_csv(fileNameTemperature,header=None)
valuesTemperature = valuesTemperature[0].to_numpy()
valuesDensity2 = pd.read_csv(fileNameDensity2,header=None)
valuesDensity2 = valuesDensity2[0].to_numpy()
valuesVelocityX2 = pd.read_csv(fileNameVelocityX2,header=None)
valuesVelocityX2 = valuesVelocityX2[0].to_numpy()
valuesVelocityY2 = pd.read_csv(fileNameVelocityY2,header=None)
valuesVelocityY2 = valuesVelocityY2[0].to_numpy()
valuesVelocityZ2 = pd.read_csv(fileNameVelocityZ2,header=None)
valuesVelocityZ2 = valuesVelocityZ2[0].to_numpy()
valuesTemperature2 = pd.read_csv(fileNameTemperature2,header=None)
valuesTemperature2 = valuesTemperature2[0].to_numpy()
k = 0
dx = length/jMax
rho = np.zeros((lMax,jMax,nout))
rhouX = np.zeros((lMax,jMax,nout))
rhouY = np.zeros((lMax,jMax,nout))
rhouZ = np.zeros((lMax,jMax,nout))
rt = np.zeros((lMax,jMax,nout))
rho2 = np.zeros((lMax,jMax,nout))
rhouX2 = np.zeros((lMax,jMax,nout))
rhouY2 = np.zeros((lMax,jMax,nout))
rhouZ2 = np.zeros((lMax,jMax,nout))
rt2 = np.zeros((lMax,jMax,nout))
for t in range(nout):
    for j in range(jMax):
        for l in range(lMax):
            rho[l][j][t] = valuesDensity[k]
            rhouX[l][j][t] = valuesVelocityX[k]
            rhouY[l][j][t] = valuesVelocityY[k]
            rhouZ[l][j][t] = valuesVelocityZ[k]
            rt[l][j][t] = valuesTemperature[k]
            rho2[l][j][t] = valuesDensity2[k]
            rhouX2[l][j][t] = valuesVelocityX2[k]
            rhouY2[l][j][t] = valuesVelocityY2[k]
            rhouZ2[l][j][t] = valuesVelocityZ2[k]
            rt2[l][j][t] = valuesTemperature2[k]
            k=k+1

t = 50 #Output step
y = np.zeros(10)
y2 = np.zeros(10)
x = np.zeros(10)
for j in range(jMax):
    for i in range(10):
        x[i] = j*dx+i*dx/9.0
        density = 0
        velocityX = 0
        velocityY = 0
        velocityZ = 0
        temperature = 0
        density2 = 0
        velocityX2 = 0
        velocityY2 = 0
        velocityZ2 = 0
        temperature2 = 0
        for l in range(lMax):
            density += rho[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityX += rhouX[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityY += rhouY[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityZ += rhouZ[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            temperature += rt[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            density2 += rho2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityX2 += rhouX2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityY2 += rhouY2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            velocityZ2 += rhouZ2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
            temperature2 += rt2[l][j][t]*getFunction(basis,l,(2/dx)*(x[i]-(j*dx+dx/2)))
        velocityX/=density
        velocityY/=density
        velocityZ/=density
        temperature = (temperature-density*(velocityX**2+velocityY**2+velocityZ**2))/(3*density)
        velocityX2/=density2
        velocityY2/=density2
        velocityZ2/=density2
        temperature2 = (temperature2-density2*(velocityX2**2+velocityY2**2+velocityZ2**2))/(3*density2)
        y[i] = density*1e18
        # y[i] = temperature
        y2[i] = density2*1e18
        # y2[i] = temperature2
    plt.plot(x-20,y,color='red')
    plt.plot(x-20,y2,color='k')
plt.plot(0,0,color='k',label = 'Gkeyll')
plt.plot(0,0,color='red',label='Janev-Smith approximation')

ax.set_yscale('log')
plt.ylim(5e13,2e19)
# plt.ylim(28,70)

plt.xlim(-20,20)
plt.legend()
plt.show()