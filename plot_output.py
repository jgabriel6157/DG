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

# Define vectors for parameters and file suffixes
# suffixes = ['CXapprox','Gkeyll','JSFull','k','e','c']
# labels = ["Janev-Smith approximation","Meier (Gkeyll)","Janev-Smith","Krstic and Schultz","Krstic and Schultz","approx"]
# jMaxVector = [240, 240, 60,60,96,240]
# lengthVector = [40.0, 40.0, 40.0,40,40,40]
# basisVector = ['legendre', 'legendre', 'legendre','legendre','legendre','legendre']
# noutVector = [156, 400, 47,80,116,150]
# tVector = [-1,-1,40,40,70,-1]
# dxVector = [l / j for l, j in zip(lengthVector, jMaxVector)]
# colorVector = ['red','black','green','blue','orange','violet']
# lMaxVector = [2,2,2,2,2,2]


data_dict = {
    "CXapprox": ["Janev-Smith approximation", 240, 40.0, "legendre", 156, -1, "cyan", 2],
    "g": ["Meier (Gkeyll)", 240, 40.0, "legendre", 440, -1, "blue", 2],
    "JSFull": ["Janev-Smith", 60, 40.0, "legendre", 47, 40, "green", 2],
    "k": ["Krstic and Schultz", 60, 40.0, "legendre", 80, 40, "blue", 2],
    "e": ["Krstic and Schultz", 96, 40.0, "legendre", 116, -1, "orange", 2],
    "c": ["approx", 240, 40.0, "legendre", 150, -1, "violet", 2],
    "jss" : ["JS approx, nvx = 31", 240, 40, "legendre", 440, -1, "red", 2],
    "jsl" : ["JS approx, nvx = 127", 240, 40, "legendre", 440, -1, "cyan", 2],
    "jsf" : ["Krstic 7", 240, 40, "legendre", 168, -1, "purple", 2],
    "iz" : ["GUERNICA", 60, 40, "legendre", 100, -1, "red", 2]
}

reaction = 0 #0 for CX, 1 for ionization

plotting = 0 #0 for density, 1 for temperature

if reaction == 0:
    desiredPlot = {"jss","e","JSFull","g","jsl","jsf"}

    if plotting == 0:
        densityDegasData = np.loadtxt('d2-ndensity-cxonly.dat')
        positionDegas = densityDegasData[:,0]
        densityDegas = densityDegasData[:,1]
    if plotting == 1:
        temperatureDegasData = np.loadtxt('d2-ntemperature-cxonly.dat')
        positionDegas = temperatureDegasData[:,0]
        temperatureDegas = temperatureDegasData[:,1]

if reaction == 1:
    desiredPlot = {"iz"}

    if plotting == 0:
        densityDegasData = np.loadtxt('d2-ndensity-ionizonly.dat')
        positionDegas = densityDegasData[:,0]
        densityDegas = densityDegasData[:,1]
    if plotting == 1:
        temperatureDegasData = np.loadtxt('d2-ntemperature-ionizonly.dat')
        positionDegas = temperatureDegasData[:,0]
        temperatureDegas = temperatureDegasData[:,1]

# Initialize data structures
results = []

# for i, suffix in enumerate(suffixes):
for suffix, info in data_dict.items():
    if suffix not in desiredPlot:
        continue
    # Set parameters for this iteration
    jMax = info[1]
    length = info[2]
    basis = info[3]
    nout = info[4]+1
    dx = info[2]/info[1]
    lMax = info[7]+1
    color = info[6]
    label = info[0]

    # Read data files dynamically
    valuesDensity = pd.read_csv(f'Density{suffix}.csv', header=None)[0].to_numpy()
    valuesVelocityX = pd.read_csv(f'VelocityX{suffix}.csv', header=None)[0].to_numpy()
    valuesVelocityY = pd.read_csv(f'VelocityY{suffix}.csv', header=None)[0].to_numpy()
    valuesVelocityZ = pd.read_csv(f'VelocityZ{suffix}.csv', header=None)[0].to_numpy()
    valuesTemperature = pd.read_csv(f'Temperature{suffix}.csv', header=None)[0].to_numpy()

    # Initialize arrays
    rho = np.zeros((lMax, jMax, nout))
    rhouX = np.zeros((lMax, jMax, nout))
    rhouY = np.zeros((lMax, jMax, nout))
    rhouZ = np.zeros((lMax, jMax, nout))
    rt = np.zeros((lMax, jMax, nout))

    # Populate arrays
    k = 0
    for t in range(nout):
        for j in range(jMax):
            for l in range(lMax):
                rho[l][j][t] = valuesDensity[k]
                rhouX[l][j][t] = valuesVelocityX[k]
                rhouY[l][j][t] = valuesVelocityY[k]
                rhouZ[l][j][t] = valuesVelocityZ[k]
                rt[l][j][t] = valuesTemperature[k]
                k += 1

    # Evaluate and plot results for the current suffix
    t = info[5]  # Output step
    if (t>nout):
        t = -1
        print(label+" output at t = "+str(nout))
    y = np.zeros(10)
    x = np.zeros(10)
    for j in range(jMax):
        for idx in range(10):
            x[idx] = j * dx + idx * dx / 9.0
            density = 0
            velocityX = 0
            velocityY = 0
            velocityZ = 0
            temperature = 0
            for l in range(lMax):
                weight = getFunction(basis, l, (2 / dx) * (x[idx] - (j * dx + dx / 2)))
                density += rho[l][j][t] * weight
                velocityX += rhouX[l][j][t] * weight
                velocityY += rhouY[l][j][t] * weight
                velocityZ += rhouZ[l][j][t] * weight
                temperature += rt[l][j][t] * weight

            velocityX /= density
            velocityY /= density
            velocityZ /= density
            temperature = (temperature - density * (velocityX**2 + velocityY**2 + velocityZ**2)) / (3 * density)
            if plotting == 1:
                y[idx] = temperature
            if plotting == 0:
                y[idx] = density*1E18

        plt.plot(x - 20, y, color=color)
    
    plt.plot(0,0,color=color,label = f'{label}')

if reaction == 0:
    if plotting == 0:
        plt.plot(positionDegas,densityDegas,'k--',label='DEGAS2')
        ax.set_yscale('log')
        plt.ylim(5e13,2e19)
    if plotting == 1:
        plt.plot(positionDegas,temperatureDegas,'k--',label='DEGAS2')
        plt.ylim(28,70)
if reaction == 1:
    if plotting == 0:
        plt.plot(positionDegas,densityDegas,'k--',label='DEGAS2')
        ax.set_yscale('log')
        plt.ylim(1e12,2e19)
    if plotting == 1:
        plt.plot(positionDegas,temperatureDegas,'k--',label='DEGAS2')
        plt.ylim(7,19)

plt.xlim(-20,20)
plt.legend()
plt.show()