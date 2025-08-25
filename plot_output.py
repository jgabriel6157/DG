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
plt.rcParams.update({'font.size': 32})
plt.rcParams['lines.linewidth'] = 4
fig = plt.figure(figsize=(12,12))
ax = fig.gca()

data_dict = {
    "CXapprox": ["Janev-Smith approximation", 240, 40.0, "legendre", 156, -1, "cyan", 2],
    "g": ["Meier", 240, 40.0, "legendre", 440, -1, "blue", 2],
    "JSFull": ["Janev-Smith", 60, 40.0, "legendre", 47, -1, "green", 2],
    "k": ["Krstic and Schultz", 60, 40.0, "legendre", 80, -1, "blue", 2],
    "e": ["GUERNICA", 96, 40.0, "legendre", 116, -1, "red", 2],
    "bad": ["LaBombard1", 60, 40.0, "legendre", 110, -1, "orange", 2],
    "badS": ["LaBombard2", 48, 40.0, "legendre", 88, -1, "orange", 2],
    "la": ["LaBombard avg", 48, 40.0, "legendre", 88, -1, "orange", 2],
    "lna": ["LaBombard not avg", 48, 40.0, "legendre", 88, -1, "orange", 2],
    "lbd": ["LaBombard", 96, 40.0, "legendre", 176, -1, "orange", 2], #nvy/nvz = 7
    "jss" : ["Averaging", 240, 40, "legendre", 440, -1, "green", 2], #nvx = 31
    "jsl" : ["Averaging approximation", 240, 40, "legendre", 440, -1, "green", 2], #nvx = 127
    "kf" : ["nvy/z=7", 240, 40, "legendre", 440, -1, "pink", 2], #nvy/z = 7
    "ks" : ["nvy/z=15", 60, 40, "legendre", 110, -1, "plum", 2], #nvy/z = 15
    "knc" : ["Low v res", 48, 40, "legendre", 88, -1, "red", 2], #nvy/z = 15
    "knc2" : ["High v res", 48, 40, "legendre", 88, -1, "deepskyblue", 2], #nvy/z = 20
    "kgh" : ["Gauss-Hermite", 48, 40, "legendre", 88, -1, "lightblue", 2], #nvy/z = 15
    "scaled" : ["KS, JS scaled", 96, 40, "legendre", 156, -1, "maroon", 2], #nvy/z = 15
    "scaledSmall" : ["KS, JS scaled", 48, 40, "legendre", 88, -1, "red", 2], #nvy/z = 7
    "ss" : ["under-resolved", 48, 40, "legendre", 88, -1, "maroon", 2], #nvy/z = 7
    "nss" : ["New", 48, 40, "legendre", 88, -1, "maroon", 2], #nvy/z = 7
    "nsl" : ["Full Operator", 192, 40, "legendre", 440, -1, "red", 2], #nvy/z = 7
    "nsv" : ["Full Operator v-res", 96, 40, "legendre", 141, -1, "maroon", 2], #nvy/z = 15
    "crec" : ["New2", 48, 40, "legendre", 88, -1, "red", 2], #nvy/z = 7, revised Crec
    "ses" : ["KS elastic, JS scaled", 48, 40, "legendre", 88, -1, "maroon", 2], #nvy/z = 7
    "a2" : ["Janev-Smith approx two", 60, 40.0, "legendre", 110, -1, "green", 2],
    "a3" : ["Non-averaging approximation", 60, 40.0, "legendre", 110, -1, "lime", 2],
    "iz" : ["GUERNICA", 60, 40, "legendre", 100, -1, "red", 2],
    "iz2" : ["GUERNICA", 60, 40, "legendre", 10, -1, "red", 2]
}

reaction = 0 #0 for CX, 1 for ionization

plotting = 1 #0 for density, 1 for temperature, 2 for scaled pressure

if reaction == 0:
    desiredPlot = {"nsl",""}

    if plotting == 0:
        densityDegasData = np.loadtxt('d2-ndensity-cxonly.dat')
        positionDegas = densityDegasData[:,0]
        densityDegas = densityDegasData[:,1]
    if plotting == 1:
        temperatureDegasData = np.loadtxt('d2-ntemperature-cxonly.dat')
        positionDegas = temperatureDegasData[:,0]
        temperatureDegas = temperatureDegasData[:,1]
    if plotting == 2:
        densityDegasData = np.loadtxt('d2-ndensity-cxonly.dat')
        temperatureDegasData = np.loadtxt('d2-ntemperature-cxonly.dat')
        positionDegas = temperatureDegasData[:,0]
        pressureDegas = densityDegasData[:,1]*temperatureDegasData[:,1]

if reaction == 1:
    desiredPlot = {"iz2"}

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
            if plotting == 2:
                y[idx] = temperature*density*1E18

        plt.plot(x - 20, y, color=color)
    
    plt.plot(0,0,color=color,label = f'{label}')

if reaction == 0:
    if plotting == 0:
        plt.plot(positionDegas,densityDegas,'k--',label='DEGAS2',dashes=(2,2))
        ax.set_yscale('log')
        # plt.ylim(5e13,2e19)
        plt.ylim(5e13,2e14)
        plt.ylabel(r'$n_n$ [m$^{-3}$]')
    if plotting == 1:
        plt.plot(positionDegas,temperatureDegas,'k--',label='DEGAS2',dashes=(2,2))
        plt.ylim(28,70)
        plt.ylabel(r'$T_n$ [eV]')
    if plotting == 2:
        plt.plot(positionDegas,pressureDegas,'k--',label='DEGAS2',dashes=(2,2))
        plt.ylim(5e14,2e20)
if reaction == 1:
    if plotting == 0:
        plt.plot(positionDegas,densityDegas,'k--',label='DEGAS2',dashes=(2,2))
        ax.set_yscale('log')
        plt.ylim(1e12,2e19)
        plt.ylabel(r'$n_n$ [m$^{-3}$]')
    if plotting == 1:
        plt.plot(positionDegas,temperatureDegas,'k--',label='DEGAS2',dashes=(2,2))
        plt.ylim(7,19)
        plt.ylabel(r'$T_n$ [eV]')

plt.xlabel(r'x [m]')
plt.xlim(-20,20)
plt.tight_layout()
# plt.legend()
plt.show()