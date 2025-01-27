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
suffixes = ['CXapprox','Gkeyll','JSFull','k']
labels = ["Janev-Smith approximation","Meier (Gkeyll)","Janev-Smith","Krstic and Schultz"]
jMaxVector = [112, 112, 60,60]
lengthVector = [40.0, 40.0, 40.0,40]
basisVector = ['legendre', 'legendre', 'legendre','legendre']
noutVector = [100, 100, 47,80]
dxVector = [l / j for l, j in zip(lengthVector, jMaxVector)]
colorVector = ['red','black', 'green','blue']
lMaxVector = [2,2,2,2]

# Initialize data structures
results = []

for i, suffix in enumerate(suffixes):
    # Set parameters for this iteration
    jMax = jMaxVector[i]
    length = lengthVector[i]
    basis = basisVector[i]
    nout = noutVector[i]+1
    dx = dxVector[i]
    lMax = lMaxVector[i]+1
    color = colorVector[i]
    label = labels[i]

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
    t = -1  # Output step
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

            y[idx] = temperature
            # y[idx] = density*1E18

        plt.plot(x - 20, y, color=color)
    
    plt.plot(0,0,color=color,label = f'{label}')

# ax.set_yscale('log')
# plt.ylim(5e13,2e19)
plt.ylim(28,70)

plt.xlim(-20,20)
plt.legend()
plt.show()