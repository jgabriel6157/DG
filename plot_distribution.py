import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

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

fileName = 'lastOutputjss.csv'
fileNameSol = 'lastOutputJS.csv'
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
    elif inputParam[0:3]=='nvy':
        nvy = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:3]=='nvz':
        nvz = int(inputParam[inputParam.index('=')+2:-1])
    elif inputParam[0:5]=='maxVX':
        domainMaxVX = assignFloat(inputParam[inputParam.index('=')+2:-1])
    if not inputParam:
        break
nout+=1
lMax+=1
inputFile.close()

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
# valuesSol = pd.read_csv(fileNameSol,header=None)
# valuesSol = valuesSol[0].to_numpy()
m = 0*20574000
dx = length/jMax
dvx = 2*domainMaxVX/(nvx-1)
# dvx = 1.0/nvx
u = np.zeros((lMax,jMax,nvx,nvy,nvz))
# uSol = np.zeros((lMax,jMax,nvx,nvy,nvz))

fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

for j in range(jMax):
    for kx in range(nvx):
        for ky in range(nvy):
            for kz in range(nvz):
                for lx in range(lMax):
                    u[lx,j,kx,ky,kz] = values[m]
                    # uSol[lx,j,kx,ky,kz] = valuesSol[m]
                    m = m+1
print(m)

# Select velocity direction
chosen_dim = "vx"  # Options: "vx", "vy", "vz"
    
x = 5
j = int(np.floor(x/dx))
xj = j*dx+dx/2

if chosen_dim == "vx":
    y = np.zeros(nvx)
    v = np.linspace(-domainMaxVX, domainMaxVX, nvx)  # Uniform grid for vx
    vy = 7
    vz = 7
    for vx in range(nvx):
        for l in range(lMax):
            y[vx] += u[l][j][vx][vy][vz]*getFunction(basis,l,(2.0/dx)*(x-xj))

elif chosen_dim == "vy":
    y = np.zeros(nvy)
    gh_points, _ = np.polynomial.hermite.hermgauss(nvy)
    gauss_hermite_points = gh_points * np.sqrt(60)
    v = gauss_hermite_points * np.sqrt(60)  # Scale GH points
    vx = 63
    vz = 7
    for vy in range(len(v)):
        for l in range(lMax):
            y[vy] += u[l][j][vx][vy][vz]*getFunction(basis,l,(2.0/dx)*(x-xj))

elif chosen_dim == "vz":
    y = np.zeros(nvz)
    gh_points, _ = np.polynomial.hermite.hermgauss(nvz)
    gauss_hermite_points = gh_points * np.sqrt(60)
    v = gauss_hermite_points * np.sqrt(60)  # Scale GH points
    vx = 63
    vy = 7
    for vz in range(len(v)):
        for l in range(lMax):
            y[vz] += u[l][j][vx][vy][vz]*getFunction(basis,l,(2.0/dx)*(x-xj))

# Maxwellian function
def maxwellian(v, rho, u, T):
    return rho / np.sqrt(2 * np.pi * T) * np.exp(- (v - u) ** 2 / (2 * T))

# Compute Moments
def fit_maxwellian():
    rho = np.trapz(y, v)
    u_mean = np.trapz(v * y, v) / rho
    T = np.trapz((v - u_mean) ** 2 * y, v) / rho

    # Fit Maxwellian to Data
    params, _ = curve_fit(lambda v, rho, u, T: maxwellian(v, rho, u, T), v, y, p0=[rho, u_mean, T])
    rho_fit, u_fit, T_fit = params
    y_maxwellian = maxwellian(v, rho_fit, u_fit, T_fit)

    return y_maxwellian, rho_fit, u_fit, T_fit

# Extract and process data
y_maxwellian, rho_fit, u_fit, T_fit = fit_maxwellian()

# Plot results
plt.plot(v, y, label="Numerical Distribution")
plt.plot(v, y_maxwellian, linestyle="--", label="Fitted Maxwellian")
plt.xlabel(f"{chosen_dim} (Velocity)")
plt.ylabel("Distribution Function")
plt.legend()
plt.show()