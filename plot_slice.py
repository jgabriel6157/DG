import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

plt.rcParams.update({'font.size': 24})
plt.rcParams['lines.linewidth'] = 4
fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(projection='3d')

data_dict = {
    "g": ["Meier", 240, 40.0, "legendre", 440, -1, "blue", 2,127,15,15,31],
    "badS": ["LaBombard", 48, 40.0, "legendre", 88, -1, "orange", 2,31,15,15,31],
    "jss" : ["Averaging", 240, 40, "legendre", 440, -1, "green", 2,31,15,15,31], #nvx = 31
    "knc2" : ["Krstic & Schultz", 48, 40, "legendre", 88, -1, "red", 2,31,20,20,31], #nvy/z = 20
    "scaled" : ["KS, JS scaled", 96, 40, "legendre", 156, -1, "red", 2,31,15,15,31], #nvy/z = 15
    "scaledSmall" : ["KS, JS scaled", 48, 40, "legendre", 88, -1, "maroon", 2,31,7,7,31], #nvy/z = 7
    "iz" : ["GUERNICA", 60, 40, "legendre", 100, -1, "red", 2,31,15,15,31]
}

desiredPlot = {"knc2","jss","badS","g"}

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
    nvx = info[8]
    nvy = info[9]
    nvz = info[10]
    domainMaxVX = info[11]

    values = pd.read_csv(f'lastOutput{suffix}.csv', header=None)[0].to_numpy()

    m = 0
    dx = length/jMax
    dvx = 2*domainMaxVX/(nvx-1)
    u = np.zeros((lMax,jMax,nvx,nvy,nvz))

    for j in range(jMax):
        for kx in range(nvx):
            for ky in range(nvy):
                for kz in range(nvz):
                    for lx in range(lMax):
                        u[lx,j,kx,ky,kz] = values[m]
                        m = m+1

    vz = 7
    vy = 7
    x = 4
    j = int(np.floor(x/dx))
    xj = j*dx+dx/2
    y = np.zeros(nvx)
    v = np.zeros(nvx)
    for vx in range(nvx):
        v[vx] = -domainMaxVX + vx*dvx

        for l in range(lMax):
            y[vx] += u[l][j][vx][vy][vz]*getFunction(basis,l,(2.0/dx)*(x-xj))
    
    plt.plot(v, y/max(y),color = color)

# Compute Moments
rho = np.trapz(y, v)  # Density
u_mean = np.trapz(v * y, v) / rho  # Mean velocity
T = np.trapz((v) ** 2 * y, v) / rho  # Temperature (assuming unit mass)

T = T-(u_mean*u_mean)

print(rho)
print(u_mean)
print(T)

# Maxwellian distribution function
def maxwellian(v, rho, u, T):
    return rho / np.sqrt(2 * np.pi * T) * np.exp(- (v - u) ** 2 / (2 * T))

# Fit Maxwellian to Data
params, _ = curve_fit(lambda v, rho, u, T: maxwellian(v, rho, u, T), v, y, p0=[rho, u_mean, T])
rho_fit, u_fit, T_fit = params

# Generate Maxwellian Fit
y_maxwellian = maxwellian(v, rho_fit, u_fit, T_fit)

ion_maxwellian = maxwellian(v, rho_fit,-1.3, 60)

# Plot Original vs Maxwellian

# plt.plot(v, y_maxwellian, color='green',linestyle="--", label="Fitted Maxwellian",dashes=(2,2))
plt.plot(v, ion_maxwellian/max(ion_maxwellian), color='gray',linestyle=":", label="T=60 eV Distribution")
plt.xlabel("Velocity")
plt.ylabel("Normalized Distribution Function")
plt.legend()
plt.show()



# plt.plot(v,y)
# plt.show()