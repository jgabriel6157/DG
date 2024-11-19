import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import fnmatch
import numpy as np
from scipy.integrate import quad
from scipy import integrate
import os
import random
import time

n0 = 5*1e18
T = 20
u = np.sqrt(20)
nu_cx = (2.2*10**-14) * np.sqrt(1.672e-27 / 1.602e-19)
Tw = 2
csL = np.sqrt(20)
csR = -np.sqrt(20)
Crec = 4.98*1e18
Lz = 40
errorVal = 1e-8

def integrand_1(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csL) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nu_cx * (z+Lz/2) / v)

def integrand_2(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csR) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nu_cx * (z - Lz/2) / v)

def n_minus(z):
    integral_1, _ = quad(integrand_1, 0, np.inf, args=(z,),epsabs=errorVal,epsrel=errorVal)
    return integral_1

def n_plus(z):
    integral_2, _ = quad(integrand_2, -np.inf, 0, args=(z,),epsabs=errorVal,epsrel=errorVal)
    return integral_2

def G(v, z): #really G*n_i(z)
    if z < 0:
        ug = -u
    else:
        ug = u
    return np.exp(-((v - ug) ** 2) / (2 * T))/np.sqrt(2 * np.pi * T)

def B_plus(z, z_prime):
    def integrand(v_prime):
        return (nu_cx*n0 / v_prime) * G(v_prime, z_prime) * np.exp(nu_cx*n0 / v_prime * (z_prime - z))
    integral, _ = quad(integrand, 0, np.inf,epsabs=errorVal,epsrel=errorVal)
    return integral

def B_minus(z, z_prime):
    def integrand(v_prime):
        return (nu_cx*n0 / v_prime) * G(v_prime, z_prime) * np.exp(nu_cx*n0 / v_prime * (z_prime - z))
    integral, _ = quad(integrand, -np.inf, 0,epsabs=errorVal,epsrel=errorVal)
    return integral

num_points = 1009 #jMax*(res-1)+1
z_vals = np.linspace(-Lz/2, Lz/2, num_points)

fig = plt.figure()
ax = fig.gca()
ax.set_yscale('log')

n = np.zeros(num_points)
for i, z in enumerate(z_vals):
    # if z <= 0:
    #     n[i] = n0 * (np.square(1 / np.cosh(-(Lz/2 + z - 2) / 2)) + 1e-6)
    # else:
    #     n[i] = n0 * (np.square(1 / np.cosh((Lz/2 - z - 2) / 2)) + 1e-6)
    n[i] = n_minus(z)+n_plus(z)

# i = num_points-1
# z = 20.0
# integral_plus = np.trapz([n[j] * B_plus(z, z_vals[j]) for j in range(i)], x=z_vals[:i])
# integral_minus = np.trapz([n[j] * B_minus(z, z_vals[j]) for j in range(i, num_points)], z_vals[i:])
# print(B_plus(z,z))
# print(integral_plus)

# tol = 1e-6
# max_iter = 100
# for iter in range(max_iter):
#     print(iter)
#     n_old = n.copy()

#     # Update each n(z) value
#     for i, z in enumerate(z_vals):
#         # Calculate integrals
#         # if i==0:
#         integral_plus = np.trapz([n_old[j] * B_plus(z, z_vals[j]) for j in range(i)], z_vals[:i])
#         integral_minus = np.trapz([n_old[j] * B_minus(z, z_vals[j]) for j in range(i, num_points)], z_vals[i:])
        
#         # Update n(z)
#         n[i] = n_minus(z) + n_plus(z) + integral_plus - integral_minus
#         # print(i)
#         # print(z)
#         # print(integral_plus)
#         # print(integral_minus)
#         print(n[i])

#     # Check for convergence
#     print(np.linalg.norm(n - n_old))
#     if np.linalg.norm(n - n_old) < tol:
#         print(f"Converged in {iter} iterations.")
#         break
# else:
#     print("Did not converge within the maximum number of iterations.")

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

fileName = 'Density.csv'
fileNameSol = 'Density1.csv'
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

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
# valuesSol = pd.read_csv(fileNameSol,header=None)
# valuesSol = valuesSol[0].to_numpy()
# m = 33600
zvals = np.linspace(0,Lz,73)
#tw = 2
sol = [9.530037918593606e+18,3.443585377795066e+18,1.1798778383970465e+18,4.230035393563359e+17,1.582902613583898e+17,6.149948822451518e+16,2.472530326865515e+16,
       1.0255050148172988e+16,4380327202107013.0,1923718648188707.2,867818305504052.6,401746227850319.7,190725403969430.44,92780866066238.86,46213707070264.734,
       23548111437137.75,12262806802297.824,6519188821336.82,3533903271155.998,1950887233760.4268,1095407474005.471,624794875291.5898,361564797707.91644,
       212038887327.73297,125878241642.13467,75570061369.60495,45835537586.36156,28062240135.22628,17327437665.918846,10780669536.719915,6751335947.825615,
       4249109178.072075,2680495390.223425,1685151768.5609396,1043217405.1694884,602606475.811542,302962114.23153406,602606475.8111047,1043217405.1681938,
       1685151768.5584273,2680495390.218846,4249109178.0642195,6751335947.81239,10780669536.698355,17327437665.88374,28062240135.16903,45835537586.2684,
       75570061369.4519,125878241641.89076,212038887327.32016,361564797707.2179,624794875290.4581,1095407474003.4694,1950887233756.92,3533903271149.823,
       6519188821325.937,12262806802278.195,23548111437110.336,46213707070202.67,92780866066119.61,190725403969180.28,401746227849867.94,867818305503060.1,
       1923718648186979.5,4380327202101790.5,1.0255050148166076e+16,2.4725303268610548e+16,6.149948822447443e+16,1.582902613581344e+17,4.230035393560502e+17,
       1.1798778383881987e+18,3.4435853777890196e+18,9.530037918590149e+18]
#tw = 10
sol = [9.271348489836261e+18,3.8754428088949514e+18,1.6932338616168102e+18,7.754379043645197e+17,3.703097076228575e+17,1.822157351491099e+17,
       9.21811376675369e+16,4.767955433672934e+16,2.5171769786945176e+16,1.3524586186453064e+16,7384923039433558.0,4090955994446276.0,2296485632385380.5,
       1304838202117868.2,749739790080419.5,435267778721258.56,255142715477517.84,150906047778989.8,90007234113567.34,54108963673047.09,32770176667558.055,
       19985768033112.9,12269481498689.963,7579465572786.979,4709903064214.832,2943114536287.073,1848788332094.2627,1167104784356.1567,740140862116.5347,
       471297430651.7465,301123988646.1213,192810267953.17917,123423571847.6819,78543255380.94543,49099902858.73616,28561938896.348938,14404363732.721174,
       28561938896.348713,49099902858.74805,78543255380.95578,123423571847.68909,192810267953.1872,301123988646.12506,471297430651.776,740140862116.5443,
       1167104784356.1626,1848788332094.2788,2943114536287.0835,4709903064215.169,7579465572787.246,12269481498690.16,19985768033114.184,32770176667558.723,
       54108963673046.71,90007234113566.84,150906047778987.78,255142715477515.3,435267778721187.4,749739790080367.8,1304838202117811.5,2296485632385492.0,
       4090955994446296.5,7384923039433890.0,1.3524586186453012e+16,2.5171769786939908e+16,4.76795543367135e+16,9.21811376674128e+16,1.822157351490363e+17,
       3.703097076227858e+17,7.754379043650042e+17,1.6932338616167905e+18,3.875442808894982e+18,9.271348489836345e+18]

for m in [336*50]:
    dx = length/jMax
    dvx = 2*domainMaxVX/(nvx-1)
    # dvx = 1.0/nvx
    u = np.zeros((lMax,jMax))
    # uSol = np.zeros((lMax,jMax))

    # fig = plt.figure()
    # ax = fig.gca()
    # ax.set_yscale('log')

    for j in range(jMax):
        for lx in range(lMax):
            u[lx,j] = values[m]
            # uSol[lx,j] = valuesSol[m]
            m = m+1
    print(m)

    l2Norm = 0
    solution = 0

    res = 10
    for j in range(jMax):
        xj = j*dx+dx/2
        y = np.zeros(res)
        x = np.zeros(res)
        # sol = np.zeros(res)
        error = np.zeros(res)
        y2 = np.zeros(res)
        for i in range(res):
            x[i] = j*dx+i*dx/(res-1)
            # sol[i] = n[j*(res-1)+i]/(1e18)
            for l in range(lMax):
                y[i] += u[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
                # sol[i] += uSol[l][j]*getFunction(basis,l,(2.0/dx)*(x[i]-xj))
            # error[i] = (y[i]-sol[i])**2
            # y2[i] = y[i]**2
        # l2Norm+=integrate.simpson(error,x=x)
        # solution+=integrate.simpson(y2,x=x)
        plt.plot(x,y*1e18,color='red')
        # plt.plot(x,sol,'k:')

# print(np.sqrt(l2Norm/solution))
plt.ylim(1e8,1e19)
plt.plot(zvals,sol,'k:')
# plt.plot(z_vals+Lz/2,n/(1e18), 'k:')
plt.show()