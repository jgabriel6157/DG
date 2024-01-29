import pandas as pd
import matplotlib.pyplot as plt
import fnmatch
import numpy as np
import os
import random

fileName = 'Output.csv'

def LegendreP(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    else:
        return ((2.0*n-1.0)*x*LegendreP(n-1,x)-(n-1)*LegendreP(n-2,x))/n
    
def LegendrePorthonormal(n,x):
    return np.sqrt((2.0*n+1.0)/2.0)*LegendreP(n,x)

values = pd.read_csv('/home/jack/Documents/c++/Output.csv',header=None)
values = values[0].to_numpy()
k = 0
jMax = 32
dx = 2*np.pi/jMax
lMax = 2
u = np.zeros((lMax,jMax))
for j in range(jMax):
    for l in range(lMax):
        u[l][j] = values[k]
        k=k+1


for j in range(jMax):
    y = np.zeros(10)
    sol = np.zeros(10)
    x = np.zeros(10)
    for i in range(10):
        x[i] = j*dx+i*dx/9.0
        for l in range(lMax):
            y[i] = y[i] + u[l][j]*LegendrePorthonormal(l,(2/dx)*(x[i]-(j*dx+dx/2)))
        sol[i] = np.sin(x[i])
    plt.plot(x,y,color='red')
    plt.plot(x,sol,color='k',linestyle='--')

plt.show()