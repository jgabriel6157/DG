import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import fnmatch
import numpy as np
from scipy.integrate import quad
from scipy.integrate import simpson
import os
import random
import time

n0 = 5*10**18 
T = 20
u = 0 
nu_cx = 2.2*10**-14
Tw = 2  
csL = np.sqrt(20)  
csR = -np.sqrt(20)  
Crec = 4.98*10**18  
Lz = 40

def integrand_1(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csL) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nu_cx * np.sqrt(1.672e-27 / 1.602e-19) * (z+Lz/2) / v)

def integrand_2(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csR) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nu_cx * np.sqrt(1.672e-27 / 1.602e-19) * (z - Lz/2) / v)

def n_minus(z):
    integral_1, _ = quad(integrand_1, 0, np.inf, args=(z,))
    return integral_1

def n_plus(z):
    integral_2, _ = quad(integrand_2, -np.inf, 0, args=(z,))
    return integral_2

def G(v):
    return np.exp(-((v - u) ** 2) / (2 * T))/np.sqrt(2 * np.pi * Tw * 9.58134*10**7)

def B_plus(z, z_prime):
    def integrand(v_prime):
        return (nu_cx / v_prime) * G(v_prime) * np.exp(nu_cx*n0 / v_prime * (z_prime - z))
    integral, _ = quad(integrand, 0, np.inf)
    return integral

def B_minus(z, z_prime):
    def integrand(v_prime):
        return (nu_cx / v_prime) * G(v_prime) * np.exp(nu_cx*n0 / v_prime * (z_prime - z))
    integral, _ = quad(integrand, -np.inf, 0)
    return integral

num_points = 1000
z_vals = np.linspace(-Lz/2, Lz/2, num_points)

fig = plt.figure()
ax = fig.gca()
ax.set_yscale('log')

n = np.zeros(num_points)
for i, z in enumerate(z_vals):
    if z <= 0:
        n[i] = n0 #* (np.square(1 / np.cosh(-(Lz/2 + z - 2) / 2)) + 1e-6)
    else:
        n[i] = n0 #* (np.square(1 / np.cosh((Lz/2 - z - 2) / 2)) + 1e-6)

tol = 1e-6
max_iter = 100
for iter in range(max_iter):
    print(iter)
    n_old = n.copy()

    # Update each n(z) value
    for i, z in enumerate(z_vals):
        # Calculate integrals
        if i==0:
            integral_plus = np.trapz([n[j] * B_plus(z, z_vals[j]) for j in range(i)], z_vals[:i])
            integral_minus = np.trapz([n[j] * B_minus(z, z_vals[j]) for j in range(i, num_points)], z_vals[i:])
        
        # Update n(z)
        n[i] = n_minus(z) + n_plus(z) + integral_plus - integral_minus

    # Check for convergence
    print(np.linalg.norm(n - n_old))
    if np.linalg.norm(n - n_old) < tol:
        print(f"Converged in {iter} iterations.")
        break
else:
    print("Did not converge within the maximum number of iterations.")

plt.plot(z_vals+Lz/2,n)
plt.show()