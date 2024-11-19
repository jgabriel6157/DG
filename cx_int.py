import numpy as np
from scipy.integrate import dblquad
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline

# Constants
Lz = 40
T = 20
n0 = 5 * 10**18
Crec = 4.98 * 10**18
nuCX = (2.2 * 10**-14) * np.sqrt(1.672 * 10**-27 / (1.602 * 10**-19))
csL = np.sqrt(20)
csR = -np.sqrt(20)
Tw = 10
errorVal = 1e-8

def integrand_1(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csL) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nuCX * (z+Lz/2) / v)

def integrand_2(v, z):
    return (1 / np.sqrt(2 * np.pi * Tw)) * Crec * np.exp(-((v - csR) ** 2) / (2 * Tw)) * \
           np.exp(-n0 * nuCX * (z - Lz/2) / v)

def n_minus(z):
    integral_1, _ = quad(integrand_1, 0, np.inf, args=(z,),epsabs=errorVal,epsrel=errorVal)
    return integral_1

def n_plus(z):
    integral_2, _ = quad(integrand_2, -np.inf, 0, args=(z,),epsabs=errorVal,epsrel=errorVal)
    return integral_2

# Define the inner integrand function in terms of vp and zp
def inner_integrand(vp, zp):
    factor = nInterp(zp)*(nuCX * n0) / np.sqrt(2 * np.pi * T)
    exponent1 = -((vp - np.sqrt(T) * np.sign(zp))**2) / (2 * T)
    exponent2 = (n0 * nuCX * (zp - z)) / vp
    return factor * (1 / vp) * np.exp(exponent1) * np.exp(exponent2)

num_points = 73 #jMax*(res-1)+1
z_vals = np.linspace(-Lz/2, Lz/2, num_points)

nPoints = np.zeros(num_points)
for i, z in enumerate(z_vals):
    # nPoints[i] = n0 * (np.square(1 / np.cosh(-(Lz/2 - z*np.sign(z) - 2) / 2)) + 1e-6)
    nPoints[i] = n_minus(z)+n_plus(z)

# nInterp = interp1d(z_vals,nPoints, kind='cubic')
# nInterp = make_interp_spline(z_vals,nPoints,k=3)

tol = 1e6
maxIter = 100

for iter in range(maxIter):
    print(iter)
    nOld = nPoints.copy()
    nInterp = make_interp_spline(z_vals,nPoints,k=3)

    for i,z in enumerate(z_vals):
        integralPlus,_ = dblquad(inner_integrand, -Lz/2, z, lambda _: 0, lambda _: np.inf)
        integralMinus,_ = dblquad(inner_integrand, z, Lz/2, lambda _: -np.inf, lambda _: 0)

        nPoints[i] = n_minus(z)+n_plus(z)+integralPlus-integralMinus
        print(nPoints[i])
    
    print(np.linalg.norm(nPoints - nOld))
    if np.linalg.norm(nPoints - nOld) < tol:
        print(f"Converged in {iter} iterations.")
        break
else:
    print("Did not converge within the maximum number of iterations.")
    

# # Outer integral over zp with the inner integral as the function
# result, error = dblquad(inner_integrand, -Lz/2, zed, lambda _: 0, lambda _: np.inf)
# print(result) 
# print(error)

# result2, error2 = dblquad(inner_integrand, zed, Lz/2, lambda _: -np.inf, lambda _: 0)
# print(result2)
# print(error2)
