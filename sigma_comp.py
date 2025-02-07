# import numpy as np
# import matplotlib.pyplot as plt


# def janev_smith(E):
#     # Constants
#     A1 = 3.2345e-20  # m^2
#     A2 = 235.88
#     A3 = 0.038371
#     A4 = 3.8068e-6
#     A5 = 1.1832e-10
#     A6 = 2.3713

#     # Convert energy to keV
#     E_keV = E / 1000

#     E_keV = E_keV / 1.007276

#     # Cross-section calculation
#     sigma = A1 * np.log(A2 / E_keV + A6) / (
#         1 + A3 * E_keV + A4 * E_keV**3.5 + A5 * E_keV**5.4
#     )
#     return sigma

# def krstic_shultz(E):
#     # Constants
#     a0 = 0.160892e3
#     a1 = -0.156336e2
#     b1 = 0.108112e-1

#     E/=1.007276

#     # Cross-section calculation
#     sigma = ((a0 + a1 * np.log(E)) / (1 + b1 * np.log(E))) * 2.80028e-17 * 1e-4
#     return sigma

# # Example usage
# E_values = np.linspace(0.12, 2000, 100)  # Energy values in eV
# sigma_js = [janev_smith(E) for E in E_values]
# sigma_ks = [krstic_shultz(E) for E in E_values]

# plt.figure(figsize=(10, 6))

# plt.plot(E_values, sigma_js, label="Janev-Smith 1993", linewidth=2)

# plt.plot(E_values, sigma_ks, label="Krstic-Schultz 1998", linewidth=2, linestyle="--")

# plt.xscale("log")

# # plt.yscale("log")

# plt.xlabel("Energy (eV)", fontsize=14)

# plt.ylabel("Cross-section (m²)", fontsize=14)

# plt.title("Comparison of Cross-section Models", fontsize=16)

# plt.legend(fontsize=12)

# plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# plt.tight_layout()

# plt.show()

# Data for additional points

# energy_eV_amu = np.array([

#     1.20E-01, 2.00E-01, 5.00E-01, 1.00E+00, 2.00E+00, 5.00E+00, 1.00E+01,

#     2.00E+01, 5.00E+01, 1.00E+02, 2.00E+02, 5.00E+02, 1.00E+03, 2.00E+03,

#     5.00E+03, 1.00E+04, 2.00E+04, 5.00E+04, 1.00E+05, 2.00E+05, 5.00E+05,

#     1.00E+06, 2.00E+06, 5.00E+06, 1.00E+07

# ])



# cross_section_cm2 = np.array([

#     4.96E-15, 4.70E-15, 4.33E-15, 4.10E-15, 3.83E-15, 3.46E-15, 3.17E-15,

#     2.93E-15, 2.65E-15, 2.44E-15, 2.22E-15, 1.97E-15, 1.71E-15, 1.44E-15,

#     1.10E-15, 7.75E-16, 4.45E-16, 9.93E-17, 1.01E-17, 6.09E-19, 6.03E-21,

#     1.57E-22, 3.78E-24, 2.56E-26, 5.99E-28

# ])



# # Convert cross-section from cm^2 to m^2 (1 cm^2 = 1e-4 m^2)

# cross_section_m2 = cross_section_cm2 * 1e-4



# # Data for spin exchange points

# energy_cm_eV = np.array([

#     0.1000, 0.1995, 0.5012, 1.0000, 1.9950, 5.0120, 10.0000, 19.9500, 50.1200, 100.0000

# ])



# spin_exchange_au = np.array([

#     1.997208E+02, 2.026141E+02, 1.720196E+02, 1.621056E+02, 1.492323E+02,

#     1.327140E+02, 1.210393E+02, 1.096942E+02, 9.535670E+01, 8.510442E+01

# ])



# # Convert spin exchange from a.u. to m^2 (1 a.u. = 2.80028e-17 cm^2 = 2.80028e-21 m^2)

# spin_exchange_m2 = spin_exchange_au * 2.80028e-21



# # Plot

# plt.figure(figsize=(12, 8))

# plt.plot(E_values, sigma_js, label="Janev-Smith 1993", linewidth=2, color='blue')

# plt.plot(E_values, sigma_ks, label="Krstic-Schultz 1998", linewidth=2, linestyle="--", color='orange')

# plt.scatter(energy_eV_amu, cross_section_m2, label="Janev-Smith Data", color='blue', marker='o')

# plt.scatter(energy_cm_eV, spin_exchange_m2, label="Krstic Data", color='orange', marker='x')



# # Labels and formatting

# plt.xscale("log")

# # plt.yscale("log")

# plt.xlabel("Energy (eV)", fontsize=14)

# plt.ylabel("Cross-section (m²)", fontsize=14)

# plt.xlim((0.1, 2000))
# plt.ylim((1E-19,7E-19))

# plt.title("Comparison of Cross-section Models with Additional Data", fontsize=16)

# plt.legend(fontsize=12)

# plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# plt.tight_layout()

# plt.show()

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

#   Returns maxwellian averaged <sigma V) for charge exchange of atomic 
#       hydrogen. Coefficients are taken
#       from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
#       Springer-Verlag, 1987, p.272.

def sigmav_cx_h0(T,E):

    #	Input:
	#		T	- List, np.array(*) or float, ion [neutral] temperature (eV)
	#		E	- List, np.array(*) or float, neutral [ion] mono-energy (eV)

	#	Output:
	#		returns <sigma V> for 0.1 < Te < 2e4 and 0.1 < E < 2e4 
	#		Units: m^3/s

    T=np.array(T) # Converts T and E to np.array if not already
    E=np.array(E)

    if T.size!=E.size:
	    raise Exception("number of elements of T and E are different!")

    alpha=np.zeros((9,9))

    alpha[:,0:3]=[[-1.829079581680e+01,	1.640252721210e-01,	3.364564509137e-02],
      [ 2.169137615703e-01,    -1.106722014459e-01,    -1.382158680424e-03],
      [ 4.307131243894e-02,		8.948693624917e-03,    -1.209480567154e-02],
      [-5.754895093075e-04,     6.062141761233e-03,     1.075907881928e-03],
      [-1.552077120204e-03,    -1.210431587568e-03,     8.297212635856e-04],
      [-1.876800283030e-04,    -4.052878751584e-05,    -1.907025662962e-04],
      [ 1.125490270962e-04,     2.875900435985e-05,     1.338839628570e-05],
      [-1.238982763007e-05,    -2.616998139678e-06,    -1.171762874107e-07],
      [ 4.163596197181e-07,     7.558092849125e-08,    -1.328404104165e-08]]

    alpha[:,3:6]=[[ 9.530225559189e-03,    -8.519413589968e-04,    -1.247583860943e-03],
      [ 7.348786286628e-03,    -6.343059502294e-04,    -1.919569450380e-04],
      [-3.675019470470e-04,     1.039643390686e-03,    -1.553840717902e-04],
      [-8.119301728339e-04,     8.911036876068e-06,     3.175388949811e-05],
      [ 1.361661816974e-04,    -1.008928628425e-04,     1.080693990468e-05],
      [ 1.141663041636e-05,     1.775681984457e-05,    -3.149286923815e-06],
      [-4.340802793033e-06,    -7.003521917385e-07,     2.318308730487e-07],
      [ 3.517971869029e-07,    -4.928692832866e-08,     1.756388998863e-10],
      [-9.170850253981e-09,     3.208853883734e-09,    -3.952740758950e-10]]

    alpha[:,6:9]=[[ 3.014307545716e-04,    -2.499323170044e-05,     6.932627237765e-07],
      [ 4.075019351738e-05,    -2.850044983009e-06,     6.966822400446e-08],
      [ 2.670827249272e-06,     7.695300597935e-07,    -3.783302281524e-08],
      [-4.515123641755e-06,     2.187439283954e-07,    -2.911233951880e-09],
      [ 5.106059413591e-07,    -1.299275586093e-07,     5.117133050290e-09],
      [ 3.105491554749e-08,     2.274394089017e-08,    -1.130988250912e-09],
      [-6.030983538280e-09,    -1.755944926274e-09,     1.005189187279e-10],
      [-1.446756795654e-10,     7.143183138281e-11,    -3.989884105603e-12],
      [ 2.739558475782e-11,    -1.693040208927e-12,     6.388219930167e-14]]

    #   Limits values to >= 0.1 and <= 2.01e4

    E2=np.maximum(E,.1)
    E2=np.minimum(E2,2.01e4)
    T2=np.maximum(T,.1)
    T2=np.minimum(T2,2.01e4)

    alogE=np.log(E2)
    alogT=np.log(T2)

    result=np.zeros(E2.shape)
    for i in range(9):
        for j in range(9):
            result = result+alpha[j,i]*alogE**i*alogT**j
    
    return np.e**result

# Constants
m_H = 1.6735575e-27  # kg (mass of hydrogen atom)
eV_to_J = 1.60218e-19  # Conversion from eV to Joules

# Janev-Smith cross-section function
def sigma_JS(E):
    E_min, E_max = 0.12, 4e5
    A1, A2, A3, A4, A5, A6 = 3.2345e-20, 235.88, 0.038371, 3.8068e-6, 1.1832e-10, 2.3713
    
    E = np.clip(E, E_min, E_max) / 1000  # Convert eV to keV
    return A1 * np.log(A2 / E + A6) / (1 + A3 * E + A4 * E**3.5 + A5 * E**5.4)  # m^2

# Maxwellian distribution in 3D
def f_M(vp, T):
    v_th = np.sqrt(2 * T * eV_to_J / m_H)  # Thermal velocity
    coeff = (m_H / (2 * np.pi * T * eV_to_J)) ** (3/2)
    return coeff * np.exp(-m_H * np.dot(vp, vp) / (2 * T * eV_to_J))

# Integrand function
def integrand(vx, vy, vz, v, T):
    vp = np.array([vx, vy, vz])
    v_rel = np.linalg.norm(v - vp)
    E_rel = 0.5 * m_H * v_rel**2 / eV_to_J  # Convert to eV
    return f_M(vp, T) * sigma_JS(E_rel) * v_rel

# Perform numerical integration
def compute_integral(E, T):
    v = np.sqrt(2 * E * eV_to_J / m_H) * np.array([1, 0, 0])  # Convert energy to velocity
    bounds = [[-5e4, 5e4], [-5e4, 5e4], [-5e4, 5e4]]  # Integration limits in m/s
    result, error = spi.nquad(integrand, bounds, args=(v, T))
    return result

# Example usage
# E_test = 10  # eV
# T_test = 60  # eV
# result = compute_integral(E_test, T_test)
# print("Integral result:", result)
# print(sigmav_cx_h0(T_test,E_test))
# Define energy range
E_values = np.logspace(-1, 3.3, 20)  # Energy from 0.1 eV to 2000 eV (log scale)
T_test = 60  # eV

# Compute integral for each energy
integral_values = [compute_integral(E, T_test)*1e6 for E in E_values]
averaged_values = [sigmav_cx_h0(E, T_test)/10 for E in E_values]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(E_values, integral_values, marker='o', linestyle='-', label = 'integral')
plt.plot(E_values, averaged_values, marker='o', linestyle='-', label = 'averaged')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy (eV)")
plt.ylabel("<sigma v> (m^3/s)")
plt.legend()
plt.grid(True)
plt.show()

# for E_value in [0.1,1,10,100,1000,10000]:
#     T_values = np.logspace(-1,4.3,200)
#     sigmas = [sigmav_cx_h0(E_value, T) for T in T_values]
#     plt.plot(T_values, sigmas, label = str(E_value))
#     plt.xscale("log")
#     plt.yscale("log")
# plt.ylim(1e-9,2e-6)
# plt.xlim(0.1,2e4)
# plt.legend()
# plt.show()