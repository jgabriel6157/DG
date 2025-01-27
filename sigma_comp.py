import numpy as np
import matplotlib.pyplot as plt


def janev_smith(E):
    # Constants
    A1 = 3.2345e-20  # m^2
    A2 = 235.88
    A3 = 0.038371
    A4 = 3.8068e-6
    A5 = 1.1832e-10
    A6 = 2.3713

    # Convert energy to keV
    E_keV = E / 1000

    E_keV = E_keV / 1.007276

    # Cross-section calculation
    sigma = A1 * np.log(A2 / E_keV + A6) / (
        1 + A3 * E_keV + A4 * E_keV**3.5 + A5 * E_keV**5.4
    )
    return sigma

def krstic_shultz(E):
    # Constants
    a0 = 0.160892e3
    a1 = -0.156336e2
    b1 = 0.108112e-1

    E/=1.007276

    # Cross-section calculation
    sigma = ((a0 + a1 * np.log(E)) / (1 + b1 * np.log(E))) * 2.80028e-17 * 1e-4
    return sigma

# Example usage
E_values = np.linspace(0.12, 2000, 100)  # Energy values in eV
sigma_js = [janev_smith(E) for E in E_values]
sigma_ks = [krstic_shultz(E) for E in E_values]

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

energy_eV_amu = np.array([

    1.20E-01, 2.00E-01, 5.00E-01, 1.00E+00, 2.00E+00, 5.00E+00, 1.00E+01,

    2.00E+01, 5.00E+01, 1.00E+02, 2.00E+02, 5.00E+02, 1.00E+03, 2.00E+03,

    5.00E+03, 1.00E+04, 2.00E+04, 5.00E+04, 1.00E+05, 2.00E+05, 5.00E+05,

    1.00E+06, 2.00E+06, 5.00E+06, 1.00E+07

])



cross_section_cm2 = np.array([

    4.96E-15, 4.70E-15, 4.33E-15, 4.10E-15, 3.83E-15, 3.46E-15, 3.17E-15,

    2.93E-15, 2.65E-15, 2.44E-15, 2.22E-15, 1.97E-15, 1.71E-15, 1.44E-15,

    1.10E-15, 7.75E-16, 4.45E-16, 9.93E-17, 1.01E-17, 6.09E-19, 6.03E-21,

    1.57E-22, 3.78E-24, 2.56E-26, 5.99E-28

])



# Convert cross-section from cm^2 to m^2 (1 cm^2 = 1e-4 m^2)

cross_section_m2 = cross_section_cm2 * 1e-4



# Data for spin exchange points

energy_cm_eV = np.array([

    0.1000, 0.1995, 0.5012, 1.0000, 1.9950, 5.0120, 10.0000, 19.9500, 50.1200, 100.0000

])



spin_exchange_au = np.array([

    1.997208E+02, 2.026141E+02, 1.720196E+02, 1.621056E+02, 1.492323E+02,

    1.327140E+02, 1.210393E+02, 1.096942E+02, 9.535670E+01, 8.510442E+01

])



# Convert spin exchange from a.u. to m^2 (1 a.u. = 2.80028e-17 cm^2 = 2.80028e-21 m^2)

spin_exchange_m2 = spin_exchange_au * 2.80028e-21



# Plot

plt.figure(figsize=(12, 8))

plt.plot(E_values, sigma_js, label="Janev-Smith 1993", linewidth=2, color='blue')

plt.plot(E_values, sigma_ks, label="Krstic-Schultz 1998", linewidth=2, linestyle="--", color='orange')

plt.scatter(energy_eV_amu, cross_section_m2, label="Janev-Smith Data", color='blue', marker='o')

plt.scatter(energy_cm_eV, spin_exchange_m2, label="Krstic Data", color='orange', marker='x')



# Labels and formatting

plt.xscale("log")

# plt.yscale("log")

plt.xlabel("Energy (eV)", fontsize=14)

plt.ylabel("Cross-section (m²)", fontsize=14)

plt.xlim((0.1, 2000))
plt.ylim((1E-19,7E-19))

plt.title("Comparison of Cross-section Models with Additional Data", fontsize=16)

plt.legend(fontsize=12)

plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()

plt.show()
