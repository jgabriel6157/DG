import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fileName = 'Moments.csv'
fileNameCons = 'MomentsCons.csv'

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
mass = np.zeros(int(len(values)/4))
momentum = np.zeros(int(len(values)/4))
energy = np.zeros(int(len(values)/4))
entropy = np.zeros(int(len(values)/4))

values2 = pd.read_csv(fileNameCons,header=None)
values2 = values2[0].to_numpy()
mass2 = np.zeros(int(len(values2)/4))
momentum2 = np.zeros(int(len(values2)/4))
energy2 = np.zeros(int(len(values2)/4))
entropy2 = np.zeros(int(len(values2)/4))

count = 0
counter = 0
for value in values:
    if count == 0:
        mass[counter] = float(value)
    elif count == 1:
        momentum[counter] = float(value)
    elif count == 2:
        energy[counter] = float(value)
    elif count == 3:
        entropy[counter] = float(value)
        counter+=1
        count = -1
    count+=1

count = 0
counter = 0
for value2 in values2:
    if count == 0:
        mass2[counter] = float(value2)
    elif count == 1:
        momentum2[counter] = float(value2)
    elif count == 2:
        energy2[counter] = float(value2)
    elif count == 3:
        entropy2[counter] = float(value2)
        counter+=1
        count = -1
    count+=1

plt.plot(abs(mass2),'r-')
plt.plot(abs(energy2),'r--')

plt.plot(abs(mass),'b-')
plt.plot(abs(energy),'b--')
plt.yscale('log')

# plt.plot(entropy2,'r-')
# plt.plot(entropy,'b--')

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()