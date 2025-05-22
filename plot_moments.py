import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fileName = 'MomentsPeriodic.csv'

values = pd.read_csv(fileName,header=None)
values = values[0].to_numpy()
mass = np.zeros(int(len(values)/4))
momentum = np.zeros(int(len(values)/4))
energy = np.zeros(int(len(values)/4))
entropy = np.zeros(int(len(values)/4))
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

# # Create a figure and subplots
# fig, axs = plt.subplots(2, 2)

# # Plot each array in its own subplot
# axs[0, 0].plot(mass)
# axs[0, 0].set_title('(Mass(t)-Mass(0))/Mass(0)')
# axs[0, 1].plot(momentum)
# axs[0, 1].set_title('Momentum(t)')
# axs[1, 0].plot(energy)
# axs[1, 0].set_title('(Energy(t)-Energy(0))/Energy(0)')
# axs[1, 1].plot(entropy)
# axs[1, 1].set_title('(Entropy(t)-Entropy(0))/Entropy(0)')

plt.rcParams.update({'font.size': 32})
plt.rcParams['lines.linewidth'] = 4
fig = plt.figure(figsize=(12,12))
ax = plt.gca()

plt.plot(abs(mass),color='b',label=r'$\Delta$Mass')
plt.plot(abs(momentum),color='green',label=r'$\Delta$Momentum')
plt.plot(abs(energy),color='orange',label=r'$\Delta$Energy')

ax.set_yscale('log')

# Adjust layout
plt.tight_layout()

plt.legend()

# Show the plot
plt.show()