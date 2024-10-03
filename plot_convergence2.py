import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

x = [4,8,12,16,20,24,28,32]

# Data for each order
# Order 1
t1 = [0.00417163, 0.01456, 0.0322, 0.055, 0.085, 0.125, 0.1647, 0.213307]
y1 = [0.294571, 0.057, 0.0220482, 0.0115, 0.007, 0.004775, 0.00345, 0.00261685]

# Order 2
t2 = [0.0117978, 0.04304, 0.095, 0.166077, 0.258, 0.369, 0.509, 0.651782]
y2 = [0.0234314, 0.00286, 0.00085, 0.000358, 0.00018, 0.000106, 6.69e-5, 4.482e-5]

# Order 4
t4 = [0.05, 0.188867, 0.42368, 0.736451, 1.14436, 1.68, 2.24874, 3.01426]
y4 = [0.000206, 6.58286e-6, 8.75e-7, 2.09484e-7, 6.937e-8, 2.81658e-8, 1.333e-8, 7.06298e-9]

# Order 8
t8 = [0.287074, 1.07451, 2.4, 4.25398, 6.53409, 9.44362, 12.9971, 16.8101]
y8 = [2.08436e-7, 2.60434e-8, 6.945e-9, 3.25496e-9, 1.66679e-9, 9.65e-10, 6.07861e-10, 4.0747e-10]

# Plot the error data
plt.plot(x, y1, marker='s', label='Order 1', color='blue')
plt.plot(x, y2, marker='s', label='Order 2', color='green')
plt.plot(x, y4, marker='s', label='Order 4', color='orange')
plt.plot(x, y8, marker='s', label='Order 8', color='red')

# Annotate the time values
for i, xi in enumerate(x):
    plt.annotate(f'{t1[i]:.3f}s', (xi, y1[i]), textcoords="offset points", xytext=(0,5), ha='center', color='blue')
    plt.annotate(f'{t2[i]:.3f}s', (xi, y2[i]), textcoords="offset points", xytext=(0,5), ha='center', color='green')
    plt.annotate(f'{t4[i]:.3f}s', (xi, y4[i]), textcoords="offset points", xytext=(0,5), ha='center', color='orange')
    plt.annotate(f'{t8[i]:.3f}s', (xi, y8[i]), textcoords="offset points", xytext=(0,5), ha='center', color='red')

# Formatting the plot
plt.yscale("log")
# plt.xlabel("Number of Spatial Cells")
# plt.ylabel("L2 Error Norm")
# plt.title("Error and Time for Different Polynomial Orders")
# plt.legend()

# Show plot
plt.show()
