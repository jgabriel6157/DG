import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

plt.rcParams.update({'font.size': 22})

x = np.array([4, 8, 12, 16, 20, 24, 28, 32])

# Order 1 data
t1 = np.array([0.00417163, 0.01456, 0.0322, 0.055, 0.085, 0.125, 0.1647, 0.213307])
y1 = np.array([0.294571, 0.057, 0.0220482, 0.0115, 0.007, 0.004775, 0.00345, 0.00261685])

# Order 2 data
t2 = np.array([0.0117978, 0.04304, 0.095, 0.166077, 0.258, 0.369, 0.509, 0.651782])
y2 = np.array([0.0234314, 0.00286, 0.00085, 0.000358, 0.00018, 0.000106, 6.69e-5, 4.482e-5])

# Order 4 data
t4 = np.array([0.05, 0.188867, 0.42368, 0.736451, 1.14436, 1.68, 2.24874, 3.01426])
y4 = np.array([0.000206, 6.58286e-6, 8.75e-7, 2.09484e-7, 6.937e-8, 2.81658e-8, 1.333e-8, 7.06298e-9])

# Order 8 data
t8 = np.array([0.287074, 1.07451, 2.4, 4.25398, 6.53409, 9.44362, 12.9971, 16.8101])
y8 = np.array([2.08436e-7, 2.60434e-8, 6.945e-9, 3.25496e-9, 1.66679e-9, 9.65e-10, 6.07861e-10, 4.0747e-10])

# Create interpolation functions for time and error for each order
time_interp_1st_order = interp1d(t1, x, kind='linear', fill_value="extrapolate")
error_interp_1st_order = interp1d(x, y1, kind='linear', fill_value="extrapolate")

time_interp_2nd_order = interp1d(t2, x, kind='linear', fill_value="extrapolate")
error_interp_2nd_order = interp1d(x, y2, kind='linear', fill_value="extrapolate")

time_interp_4th_order = interp1d(t4, x, kind='linear', fill_value="extrapolate")
error_interp_4th_order = interp1d(x, y4, kind='linear', fill_value="extrapolate")

time_interp_8th_order = interp1d(t8, x, kind='linear', fill_value="extrapolate")
error_interp_8th_order = interp1d(x, y8, kind='linear', fill_value="extrapolate")

# Function to get grid cells and accuracy for a given time and order
def query_simulation(order, desired_time):
    if order == 1:
        grid_cells_for_time = time_interp_1st_order(desired_time)
        accuracy_for_time = error_interp_1st_order(grid_cells_for_time)
    elif order == 2:
        grid_cells_for_time = time_interp_2nd_order(desired_time)
        accuracy_for_time = error_interp_2nd_order(grid_cells_for_time)
    elif order == 4:
        grid_cells_for_time = time_interp_4th_order(desired_time)
        accuracy_for_time = error_interp_4th_order(grid_cells_for_time)
    elif order == 8:
        grid_cells_for_time = time_interp_8th_order(desired_time)
        accuracy_for_time = error_interp_8th_order(grid_cells_for_time)
    else:
        raise ValueError("Order not supported. Choose 1, 2, 4, or 8.")
    
    return grid_cells_for_time, accuracy_for_time

# Example: Querying for 1 second at different orders
desired_time = 0.8

for order in [1, 2, 4, 8]:
    grid_cells, accuracy = query_simulation(order, desired_time)
    print(f"At order {order}, for a simulation time of {desired_time:.2f} seconds:")
    print(f" - Number of grid cells used: {grid_cells:.2f}")
    print(f" - Corresponding accuracy (L2 error norm): {accuracy:.3e}\n")
