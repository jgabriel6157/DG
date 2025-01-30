import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

x = np.array([4,8,16,32])

y1 = np.array([4e-1,9.1e-2,2.3e-2,5.7e-3])
t1 = np.array([2.19,3.5,8.13,19.6])
y2 = np.array([4.3e-2,6.3e-3,8e-4,1e-4])
t2 = np.array([3.75,7.31,15.3,38.4])
y4 = np.array([3.1e-4,9.6e-6,3.2e-7,1e-8])
t4 = np.array([8.94,20,45,115])
y8 = np.array([2.5e-9,4.8e-12,2.2e-13,5e-13])
t8 = np.array([32,68.3,163,665])

c1 = [18.61,10.69,4.38]
e1 = [2.018e-2,4.448e-3,2.812e-4]
c2 = [31.16,18.56,7.64]
e2 = [6.605e-3,6.879e-4,3.676e-5]
c3 = [26.18,11.2,3.78]
e3 = [3.545e-4,5.888e-6,2.637e-9]
c4 = [31.72,13.76,4.66]
e4 = [1.121e-4,2.918e-6,2.088e-9]

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
desired_time = 38

for order in [1, 2, 4, 8]:
    grid_cells, accuracy = query_simulation(order, desired_time)
    print(f"At order {order}, for a simulation time of {desired_time:.2f} seconds:")
    print(f" - Number of grid cells used: {grid_cells:.2f}")
    print(f" - Corresponding accuracy (L2 error norm): {accuracy:.3e}\n")

plt.plot(x,y1,marker='s')
plt.plot(x,y2,marker='s')
plt.plot(x,y4,marker='s')
plt.plot(x,y8,marker='s')

plt.plot(c1,e1,linestyle='--',color='k')
plt.plot(c2,e2,linestyle='--',color='k')
plt.plot(c3,e3,linestyle='--',color='k')
plt.plot(c4,e4,linestyle='--',color='k')

plt.yscale("log")
plt.show()