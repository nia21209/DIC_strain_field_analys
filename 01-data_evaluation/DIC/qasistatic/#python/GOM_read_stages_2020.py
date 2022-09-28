# -----------------------------------------------------------------------------
# Read stages from the GOM files writen by Anton Nischler @ LLK 26.09.2022
# Load measured strain field from GOM-System and save the data for interpolation
# -----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys

# ---------------------------------------------------------
# Read the input data 
# ---------------------------------------------------------
try:

    gom_data = pd.read_csv('surface_yx_0019.csv',\
                            header=None,\
                            skiprows=[0,1,2,3,4,5],\
                            sep=';')

except:
    print('Error 1: File *.gom could not been opend!')
    sys.exit

# ---------------------------------------------------------
# Generate Meshgrid
# ---------------------------------------------------------
n = len(gom_data.iloc[:,0])

dx = 0
dy = 0

for i in range(0,n-1,1):
    
    dx_i = abs(gom_data.iloc[i,1] - gom_data.iloc[i+1,1])
    dy_i = abs(gom_data.iloc[i,2] - gom_data.iloc[i+1,2])
    
    if dx_i >= dx and dx_i < 1:
        dx = dx_i
    else:
        pass

    if dy_i >= dy and dy_i < 1:
        dy = dy_i
    else:
        pass

x_min = np.min(gom_data.iloc[:,1])
x_max = np.max(gom_data.iloc[:,1])

y_min = np.min(gom_data.iloc[:,2])
y_max = np.max(gom_data.iloc[:,2])

x_len = x_max - x_min
y_len = y_max - y_min

n = int(x_len/dx)
m = int(y_len/dy)

x = np.linspace(x_min,x_max,n*1)
y = np.linspace(y_min,y_max,m*1)

grid_x, grid_y = np.meshgrid(x,y)

# ---------------------------------------------------------
# Interpolation
# ---------------------------------------------------------
points = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,2]])

grid_z_0 = griddata(points, gom_data.iloc[:,5], (grid_x, grid_y), method='nearest')

# ---------------------------------------------------------
# Calculations
# ---------------------------------------------------------
eps_x_mean = np.mean(gom_data.iloc[:,4])
eps_y_mean = np.nanmean(grid_z_0)#(gom_data.iloc[:,5])
print('eps_x_mean: ', eps_x_mean)
print('eps_y_mean: ', eps_y_mean)

# ---------------------------------------------------------
# Plot the Data
# ---------------------------------------------------------
plt.subplot(111)
plt.imshow(grid_z_0, cmap='jet', extent=(x_min,x_max,y_min,y_max), origin='lower',vmin=0,vmax=1.)
plt.colorbar()
plt.grid()
plt.show()