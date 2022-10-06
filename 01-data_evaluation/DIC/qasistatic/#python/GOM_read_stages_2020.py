# -----------------------------------------------------------------------------
# Read stages from the GOM files writen by Anton Nischler @ LLK 26.09.2022
# Load measured strain field from GOM-System and save the data for interpolation
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
import os
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs

# ---------------------------------------------------------
# Read the input data 
# ---------------------------------------------------------

Rp = 0.2
stage_min = 1
stage_max = 101

try:
    os.remove('stages_means.csv')
except:
    pass

for i in range(stage_min,stage_max+1,2):
    
    if i < 10:
        i_file_id = 'surface_yx_000' + str(i)
    elif i < 100:
        i_file_id = 'surface_yx_00' + str(i)
    else:
        i_file_id = 'surface_yx_0' + str(i)

    try:

        gom_data = pd.read_csv('../GM-26_TD_tens/field_quantities/surface_xy/%s.csv'% i_file_id,\
                                header=None,\
                                skiprows=[0,1,2,3,4,5],\
                                sep=';')

    except:
        print('Error 1: File *.gom could not been opend!')
        sys.exit

# ---------------------------------------------------------
# Calculations
# ---------------------------------------------------------
    eps_y_mean = np.mean(gom_data.iloc[:,5])
    
    if eps_y_mean > Rp:

        eps_x_mean = np.mean(gom_data.iloc[:,4])
        nu_yx = -eps_x_mean / eps_y_mean

        print('-------------------------------------')
        print(i)
        print('eps_x_mean: ', eps_x_mean)
        print('eps_y_mean: ', eps_y_mean)
        print('nu_yx: ', nu_yx)
        print('-------------------------------------')

# ---------------------------------------------------------
# Write results to ASCII file
# ---------------------------------------------------------
        ascii_file = open('stages_means.csv', 'a')
        ascii_file.write('%s,%0.3f,%0.3f,%0.3f\n' % (i_file_id,eps_x_mean,eps_y_mean,nu_yx))
        ascii_file.close()

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

        x = np.linspace(x_min,x_max,n*10)
        y = np.linspace(y_min,y_max,m*10)

        grid_x, grid_y = np.meshgrid(x,y)

# ---------------------------------------------------------
# Interpolation
# ---------------------------------------------------------
        points = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,2]])

        grid_z_0 = griddata(points, gom_data.iloc[:,5], (grid_x, grid_y), method='cubic')

# ---------------------------------------------------------
# Plot the Data
# ---------------------------------------------------------
        plt.figure(figsize=(8,8))
        plt.imshow(grid_z_0, cmap='jet', extent=(x_min,x_max,y_min,y_max), origin='lower',vmin=0,vmax=eps_y_mean*1.5)
        plt.colorbar()
        #plt.grid()

        try:
            os.remove('%s.png' % i_file_id)
        except:
            pass

        plt.savefig('%s.png' % i_file_id)
        plt.title('%s;e_y=%0.2f;e_x=%0.2f;nu_yx=%0.2f' % (i_file_id,eps_y_mean,eps_x_mean,nu_yx))
        #plt.show()


        X = np.column_stack([gom_data.iloc[:,4], gom_data.iloc[:,5]])
        kmeans = KMeans(n_clusters=1).fit(X)
        y_pred = KMeans(n_clusters=1).fit_predict(X)
        centroids = kmeans.cluster_centers_
        print(centroids)

        plt.figure(figsize=(6,6))
        plt.scatter(X[:, 0], X[:, 1], c=y_pred)
        plt.scatter(centroids[:,0],centroids[:,1],marker='x',s=169,linewidths=3)
        plt.plot((0,-5),(0,10))
        plt.title("kmeans")
        plt.xlabel('eps_xx')
        plt.ylabel('eps_yy')
        plt.xlim((eps_x_mean*1.3,0))
        plt.ylim((0,eps_y_mean*1.3))

        try:
            os.remove('%s_kmeans.png' % i_file_id)
        except:
            pass

        plt.savefig('%s_kmeans.png' % i_file_id)

        plt.grid()
        #plt.show()
    
    else:
        pass





# plt.figure()
# n, bins, patches = plt.hist(gom_data.iloc[:,5],50,density=True,facecolor='g',alpha=0.75)
# plt.grid()
# #plt.show()



