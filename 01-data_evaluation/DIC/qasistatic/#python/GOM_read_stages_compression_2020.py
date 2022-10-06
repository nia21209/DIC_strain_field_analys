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

n_ROI = 2
Rp = 0.2
stage_min = 9
stage_max = 9

def grid(d_ROI):

    n = len(d_ROI[:,0])

    dx = 0
    dy = 0

    for i in range(0,n-1,1):

        dx_i = abs(d_ROI[i,1] - d_ROI[i+1,1])
        dy_i = abs(d_ROI[i,2] - d_ROI[i+1,2])
            
        if dx_i >= dx and dx_i < 1:
            dx = dx_i
        else:
            pass

        if dy_i >= dy and dy_i < 1:
            dy = dy_i
        else:
            pass
        
    print(dx,dy)

    x_min = np.min(d_ROI[:,1])
    x_max = np.max(d_ROI[:,1])

    y_min = np.min(d_ROI[:,2])
    y_max = np.max(d_ROI[:,2])

    x_len = x_max - x_min
    y_len = y_max - y_min

    n = int(x_len/dx)
    m = int(y_len/dy)

    x = np.linspace(x_min,x_max,n*10)
    y = np.linspace(y_min,y_max,m*10)

    grid_x, grid_y = np.meshgrid(x,y)

    return (grid_x,grid_y)

# ----------------------------------------------------------------------

try:
    os.remove('stages_means.csv')
except:
    pass

for i in range(stage_min,stage_max+1,2):
    
    if i < 10:
        i_file_id = 'surface_yx_20_10_000' + str(i)
    elif i < 100:
        i_file_id = 'surface_yx_20_10_00' + str(i)
    else:
        i_file_id = 'surface_yx_20_10_0' + str(i)

    try:
        gom_data = pd.read_csv('../GM-28_RD_comp_tens/field_quantities/surface_yx/%s.csv'% i_file_id,\
                                header=None,\
                                skiprows=[0,1,2,3,4,5],\
                                sep=';')

    except:
        print('Error 1: File *.gom could not been opend!')
        sys.exit

    X = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,1]])
    y_pred = KMeans(n_clusters=n_ROI).fit_predict(X)
    
    d_ROI = {}
    for m in range(1,n_ROI+1,1):

        d_ROI['ROI_{}'.format(m)] = np.empty((len(gom_data),len(gom_data.columns)))
        d_ROI['ROI_{}'.format(m)][:] = np.nan
    
        for j in range(0,len(gom_data.iloc[:,0]),1):
            y_pred_i = y_pred[j]
            
            if y_pred_i == m-1:
                d_ROI['ROI_{}'.format(m)][j,:] = gom_data.iloc[j,:]
            else:
                pass
        
        d_ROI['ROI_{}'.format(m)] = d_ROI['ROI_{}'.format(m)][~np.isnan(d_ROI['ROI_{}'.format(m)]).any(axis=1)]

    grid_x_1, grid_y_1 = grid(d_ROI['ROI_1'])
    grid_x_2, grid_y_2 = grid(d_ROI['ROI_2'])
    
    # plt.figure()
    # plt.scatter(d_ROI['ROI_1'][:,1],d_ROI['ROI_1'][:,2])
    # plt.scatter(d_ROI['ROI_2'][:,1],d_ROI['ROI_2'][:,2])
    # plt.grid()
    # plt.show()

# ---------------------------------------------------------
# Interpolation
# ---------------------------------------------------------
    points1 = np.column_stack([d_ROI['ROI_1'][:,1], d_ROI['ROI_1'][:,2]])
    grid_z_0 = griddata(points1, d_ROI['ROI_1'][:,5], (grid_x_1, grid_y_1), method='cubic')

    points2 = np.column_stack([d_ROI['ROI_2'][:,1], d_ROI['ROI_2'][:,2]])
    grid_z_1 = griddata(points2, d_ROI['ROI_2'][:,5], (grid_x_2, grid_y_2), method='cubic')

# ---------------------------------------------------------
# Plot the Data
# ---------------------------------------------------------
    plt.figure(figsize=(6,8))
    plt.pcolormesh(grid_x_1,grid_y_1,grid_z_0,cmap='jet', vmin=-3, vmax=1)
    plt.pcolormesh(grid_x_2,grid_y_2,grid_z_1,cmap='jet', vmin=-3, vmax=1)
    plt.colorbar()
    plt.grid()
    plt.show()

    # plt.figure(figsize=(8,8))
    # plt.imshow(grid_z_0, cmap='jet', extent=(-2,-1,-1,1), origin='lower',vmin=-3,vmax=3) 
    # plt.imshow(grid_z_1, cmap='jet', extent=(1,2,-1,1), origin='lower',vmin=-3,vmax=3)
    # plt.colorbar()
    # plt.grid()

    # try:
    #     os.remove('%s.png' % i_file_id)
    # except:
    #     pass

    # #plt.savefig('%s.png' % i_file_id)
    # #plt.title('%s;e_y=%0.2f;e_x=%0.2f;nu_yx=%0.2f' % (i_file_id,eps_y_mean,eps_x_mean,nu_yx))
    # plt.show()

# # ---------------------------------------------------------
# # Calculations
# # ---------------------------------------------------------
#     eps_y_mean = abs(np.mean(gom_data.iloc[:,5]))
    
#     if eps_y_mean > Rp:

#         eps_x_mean = np.mean(gom_data.iloc[:,4])
#         nu_yx = -eps_x_mean / eps_y_mean

#         print('-------------------------------------')
#         print(i)
#         print('eps_x_mean: ', eps_x_mean)
#         print('eps_y_mean: ', eps_y_mean)
#         print('nu_yx: ', nu_yx)
#         print('-------------------------------------')

# # ---------------------------------------------------------
# # Write results to ASCII file
# # ---------------------------------------------------------
#         ascii_file = open('stages_means.csv', 'a')
#         ascii_file.write('%s,%0.3f,%0.3f,%0.3f\n' % (i_file_id,eps_x_mean,eps_y_mean,nu_yx))
#         ascii_file.close()

# # ---------------------------------------------------------
# # Generate Meshgrid
# # ---------------------------------------------------------
#         n = len(gom_data.iloc[:,0])

#         dx = 0
#         dy = 0

#         for i in range(0,n-1,1):
            
#             dx_i = abs(gom_data.iloc[i,1] - gom_data.iloc[i+1,1])
#             dy_i = abs(gom_data.iloc[i,2] - gom_data.iloc[i+1,2])
            
#             if dx_i >= dx and dx_i < 1:
#                 dx = dx_i
#             else:
#                 pass

#             if dy_i >= dy and dy_i < 1:
#                 dy = dy_i
#             else:
#                 pass

#         x_min = np.min(gom_data.iloc[:,1])
#         x_max = np.max(gom_data.iloc[:,1])

#         y_min = np.min(gom_data.iloc[:,2])        X = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,1]])
#         kmeans = KMeans(n_clusters=2,n_init=10).fit(X)
#         y_pred = KMeans(n_clusters=2).fit_predict(X)
#         centroids = kmeans.cluster_centers_
#         print(centroids)
#         for i in range(0,len(y_pred),1):
#             print(y_pred[i])

#         x = np.linspace(x_min,0,n*10)
#         y = np.linspace(y_min,y_max,m*10)

#         grid_x, grid_y = np.meshgrid(x,y)

# # ---------------------------------------------------------
# # Interpolation
# # ---------------------------------------------------------
#         points = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,2]])

#         grid_z_0 = griddata(points, gom_data.iloc[:,5], (grid_x, grid_y), method='cubic')

# # ---------------------------------------------------------
# # Plot the Data
# # ---------------------------------------------------------
#         # plt.figure(figsize=(8,8))
#         # plt.imshow(grid_z_0, cmap='jet', extent=(x_min,x_max,y_min,y_max), origin='lower',vmin=-3,vmax=3)
#         # plt.colorbar()
#         # #plt.grid()

#         # try:
#         #     os.remove('%s.png' % i_file_id)
#         # except:
#         #     pass

#         # #plt.savefig('%s.png' % i_file_id)
#         # plt.title('%s;e_y=%0.2f;e_x=%0.2f;nu_yx=%0.2f' % (i_file_id,eps_y_mean,eps_x_mean,nu_yx))
#         # #plt.show()


#         X = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,1]])
#         kmeans = KMeans(n_clusters=2,n_init=10).fit(X)
#         y_pred = KMeans(n_clusters=2).fit_predict(X)
#         centroids = kmeans.cluster_centers_
#         print(centroids)
#         for i in range(0,len(y_pred),1):
#             print(y_pred[i])
      
#         plt.figure(figsize=(6,6))
#         plt.scatter(X[:, 0], X[:, 1], c=y_pred)
#         plt.scatter(centroids[:,0],centroids[:,1],marker='x',s=169,linewidths=3)
#         #plt.plot((0,-5),(0,10))
#         plt.title("kmeans")
#         plt.xlabel('eps_xx')
#         plt.ylabel('eps_yy')
#         #plt.xlim((eps_x_mean*1.3,0))
#         #plt.ylim((0,eps_y_mean*1.3))

#         try:
#             os.remove('%s_kmeans.png' % i_file_id)
#         except:
#             pass

#         #plt.savefig('%s_kmeans.png' % i_file_id)

#         plt.grid()
#         plt.show()
    
#     else:
#         pass

#     # plt.figure()
#     # plt.plot(gom_data.iloc[:,1],gom_data.iloc[:,2],'.')
#     # plt.grid()
#     # plt.show()