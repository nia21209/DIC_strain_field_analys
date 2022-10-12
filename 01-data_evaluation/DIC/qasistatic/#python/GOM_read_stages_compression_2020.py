# -----------------------------------------------------------------------------
# Read stages from the GOM files writen by Anton Nischler @ LLK 26.09.2022
# Load measured strain field from GOM-System and save the data for interpolation
# -----------------------------------------------------------------------------
from cProfile import label
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
import os
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_blobs

specimen_id = str('GM-28_RD_comp_tens')
surface_id = str('surface_yz')
xy_surface = False
file_name = str('surface_yz_40_20')

stage_min = 1
stage_max = 100

threshold_quality = 3
quality_filter = False

Rp = 0.1

try:
    os.remove('%s_mean_plastic.dat' % file_name)
except:
    pass

ascii_file = open('%s_mean_plastic.dat' % file_name, 'a')

if xy_surface == True:
    ascii_file.write('stage,eps_xx_mean_in_mBtG,eps_yy_mean_in_mBtG,nu_yx_mean_in_mBtG,eps_xx_mean_out_mBtG,eps_yy_mean_out_mBtG,nu_yx_mean_out_mBtG\n')
else:
    ascii_file.write('stage,eps_zz_mean_in_mBtG,eps_yy_mean_in_mBtG,nu_yz_mean_in_mBtG,eps_zz_mean_out_mBtG,eps_yy_mean_out_mBtG,nu_yz_mean_out_mBtG\n')

ascii_file.close()

# ----------------------------------------------------------------------
# Function for generating the Meshgrids
# ----------------------------------------------------------------------
def grid(d_ROI):

    n = len(d_ROI[:,0])

    dx = 0
    dy = 0

    for i in range(0,n-1,1):

        dx_i = abs(d_ROI[i,1] - d_ROI[i+1,1])
        dy_i = abs(d_ROI[i,2] - d_ROI[i+1,2])
 
        if dx_i >= dx and dx_i < .2:
            dx = dx_i
        else:
            pass

        if dy_i >= dy and dy_i < .2:
            dy = dy_i
        else:
            pass

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
# Function for separation of the mBtG
# ----------------------------------------------------------------------
def calculation_mBtG(DIC_data):

    X = np.column_stack([DIC_data.iloc[:,4], DIC_data.iloc[:,5]])
    
    kmeans = KMeans(n_clusters=2).fit(X)
    centroids = kmeans.cluster_centers_
    #threshold_one_sigma = abs(0.674490 * np.min(centroids[:,1]))
    
    spread = abs(np.max(centroids[:,1]) - np.min(centroids[:,1]))

    if spread < 1 and Rp >= abs(np.min(centroids[:,1])):

        kmeans = KMeans(n_clusters=1).fit(X)
        y_pred = KMeans(n_clusters=1).fit_predict(X)
        centroids = kmeans.cluster_centers_

        plt.figure(figsize=(6,6))
        plt.scatter(X[:, 0], X[:, 1], c=y_pred)
        plt.scatter(centroids[:,0],centroids[:,1],marker='x',s=169,linewidths=3)
        plt.title("kmeans - %s" % i_file_id)
        plt.xlabel('eps_xx')
        plt.ylabel('eps_yy')
        plt.grid()

        try:
            os.remove('%s_kmeans.png' % i_file_id)
        except:
            pass

        plt.savefig('%s_kmeans.png' % i_file_id)
        #plt.show()

        ascii_file = open('%s_mean_plastic.dat' % file_name, 'a')
        ascii_file.write('%i,nan,nan,nan,%0.6f,%0.6f,%0.6f\n' % (i,\
                                                     centroids[0,0],\
                                                     centroids[0,1],\
                                                     -centroids[0,0]/centroids[0,1]))
        ascii_file.close()

    elif spread < 1 and Rp < abs(np.min(centroids[:,1])):

        kmeans = KMeans(n_clusters=1).fit(X)
        y_pred = KMeans(n_clusters=1).fit_predict(X)
        centroids = kmeans.cluster_centers_

        plt.figure(figsize=(6,6))
        plt.scatter(X[:, 0], X[:, 1], c=y_pred)
        plt.scatter(centroids[:,0],centroids[:,1],marker='x',s=169,linewidths=3)
        plt.title("kmeans - %s" % i_file_id)
        plt.xlabel('eps_xx')
        plt.ylabel('eps_yy')
        plt.grid()

        try:
            os.remove('%s_kmeans.png' % i_file_id)
        except:
            pass

        plt.savefig('%s_kmeans.png' % i_file_id)
        #plt.show()

        ascii_file = open('%s_mean_plastic.dat' % file_name, 'a')
        ascii_file.write('%i,%0.6f,%0.6f,%0.6f,nan,nan,nan\n' % (i,\
                                                     centroids[0,0],\
                                                     centroids[0,1],\
                                                     -centroids[0,0]/centroids[0,1]))
        ascii_file.close()
    
    else:
        kmeans = KMeans(n_clusters=2).fit(X)
        y_pred = KMeans(n_clusters=2).fit_predict(X)
        centroids = kmeans.cluster_centers_

        plt.figure(figsize=(6,6))
        plt.scatter(X[:, 0], X[:, 1], c=y_pred)
        plt.scatter(centroids[:,0],centroids[:,1],marker='x',s=169,linewidths=3)
        plt.title("kmeans - %s" % i_file_id)
        plt.xlabel('eps_xx')
        plt.ylabel('eps_yy')
        plt.grid()

        try:
            os.remove('%s_kmeans.png' % i_file_id)
        except:
            pass

        plt.savefig('%s_kmeans.png' % i_file_id)
        #plt.show()

        # Calculation of the mean values in and out of the mBtG
        index_eps_yy_min = np.argwhere(centroids == np.min(centroids[:,1]))
        eps_xx_mean_in_mBtG = centroids[index_eps_yy_min[0,0],0]
        eps_yy_mean_in_mBtG = np.min(centroids[:,1])
        nu_yx_mean_in_mBtG = -eps_xx_mean_in_mBtG / eps_yy_mean_in_mBtG

        index_eps_yy_max = np.argwhere(centroids == np.max(centroids[:,1]))
        eps_xx_mean_out_mBtG = centroids[index_eps_yy_max[0,0],0]
        eps_yy_mean_out_mBtG = np.max(centroids[:,1])
        nu_yx_mean_out_mBtG = -eps_xx_mean_out_mBtG / eps_yy_mean_out_mBtG

        print('----- %s -----' % i_file_id)
        print('eps_xx_mean_in_mBtG = %0.3f' % eps_xx_mean_in_mBtG)
        print('eps_yy_mean_in_mBtG = %0.3f' % eps_yy_mean_in_mBtG)
        print('nu_yx_mean_in_mBtG = %0.3f' % nu_yx_mean_in_mBtG)

        print('-------------------')
        print('eps_xx_mean_out_mBtG = %0.3f' % eps_xx_mean_out_mBtG)
        print('eps_yy_mean_out_mBtG = %0.3f' % eps_yy_mean_out_mBtG)
        print('nu_yx_mean_out_mBtG = %0.3f' % nu_yx_mean_out_mBtG)

        ascii_file = open('%s_mean_plastic.dat' % file_name, 'a')
        ascii_file.write('%i,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f\n' % (i,\
                                                                       eps_xx_mean_in_mBtG,\
                                                                       eps_yy_mean_in_mBtG,\
                                                                       nu_yx_mean_in_mBtG,\
                                                                       eps_xx_mean_out_mBtG,\
                                                                       eps_yy_mean_out_mBtG,\
                                                                       nu_yx_mean_out_mBtG))
        ascii_file.close()

    return (centroids)

# ----------------------------------------------------------------------
# Function for interpolation
# ----------------------------------------------------------------------
def field_interpolation(ROI_, grid_x_, grid_y_):

    points_ = np.column_stack([ROI_[:,1], ROI_[:,2]])  
    grid_z_ = griddata(points_, ROI_[:,5], (grid_x_, grid_y_), method='cubic')

    return grid_z_

# ----------------------------------------------------------------------
# Read the input data; ROI generation; Plot strain field
# ----------------------------------------------------------------------

try:
    os.remove('stages_means.csv')
except:
    pass

for i in range(stage_min,stage_max+1,2):
    
    if i < 10:
        i_file_id = file_name + '_000' + str(i)
    elif i < 100:
        i_file_id = file_name + '_00' + str(i)
    else:
        i_file_id = file_name + '_0' + str(i)

    try:
        gom_data = pd.read_csv('../%s/field_quantities/%s/%s.csv'% (specimen_id, surface_id, i_file_id),\
                                header=None,\
                                skiprows=[0,1,2,3,4,5],\
                                sep=';')

    except:
        print('Error 1: File *.gom could not been opend!')
        sys.exit

# ----------------------------------------------------------------------
# Quality filter
# ----------------------------------------------------------------------

    if quality_filter is True:

        length_gom_data = len(gom_data.iloc[:,0])

        for j in range(0,length_gom_data,1):

            if gom_data.iloc[j,6] < threshold_quality:
                gom_data.iloc[j,:] = np.nan
            
            else:
                pass
    elif quality_filter is False:
        pass

    gom_data = gom_data.dropna()

# ----------------------------------------------------------------------
# Finde Region of Interests (ROI) an separate data set
# ----------------------------------------------------------------------

    if xy_surface == True:
        X = np.column_stack([gom_data.iloc[:,1], gom_data.iloc[:,1]])
    else:
        X = np.column_stack([gom_data.iloc[:,2], gom_data.iloc[:,2]])

    db = DBSCAN(eps=0.3, min_samples=10).fit(X)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    
    d_ROI = {}
    d_grid_x = {}
    d_grid_y = {}
    d_grid_z = {}

    for m in range(1,n_clusters_+1,1):

        d_ROI['ROI_{}'.format(m)] = np.empty((len(gom_data),len(gom_data.columns)))
        d_ROI['ROI_{}'.format(m)][:] = np.nan
    
        for j in range(0,len(gom_data.iloc[:,0]),1):

            y_pred_i = labels[j]
            
            if y_pred_i == m-1:
                d_ROI['ROI_{}'.format(m)][j,:] = gom_data.iloc[j,:]
            else:
                pass
        
        d_ROI['ROI_{}'.format(m)] = d_ROI['ROI_{}'.format(m)][~np.isnan(d_ROI['ROI_{}'.format(m)]).any(axis=1)]

        d_grid_x['grid_x_{}'.format(m)], d_grid_y['grid_y_{}'.format(m)] = grid(d_ROI['ROI_{}'.format(m)])  

# ---------------------------------------------------------
# Interpolation
# ---------------------------------------------------------
        d_grid_z['grid_z_{}'.format(m)] = field_interpolation(d_ROI['ROI_{}'.format(m)], d_grid_x['grid_x_{}'.format(m)], d_grid_y['grid_y_{}'.format(m)])

# ---------------------------------------------------------
# Plot the Data
# ---------------------------------------------------------
    plt.figure(figsize=(6,8))

    for m in range(1,n_clusters_+1,1):
        plt.pcolormesh(d_grid_x['grid_x_{}'.format(m)], d_grid_y['grid_y_{}'.format(m)], d_grid_z['grid_z_{}'.format(m)], cmap='jet', vmin=-3, vmax=1)

    plt.colorbar()
    plt.grid()
    plt.title('%s' % i_file_id)

    try:
        os.remove('%s.png' % i_file_id)
    except:
        pass

    plt.savefig('%s.png' % i_file_id)
    #plt.show()

# ----------------------------------------------------------------------
# Separate mBtG 
# ----------------------------------------------------------------------
    data = calculation_mBtG(gom_data)

# # --------------------------------------------------------------------
# # Calculations
# # --------------------------------------------------------------------
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

# # --------------------------------------------------------------------
# # Write results to ASCII file
# # ---------------------------------------------------------
#         ascii_file = open('stages_means.csv', 'a')
#         ascii_file.write('%s,%0.3f,%0.3f,%0.3f\n' % (i_file_id,eps_x_mean,eps_y_mean,nu_yx))
#         ascii_file.close()


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

# for later

    # count_BtG_in = 0
    # count_BtG_out = 0

    # for i in range(0,len(DIC_data.iloc[:,0]),1):
        
    #     delta = abs(centroids[1,1] - DIC_data.iloc[i,5])

    #     if delta <= threshold_one_sigma:
    #         count_BtG_in = count_BtG_in + 1

    #     else:
    #         count_BtG_out = count_BtG_out + 1
    
    # data_BtG_in = np.zeros((count_BtG_in,len(DIC_data.iloc[0,:])))
    # data_BtG_out = np.zeros((count_BtG_out,len(DIC_data.iloc[0,:])))

    # count_BtG_in = 0
    # count_BtG_out = 0

    # for i in range(0,len(DIC_data.iloc[:,0]),1):
        
    #     delta = abs(centroids[1,1] - DIC_data.iloc[i,5])

    #     if delta <= threshold_one_sigma:
    #         data_BtG_in[count_BtG_in,:] = DIC_data.iloc[i,:]
    #         count_BtG_in = count_BtG_in + 1

    #     else:
    #         data_BtG_out[count_BtG_out,:] = DIC_data.iloc[i,:]
    #         count_BtG_out = count_BtG_out + 1