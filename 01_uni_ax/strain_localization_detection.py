# -----------------------------------------------------------------------------
# Filter HeV from the GOM files writen by Anton Nischler @ LLK 15.11.2020
# Filter measured strain field from GOM-System to calculate the highly strained
# volume
# -----------------------------------------------------------------------------

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
#from matplotlib import cm

def strain_localization_detection(id_field,first_stage,last_stage,thr_1,thr_2):
    
    plt.figure('Histogram')
    n, bins, patches = plt.hist((np.nan_to_num(id_field['%s' % last_stage][:,9])*100),100, density=False, facecolor='g')#alpha=0.75
    plt.show()
    plt.grid(which='minor')
    plt.close()
    
    def filter_HeV(Index_x_max,Index_y_max):
    
        print('Calculate the highly_strained_region with high strained amplitudes')
        
        # ---------------------------------------------------------------------
        # ------------------------ Buid xy meshgrid ---------------------------
        # ---------------------------------------------------------------------
    
        # Initialize the array for the field variables
        x_1 = np.zeros((Index_x_max, Index_y_max))
        x_2 = np.zeros((Index_x_max, Index_y_max))
    
        # Evaluate dx
        n = 1
        
        dIndex_x = abs(abs(id_field['%s' % last_stage][0, 0]) - abs(id_field['%s' % last_stage][1, 0]))
        dx = abs(abs(id_field['%s' % last_stage][0, 2]) - abs(id_field['%s' % last_stage][1, 2]))

        while dIndex_x != 1:
            dIndex_x = abs(abs(id_field['%s' % last_stage][n, 0]) - abs(id_field['%s' % last_stage][n+1, 0]))
            dx = abs(abs(id_field['%s' % last_stage][n, 2]) - abs(id_field['%s' % last_stage][n+1, 2]))
            if n == length_id_field:
                print('Fatal error: Facette length in x direction can not be calculated')
                break
            n = n + 1
            
        # Evaluate dy
        n = 1

        dIndex_y = abs(abs(id_field['%s' % last_stage][0, 1]) - abs(id_field['%s' % last_stage][1, 1]))
        dy = abs(abs(id_field['%s' % last_stage][0, 3]) - abs(id_field['%s' % last_stage][1, 3]))

        while dIndex_y != 1:
            dIndex_y = abs(abs(id_field['%s' % last_stage][n, 1]) - abs(id_field['%s' % last_stage][n+1, 1]))
            dy = abs(abs(id_field['%s' % last_stage][n, 3]) - abs(id_field['%s' % last_stage][n+1, 3]))
            if n == length_id_field:
                print('Fatal error: Facette length in x direction can not be calculated')
                break
            n = n + 1
        
        scale_factor = dx/dy
    
        # Array for x coordinates
        for k in range(0, Index_y_max, 1):
            x = 0
            for l in range(0, Index_x_max, 1):
                x_1[l, k] = x
                x = x + dx
    
        # Array for y coordinates
        y = 0
        for k in range(0, Index_y_max, 1):
            for l in range(0, Index_x_max, 1):
                x_2[l, k] = y
            y = y + dy
           
        # ---------------------------------------------------------------------
        # -------------- Build grid for highly strained region ----------------
        # ---------------------------------------------------------------------

        id_fild_lll_error=([])
        id_fild_ull_error=([])
        
        print('-------------------------------------------------------------------')
        Error = float(input('Define the minimum strain_22 in % to correct measurment errors: '))/100
        print('-------------------------------------------------------------------')
        
        # filter Errors from strain field
        for i in range(0,length_id_field,1):
            
            if abs(id_field['%s' % last_stage][i,9]) > abs(Error):
                id_fild_lll_error = np.append(id_fild_lll_error,0.)
                id_fild_ull_error = np.append(id_fild_ull_error,0.)
            
            else:
                id_fild_lll_error = np.append(id_fild_lll_error,id_field['%s' % last_stage][i,9])
                id_fild_ull_error = np.append(id_fild_ull_error,id_field['%s' % first_stage][i,9])
                   
        strain_22_lll = np.zeros((Index_x_max, Index_y_max))
        strain_22_ull = np.zeros((Index_x_max, Index_y_max))
            
        strain_22_thr_lll = thr_1 * np.min(id_fild_lll_error)
#        print('threshold strain_22 at the lower load level: ', strain_22_thr_lll, ' (%)')
        
        for j in range(0, length_id_field, 1):
            
            if id_field['%s' % last_stage][j,9] <= strain_22_thr_lll:
                strain_22_lll[int(id_field['%s' % last_stage][j,0]), int(id_field['%s' % last_stage][j,1])] = id_fild_lll_error[j]#id_field['%s' % last_stage][j,9]
                strain_22_ull[int(id_field['%s' % first_stage][j,0]), int(id_field['%s' % first_stage][j,1])] = id_fild_ull_error[j]#id_field['%s' % first_stage][j,9]

            else:
                pass
        
        # ---------------------------------------------------------------------
        #  Build grid for highly strained region with high strained amplitudes
        # ---------------------------------------------------------------------
        
        strain_22_a = (strain_22_ull - strain_22_lll) / 2
        strain_22_a_thr = thr_2 * np.max(strain_22_a)
#        print('threshold strain_22 amplitude: ', strain_22_a_thr, ' (%)')
        counter = 0
                
        for j in range(0,Index_y_max,1):
            for i in range(0,Index_x_max,1):
                
                strain_22_a_i = strain_22_a[i,j]
                
                if strain_22_a_i >= strain_22_a_thr:
                    counter = counter + 1
                    pass
                else:
                    strain_22_a[i,j] = 0.
        
        strain_22_a_av = np.sum(strain_22_a) / counter
        A_HeV = dx * dy * counter
                    
        # ---------------------------------------------------------------------
        # ------------------------- Plot scalar field -------------------------
        # ---------------------------------------------------------------------
        
        print('Plot the highly_strained_region with high strained amplitudes')
    
        # Mask 0.0 to NaN for plotting
        strain_22_lll[strain_22_lll == 0.0] = np.nan
        strain_22_ull[strain_22_ull == 0.0] = np.nan
        strain_22_a[strain_22_a == 0.0] = np.nan
    
        try:
            # Create target Directory
            os.mkdir('plot_HeV')
        except:
            print("Directory plot_HeV allready exists")
    
        # Plot in LaTex-Style
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
        fig, ax = plt.subplots(figsize=(5,5/scale_factor))
        minorLocator_x = MultipleLocator(1)
        minorLocator_y = MultipleLocator(1)
        field_plot = ax.pcolormesh(x_1, x_2,\
                                   strain_22_lll * 100, cmap='jet',\
                                   vmin=-1, vmax=0)
        plt.tick_params(axis='both', which='both', top='true', right='true',\
                        direction='in', grid_color='lightgray', grid_linewidth='0.5')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(r'strain field $\varepsilon_{22}|_{\mathrm{LLL}} (\%)$')
        fig.colorbar(field_plot, ax=ax)
        ax.xaxis.set_minor_locator(minorLocator_x)
        ax.yaxis.set_minor_locator(minorLocator_y)
        plt.tight_layout()
        plt.savefig('plot_HeV/strain_22_lll.pdf', format = 'pdf')
#        plt.show()
        plt.close()
           
        # Plot in LaTex-Style
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
        fig, ax = plt.subplots(figsize=(5,5/scale_factor))
        minorLocator_x = MultipleLocator(1)
        minorLocator_y = MultipleLocator(1)
        field_plot = ax.pcolormesh(x_1, x_2,\
                                   strain_22_ull * 100, cmap='jet',\
                                   vmin=-1, vmax=1)
        plt.tick_params(axis='both', which='both', top='true', right='true',\
                        direction='in', grid_color='lightgray', grid_linewidth='0.5')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(r'strain field $\varepsilon_{22}|_{\mathrm{ULL}} (\%)$')
        fig.colorbar(field_plot, ax=ax)
        ax.xaxis.set_minor_locator(minorLocator_x)
        ax.yaxis.set_minor_locator(minorLocator_y)
        plt.tight_layout()
        plt.savefig('plot_HeV/strain_22_ull.pdf', format = 'pdf')
#        plt.show()
        plt.close()
        
        # Plot in LaTex-Style
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
        fig, ax = plt.subplots(figsize=(5,5/scale_factor))
        minorLocator_x = MultipleLocator(1)
        minorLocator_y = MultipleLocator(1)
        field_plot = ax.pcolormesh(x_1, x_2,\
                                   strain_22_a * 100, cmap='jet',\
                                   vmin=0, vmax=strain_22_a_av)
        plt.tick_params(axis='both', which='both', top='true', right='true',\
                        direction='in', grid_color='lightgray', grid_linewidth='0.5')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(r'amplitude strain field $\varepsilon_{22}|_{\mathrm{a}} (\%)$')
        fig.colorbar(field_plot, ax=ax)
        ax.xaxis.set_minor_locator(minorLocator_x)
        ax.yaxis.set_minor_locator(minorLocator_y)
        plt.tight_layout()
        plt.savefig('plot_HeV/strain_22_a.pdf', format = 'pdf')
        plt.show()
        plt.close()
        
        return (strain_22_thr_lll,strain_22_a_thr,strain_22_a_av,A_HeV,Error)
        
    
    
    ull_rows = len(id_field['%s' % first_stage][:,9])
    lll_rows = len(id_field['%s' % last_stage][:,9])
    
    if ull_rows <= lll_rows:
        
        # ---------------------------------------------------------------------
        # --------------- Buid scalar field for plots etc. --------------------
        # ---------------------------------------------------------------------
        Index_x_max = int(np.max(id_field['%s' % first_stage][:,0])) + 1
        Index_y_max = int(np.max(id_field['%s' % first_stage][:,1])) + 1
        
        length_id_field = int(id_field['%s' % first_stage].shape[0])  # max rows of id_field
        
        strain_22_thr_lll,strain_22_a_thr,strain_22_a_av,A_HeV,Error = filter_HeV(Index_x_max,Index_y_max)
            
    else:
        
        # ---------------------------------------------------------------------
        # --------------- Buid scalar field for plots etc. --------------------
        # ---------------------------------------------------------------------
        Index_x_max = int(np.max(id_field['%s' % last_stage][:,0])) + 1
        Index_y_max = int(np.max(id_field['%s' % last_stage][:,1])) + 1
        
        length_id_field = int(id_field['%s' % last_stage].shape[0])  # max rows of id_field
        
        strain_22_thr_lll,strain_22_a_thr,strain_22_a_av,A_HeV,Error = filter_HeV(Index_x_max,Index_y_max)
    
    return (strain_22_thr_lll,strain_22_a_thr,strain_22_a_av,A_HeV,Error)
