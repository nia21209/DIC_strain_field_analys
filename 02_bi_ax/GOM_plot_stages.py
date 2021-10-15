# -----------------------------------------------------------------------------
# Plot stages from the GOM files writen by Anton Nischler @ LLK 15.11.2020
# Plot measured strain field from GOM-System and save the plots in the working
# directory
# -----------------------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import sys

def strain_field_plot(id_field,first_stage,last_stage):
    
    eps_min_22 = -0.5
    eps_max_22 = 0.5
    
    eps_min_2 = -.5
    eps_max_2 = .5
    
    id_x_1 = {}
    id_x_2 = {}
    id_strain_22 = {}
    id_strain_2 = {}
    
    for i in range(first_stage,last_stage+1,1):
        
    # --------------------------------------------------------------------------
    # --------------- Buid scalar field for plots etc. -------------------------
    # --------------------------------------------------------------------------
        name = 'stage-' + str(i)

        print('Plot GOM data from %s' % name)
    
        Index_x_max = int(np.max(id_field['%s' % i][:, 0])) + 1  # Find biggest x Index
        Index_y_max = int(np.max(id_field['%s' % i][:, 1])) + 1  # Find biggest y Index
    
        length_id_field = int(id_field['%s' % i].shape[0])  # max rows of id_field
    
    # --------------------------------------------------------------------------
    # ------------------------ Buid xy meshgrid --------------------------------
    # --------------------------------------------------------------------------
    
        # Initialize the array for the field variables
        x_1 = np.zeros((Index_x_max, Index_y_max))
        x_2 = np.zeros((Index_x_max, Index_y_max))
    
        # Evaluate dx
        n = 1
        
        dIndex_x = abs(abs(id_field['%s' % i][0, 0]) - abs(id_field['%s' % i][1, 0]))
        dx = abs(abs(id_field['%s' % i][0, 2]) - abs(id_field['%s' % i][1, 2]))

        while dIndex_x != 1:
            dIndex_x = abs(abs(id_field['%s' % i][n, 0]) - abs(id_field['%s' % i][n+1, 0]))
            dx = abs(abs(id_field['%s' % i][n, 2]) - abs(id_field['%s' % i][n+1, 2]))
            if n == length_id_field:
                print('Fatal error: Facette length in x direction can not be calculated')
                break
            n = n + 1
            
        # Evaluate dy
        n = 1

        dIndex_y = abs(abs(id_field['%s' % i][0, 1]) - abs(id_field['%s' % i][1, 1]))
        dy = abs(abs(id_field['%s' % i][0, 3]) - abs(id_field['%s' % i][1, 3]))

        while dIndex_y != 1:
            dIndex_y = abs(abs(id_field['%s' % i][n, 1]) - abs(id_field['%s' % i][n+1, 1]))
            dy = abs(abs(id_field['%s' % i][n, 3]) - abs(id_field['%s' % i][n+1, 3]))
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
    
        # Update coordinates list
        id_x_1.update({'%s' % name:x_1})
        id_x_2.update({'%s' % name:x_2})
    
    # --------------------------------------------------------------------------
    # ------------------ Build grid for scalars (strain) -----------------------
    # --------------------------------------------------------------------------
    
        strain_22 = np.zeros((Index_x_max, Index_y_max))
        strain_2 = np.zeros((Index_x_max, Index_y_max))
    
        for j in range(0, length_id_field, 1):
            strain_22[int(id_field['%s' % str(i)][j,0]), int(id_field['%s' % str(i)][j,1])] = id_field['%s' % str(i)][j,9]
            strain_2[int(id_field['%s' % str(i)][j,0]), int(id_field['%s' % str(i)][j,1])] = id_field['%s' % str(i)][j,12]
    
        # Mask 0.0 to NaN for plotting
        strain_22[strain_22 == 0.0] = np.nan
        strain_2[strain_2 == 0.0] = np.nan
    
        # Update the strain field list
        id_strain_22.update({'%s' % name: strain_22})
        id_strain_2.update({'%s' % name: strain_2})
    
    # --------------------------------------------------------------------------
    # ------------------------- Plot scalar field ------------------------------
    # --------------------------------------------------------------------------
    
        try:
        # Create target Directory
            os.mkdir(name)
        except:
            print("Directory ", name, " allready exists")
            sys.exit
    
        # Plot in LaTex-Style
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
        fig, ax = plt.subplots(figsize=(5,5/scale_factor))
        minorLocator_x = MultipleLocator(1)
        minorLocator_y = MultipleLocator(1)
        if i == first_stage:
            field_plot = ax.pcolormesh(id_x_1['%s' % name], id_x_2['%s' % name],\
                                   id_strain_22['%s' % name] * 100, cmap='jet',\
                                   vmin=eps_min_22, vmax=eps_max_22)
        else:
            field_plot = ax.pcolormesh(id_x_1['%s' % name], id_x_2['%s' % name],\
                                   id_strain_22['%s' % name] * 100, cmap='jet',\
                                   vmin=eps_min_22, vmax=eps_max_22)
        plt.tick_params(axis='both', which='both', top='true', right='true',\
                        direction='in', grid_color='lightgray', grid_linewidth='0.5')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(r'strain field $\varepsilon_{22}$ (\%)')
        fig.colorbar(field_plot, ax=ax)
        ax.xaxis.set_minor_locator(minorLocator_x)
        ax.yaxis.set_minor_locator(minorLocator_y)
        fig.tight_layout()
        plt.savefig('%s/strain_22.pdf' % name, format = 'pdf')
        
        if i == last_stage:
            plt.show()
        else:
            pass
            
        plt.close()
    
    
        # Plot in LaTex-Style
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
        fig, ax = plt.subplots(figsize=(5,5/scale_factor))
        minorLocator_x = MultipleLocator(1)
        minorLocator_y = MultipleLocator(1)
        field_plot = ax.pcolormesh(id_x_1['%s' % name], id_x_2['%s' % name],\
                                   id_strain_2['%s' % name] * 100, cmap='jet',\
                                   vmin=eps_min_2, vmax=eps_max_2)#
        plt.tick_params(axis='both', which='both', top='true', right='true',\
                        direction='in', grid_color='lightgray', grid_linewidth='0.5')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(r'strain field $\varepsilon_{2}$ (\%)')
        fig.colorbar(field_plot, ax=ax)
        ax.xaxis.set_minor_locator(minorLocator_x)
        ax.yaxis.set_minor_locator(minorLocator_y)
        fig.tight_layout()
        plt.savefig('%s/strain_2.pdf' % name, format = 'pdf')
#        plt.show()
        plt.close()
