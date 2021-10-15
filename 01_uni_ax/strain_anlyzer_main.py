# -----------------------------------------------------------------------------
# Anton Nischler @ LLK 15.10.2021
# Main program to analyse a strain field from GOM
# ASCI file needs to be exported with "Nischler"-formate
# -----------------------------------------------------------------------------

import os

# -----------------------------------------------------------------------------
# Call the GOM_read_stages function to get a list from the upper and lower
# load level and return

import GOM_read_stages as GRS

print('---------------------------------------------------------')
specimen_id = str(input('Enter the specimen_id (string): '))
first_stage = int(input('Enter the number of the stage from the upper load level (integer): '))
last_stage = int(input('Enter the number of the stage from the lower load level (integer): '))
thr_1 = float(input('Define the first threshold for the lower load level (float): '))
thr_2 = float(input('Define the second threshold for the strain amplitude (float): '))
print('---------------------------------------------------------')


#os.chdir('%s/GOM/' % specimen_id)
id_field = GRS.strain_field_read(specimen_id,first_stage,last_stage)
print(id_field)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Call the GOM_plot_stages function to plot the upper and lower strain fields

import GOM_plot_stages as GPS

plot_switch = 1

if plot_switch == 1:
    
    GPS.strain_field_plot(id_field,first_stage,last_stage)
    
else:
    pass
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Call the strain_localization_detection function to filter the highly strained 
# reagion

# import strain_localization_detection as HeV

# t = 2                               # sheet thickness
# surface_correction = 1.452          # correction for curved surface and anti buckling device

# strain_22_thr_lll,strain_22_a_thr,strain_22_a_av,A_HeV,Error = HeV.strain_localization_detection(id_field,first_stage,last_stage,thr_1,thr_2)
# print('---------------------------------------------------------')
# print('strain_22_thr_lll: ', strain_22_thr_lll*100, ' (%)')
# print('strain_22_a_thr: ', strain_22_a_thr*100, ' (%)')
# print('strain_22_a_av :', strain_22_a_av*100, ' (%)')
# print('A_HeV: ', A_HeV*surface_correction, ' (mm^2)')
# print('V_HeV: ', A_HeV*t, ' (mm^3)')
# print('---------------------------------------------------------')

#-----------------------------------------------------------------------
#------------------------ Write Results --------------------------------
#-----------------------------------------------------------------------
# file = open('%s.dat' % specimen_id, 'w')

# file.write('---- Input ----\n')
# file.write('threshold HeV: %0.2f\n' % thr_1)
# file.write('threshold HeV amplitude: %0.2f\n' % thr_2)
# file.write('Error correction: %0.5f\n' % Error)
# file.write('\n')
# file.write('---- Output ----\n')
# file.write('strain_22 threshold at the lower load level: %0.3f ' % (strain_22_thr_lll*100))
# file.write('(%)\n')
# file.write('strain_22 amplitude threshold: %0.3f' % (strain_22_a_thr*100))
# file.write('(%)\n')
# file.write('average strain_22 amplitude in BTG: %0.3f' % (strain_22_a_av*100))
# file.write('(%)\n')
# file.write('highly strained region: %0.3f (mm^2)\n' % (A_HeV*surface_correction))
# file.write('highly strained volume: %0.3f (mm^3)\n' % (A_HeV*t))
      
# file.close()
