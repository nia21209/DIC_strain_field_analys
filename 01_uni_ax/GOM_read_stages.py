# -----------------------------------------------------------------------------
# Read stages from the GOM files writen by Anton Nischler @ LLK 15.11.2020
# Load measured strain field from GOM-System and save the data
# -----------------------------------------------------------------------------

# Import Moduls
import numpy as np
import pandas as pd
import sys


def strain_field_read(specimen_id,first_stage,last_stage):
    
    id_field = {}  # initialize list for field variables

    
    # read raw data from *.gom file with pandas
    
    for i in range(first_stage, last_stage+1,1):
        i_file_id = specimen_id + '-Stufe-0-' + str(i)
        
        print('Read GOM data from stage %s' % i_file_id)
        
        try:
            
            gom_stage_i = pd.read_csv('%s.gom' % (i_file_id), \
                                      header = None, \
                                      skiprows = [0,1,2,3,4,5,6,7,8,9,10,11,12], \
                                      sep = ';')

        except:
            print('Error 1: File %s.gom could not opend!' % (i_file_id))
            sys.exit
            
        # find maximum rows from each column for the nxm data_matrix
        rows = 0
    
        for n in range(0, 13, 1):
            data = gom_stage_i.iloc[:, n]
            data = data.dropna()
    
            rows_max = len(data)
    
            if rows <= rows_max:
                    rows = rows_max
            else:
                pass
        
        # initialyse the data_matrix
        data_matrix = np.zeros((rows, 13))
        
        # replacing missing values in rows with 0.0
        for k in range(0, 13, 1):
            for l in range(0, rows, 1):
                try:
                    float(gom_stage_i.iloc[l, k])
                    data_matrix[l, k] = gom_stage_i.iloc[l, k]  # if number in GOM_ij exists write to data_matrix
                except:
                    data_matrix[l, k] = 0.0  # if no number exists write zero in data_matrix
        
        # replace NaN's with 0.0
        data_matrix = np.nan_to_num(data_matrix, copy=True)  
        
        # append data_matrix (GOM stage data) to list field
        id_field.update({'%s' % str(i):data_matrix})  
                 
    return id_field
