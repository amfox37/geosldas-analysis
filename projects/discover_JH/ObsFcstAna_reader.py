#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os


# Precision of fortran tag
int_precision = 'int32'
# Precision of data in input file
float_precision = 'float32'
# Precision of data in input file
logical_precision = 'int32'


# Initialize outputs in case file does not exist or is empty
nodata = -9999
date_time = {
    'year': nodata,
    'month': nodata,
    'day': nodata,
    'hour': nodata,
    'min': nodata,
    'sec': nodata,
    'dofyr': nodata,
    'pentad': nodata
}

obs_assim = []
obs_species = []
obs_tilenum = []
obs_lon = []
obs_lat = []
obs_obs = []
obs_obsvar = []
obs_fcst = []
obs_fcstvar = []
obs_ana = []
obs_anavar = []

# Determine machine format
machfmt = 'b'


# Get a list of files with a similar name in a directory
test_data = '/Users/amfox/Desktop/GEOSldas_diagnostics/test_data'
path = test_data
file_name = 'OLv7_M36.ens_avg'
file_ext = '.bin'
files = [file for file in os.listdir(path) if file.startswith(file_name) and file.endswith(file_ext)]

print(files)


# Open each file in turn
mode = 'rb' if machfmt == 'b' else 'rl'

for file in files:
    with open(os.path.join(path, file), mode) as ifp:
        print ('Reading file ', file, '...')
        
        # Read N_obs and time stamp entry
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        N_obs = np.fromfile(ifp, int_precision, 1)
        N_obs = int(N_obs)
        year = np.fromfile(ifp, int_precision, 1)
        month = np.fromfile(ifp, int_precision, 1)
        day = np.fromfile(ifp, int_precision, 1)
        hour = np.fromfile(ifp, int_precision, 1)
        minute = np.fromfile(ifp, int_precision, 1)
        second = np.fromfile(ifp, int_precision, 1)
        dofyr = np.fromfile(ifp, int_precision, 1)
        pentad = np.fromfile(ifp, int_precision, 1)
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        # Populate date_time structure
        date_time = {
            'year': year,
            'month': month,
            'day': day,
            'hour': hour,
            'min': minute,
            'sec': second,
            'dofyr': dofyr,
            'pentad': pentad
        }
        
        # Read observation assim flag
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        tmp_data = np.fromfile(ifp, logical_precision, N_obs)
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        tmp_data2 = np.zeros((N_obs, 1))
        indices = np.where(tmp_data != 0)[0]
        tmp_data2[indices] = 1
        obs_assim = np.append(obs_assim, tmp_data2)
        
        # Read species information
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_species = np.append(obs_species, np.fromfile(ifp, int_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        
        # Read tile number information
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_tilenum = np.append(obs_tilenum, np.fromfile(ifp, int_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read longitude
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_lon = np.append(obs_lon, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read latitude
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_lat = np.append(obs_lat, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        
        # Read observation value
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_obs = np.append(obs_obs, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read observation variance
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_obsvar = np.append(obs_obsvar, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read observation-space model forecast value
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_fcst = np.append(obs_fcst, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read observation-space model forecast variance
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_fcstvar = np.append(obs_fcstvar, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read observation-space analysis value
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_ana = np.append(obs_ana, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)

        # Read observation-space analysis variance
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        obs_anavar = np.append(obs_anavar, np.fromfile(ifp, float_precision, N_obs))
        fortran_tag = np.fromfile(ifp, int_precision, 1)
        
        
# Close file
    ifp.close()
    
print('Total number of obs = ',len(obs_assim))

