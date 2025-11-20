import os
import numpy as np
from datetime import datetime

from dateutil.relativedelta import relativedelta

###############################

expt_name = 'LS_DAv8_M36'

# Define the path directory
path = f'/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/{expt_name}/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg'

# Define the common file name start
file_name_start = f'{expt_name}.ens_avg.ldas_ObsFcstAna.'

start_date = datetime(2020, 1, 1)
end_date = datetime(2020, 3, 1)

start_date_str = start_date.strftime('%Y%m%d')
end_date_str = end_date.strftime('%Y%m%d')

max_tilenum = 112573
max_speciesnum = 13

###############################

def read_obsfcstana_no_datetime(path, file_name, printflag=False):
    # Define precisions
    int_precision = 'int32'
    float_precision = 'float32'
    logical_precision = 'int32'

    # Initialize lists for outputs
    date_time_list = []
    obs_assim_list = []
    obs_species_list = []
    obs_tilenum_list = []
    obs_lon_list = []
    obs_lat_list = []
    obs_obs_list = []
    obs_obsvar_list = []
    obs_fcst_list = []
    obs_fcstvar_list = []
    obs_ana_list = []
    obs_anavar_list = []

    machfmt = 'b'
    file_ext = '.bin'
    # Build full file paths (note: file already includes the path)
    files = [os.path.join(root, file) 
            for root, dirs, files in os.walk(path) 
            for file in files if file.startswith(file_name) and file.endswith(file_ext)]

    if printflag:
        print(files)

    mode = 'rb' if machfmt == 'b' else 'rl'

    for file in files:
        with open(file, mode) as ifp:  # file already includes the path
            if printflag:
                print('Reading file', file, '...')
            
            # Read header and time stamp data
            _ = np.fromfile(ifp, int_precision, 1)  # fortran_tag
            N_obs = int(np.fromfile(ifp, int_precision, 1)[0])
            # Read time components
            year    = np.fromfile(ifp, int_precision, 1)
            month   = np.fromfile(ifp, int_precision, 1)
            day     = np.fromfile(ifp, int_precision, 1)
            hour    = np.fromfile(ifp, int_precision, 1)
            minute  = np.fromfile(ifp, int_precision, 1)
            second  = np.fromfile(ifp, int_precision, 1)
            dofyr   = np.fromfile(ifp, int_precision, 1)
            pentad  = np.fromfile(ifp, int_precision, 1)
            _ = np.fromfile(ifp, int_precision, 1)  # fortran_tag

            # Create a single dictionary for the timestamp info and extend the list
            date_time_tmp = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
                'pentad': pentad
            }
            date_time_list.extend([date_time_tmp] * N_obs) 

            # Read observation assimilation flag
            _ = np.fromfile(ifp, int_precision, 1)
            tmp_data = np.fromfile(ifp, logical_precision, N_obs)
            _ = np.fromfile(ifp, int_precision, 1)
            # Vectorized conversion: nonzero becomes 1, else 0.
            tmp_data2 = (tmp_data != 0).astype(np.int32).reshape(-1, 1)
            obs_assim_list.append(tmp_data2)

            # Read species information
            _ = np.fromfile(ifp, int_precision, 1)
            obs_species_list.append(np.fromfile(ifp, int_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)
            
            # Read tile number information
            _ = np.fromfile(ifp, int_precision, 1)
            obs_tilenum_list.append(np.fromfile(ifp, int_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read longitude
            _ = np.fromfile(ifp, int_precision, 1)
            obs_lon_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read latitude
            _ = np.fromfile(ifp, int_precision, 1)
            obs_lat_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)
            
            # Read observation value
            _ = np.fromfile(ifp, int_precision, 1)
            obs_obs_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read observation variance
            _ = np.fromfile(ifp, int_precision, 1)
            obs_obsvar_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read forecast value
            _ = np.fromfile(ifp, int_precision, 1)
            obs_fcst_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read forecast variance
            _ = np.fromfile(ifp, int_precision, 1)
            obs_fcstvar_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read analysis value
            _ = np.fromfile(ifp, int_precision, 1)
            obs_ana_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

            # Read analysis variance
            _ = np.fromfile(ifp, int_precision, 1)
            obs_anavar_list.append(np.fromfile(ifp, float_precision, N_obs))
            _ = np.fromfile(ifp, int_precision, 1)

    # After processing all files, concatenate lists into numpy arrays
    obs_assim = np.concatenate(obs_assim_list) if obs_assim_list else np.array([])
    obs_species = np.concatenate(obs_species_list) if obs_species_list else np.array([])
    obs_tilenum = np.concatenate(obs_tilenum_list) if obs_tilenum_list else np.array([])
    obs_lon = np.concatenate(obs_lon_list) if obs_lon_list else np.array([])
    obs_lat = np.concatenate(obs_lat_list) if obs_lat_list else np.array([])
    obs_obs = np.concatenate(obs_obs_list) if obs_obs_list else np.array([])
    obs_obsvar = np.concatenate(obs_obsvar_list) if obs_obsvar_list else np.array([])
    obs_fcst = np.concatenate(obs_fcst_list) if obs_fcst_list else np.array([])
    obs_fcstvar = np.concatenate(obs_fcstvar_list) if obs_fcstvar_list else np.array([])
    obs_ana = np.concatenate(obs_ana_list) if obs_ana_list else np.array([])
    obs_anavar = np.concatenate(obs_anavar_list) if obs_anavar_list else np.array([])

    return (obs_assim, obs_species, obs_tilenum, obs_lon, obs_lat, 
            obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar)


# Define the print flag
printflag = False

# Loop over the dates
current_date = start_date

while current_date <= end_date:
    # Define the file name for the current date
    file_name = file_name_start + current_date.strftime('%Y%m')
    
    # Call the read_obsfcstana function for the current file
    try:
        obs_assim, species, tilenum, lon, lat, obs, obsvar, fcst, fcstvar, ana, anavar = read_obsfcstana_no_datetime(path, file_name, printflag)
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        current_date += relativedelta(months=1)
        continue
    except Exception as e:
        print(f"Error processing file {file_name}: {e}")
        current_date += relativedelta(months=1)
        continue
    # Initialize arrays

    # Use obs_assim as a mask, obs_assim = 1 means assimilated
    obs_assim = obs_assim.ravel()
    obs = obs[obs_assim == 1]
    species = species[obs_assim == 1]
    tilenum = tilenum[obs_assim == 1]
    lon = lon[obs_assim == 1]
    lat = lat[obs_assim == 1]
    fcst = fcst[obs_assim == 1]
    ana = ana[obs_assim == 1]    

    obs_cnt  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    obs_sum  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    obs2_sum = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    fcst_sum  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    fcst2_sum = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    ana_sum  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    ana2_sum = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    omf_sum  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    omf2_sum = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    oma_sum  = np.zeros((max_tilenum + 1, max_speciesnum + 1))
    oma2_sum = np.zeros((max_tilenum + 1, max_speciesnum + 1))

    # Calculate the difference between the observation and forecast and observation and analysis
    omf = obs - fcst
    oma = obs - ana 

    # Find unique species values and their number
    unique_species, counts = np.unique(species, return_counts=True)
    num_unique_species = len(unique_species)

    # Find unique tilenum values
    unique_tilenum = np.unique(tilenum)

    # Find the number of unique tilenum values
    num_unique_tilenum = len(unique_tilenum)

    # Print the number of unique tilenum values
    print(f"Read obs for: {current_date.strftime('%Y%m')}")

    # Sort the arrays based on tilenum
    sort_indices = np.argsort(tilenum)
    sorted_tilenum = tilenum[sort_indices]
    sorted_species = species[sort_indices]
    sorted_obs = obs[sort_indices]
    sorted_fcst = fcst[sort_indices]
    sorted_ana = ana[sort_indices]
    sorted_omf = omf[sort_indices]
    sorted_oma = oma[sort_indices]

    # Find the unique tilenum values and their counts
    unique_tilenum, counts = np.unique(sorted_tilenum, return_counts=True)

    # Calculate the indices where the groups should be split
    split_indices = np.cumsum(counts)[:-1]

    # Split the sorted arrays based on the split indices
    tilenum_tile = np.split(sorted_tilenum, split_indices)
    species_tile = np.split(sorted_species, split_indices)
    obs_tile = np.split(sorted_obs, split_indices)
    fcst_tile = np.split(sorted_fcst, split_indices)
    ana_tile = np.split(sorted_ana, split_indices)
    omf_tile = np.split(sorted_omf, split_indices)
    oma_tile = np.split(sorted_oma, split_indices)

    # Loop over the unique tiles

    for i in range(num_unique_tilenum):
        tc = int(tilenum_tile[i][0])  # Current tile number

        # Create a dictionary to store indices for each species in the current tile
        species_indices_dict = {sc: np.where(species_tile[i] == sc)[0] for sc in unique_species}

        for sc in unique_species:
            species_indices = species_indices_dict[sc]

            if len(species_indices) > 0:
                sc = int(sc)  # Current species number
                obs_cnt[tc, sc] += len(species_indices)
                obs_sum[tc, sc] += np.sum(obs_tile[i][species_indices])
                obs2_sum[tc, sc] += np.sum(obs_tile[i][species_indices]**2)
                fcst_sum[tc, sc] += np.sum(fcst_tile[i][species_indices])
                fcst2_sum[tc, sc] += np.sum(fcst_tile[i][species_indices]**2)
                ana_sum[tc, sc] += np.sum(ana_tile[i][species_indices])
                ana2_sum[tc, sc] += np.sum(ana_tile[i][species_indices]**2)
                omf_sum[tc, sc] += np.sum(omf_tile[i][species_indices])
                omf2_sum[tc, sc] += np.sum(omf_tile[i][species_indices]**2)
                oma_sum[tc, sc] += np.sum(oma_tile[i][species_indices])
                oma2_sum[tc, sc] += np.sum(oma_tile[i][species_indices]**2)


    # Write this values out to a npz file
    np.savez(f'{path}/{expt_name}.ens_avg.ldas_ObsFcstAna.summed.{current_date.strftime("%Y%m")}.npz',
             obs_cnt=obs_cnt, obs_sum=obs_sum, obs2_sum=obs2_sum, fcst_sum=fcst_sum, fcst2_sum=fcst2_sum, ana_sum=ana_sum, ana2_sum=ana2_sum, omf_sum=omf_sum, omf2_sum=omf2_sum, oma_sum=oma_sum, oma2_sum=oma2_sum,
             num_unique_tilenum=num_unique_tilenum, num_unique_species=num_unique_species)

    current_date += relativedelta(months=1)