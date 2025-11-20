import os
from read_obsfcstana import read_obsfcstana_extend_datetime

# Define the path (no relative paths!) and start of file name
path = '/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/fp_scaled/data/MLT_DA/Y2020/M07/D02'
file_name = 'MLT_DA.ens_avg.ldas_ObsFcstAna.20200701_1200z'

# Call the function
results = read_obsfcstana_extend_datetime(path, file_name, printflag=True)

# Unpack the results
(date_time, obs_species, obs_tilenum, obs_lon, obs_lat, obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar) = results

# Print the first few values from each variable
print("Date Time:", date_time[:5])
print("Obs Species:", obs_species[:5])
print("Obs Tile Number:", obs_tilenum[:5])
print("Obs Longitude:", obs_lon[:5])
print("Obs Latitude:", obs_lat[:5])
print("Obs Observations:", obs_obs[:5])
print("Obs Observation Variance:", obs_obsvar[:5])
print("Obs Forecast:", obs_fcst[:5])
print("Obs Forecast Variance:", obs_fcstvar[:5])
print("Obs Analysis:", obs_ana[:5])
print("Obs Analysis Variance:", obs_anavar[:5])