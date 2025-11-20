#!/usr/bin/env python3

"""
GEOSldas monthly OmF Script

This script aggregates precomputed monthly statistics from a GEOSldas experiment. 
The statistics are originally calculated using the `compute_monthly_stats_OL` 
function, which cross-masks OL observations with those used in a Data Assimilation (DA) 
run and replaces them with _SCALED_ DA observations to ensure consistent evaluation.

This script does NOT compute raw observation statistics directly. Instead, it:

1. Loads precomputed tile-level monthly statistics from NetCDF files
2. Aggregates these statistics globally by user-defined observation species groups
3. Computes monthly time series of number of obs, OmF (Observation-minus-Forecast) mean and standard deviation
4. Saves the aggregated species-group time series to a new NetCDF file for further diagnostics

Usage:
    On NASA's Discover:
    $ module load python/GEOSpyD
    $ ./run_compute_monthly_stats_OL.py
    Or to run in the background:
    $ nohup ./run_compute_monthly_stats_OL.py > out.log &
    
Requirements:
    - Python 3.x
    - Modules: numpy, netCDF4, dateutil (included in GEOSpyD)

Author: A. M. Fox  
Last Modified: Apr 8, 2025
"""

import numpy as np
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset, stringtochar

from helper.read_GEOSldas import read_tilecoord, read_obs_param
from helper.util import make_folder

import warnings; warnings.filterwarnings("ignore")
import sys 
import io

sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

expdir = '/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/'
expid = 'LS_OLv8_M36'
domain = 'SMAP_EASEv2_M36_GLOBAL'

crossmasked_obs = True  # Set to True if using crossmasked observations

start_time = datetime(2000,6,1)
end_time = datetime(2024,4,1)

# Set the input stats file name based on whether crossmasked observations are used
if crossmasked_obs:
    if "OL" not in expid:
        print("[WARNING] Are you sure you want to use crossmasked observations, not an OL experiment ID?")
        sys.exit(1)
    # Crossmasked observations can be used in the OL experiment
    stats_file_end = '_stats_CROSSMASKED.nc4'
else:
    stats_file_end = '_stats.nc4'


# Base directory for storing monthly files
# This can be the same as the experiment directory (expdir) or a different location
out_path = os.path.join(expdir, expid, 'output', domain, 'ana', 'ens_avg')

# Variable list for computing sum and sum of squared
var_list = ['obs_obs', 'obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']

# Read tilecoord and obsparam for tile and obs species information
ftc = os.path.join(expdir, expid, 'output', domain, 'rc_out', f'{expid}.ldas_tilecoord.bin')
tc = read_tilecoord(ftc)
n_tile = tc['N_tile']

# Construct the file path dynamically using start_time
fop = os.path.join(
    expdir, expid, 'output', domain, 'rc_out',
    'Y' + start_time.strftime('%Y'),
    'M' + start_time.strftime('%m'),
    f"{expid}.ldas_obsparam.{start_time.strftime('%Y%m%d')}_0000z.txt"
)
obs_param = read_obs_param(fop)
n_spec = len(obs_param)

# Define the observation species groups based on the indices
# These indices should correspond to the species in the obs_param file
# TODO: Make this based on obs_param instead of hardcoding

species_groups = {
    "SMOS": [0, 1, 2, 3],
    "SMAP": [4, 5, 6, 7],
    "ASCAT": [8, 9, 10],
    "MODIS": [11, 12]
}
# Number of species groups
n_spec_groups = len(species_groups)

# Initialize storage for time series data
N_data_group_all_months = {group: [] for group in species_groups}
OmF_mean_all_months = {group: [] for group in species_groups}
OmF_stdv_all_months = {group: [] for group in species_groups}
monthly_timestamps = []  # To store the timestamps for each month

# Time loop: processing data at monthly time step
date_time = start_time
while date_time < end_time:

    # File to store monthly statistics    
    fout_path = out_path + '/Y' + date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'
    make_folder(fout_path)
    
    fout = fout_path + expid + '.ens_avg.ldas_ObsFcstAna.' + date_time.strftime('%Y%m') + stats_file_end

    # Read monthly data if file exists
    if os.path.isfile(fout):
        print('Reading sums from monthly file: ' + fout)
        with Dataset(fout, 'r') as nc:
            mN_data = nc.variables['N_data'][:]  # Shape: (tiles, species)
            moxf_sum = nc.variables['obsxfcst_sum'][:]  # Shape: (tiles, species)
            mdata_sum = {var: nc.variables[var + '_sum'][:] for var in var_list}
            mdata2_sum = {var: nc.variables[var + '2_sum'][:] for var in var_list}
    else:
        print('[WARNING] Cannot find monthly file: ' + fout)
        date_time += relativedelta(months=1)
        continue

    # Compute metrics for each species group
    for group, indices in species_groups.items():
        # Sum across tiles for the species in this group
        N_data_group = np.sum(mN_data[:, indices], axis=1)  # Total count for the group
        OmF_sum_group = np.sum(moxf_sum[:, indices], axis=1)  # Total OmF sum for the group

        # Convert N_data_group to float to allow NaN values
        N_data_group = N_data_group.astype(float)

        # Filter out invalid data
        N_data_group[N_data_group == 0] = np.nan

        # Compute group-level means and variances
        data_mean = {var: np.sum(mdata_sum[var][:, indices], axis=1) / N_data_group for var in var_list}
        data2_mean = {var: np.sum(mdata2_sum[var][:, indices], axis=1) / N_data_group for var in var_list}
        data_var = {var: data2_mean[var] - data_mean[var]**2 for var in var_list}

        # Compute OmF_mean and OmF_stdv
        OmF_mean = data_mean['obs_obs'] - data_mean['obs_fcst']
        OmF_stdv = np.sqrt(
            data_var['obs_obs'] + data_var['obs_fcst'] -
            2 * (OmF_sum_group / N_data_group - data_mean['obs_obs'] * data_mean['obs_fcst'])
        )

        # Ensure arrays are writable
        OmF_mean = np.copy(OmF_mean)
        OmF_stdv = np.copy(OmF_stdv)

        # Aggregate results to single values
        N_data_group_total = np.nansum(N_data_group)  # Total N_data for the group
        OmF_mean_avg = np.nanmean(OmF_mean)  # Mean OmF for the group
        OmF_stdv_avg = np.nanmean(OmF_stdv)  # Stdv OmF for the group

        # Store the results for the current month
        N_data_group_all_months[group].append(N_data_group_total)
        OmF_mean_all_months[group].append(OmF_mean_avg)
        OmF_stdv_all_months[group].append(OmF_stdv_avg)

    # Store the current timestamp for the month
    monthly_timestamps.append(date_time)

    # Increment to the next month
    date_time += relativedelta(months=1)

# Save the aggregated results to a new NetCDF file
out_nc_file = out_path + expid + '_monthly_OmF.nc4'

with Dataset(out_nc_file, 'w', format='NETCDF4') as nc:
    # Create dimensions
    nc.createDimension('time', len(monthly_timestamps))
    nc.createDimension('species_group', len(species_groups))

    # Create dimension for string lengths of group names
    max_str_len = max(len(name) for name in species_groups.keys())
    nc.createDimension('name_strlen', max_str_len)

    # Create variables
    time_var = nc.createVariable('time', 'f8', ('time',))
    group_names_var = nc.createVariable('species_group_name', 'S1', ('species_group', 'name_strlen'))
    N_data_var = nc.createVariable('N_data', 'i4', ('time', 'species_group'))
    OmF_mean_var = nc.createVariable('OmF_mean', 'f4', ('time', 'species_group'))
    OmF_stdv_var = nc.createVariable('OmF_stdv', 'f4', ('time', 'species_group'))

    # Assign time data
    time_var[:] = np.array([dt.timestamp() for dt in monthly_timestamps])

    # Assign species group names as string char arrays
    group_name_array = stringtochar(np.array(list(species_groups.keys()), f'S{max_str_len}'))
    group_names_var[:, :] = group_name_array

    # Assign the main data arrays
    group_order = list(species_groups.keys())
    N_data_var[:, :] = np.array([N_data_group_all_months[group] for group in group_order]).T
    OmF_mean_var[:, :] = np.array([OmF_mean_all_months[group] for group in group_order]).T
    OmF_stdv_var[:, :] = np.array([OmF_stdv_all_months[group] for group in group_order]).T

    # Optional metadata
    N_data_var.setncattr('description', 'Total observations per species group per month')
    OmF_mean_var.setncattr('description', 'Monthly global mean of OmF per species group')
    OmF_stdv_var.setncattr('description', 'Monthly global stdv of OmF per species group')
    nc.setncattr('species_group_mapping', str(species_groups))  