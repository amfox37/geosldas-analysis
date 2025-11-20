#!/usr/bin/env python3

"""
GEOSldas Openloop Diagnostics Script

This script computes and aggregates monthly statistics for GEOSldas offline 
simulations (OL mode). It processes observation, forecast, and analysis data 
to generate diagnostic metrics for further analysis.

Calls compute_monthly_stats_OL function that reads in ObsFcstAna files from both 
OL and DA experiments, and crossmasks the OL experiment obervationss with those 
used in the DA experiment. It replaces the obs values in the OL experiment with 
the corresponding _SCALED_ observation values from the DA experiment which can 
then be used in calculating OmF etc.

The script performs the following main tasks:
1. Reads tile and observation parameter information
2. Iterates through monthly data files, computing or loading precomputed statistics
3. Aggregates monthly statistics over the specified time period
4. Saves aggregated metrics for further use

Usage:
    On NASA's Discover:
    $ module load python/GEOSpyD
    $ ./run_compute_monthly_stats_OL.py
     or to run in the background,
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
from netCDF4 import Dataset

from helper.read_GEOSldas import read_tilecoord, read_obs_param
from helper.util import make_folder
from helper.compute_monthly_stats_OL import compute_monthly_stats_OL
from helper.write_nc4 import write_sums_nc4

import warnings; warnings.filterwarnings("ignore")
import sys 
import io

sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

expdir = '/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/'
expid = 'LS_OLv8_M36'
domain = 'SMAP_EASEv2_M36_GLOBAL'

da_expdir='/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_DAv8_M36/'
da_expid='LS_DAv8_M36'

start_time = datetime(2000,6,1)
end_time = datetime(2024,4,1)

# Define a minimum threshold for the temporal data points to ensure statistical reliability
# of the computed metrics. 
Nmin = 20

# Base directory for storing monthly files
# This can be the same as the experiment directory (expdir) or a different location
out_path_mo = '/discover/nobackup/qliu/SMAP_diag/' +expid+'/output/'+domain+'/ana/ens_avg/'

# Directory for diagnostic plots
out_path = '/discover/nobackup/qliu/SMAP_diag/'
make_folder(out_path)

# Variable list for computing sum and sum of squared
var_list = ['obs_obs', 'obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']

# Read tilecoord and obsparam for tile and obs species information
ftc = expdir+expid+'/output/'+domain+'/rc_out/'+expid+'.ldas_tilecoord.bin'
tc = read_tilecoord(ftc)
n_tile = tc['N_tile']

fop = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
obs_param = read_obs_param(fop)
n_spec = len(obs_param)

# Initialize statistical metrics 
data_sum = {}
data2_sum = {}
N_data = np.zeros((n_tile, n_spec))
oxf_sum = np.zeros((n_tile, n_spec))
oxa_sum = np.zeros((n_tile, n_spec))
fxa_sum = np.zeros((n_tile, n_spec))

for var in var_list:
    data_sum[var] = np.zeros((n_tile, n_spec))
    data2_sum[var] = np.zeros((n_tile, n_spec))

# Time loop: processing data at monthly time step
date_time = start_time
while date_time < end_time:
    # File to store monthly statistics    
    fout_path = out_path_mo + '/Y'+ date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'
    make_folder(fout_path)
    
    fout = fout_path + expid+'.ens_avg.ldas_ObsFcstAna.' + date_time.strftime('%Y%m') +'_stats.nc4'

    # Read monthly data if file exists, otherwise compute monthly statistics first   
    if os.path.isfile(fout):
        print('read sums from  monthly file: '+fout)
        mdata_sum = {}
        mdata2_sum = {}
        with Dataset(fout,'r') as nc:
            mN_data = nc.variables['N_data'][:]
            moxf_sum = nc.variables['obsxfcst_sum'][:]
            moxa_sum = nc.variables['obsxana_sum'][:]
            mfxa_sum = nc.variables['fcstxana_sum'][:]
            for var in var_list:
                mdata_sum[var] = nc.variables[var+'_sum'][:]
                mdata2_sum[var] = nc.variables[var+'2_sum'][:]
    else:
        print('compute monthly sums for '+date_time.strftime('%Y%m'))
        mN_data, mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum = \
                 compute_monthly_stats_OL(expdir,expid,domain,date_time,tc,obs_param,var_list, da_expdir, da_expid)
        print('save to monthly file: '+fout)
        write_sums_nc4(fout, mN_data,mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum, obs_param)

    # Aggregate monthly data
    N_data += mN_data
    oxf_sum += moxf_sum
    oxa_sum += moxa_sum
    fxa_sum += mfxa_sum
   
    for var in var_list:
        data_sum[var] += mdata_sum[var] 
        data2_sum[var] += mdata2_sum[var]  
        
    date_time =date_time + relativedelta(months=1)
