"""
This script computes monthly statistics for observation, forecast, and analysis data
from the GEOSldas system. It processes data for a given month, aggregates statistics
over tiles and observation species, and calculates metrics such as observation-forecast,
observation-analysis, and forecast-analysis correlations. The script also cross-masks
OPENLOOP (OL) and Data Assimilation (DA) observations and replaces the OL observations 
with _SCALED_ observations from the DA experiment so comparable OmF can be calculated later.

Functions:
- compute_monthly_stats_OL: Main function to compute monthly statistics.

Dependencies:
- numpy
- datetime
- dateutil.relativedelta
- helper.read_GEOSldas (custom module for reading GEOSldas data)

Usage:
Called from or run the script directly to compute statistics for a specific experiment setup.

Author: A. M. Fox
Last Modified: Apr 8, 2025
"""


import numpy as np
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from helper.read_GEOSldas import read_ObsFcstAna, read_tilecoord, read_obs_param

def compute_monthly_stats_OL(expdir, expid, domain, this_month, tc, obs_param, var_list, da_expdir, da_expid):

    n_tile = tc['N_tile']
    n_spec = len(obs_param)

    start_time = this_month.replace(day=1,hour=3) 
    end_time = start_time + relativedelta(months=1)

    data_sum = {}
    data2_sum = {}

    N_data = np.zeros((n_tile, n_spec))
    oxf_sum = np.zeros((n_tile, n_spec))
    oxa_sum = np.zeros((n_tile, n_spec))
    fxa_sum = np.zeros((n_tile, n_spec))

    for var in var_list:
        data_sum[var] = np.zeros((n_tile, n_spec))
        data2_sum[var] = np.zeros((n_tile, n_spec))

    date_time = start_time
    while date_time < end_time:

        fname = expdir + expid + '/output/' + domain + '/ana/ens_avg/Y' + \
                date_time.strftime('%Y') + '/M' + \
                date_time.strftime('%m') + '/' + \
                expid + '.ens_avg.ldas_ObsFcstAna.' + \
                date_time.strftime('%Y%m%d_%H%M') + 'z.bin'

        sub_fname = da_expdir + da_expid + '/output/' + domain + '/ana/ens_avg/Y' + \
                    date_time.strftime('%Y') + '/M' + \
                    date_time.strftime('%m') + '/' + \
                    da_expid + '.ens_avg.ldas_ObsFcstAna.' + \
                    date_time.strftime('%Y%m%d_%H%M') + 'z.bin'

        print("fname: ", fname)

        OFA = read_ObsFcstAna(fname)
        OFA_sub = read_ObsFcstAna(sub_fname)

        if len(OFA['obs_tilenum']) > 0 and len(OFA_sub['obs_tilenum']) > 0:
            # Build lookup array for subset obs using structured dtype for fast matching
            subset_keys = np.core.defchararray.add(
                np.core.defchararray.add(OFA_sub['obs_species'].astype(str),
                                         np.round(OFA_sub['obs_lat'], 6).astype(str)),
                np.round(OFA_sub['obs_lon'], 6).astype(str)
            )
            subset_lookup = dict(zip(subset_keys, OFA_sub['obs_obs']))

            # Initialize full size variable to keep values
            data_tile = {var: np.full((n_tile, n_spec), np.nan) for var in var_list}

            for ispec in np.arange(n_spec):
                this_species = int(obs_param[ispec]['species'])
                masked_data = {}
                if obs_param[ispec]['assim'] == 'T':
                    mask = np.logical_and(OFA['obs_species'] == this_species, OFA['obs_assim'] == 1)
                else:
                    mask = OFA['obs_species'] == this_species

                lat = OFA['obs_lat'][mask]
                lon = OFA['obs_lon'][mask]
                species = OFA['obs_species'][mask]
                tilenum = OFA['obs_tilenum'][mask]

                keys = np.core.defchararray.add(
                    np.core.defchararray.add(species.astype(str),
                                             np.round(lat, 6).astype(str)),
                    np.round(lon, 6).astype(str)
                )

                matched_idx = np.array([key in subset_lookup for key in keys])
                if not np.any(matched_idx):
                    continue

                masked_tilenum = tilenum[matched_idx]
                masked_data['obs_obs'] = np.array([subset_lookup[key] for key in keys[matched_idx]])

                for var in var_list:
                    if var != 'obs_obs':
                        masked_data[var] = OFA[var][mask][matched_idx]

                tile_idx = np.where(np.isin(tc['tile_id'], masked_tilenum))[0]

                for var in var_list:
                    data_tile[var][tile_idx, ispec] = masked_data[var]

            is_valid = ~np.isnan(data_tile['obs_obs'])
            N_data[is_valid] += 1
            oxf_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_fcst'][is_valid]
            oxa_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_ana'][is_valid]
            fxa_sum[is_valid] += data_tile['obs_fcst'][is_valid] * data_tile['obs_ana'][is_valid]
            for var in var_list:
                data_sum[var][is_valid] += data_tile[var][is_valid]
                data2_sum[var][is_valid] += data_tile[var][is_valid] ** 2

        date_time = date_time + timedelta(seconds=10800)

    return N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum

if __name__ == '__main__':
    date_time = datetime(2015,5,1)
    expdir = '/gpfsm/dnb05/projects/p51/SMAP_Nature/'
    expid = 'SPL4SM_Vv8010_OPENLOOP'
    da_expdir = '/gpfsm/dnb05/projects/p51/SMAP_Nature/'
    da_expid = 'SPL4SM_Vv8010_DA_RUN'
    domain = 'SMAP_EASEv2_M09_GLOBAL'
    var_list = ['obs_obs', 'obs_obsvar', 'obs_fcst', 'obs_fcstvar', 'obs_ana', 'obs_anavar']
    ftc = expdir+expid+'/output/'+domain+'/rc_out/'+expid+'.ldas_tilecoord.bin'
    tc = read_tilecoord(ftc)

    fop = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
    obs_param = read_obs_param(fop)

    N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum = \
           compute_monthly_stats_OL(expdir, expid, domain, date_time, tc, obs_param, var_list, da_expdir, da_expid)
