#!/usr/bin/env python3
"""
Latitude-masked (tile center lat < 37.5N) monthly spatial-mean OmF/OmA pkl
files, for the 6 gated-dense/coherency-A/coherency-B (DA + OL-cross-mask)
datasets already built by run_cygl1_gated_dense.py / run_cygl1_gated_dense_OL_xmask.py
/ run_cygl1_coherency.py / run_cygl1_coherency_OL_xmask.py.

Restricts the spatial average to tiles where CYGNSS L1 can actually observe
(the AZ-box domain spans ~29.3-39.6N; CYGNSS L1 has no obs north of ~37.5N),
by re-reading the same monthly per-tile sums nc4 files those scripts already
wrote and re-aggregating with a tile mask, instead of the full domain.

Reuses calc_spatial_stats_from_sums's exact logic (postproc_ObsFcstAna.py),
with one inserted masking step. Tile axis in the sums files is LOCAL tile
index (position = OFA obs_tilenum - 1), NOT the global tile_id field in the
tilecoord file -- confirmed empirically (tile_id range 4384-7651 for this
909-tile AZ domain, not sequential/positional). Mask is built directly from
tile_coord['com_lat'] array position, matching that local-index convention.

Output: one pkl per dataset containing 'monthly' (same schema as the
unmasked spatial_stats_*.pkl: O_mean/O_stdv/F_mean/F_stdv/A_mean/A_stdv/
OmF_mean/OmF_stdv/OmA_mean/OmA_stdv/N_data/date_vec, shape (12, n_species))
plus 'full_period' (N-weighted pooled OmF/OmA mean+stdv per species AND per
species group [Tb 1-8, SM 9-12, CygL1 13], pooled across the 12 months).
"""
import sys
sys.path.append('/discover/nobackup/projects/land_da/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/GEOSldas_App/util/shared/python/')
import warnings; warnings.filterwarnings("ignore")
import os
import pickle
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from dateutil.relativedelta import relativedelta

from read_GEOSldas import read_tilecoord

SUM_ROOT = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/postproc_paired_density/'
OUT_PATH = SUM_ROOT + 'stats_output/'

LAT_THRESH = 37.5

TILECOORD_FILE = ('/discover/nobackup/projects/land_da/cygl1_operator_test/'
                   'DAv8_M36_AZ_paired_cygl1_dense_gated/output/SMAP_EASEv2_M36_GLOBAL/'
                   'rc_out/DAv8_M36_AZ_paired_cygl1_dense_gated.ldas_tilecoord.bin')

SPECIES_NAMES = [
    "SMOS_Tbh_A", "SMOS_Tbh_D", "SMOS_Tbv_A", "SMOS_Tbv_D",
    "SMAP_Tbh_A", "SMAP_Tbh_D", "SMAP_Tbv_A", "SMAP_Tbv_D",
    "ASCAT_A", "ASCAT_B", "ASCAT_C", "CYGNSS_SM_6hr",
    "CYGNSS_L1",
]
GROUPS = {"Tb (K)": list(range(0, 8)), "SM (m3/m3)": list(range(8, 12)), "CygL1 (dB)": [12]}

DATASETS = {
    'DA_gated_dense':                     (datetime(2020, 1, 1), datetime(2021, 1, 1)),
    'OL_paired_monitor_xmask_gated_dense': (datetime(2020, 1, 1), datetime(2021, 1, 1)),
    'DA_coherency_screened':               (datetime(2020, 1, 1), datetime(2021, 1, 1)),
    'OL_paired_monitor_xmask_coherency_A': (datetime(2020, 1, 1), datetime(2021, 1, 1)),
    'DA_coherency_randmatch':              (datetime(2020, 1, 1), datetime(2021, 1, 1)),
    'OL_paired_monitor_xmask_coherency_B': (datetime(2020, 1, 1), datetime(2021, 1, 1)),
}


def build_tile_mask():
    tc = read_tilecoord(TILECOORD_FILE)
    mask = tc['com_lat'] < LAT_THRESH   # positional -- axis-0 index IS local tile number - 1
    return mask, tc['N_tile']


def calc_spatial_stats_masked(sum_path, exptag, start_time, end_time, tile_mask):
    """Exact reproduction of postproc_ObsFcstAna.calc_spatial_stats_from_sums,
    with masked-out tiles (tile_mask False) forced to N_data=0 before
    aggregation -- so they contribute nothing to the monthly spatial mean,
    same as a tile with genuinely zero obs that month."""
    var_list = ['obs_obs', 'obs_fcst', 'obs_ana']

    O_mean=[]; O_stdv=[]; F_mean=[]; F_stdv=[]; A_mean=[]; A_stdv=[]
    OmF_mean=[]; OmF_stdv=[]; OmA_mean=[]; OmA_stdv=[]; N_data=[]; date_vec=[]

    current_time = start_time
    while current_time < end_time:
        mo_path = sum_path + '/Y' + current_time.strftime('%Y') + '/M' + current_time.strftime('%m') + '/'
        fnc4_sums = mo_path + exptag + '.ens_avg.ldas_ObsFcstAna_sums.' + current_time.strftime('%Y%m') + '.nc4'

        mdata_sum = {}
        mdata2_sum = {}
        with Dataset(fnc4_sums, 'r') as nc:
            mN_data  = nc.variables['N_data'][:].astype(float)
            moxf_sum = nc.variables['obsxfcst_sum'][:]
            moxa_sum = nc.variables['obsxana_sum'][:]
            for var in var_list:
                mdata_sum[var]  = nc.variables[var + '_sum'][:]
                mdata2_sum[var] = nc.variables[var + '2_sum'][:]

        # --- inserted masking step: zero out excluded tiles (lat >= 37.5N) ---
        mN_data[~tile_mask, :] = 0
        # -----------------------------------------------------------------

        moxf_sum[mN_data == 0] = np.nan
        moxa_sum[mN_data == 0] = np.nan
        for var in var_list:
            mdata_sum[var][mN_data == 0]  = np.nan
            mdata2_sum[var][mN_data == 0] = np.nan

        for var in var_list:
            mN_data[np.isnan(mdata_sum[var])] = np.nan
            mN_data[mN_data == 0] = np.nan

        for var in var_list:
            mdata_sum[var][np.isnan(mN_data)]  = np.nan
            mdata2_sum[var][np.isnan(mN_data)] = np.nan
            moxf_sum[np.isnan(mN_data)] = np.nan
            moxa_sum[np.isnan(mN_data)] = np.nan

        N_data_mo = np.nansum(mN_data, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            OxF_mean = np.nansum(moxf_sum, axis=0) / np.where(N_data_mo > 0, N_data_mo, np.nan)
            OxA_mean = np.nansum(moxa_sum, axis=0) / np.where(N_data_mo > 0, N_data_mo, np.nan)

        data_mean = {}; data2_mean = {}; data_var = {}
        for var in var_list:
            with np.errstate(divide='ignore', invalid='ignore'):
                data_mean[var]  = np.nansum(mdata_sum[var],  axis=0) / np.where(N_data_mo > 0, N_data_mo, np.nan)
                data2_mean[var] = np.nansum(mdata2_sum[var], axis=0) / np.where(N_data_mo > 0, N_data_mo, np.nan)
            data_var[var] = data2_mean[var] - data_mean[var] ** 2

        O_mean_mo = data_mean['obs_obs']
        F_mean_mo = data_mean['obs_fcst']
        A_mean_mo = data_mean['obs_ana']
        O_var = data_var['obs_obs']; F_var = data_var['obs_fcst']; A_var = data_var['obs_ana']

        OmF_mean_mo = O_mean_mo - F_mean_mo
        OmA_mean_mo = O_mean_mo - A_mean_mo
        OmF_stdv_mo = np.sqrt(O_var + F_var - 2 * (OxF_mean - O_mean_mo * F_mean_mo))
        OmA_stdv_mo = np.sqrt(O_var + A_var - 2 * (OxA_mean - O_mean_mo * A_mean_mo))

        N_data.append(N_data_mo)
        O_mean.append(O_mean_mo); O_stdv.append(np.sqrt(O_var))
        F_mean.append(F_mean_mo); F_stdv.append(np.sqrt(F_var))
        A_mean.append(A_mean_mo); A_stdv.append(np.sqrt(A_var))
        OmF_mean.append(OmF_mean_mo); OmF_stdv.append(OmF_stdv_mo)
        OmA_mean.append(OmA_mean_mo); OmA_stdv.append(OmA_stdv_mo)
        date_vec.append(current_time.strftime('%Y%m'))
        current_time = current_time + relativedelta(months=1)

    return {
        'O_mean': np.array(O_mean),   'O_stdv': np.array(O_stdv),
        'F_mean': np.array(F_mean),   'F_stdv': np.array(F_stdv),
        'A_mean': np.array(A_mean),   'A_stdv': np.array(A_stdv),
        'OmF_mean': np.array(OmF_mean), 'OmF_stdv': np.array(OmF_stdv),
        'OmA_mean': np.array(OmA_mean), 'OmA_stdv': np.array(OmA_stdv),
        'N_data': np.array(N_data), 'date_vec': date_vec,
    }


def pooled(N, mean, stdv):
    valid = (N > 0) & ~np.isnan(mean)
    if not valid.any():
        return np.nan, np.nan, 0
    N, mean, stdv = N[valid], mean[valid], stdv[valid]
    Ntot = N.sum()
    pooled_mean = (N * mean).sum() / Ntot
    pooled_var = (N * (stdv ** 2 + mean ** 2)).sum() / Ntot - pooled_mean ** 2
    return pooled_mean, np.sqrt(max(pooled_var, 0)), int(Ntot)


def build_full_period(monthly):
    per_species = {}
    for i, name in enumerate(SPECIES_NAMES):
        mF, sF, nF = pooled(monthly['N_data'][:, i], monthly['OmF_mean'][:, i], monthly['OmF_stdv'][:, i])
        mA, sA, nA = pooled(monthly['N_data'][:, i], monthly['OmA_mean'][:, i], monthly['OmA_stdv'][:, i])
        per_species[name] = {'OmF_mean': mF, 'OmF_stdv': sF, 'OmA_mean': mA, 'OmA_stdv': sA, 'N': nF}

    per_group = {}
    for gname, idxs in GROUPS.items():
        N = monthly['N_data'][:, idxs]; mF = monthly['OmF_mean'][:, idxs]; sF = monthly['OmF_stdv'][:, idxs]
        mA = monthly['OmA_mean'][:, idxs]; sA = monthly['OmA_stdv'][:, idxs]
        valid = (N > 0) & ~np.isnan(mF)
        if valid.any():
            Nf, mFf, sFf = N[valid], mF[valid], sF[valid]
            Ntot = Nf.sum()
            pmF = (Nf * mFf).sum() / Ntot
            psF = np.sqrt(max((Nf * (sFf ** 2 + mFf ** 2)).sum() / Ntot - pmF ** 2, 0))
            mAf, sAf = mA[valid], sA[valid]
            pmA = (Nf * mAf).sum() / Ntot
            psA = np.sqrt(max((Nf * (sAf ** 2 + mAf ** 2)).sum() / Ntot - pmA ** 2, 0))
            per_group[gname] = {'OmF_mean': pmF, 'OmF_stdv': psF, 'OmA_mean': pmA, 'OmA_stdv': psA, 'N': int(Ntot)}
        else:
            per_group[gname] = {'OmF_mean': np.nan, 'OmF_stdv': np.nan, 'OmA_mean': np.nan, 'OmA_stdv': np.nan, 'N': 0}

    return {'per_species': per_species, 'per_group': per_group}


def main():
    tile_mask, n_tile = build_tile_mask()
    n_kept = int(tile_mask.sum())
    print(f'tile mask: {n_kept}/{n_tile} tiles kept (lat<{LAT_THRESH}), '
          f'{n_tile - n_kept} excluded ({100*(n_tile-n_kept)/n_tile:.1f}%)')

    results = {}
    for exptag, (start_time, end_time) in DATASETS.items():
        sum_path = SUM_ROOT + exptag + '/'
        monthly = calc_spatial_stats_masked(sum_path, exptag, start_time, end_time, tile_mask)
        full_period = build_full_period(monthly)
        out = {'monthly': monthly, 'full_period': full_period,
               'mask_info': {'lat_thresh': LAT_THRESH, 'n_tile': n_tile, 'n_kept': n_kept}}

        pkl_file = OUT_PATH + f'spatial_stats_{exptag}_latmask37p5_{start_time.strftime("%Y%m")}_' \
                   f'{(end_time - relativedelta(months=1)).strftime("%Y%m")}.pkl'
        with open(pkl_file, 'wb') as f:
            pickle.dump(out, f)
        print(f'wrote {pkl_file}')
        results[exptag] = out

    # spot-check comparison for gated-dense: full-domain vs lat-masked
    print()
    print('=== spot check: gated-dense full-domain (unmasked, existing pkl) vs lat-masked ===')
    with open(OUT_PATH + 'spatial_stats_DA_gated_dense_202001_202012.pkl', 'rb') as f:
        unmasked = pickle.load(f)
    unmasked_full = build_full_period(unmasked)
    masked_full = results['DA_gated_dense']['full_period']
    for gname in GROUPS:
        u = unmasked_full['per_group'][gname]
        m = masked_full['per_group'][gname]
        print(f'{gname:14s} unmasked OmF_stdv={u["OmF_stdv"]:.4f}(N={u["N"]})   '
              f'latmask OmF_stdv={m["OmF_stdv"]:.4f}(N={m["N"]})')


if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
