#!/usr/bin/env python3
"""
obs_level (task D), step 1: observation-level join of the CYGNSS L1 operator table to what the
DA arm did with each L1 observation, plus the matched DA/OL monitor obs needed for the
monitor-attribution step.

For each month:
  DA ens_avg ObsFcstAna (default DA_L1_full_xc015_coh040216_fixedop), all species:
      obs (scaled), obsvar, fcst, fcstvar, ana, assim_flag
  OL ens_avg ObsFcstAna (OLv8_M36_AZ_fixedop): fcst (+ fcstvar) on the same obs
      matched on (T, species, tilenum, lon, lat) -- cross-mask as in cygl1_noise_gain_split.py;
      readers reused from extract_obs_table.py (load_ofa: fill screening, duplicate-key drop)
  L1 rows are then joined to the per-obs OL operator table ol_l1_obs_<yyyymm>.parquet on
  (T, tilenum, lon, lat) -- the same key extract_obs_table.py used to fill obs_da.

Outputs (one file per year/range, not per day):
  obs_level_l1_<start>_<end>.parquet   one row per L1 obs in the DA arm (with OL fcst + table attrs);
                                       plus OL-table obs NOT in the DA arm (in_da=False) so the
                                       "which obs reach the analysis" question can be answered
  obs_level_mon_<start>_<end>.parquet  one row per monitor obs present in DA and OL
                                       (SMAP/SMOS Tb, L3, ASCAT): T, grp, species, tilenum,
                                       obs_da, obsvar, fcst_da, fcst_ol, fcstvar_da

Usage:
  obs_level_extract.py --start 202006 --end 202006
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from extract_obs_table import load_ofa  # noqa: E402  (reuse the existing OFA reader)

HERE = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/operator_skill/'
KEYS = ['T', 'species', 'tilenum', 'lon', 'lat']
TABLE_COLS = ['T', 'tilenum', 'lon', 'lat', 'obs', 'fcst', 'fcstvar', 'obs_sc_full', 'obs_sc_coh', 'obs_da',
              'inn_full', 'inn_coh', 'inn_raw', 'geom_db', 'SFMC_FCST', 'sfmc_fp', 'sp_inc_angle', 'coherency_ratio',
              'ddm_snr_db', 'ddm_snr', 'sc_num', 'srtm_slope', 'coefficient_weighted_opacity', 'n_eff', 'w_max',
              'w_owner', 'n_support', 'sp_nearest_tile_ig', 'sp_nearest_tile_jg', 'm_std_coh', 'o_std_coh']


def do_month(month, a):
    t0 = time.time()
    da = load_ofa(a.da_expid, month, ['obs', 'obsvar', 'fcst', 'fcstvar', 'ana', 'assim_flag'])
    ol = load_ofa(a.ol_expid, month, ['obs', 'fcst', 'fcstvar'])
    m = da.merge(ol[KEYS + ['obs', 'fcst', 'fcstvar']].rename(
        columns={'obs': 'obs_ol', 'fcst': 'fcst_ol', 'fcstvar': 'fcstvar_ol'}), on=KEYS, how='left')
    m = m.rename(columns={'obs': 'obs_da', 'fcst': 'fcst_da', 'fcstvar': 'fcstvar_da', 'ana': 'ana_da'})

    # ---- L1
    l1 = m[m.grp == 'L1'].drop(columns=['grp'])
    tab = pd.read_parquet(f'{HERE}obs_table/ol_l1_obs_{month:%Y%m}.parquet', columns=TABLE_COLS)
    tab = tab.rename(columns={'obs': 'obs_raw', 'fcst': 'fcst_ol_tab', 'fcstvar': 'fcstvar_ol_tab',
                              'obs_da': 'obs_da_tab'})
    n_da = len(l1)
    j = tab.merge(l1, on=['T', 'tilenum', 'lon', 'lat'], how='outer', indicator=True)
    j['in_da'] = j._merge != 'left_only'
    j['in_tab'] = j._merge != 'right_only'
    j = j.drop(columns='_merge')
    both = j.in_da & j.in_tab
    chk_o = (j.obs_da[both] - j.obs_da_tab[both]).abs().max()
    chk_f = (j.fcst_ol[both] - j.fcst_ol_tab[both]).abs().max()
    print(f'{month:%Y-%m}: DA L1 {n_da}, with OL fcst {l1.fcst_ol.notna().sum()}, OL table {len(tab)}, '
          f'matched {both.sum()}, DA-only {(~j.in_tab).sum()}, assim_flag=1 {int((l1.assim_flag == 1).sum())}; '
          f'max|obs_da-table obs_da| {chk_o:.2e}, max|F_OL-table F| {chk_f:.2e}; '
          f'DA {len(da)} OL {len(ol)} rows, {time.time() - t0:.0f} s', flush=True)
    j['month_dir'] = f'{month:%Y%m}'

    # ---- monitors (present in both DA and OL)
    mon = m[(m.grp != 'L1') & m.fcst_ol.notna()][['T', 'grp', 'species', 'tilenum', 'lon', 'lat', 'obs_da',
                                                    'obsvar', 'fcst_da', 'fcst_ol', 'fcstvar_da', 'assim_flag']]
    for c in ['obs_da', 'obsvar', 'fcst_da', 'fcst_ol', 'fcstvar_da']:
        mon[c] = mon[c].astype(np.float32)
    return j, mon


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--start', required=True, help='yyyymm')
    p.add_argument('--end', required=True, help='yyyymm (inclusive)')
    p.add_argument('--da-expid', default='DA_L1_full_xc015_coh040216_fixedop')
    p.add_argument('--ol-expid', default='OLv8_M36_AZ_fixedop')
    p.add_argument('--out-dir', default=HERE + 'da_link/obs_level/data/')
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    month, last = datetime.strptime(a.start, '%Y%m'), datetime.strptime(a.end, '%Y%m')
    L, M = [], []
    while month <= last:
        j, mon = do_month(month, a)
        L.append(j); M.append(mon)
        month += relativedelta(months=1)
    L, M = pd.concat(L, ignore_index=True), pd.concat(M, ignore_index=True)
    L.to_parquet(f'{a.out_dir}obs_level_l1_{a.start}_{a.end}.parquet')
    M.to_parquet(f'{a.out_dir}obs_level_mon_{a.start}_{a.end}.parquet')
    print(f'wrote {len(L)} L1 rows, {len(M)} monitor rows to {a.out_dir}')


if __name__ == '__main__':
    main()
