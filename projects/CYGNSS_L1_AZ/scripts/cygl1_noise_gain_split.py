#!/usr/bin/env python3
"""
Noise-vs-gain split of a DA arm's monitor O-F MSE change relative to an OL.

For every obs present in both the DA and OL ens_avg ObsFcstAna files (matched
on time, species, tile, lon, lat), using the DA run's (scaled) obs O:
    inn  = O - F_OL           (OL innovation on the DA obs)
    dF   = F_DA - F_OL        (what the DA did to the forecast)
    dMSE = (O-F_DA)^2 - (O-F_OL)^2 = dF^2 - 2 dF inn
All terms are reported as % of the OL MSE sum(inn^2):
    noise = sum(dF^2)            (cost of moving the forecast at all)
    gain  = -2 sum(dF inn)       (benefit of moving it toward the obs)
    alpha_opt = sum(dF inn) / sum(dF^2) = -gain / (2 noise)
        -> scale factor on the increments that would minimise MSE
    dMSE_at_alpha_opt = -(sum dF inn)^2 / (sum dF^2 sum inn^2)
Note: MSE here includes the O-F mean, unlike the OmF_stdv scores.

Usage:
  cygl1_noise_gain_split.py --da-expid DA_L1_full_xc015_fixedop \\
      --ol-expid OLv8_M36_AZ_fixedop --start 202004 --end 202006
"""
import argparse
import glob
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta
from datetime import datetime

EXPDIR = '/discover/nobackup/projects/land_da/cygl1_operator_test/'
DOMAIN = 'SMAP_EASEv2_M36_GLOBAL'
FILL = 1.e14

GROUPS = [('SMOS',  lambda d: d.startswith('SMOS')),
          ('SMAP',  lambda d: d.startswith('SMAP')),
          ('ASCAT', lambda d: d.startswith('ASCAT')),
          ('L3',    lambda d: d == 'CYGNSS_SM_6hr'),
          ('L1',    lambda d: d.startswith('CYGNSS_L1'))]


def read_ofa(fname):
    with Dataset(fname) as nc:
        if len(nc.dimensions['n_obs']) == 0:
            return None
        ids = nc['obsparam_species_id'][:]
        descr = nc['obsparam_descr'][:]
        id2grp = {}
        for i, d in zip(ids, descr):
            for g, f in GROUPS:
                if f(str(d)):
                    id2grp[int(i)] = g
        df = pd.DataFrame({k: np.asarray(nc[k][:]).astype(np.float64 if k in
                           ('obs', 'fcst', 'lon', 'lat') else np.int64)
                           for k in ('species', 'tilenum', 'lon', 'lat', 'obs', 'fcst')})
    df['grp'] = df['species'].map(id2grp)
    df = df[(df.obs.abs() < FILL) & (df.fcst.abs() < FILL) & df.grp.notna()]
    df['lon'] = df.lon.round(4)
    df['lat'] = df.lat.round(4)
    return df


def load_month(expid, month):
    d = f"{EXPDIR}{expid}/output/{DOMAIN}/ana/ens_avg/{month:Y%Y/M%m}/"
    files = sorted(glob.glob(d + '*ldas_ObsFcstAna.*.nc4'))
    out = []
    for f in files:
        df = read_ofa(f)
        if df is not None:
            df['t'] = os.path.basename(f).split('.')[-2]
            out.append(df)
    return pd.concat(out, ignore_index=True), len(files)


def split(m):
    inn = m.obs_da - m.fcst_ol
    dF = m.fcst_da - m.fcst_ol
    s_ii, s_dd, s_di = (inn**2).sum(), (dF**2).sum(), (dF * inn).sum()
    return dict(N=len(m),
                dMSE=100 * (s_dd - 2 * s_di) / s_ii,
                noise=100 * s_dd / s_ii,
                gain=-200 * s_di / s_ii,
                rms_dF=np.sqrt(s_dd / len(m)),
                corr=np.corrcoef(dF, inn)[0, 1],
                alpha_opt=s_di / s_dd,
                dMSE_alpha=-100 * s_di**2 / (s_dd * s_ii))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--da-expid', required=True)
    p.add_argument('--ol-expid', default='OLv8_M36_AZ_fixedop')
    p.add_argument('--start', required=True, help='yyyymm')
    p.add_argument('--end', required=True, help='yyyymm (inclusive)')
    p.add_argument('--csv', help='optional output csv path')
    a = p.parse_args()

    keys = ['t', 'species', 'tilenum', 'lon', 'lat']
    month = datetime.strptime(a.start, '%Y%m')
    last = datetime.strptime(a.end, '%Y%m')
    rows, pooled = [], []
    while month <= last:
        da, nda = load_month(a.da_expid, month)
        ol, nol = load_month(a.ol_expid, month)
        if nda != nol:
            print(f'WARNING {month:%Y%m}: {nda} DA files vs {nol} OL files')
        da = da.drop_duplicates(keys, keep=False)
        ol = ol.drop_duplicates(keys, keep=False)
        m = da.merge(ol[keys + ['obs', 'fcst']], on=keys, suffixes=('_da', '_ol'))
        pooled.append(m)
        for g, _ in GROUPS:
            mg = m[m.grp == g]
            if len(mg) > 1:
                rows.append(dict(month=f'{month:%Y-%m}', grp=g, **split(mg)))
        month += relativedelta(months=1)
    m = pd.concat(pooled)
    for g, _ in GROUPS:
        mg = m[m.grp == g]
        if len(mg) > 1:
            rows.append(dict(month='all', grp=g, **split(mg)))

    res = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(f'\n{a.da_expid} vs {a.ol_expid}: % of OL MSE (dMSE = noise + gain)\n')
    print(res.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    if a.csv:
        res.to_csv(a.csv, index=False)
        print('wrote', a.csv)


if __name__ == '__main__':
    main()
