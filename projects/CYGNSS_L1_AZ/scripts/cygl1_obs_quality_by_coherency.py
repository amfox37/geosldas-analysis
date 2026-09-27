#!/usr/bin/env python3
"""
CYGNSS L1 obs quality vs raw-granule QC fields (coherency_ratio etc.) for an
assimilating L1 arm on the FIXED operator build.

Join (all exact keys, no value matching):
  DA OFA L1 row (window T, lon, lat)  ==  obs-file row (window of its own
      timestamp, sp_lon, sp_lat as float32)  -- valid since reader fix 11dfdb1,
      which writes the selected obs' exact sp lon/lat into the OFA
  obs-file row (year, day, sc_num, sample_id, ch_id)  ==  QC-pass CSV row
      (CYGNSS_operator/artifacts/out_images/cygnss_qc_m36_window_counts_*)
  DA OFA row (T, species, tilenum, lon, lat)  ==  OL OFA row  -> F_OL

Per obs (DA's scaled obs O throughout):
  inn    = O - F_OL          (DA-independent innovation)
  omf/oma/amf from the DA run itself (Desroziers: R ~ E[(O-A)(O-F)])
  monitor pairing: nearest-in-time SMAP Tbh / CYGNSS L3 innovation (O_DA - F_OL)
      on the same tile within +-12 h -> "does this obs see the same SM error"
      (expected sign: L1 vs SMAP Tb negative, L1 vs L3 positive)

Usage:
  cygl1_obs_quality_by_coherency.py --da-expid DA_L1_full_xc015_fixedop \\
      --start 202001 --end 202006
"""
import argparse
import glob
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

EXPDIR = '/discover/nobackup/projects/land_da/cygl1_operator_test/'
DOMAIN = 'SMAP_EASEv2_M36_GLOBAL'
OBS_ROOT = EXPDIR + 'CYGNSS_L1/'
QC_ROOT = '/discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images/'
OUT_DIR = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/obs_quality_by_coherency/'
EPOCH = datetime(2020, 1, 1)
FILL = 1.e14
QC_COLS = ['sample_id', 'ch_id', 'coherency_state', 'coherency_ratio', 'ddm_snr',
           'sp_rx_gain', 'srtm_slope', 'modis_land_cover', 'pekel_sp_water_percentage_5km',
           'brcs_peak', 'brcs_crop_sum']


def species_ids(nc):
    ids = nc['obsparam_species_id'][:]
    de = [str(d) for d in nc['obsparam_descr'][:]]
    grp = {}
    for i, d in zip(ids, de):
        if d.startswith('CYGNSS_L1'):
            grp[int(i)] = 'L1'
        elif d in ('SMAP_L1C_Tbh_A', 'SMAP_L1C_Tbh_D'):
            grp[int(i)] = 'SMAPh'
        elif d == 'CYGNSS_SM_6hr':
            grp[int(i)] = 'L3'
    return grp


def read_ofa(fn, cols):
    with Dataset(fn) as nc:
        if len(nc.dimensions['n_obs']) == 0:
            return None
        grp = species_ids(nc)
        d = pd.DataFrame({k: np.asarray(nc[k][:]) for k in ['species', 'tilenum', 'lon', 'lat'] + cols})
    d['grp'] = d.species.map(grp)
    d = d[d.grp.notna()]
    for c in cols:
        d = d[np.abs(d[c]) < FILL]
    t = datetime.strptime(fn.split('.')[-2], '%Y%m%d_%H%Mz')
    d['T'] = (t - EPOCH).total_seconds()
    return d


def load_ofa(expid, month, cols):
    d = f'{EXPDIR}{expid}/output/{DOMAIN}/ana/ens_avg/{month:Y%Y/M%m}/'
    out = [read_ofa(f, cols) for f in sorted(glob.glob(d + '*ldas_ObsFcstAna.*.nc4'))]
    return pd.concat([o for o in out if o is not None], ignore_index=True)


def load_obs_files(month):
    rows = []
    day = month - timedelta(days=1)        # previous day: 00z window pulls its tail
    while day < month + relativedelta(months=1):
        fn = f'{OBS_ROOT}{day:Y%Y/M%m}/cygnss_l1_ddm3x5_crop_scalar_m36_{day:%Y%m%d}_all_cyg.nc4'
        if os.path.exists(fn):
            with Dataset(fn) as nc:
                rows.append(pd.DataFrame({k: np.asarray(nc[k][:]) for k in
                            ['status', 'year', 'day', 'sc_num', 'sample_id', 'ch_id', 'sp_lon', 'sp_lat',
                             'ddm_timestamp_utc_sec', 'sp_inc_angle', 'sp_nearest_tile_distance_km']}))
        day += timedelta(days=1)
    r = pd.concat(rows, ignore_index=True)
    r = r[r.status == 1]
    yr0 = r.year.map(lambda y: (datetime(int(y), 1, 1) - EPOCH).total_seconds())
    t = yr0 + (r.day - 1) * 86400 + r.ddm_timestamp_utc_sec
    r['T'] = np.ceil((t - 5400) / 10800) * 10800
    r['lon'] = r.sp_lon.astype(np.float32)
    r['lat'] = r.sp_lat.astype(np.float32)
    return r.drop_duplicates(['T', 'lon', 'lat'], keep=False)


_qc_cache = {}
def qc_lookup(r):
    out = []
    for (y, dy, sc), g in r.groupby(['year', 'day', 'sc_num']):
        ymd = (datetime(int(y), 1, 1) + timedelta(days=int(dy) - 1)).strftime('%Y%m%d')
        fn = f'{QC_ROOT}cygnss_qc_m36_window_counts_{ymd}_cyg{int(sc):02d}/cygnss_l1_qc_pass_{ymd}_cyg{int(sc):02d}.csv'
        if not os.path.exists(fn):
            continue
        q = pd.read_csv(fn, usecols=QC_COLS).drop_duplicates(['sample_id', 'ch_id'])
        out.append(g.merge(q, on=['sample_id', 'ch_id'], how='inner'))
    return pd.concat(out, ignore_index=True)


def build(a):
    keys = ['T', 'species', 'tilenum', 'lon', 'lat']
    month, last = datetime.strptime(a.start, '%Y%m'), datetime.strptime(a.end, '%Y%m')
    parts = []
    while month <= last:
        da = load_ofa(a.da_expid, month, ['obs', 'obsvar', 'fcst', 'ana'])
        ol = load_ofa(a.ol_expid, month, ['fcst'])
        da = da.drop_duplicates(keys, keep=False)
        ol = ol.drop_duplicates(keys, keep=False)
        m = da.merge(ol[keys + ['fcst']].rename(columns={'fcst': 'fcst_ol'}), on=keys)
        m['inn'] = m.obs - m.fcst_ol
        l1 = m[m.grp == 'L1'].copy()
        n0 = len(l1)
        r = load_obs_files(month)
        l1 = l1.merge(r.drop(columns=['sp_lon', 'sp_lat', 'status']), on=['T', 'lon', 'lat'])
        n1 = len(l1)
        l1 = qc_lookup(l1)
        print(f'{month:%Y-%m}: L1 OFA {n0}, matched obs file {n1}, matched QC {len(l1)}')
        # nearest-in-time monitor innovation on the same tile within +-12 h
        for g in ['SMAPh', 'L3']:
            mon = m[m.grp == g][['tilenum', 'T', 'inn']].sort_values('T')
            l1 = pd.merge_asof(l1.sort_values('T'), mon.rename(columns={'inn': f'inn_{g}'}),
                               on='T', by='tilenum', direction='nearest', tolerance=43200)
        l1['month'] = month.month
        parts.append(l1)
        month += relativedelta(months=1)
    return pd.concat(parts, ignore_index=True)


def summarize(d, col, bins, label):
    d = d.copy()
    d['bin'] = pd.cut(d[col], bins, include_lowest=True)
    rows = []
    for b, g in d.groupby('bin', observed=True):
        omf, oma, amf = g.obs - g.fcst, g.obs - g.ana, g.ana - g.fcst
        s, l = g.dropna(subset=['inn_SMAPh']), g.dropna(subset=['inn_L3'])
        rows.append(dict(bin=str(b), N=len(g), frac=len(g) / len(d),
                         inn_mean=g.inn.mean(), inn_std=g.inn.std(),
                         R_desroz=(oma * omf).mean(), R_assumed=g.obsvar.mean(),
                         HPH_desroz=(amf * omf).mean(),
                         r_SMAPh=np.corrcoef(s.inn, s.inn_SMAPh)[0, 1] if len(s) > 30 else np.nan,
                         N_SMAPh=len(s),
                         r_L3=np.corrcoef(l.inn, l.inn_L3)[0, 1] if len(l) > 30 else np.nan,
                         N_L3=len(l)))
    t = pd.DataFrame(rows)
    print(f'\n--- {label} ---')
    print(t.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    return t


def extra_checks(d, lo, hi):
    """Tail clustering by tile, land cover, SNR x coherency, keep-band summary."""
    out = {}
    d = d.copy()
    d['cls'] = np.where(d.coherency_ratio < lo, 'low', np.where(d.coherency_ratio > hi, 'high', 'mid'))
    print(f'\ncorr(coherency_ratio, ddm_snr) = {d[["coherency_ratio", "ddm_snr"]].corr().iloc[0, 1]:.3f}')
    rows = []
    ntile = d.tilenum.nunique()
    for c in ['low', 'high']:
        share = d.groupby('tilenum').cls.apply(lambda x: (x == c).mean())
        maj = share[share > 0.5].index
        rows.append(dict(tail=c, n_obs=int((d.cls == c).sum()), n_tiles=ntile,
                         tiles_majority_tail=len(maj),
                         share_of_tail_on_majority_tiles=(d[d.tilenum.isin(maj)].cls == c).sum() / (d.cls == c).sum(),
                         share_of_tail_on_top10pct_tiles=d[d.cls == c].tilenum.value_counts().head(ntile // 10).sum()
                         / (d.cls == c).sum()))
    out['tail_tile_clustering'] = pd.DataFrame(rows)
    print(out['tail_tile_clustering'].to_string(index=False))
    lc = pd.crosstab(d.modis_land_cover, d.cls, normalize='columns')
    out['landcover_share_by_class'] = lc.reset_index()
    print(lc[lc.max(axis=1) > 0.02].round(3).to_string())

    def st(g):
        s = g.dropna(subset=['inn_SMAPh'])
        omf, oma = g.obs - g.fcst, g.obs - g.ana
        return pd.Series(dict(N=len(g), inn_mean=g.inn.mean(), inn_std=g.inn.std(),
                              R_desroz=(omf * oma).mean(), R_assumed=g.obsvar.mean(),
                              r_SMAPh=np.corrcoef(s.inn, s.inn_SMAPh)[0, 1],
                              r_L3=g[['inn', 'inn_L3']].dropna().corr().iloc[0, 1]))
    d['snr_tercile'] = pd.qcut(d.ddm_snr, 3, labels=['low', 'mid', 'high'])
    out['snr_tercile_edges'] = pd.DataFrame({'edge': d.ddm_snr.quantile([0, 1 / 3, 2 / 3, 1]).values})
    t = d.groupby(['cls', 'snr_tercile'], observed=True).apply(st).reset_index()
    out['coherency_class_x_snr_tercile'] = t
    print(t.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    kb = pd.DataFrame([dict(subset='all', **st(d)), dict(subset=f'keep {lo}-{hi}', **st(d[d.cls == 'mid']))])
    out['keep_band_summary'] = kb
    print(kb.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--da-expid', required=True)
    p.add_argument('--ol-expid', default='OLv8_M36_AZ_fixedop')
    p.add_argument('--start', required=True, help='yyyymm')
    p.add_argument('--end', required=True, help='yyyymm (inclusive)')
    p.add_argument('--rebuild', action='store_true')
    p.add_argument('--keep-lo', type=float, default=0.40)
    p.add_argument('--keep-hi', type=float, default=2.16)
    a = p.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    pq = f'{OUT_DIR}{a.da_expid}_{a.start}_{a.end}_l1_qc_joined.parquet'
    if os.path.exists(pq) and not a.rebuild:
        d = pd.read_parquet(pq)
    else:
        d = build(a)
        d.to_parquet(pq)
        print('wrote', pq)

    pd.set_option('display.width', 250)
    tabs = {}
    qs = np.unique(np.nanquantile(d.coherency_ratio, np.linspace(0, 1, 11)))
    tabs['coherency_deciles_all'] = summarize(d, 'coherency_ratio', qs, 'coherency_ratio deciles, all months')
    tabs['coherency_deciles_JFM'] = summarize(d[d.month <= 3], 'coherency_ratio', qs, 'coherency_ratio deciles, Jan-Mar')
    tabs['coherency_deciles_AMJ'] = summarize(d[d.month >= 4], 'coherency_ratio', qs, 'coherency_ratio deciles, Apr-Jun')
    tabs['coherency_state'] = summarize(d, 'coherency_state', [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 'coherency_state')
    for c in ['ddm_snr', 'srtm_slope', 'sp_inc_angle']:
        qc = np.unique(np.nanquantile(d[c], np.linspace(0, 1, 6)))
        tabs[f'{c}_quintiles'] = summarize(d, c, qc, f'{c} quintiles, all months')
    print('\npekel_sp_water_percentage_5km: fraction nonzero', f'{(d.pekel_sp_water_percentage_5km > 0).mean():.4f}')
    tabs.update(extra_checks(d, a.keep_lo, a.keep_hi))
    tag = f'{a.da_expid}_{a.start}_{a.end}'
    for k, t in tabs.items():
        t.to_csv(f'{OUT_DIR}{tag}_{k}.csv', index=False)
    print(f'\nwrote {len(tabs)} tables to {OUT_DIR}{tag}_*.csv')

if __name__ == '__main__':
    main()
