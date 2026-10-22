#!/usr/bin/env python3
"""
Per-observation table of CYGNSS L1 operator fit from the fixed-operator open loop.

One row per L1 observation that the OL reader selected (one per owner tile per
3-h window). Built to ask *where the operator does well and where it does not*,
so it keeps everything we may want to stratify by.

Columns, by source:
  OL ObsFcstAna (species CYGNSS_L1_DDM3X5_CROP_SCALAR, unscaled in the OL)
      obs (raw dB), fcst = H(x) on the OL ens-mean state... (ens_avg file),
      fcstvar (ensemble spread of H(x)), T (window centre, s since 2020-01-01)
  L1 obs file (exact join on T, sp lon/lat as float32; valid since reader fix 11dfdb1)
      sample/channel/spacecraft ids, incidence angle, DDM SNR, owner tile ig/jg
      and distance, opacity/attenuation, unmatched coefficient fraction, and
      operator-structure metrics from the sparse support:
        n_support, w_max (largest C_t/sum C_t), n_eff (1/sum w^2),
        w_owner (owner tile's share), owner_is_max
  Raw-granule QC CSV (read straight from the packed tarballs, no extraction)
      coherency_ratio/state, srtm_slope, srtm_dem_alt, land cover, water %, ...
  Scaled obs: GEOSldas z-score scaling re-applied offline at the owner tile and
      pentad of T, with both L1 climatologies (full stream and coherency-screened)
      -> obs_sc_full, obs_sc_coh (+ the clim mean/std so the dynamic-range
      mismatch m_std/o_std is available per obs)
  OL inst3_1d_lndfcstana (forecast at T)
      owner-tile SFMC/RZMC/TSOIL1/TSURF (+ SFMC ens std), and the footprint
      SFMC: sum over support tiles of w_t*SFMC_t / sum of in-domain w_t
      (sfmc_fp, w_in_domain)
  geom_db = 10 log10(sum_t C_t): the state-independent (geometry/terrain/vegetation) part of H;
      fcst - geom_db is the soil-moisture (reflectivity) part
  Monitor innovations (scaled monitor obs from a DA arm that monitors them,
      minus OL fcst; same tile, nearest in time within +-12 h; dt kept)
      inn_SMAPh, inn_SMAPv, inn_SMOSh, inn_SMOSv, inn_ASCAT, inn_L3
  The DA arm's own scaled L1 obs where that arm kept the same obs (obs_da),
      used to check the offline scaling against GEOSldas.

Usage (one year per call; months run sequentially):
  extract_obs_table.py --year 2020
Outputs (gitignored): output/operator_skill/obs_table/ol_l1_obs_<yyyymm>.parquet
                      output/operator_skill/obs_table/qc_slim_<yyyy>.parquet (cache)
"""
import argparse
import glob
import os
import tarfile
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

EXPDIR = '/discover/nobackup/projects/land_da/cygl1_operator_test/'
DOMAIN = 'SMAP_EASEv2_M36_GLOBAL'
OBS_ROOT = EXPDIR + 'CYGNSS_L1/'
CLIMS = {'full': EXPDIR + 'scaling_params/cygnss_l1_z_score_clim/AZ_CYGNSS_L1_zscore_all_pentads.nc4',
         'coh': EXPDIR + 'scaling_params/cygnss_l1_z_score_clim_coh040216/AZ_CYGNSS_L1_zscore_all_pentads.nc4'}
QC_TARS = '/discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images/_tars/'
OUT_DIR = ('/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/'
           'output/operator_skill/obs_table/')
EPOCH = datetime(2020, 1, 1)
FILL = 1.e14
L1_DESCR = 'CYGNSS_L1_DDM3X5_CROP_SCALAR'
MON_GROUPS = {'SMAP_L1C_Tbh_A': 'SMAPh', 'SMAP_L1C_Tbh_D': 'SMAPh',
              'SMAP_L1C_Tbv_A': 'SMAPv', 'SMAP_L1C_Tbv_D': 'SMAPv',
              'SMOS_fit_Tbh_A': 'SMOSh', 'SMOS_fit_Tbh_D': 'SMOSh',
              'SMOS_fit_Tbv_A': 'SMOSv', 'SMOS_fit_Tbv_D': 'SMOSv',
              'ASCAT_HSAF_META_SM': 'ASCAT', 'ASCAT_HSAF_METB_SM': 'ASCAT',
              'ASCAT_HSAF_METC_SM': 'ASCAT', 'CYGNSS_SM_6hr': 'L3', L1_DESCR: 'L1'}
QC_COLS = ['sample_id', 'ch_id', 'coherency_state', 'coherency_ratio', 'ddm_snr', 'sp_rx_gain',
           'srtm_slope', 'srtm_dem_alt', 'modis_land_cover', 'pekel_sp_water_percentage_5km',
           'sp_land_confidence', 'brcs_peak', 'brcs_crop_sum']
OBS_VARS = ['status', 'year', 'day', 'sc_num', 'sample_id', 'ch_id', 'sp_lon', 'sp_lat',
            'ddm_timestamp_utc_sec', 'sp_inc_angle', 'incidence_angle_deg', 'ddm_snr_db',
            'sp_nearest_tile_ig', 'sp_nearest_tile_jg', 'sp_nearest_tile_distance_km',
            'observed_y_db', 'tile_count', 'tile_start', 'unmatched_coefficient_fraction',
            'vegetation_geometry_factor', 'coefficient_weighted_opacity',
            'coefficient_weighted_attenuation', 'active_pixel_count']
STATE_VARS = ['SFMC_FCST', 'SFMC_FCST_ENSSTD', 'RZMC_FCST', 'TSOIL1_FCST', 'TSURF_FCST']


# Obs files that the existing runs read but that have since been replaced in CYGNSS_L1/.
# July 2020: the pre-QC-gate files were swapped for the QC-gated rerun on 2026-10-12 (see the README
# in the backup dir); every run made before then (OL + all fixedop arms) read the backup version.
OBS_OVERRIDES = {'202007': EXPDIR + 'CYGNSS_L1_backup/Y2020_M07_pregate_20260629/'}


def obs_file(day):
    fn = f'cygnss_l1_ddm3x5_crop_scalar_m36_{day:%Y%m%d}_all_cyg.nc4'
    if f'{day:%Y%m}' in OBS_OVERRIDES:
        return OBS_OVERRIDES[f'{day:%Y%m}'] + fn
    return f'{OBS_ROOT}{day:Y%Y/M%m}/' + fn


def window_T(t_utc_since_epoch):
    """Centre of the 3-h analysis window containing t (windows are centred on 00,03,...z)."""
    return np.ceil((t_utc_since_epoch - 5400) / 10800) * 10800


def pentad_of_year(t):
    """GEOSldas pentad_of_year (LDAS_DateTimeMod): leap years drop Feb 29 from the count."""
    doy = t.timetuple().tm_yday
    leap = t.year % 4 == 0 and (t.year % 100 != 0 or t.year % 400 == 0)
    return (doy - 2) // 5 + 1 if (leap and doy >= 59) else (doy - 1) // 5 + 1


# ---------------------------------------------------------------- obs files
def load_obs_day(day):
    """Per-obs attributes + operator-structure metrics from one daily coefficient file."""
    fn = obs_file(day)
    if not os.path.exists(fn):
        return None
    with Dataset(fn) as nc:
        r = pd.DataFrame({k: np.asarray(nc[k][:]) for k in OBS_VARS})
        sig = (np.asarray(nc['tile_ig'][:]), np.asarray(nc['tile_jg'][:]), np.asarray(nc['coefficient_weight'][:]),
               np.asarray(nc['coefficient'][:]))
    sup = pd.DataFrame({'row': np.repeat(np.arange(len(r)), r.tile_count.values),
                        'ig': sig[0], 'jg': sig[1], 'w': sig[2], 'C': sig[3]})
    sup = sup.merge(r[['sp_nearest_tile_ig', 'sp_nearest_tile_jg']].rename_axis('row').reset_index(), on='row')
    sup['is_owner'] = (sup.ig == sup.sp_nearest_tile_ig) & (sup.jg == sup.sp_nearest_tile_jg)
    g = sup.groupby('row')
    r['n_support'] = r.tile_count
    # geometry term of H = 10 log10(sum C_t R_t) = 10 log10(sum C_t) + 10 log10(sum w_t R_t)
    r['geom_db'] = 10. * np.log10(g.C.sum().reindex(r.index))
    r['w_max'] = g.w.max().reindex(r.index)
    r['n_eff'] = (1. / g.w.apply(lambda w: (w ** 2).sum())).reindex(r.index)
    r['w_owner'] = sup[sup.is_owner].groupby('row').w.sum().reindex(r.index).fillna(0.)
    r['owner_is_max'] = (r.w_owner >= r.w_max - 1e-12)
    r['file_day'] = day.strftime('%Y%m%d')
    r['row'] = np.arange(len(r))
    sup['file_day'] = r.file_day.iloc[0]
    return r, sup[['file_day', 'row', 'ig', 'jg', 'w']]


# ------------------------------------------------------------- QC from tars
def qc_slim_for_year(year, keys, cache):
    """Pull the QC columns for exactly the (date, sc, sample, ch) in `keys`
    by streaming once through each needed yearly tarball."""
    if os.path.exists(cache):
        return pd.read_parquet(cache)
    keys = keys.drop_duplicates()
    want = {(d, int(s)): g[['sample_id', 'ch_id']] for (d, s), g in keys.groupby(['ymd', 'sc_num'])}
    parts = []
    for ty in sorted({int(d[:4]) for d, _ in want}):
        tf_path = f'{QC_TARS}cygnss_qc_m36_window_counts_{ty}.tar'
        t0, n = time.time(), 0
        with tarfile.open(tf_path) as tf:
            for m in tf:
                if not (m.isfile() and 'cygnss_l1_qc_pass_' in m.name and m.name.endswith('.csv')):
                    continue
                base = os.path.basename(m.name)[len('cygnss_l1_qc_pass_'):-4]   # yyyymmdd_cygNN
                ymd, sc = base.split('_cyg')
                k = (ymd, int(sc))
                if k not in want:
                    continue
                q = pd.read_csv(tf.extractfile(m), usecols=QC_COLS).drop_duplicates(['sample_id', 'ch_id'])
                q = q.merge(want[k], on=['sample_id', 'ch_id'])
                q['ymd'], q['sc_num'] = ymd, int(sc)
                parts.append(q)
                n += 1
        print(f'  QC tar {ty}: read {n} CSVs in {time.time() - t0:.0f} s', flush=True)
    out = pd.concat(parts, ignore_index=True)
    out.to_parquet(cache)
    return out


# ------------------------------------------------------------- GEOSldas output
def read_ofa(fn, cols):
    with Dataset(fn) as nc:
        if len(nc.dimensions['n_obs']) == 0:
            return None
        grp = {int(i): MON_GROUPS.get(str(d)) for i, d in
               zip(nc['obsparam_species_id'][:], nc['obsparam_descr'][:])}
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
    out = pd.concat([o for o in out if o is not None], ignore_index=True)
    return out.drop_duplicates(['T', 'species', 'tilenum', 'lon', 'lat'], keep=False)


def load_state(expid, month):
    d = f'{EXPDIR}{expid}/output/{DOMAIN}/cat/ens_avg/{month:Y%Y/M%m}/'
    parts = []
    for f in sorted(glob.glob(d + '*inst3_1d_lndfcstana_Nt.*.nc4')):
        t = datetime.strptime(f.split('.')[-2], '%Y%m%d_%H%Mz')
        with Dataset(f) as nc:
            s = pd.DataFrame({k: np.asarray(nc[k][0, :]) for k in STATE_VARS})
            s['ig'], s['jg'] = np.asarray(nc['IG'][:]), np.asarray(nc['JG'][:])
        s['T'] = (t - EPOCH).total_seconds()
        parts.append(s)
    s = pd.concat(parts, ignore_index=True)
    for c in STATE_VARS:
        s.loc[np.abs(s[c]) > FILL, c] = np.nan
    return s


# ------------------------------------------------------------- scaling
class Clim:
    """GEOSldas scale_obs_cygl1scal_zscore, re-implemented (owner tile, pentad of T)."""
    def __init__(self, fn, ig, jg):
        self.i0, self.j0 = int(ig.min()), int(jg.min())
        xs, ys = slice(self.i0, int(ig.max()) + 1), slice(self.j0, int(jg.max()) + 1)
        with Dataset(fn) as nc:
            assert list(nc['x'][:3]) == [0, 1, 2] and list(nc['y'][:3]) == [0, 1, 2]
            self.ndata_min, self.eps = float(nc.ndata_min), float(nc.std_epsilon)
            g = lambda v: np.ma.filled(nc[v][:, ys, xs].astype(float), np.nan)
            self.v = {k: g(k) for k in ['o_mean', 'o_std', 'm_mean', 'm_std', 'n_data']}
            self.mmin = np.ma.filled(nc['m_min'][ys, xs].astype(float), np.nan)
            self.mmax = np.ma.filled(nc['m_max'][ys, xs].astype(float), np.nan)

    def apply(self, obs, ig, jg, pentad):
        p, j, i = pentad - 1, jg - self.j0, ig - self.i0
        v = {k: a[p, j, i] for k, a in self.v.items()}
        lo, hi = self.mmin[j, i], self.mmax[j, i]
        ok = (np.isfinite(pd.DataFrame(v)).all(axis=1).values & np.isfinite(lo) & np.isfinite(hi)
              & (v['o_std'] > self.eps) & (v['m_std'] > self.eps)
              & (v['n_data'] >= self.ndata_min) & (lo <= hi))
        sc = v['m_mean'] + v['m_std'] / v['o_std'] * (obs - v['o_mean'])
        sc = np.where(ok, np.clip(sc, lo, hi), np.nan)
        return sc, v


# ------------------------------------------------------------- build
def build_month(month, a, qc, clims):
    t0 = time.time()
    keys = ['T', 'species', 'tilenum', 'lon', 'lat']
    ol = load_ofa(a.ol_expid, month, ['obs', 'obsvar', 'fcst', 'fcstvar'])
    mon = load_ofa(a.mon_expid, month, ['obs'])
    l1 = ol[ol.grp == 'L1'].drop(columns=['grp', 'species']).copy()
    n0 = len(l1)

    # obs-file attributes (previous day too: the 00z window pulls its tail)
    days, rows, sups = [], [], []
    day = month - timedelta(days=1)
    while day < month + relativedelta(months=1):
        x = load_obs_day(day)
        if x is not None:
            rows.append(x[0]); sups.append(x[1])
        day += timedelta(days=1)
    r, sup = pd.concat(rows, ignore_index=True), pd.concat(sups, ignore_index=True)
    r = r[r.status == 1]
    yr0 = r.year.map(lambda y: (datetime(int(y), 1, 1) - EPOCH).total_seconds())
    r['T'] = window_T(yr0 + (r.day - 1) * 86400 + r.ddm_timestamp_utc_sec)
    r['lon'], r['lat'] = r.sp_lon.astype(np.float32), r.sp_lat.astype(np.float32)
    r = r.drop_duplicates(['T', 'lon', 'lat'], keep=False)
    l1 = l1.merge(r.drop(columns=['status', 'sp_lon', 'sp_lat', 'tile_start']), on=['T', 'lon', 'lat'])
    n1 = len(l1)
    l1['raw_mismatch'] = (l1.obs - l1.observed_y_db).abs()

    # QC
    l1['ymd'] = [(datetime(int(y), 1, 1) + timedelta(days=int(d) - 1)).strftime('%Y%m%d')
                 for y, d in zip(l1.year, l1.day)]
    l1 = l1.merge(qc, on=['ymd', 'sc_num', 'sample_id', 'ch_id'], how='left')
    n2 = int(l1.coherency_ratio.notna().sum())

    # scaling with both climatologies
    pent = np.array([pentad_of_year(EPOCH + timedelta(seconds=float(t))) for t in l1['T']])
    l1['pentad'] = pent
    ig, jg = l1.sp_nearest_tile_ig.values.astype(int), l1.sp_nearest_tile_jg.values.astype(int)
    for name, c in clims.items():
        sc, v = c.apply(l1.obs.values, ig, jg, pent)
        l1[f'obs_sc_{name}'] = sc
        for k in ['o_mean', 'o_std', 'm_mean', 'm_std']:
            l1[f'{k}_{name}'] = v[k]

    # state at T: owner tile and footprint-weighted SFMC
    st = load_state(a.ol_expid, month)      # window T files live in T's own month dir
    l1 = l1.merge(st.rename(columns={'ig': 'sp_nearest_tile_ig', 'jg': 'sp_nearest_tile_jg'}),
                  on=['T', 'sp_nearest_tile_ig', 'sp_nearest_tile_jg'], how='left')
    fp = sup.merge(l1[['file_day', 'row', 'T']], on=['file_day', 'row'])
    fp = fp.merge(st[['T', 'ig', 'jg', 'SFMC_FCST']], on=['T', 'ig', 'jg'], how='left')
    fp['w_in'] = np.where(fp.SFMC_FCST.notna(), fp.w, 0.)
    fp['wsm'] = fp.w_in * fp.SFMC_FCST.fillna(0.)
    agg = fp.groupby(['file_day', 'row'])[['w_in', 'wsm']].sum().reset_index()
    agg['sfmc_fp'] = agg.wsm / agg.w_in.replace(0, np.nan)
    l1 = l1.merge(agg.rename(columns={'w_in': 'w_in_domain'})[['file_day', 'row', 'w_in_domain', 'sfmc_fp']],
                  on=['file_day', 'row'], how='left')

    # monitor innovations: scaled monitor obs (DA arm) - OL fcst, same tile, nearest in time
    m = mon.merge(ol[keys + ['fcst']].rename(columns={'fcst': 'fcst_ol'}), on=keys)
    m['inn'] = m.obs - m.fcst_ol
    da_l1 = m[m.grp == 'L1'][['T', 'tilenum', 'lon', 'lat', 'obs']].rename(columns={'obs': 'obs_da'})
    l1 = l1.merge(da_l1, on=['T', 'tilenum', 'lon', 'lat'], how='left')
    l1 = l1.sort_values('T')
    for g in ['SMAPh', 'SMAPv', 'SMOSh', 'SMOSv', 'ASCAT', 'L3']:
        mm = m[m.grp == g][['tilenum', 'T', 'inn']].rename(columns={'inn': f'inn_{g}'})
        mm[f'Tmon_{g}'] = mm['T']
        mm = mm.sort_values('T').drop_duplicates(['tilenum', 'T'])
        l1 = pd.merge_asof(l1, mm, on='T', by='tilenum', direction='nearest', tolerance=43200)
        l1[f'dt_{g}'] = (l1[f'Tmon_{g}'] - l1['T']) / 3600.
        l1 = l1.drop(columns=f'Tmon_{g}')

    l1['inn_raw'] = l1.obs - l1.fcst
    l1['inn_full'] = l1.obs_sc_full - l1.fcst
    l1['inn_coh'] = l1.obs_sc_coh - l1.fcst
    l1['time'] = [EPOCH + timedelta(seconds=float(t)) for t in l1['T']]
    print(f'{month:%Y-%m}: OL L1 {n0}, obs-file {n1}, QC {n2}, scaled full {l1.obs_sc_full.notna().sum()}, '
          f'in DA arm {l1.obs_da.notna().sum()}, max|O-raw| {l1.raw_mismatch.max():.2e}, '
          f'{time.time() - t0:.0f} s', flush=True)
    return l1.drop(columns=['row'])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--year', type=int, required=True)
    p.add_argument('--months', default='1-12', help='e.g. 1-12 or 6-6')
    p.add_argument('--ol-expid', default='OLv8_M36_AZ_fixedop')
    p.add_argument('--mon-expid', default='DA_L1_full_xc015_coh040216_fixedop',
                   help='run whose ObsFcstAna holds the scaled monitor obs')
    p.add_argument('--out-dir', default=OUT_DIR)
    p.add_argument('--overwrite', action='store_true')
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    m0, m1 = (int(x) for x in a.months.split('-'))
    months = [datetime(a.year, m, 1) for m in range(m0, m1 + 1)]

    # QC keys for every obs-file row in the year (+ the day before the first month)
    t0 = time.time()
    keys, day = [], months[0] - timedelta(days=1)
    while day < months[-1] + relativedelta(months=1):
        if os.path.exists(obs_file(day)):
            with Dataset(obs_file(day)) as nc:
                k = pd.DataFrame({v: np.asarray(nc[v][:]) for v in ['year', 'day', 'sc_num', 'sample_id', 'ch_id']})
            keys.append(k)
        day += timedelta(days=1)
    keys = pd.concat(keys, ignore_index=True)
    keys['ymd'] = [(datetime(int(y), 1, 1) + timedelta(days=int(d) - 1)).strftime('%Y%m%d')
                   for y, d in zip(keys.year, keys.day)]
    qc = qc_slim_for_year(a.year, keys[['ymd', 'sc_num', 'sample_id', 'ch_id']],
                          f'{a.out_dir}qc_slim_{a.year}_{m0:02d}{m1:02d}.parquet')
    print(f'QC rows {len(qc)} for {len(keys)} obs-file rows ({time.time() - t0:.0f} s)', flush=True)

    first_state = sorted(glob.glob(f'{EXPDIR}{a.ol_expid}/output/{DOMAIN}/cat/ens_avg/{months[0]:Y%Y/M%m}/'
                                   f'{a.ol_expid}.inst3_1d_lndfcstana_Nt.*.nc4'))[0]   # run start has no 00z file
    with Dataset(first_state) as nc:
        ig, jg = np.asarray(nc['IG'][:]), np.asarray(nc['JG'][:])
    clims = {k: Clim(f, ig, jg) for k, f in CLIMS.items()}

    for month in months:
        out = f'{a.out_dir}ol_l1_obs_{month:%Y%m}.parquet'
        if os.path.exists(out) and not a.overwrite:
            print('exists, skip', out)
            continue
        d = build_month(month, a, qc, clims)
        d.to_parquet(out)


if __name__ == '__main__':
    main()
