#!/usr/bin/env python3
"""
Two-sided coherency_ratio filter for the CYGNSS L1 preprocessed obs stream.

The fixed GEOSldas reader has no QC hook (it never reads `status` or any
coherency field), so filtering = writing physically subset copies of the daily
obs files, with each kept obs' coefficient support slice carried over and
tile_start re-indexed (same writer logic as thin_cygl1_nested_density_6mo.py).

Steps (subcommands):
  annotate  every obs-file row in [--beg, --end] gets its raw-granule
            coherency_ratio via the exact (year, day, sc_num, sample_id, ch_id)
            key into the QC-pass CSVs; rows with no QC match get NaN (and are
            treated as failing the filter). Saved as parquet, also used as the
            keep-list for the scaling clim (keys: window T, float32 sp lon/lat).
  report    per (window, owner tile): does the reader's selection (nearest
            specular point) stay the same, get REPLACED by an in-band
            next-nearest obs, or get LOST -- i.e. how many assimilated obs the
            filter really costs, and how many OL clim rows it cannot represent.
  write     filtered daily files for [--beg, --end] into --dst.

Usage:
  cygl1_coherency_filter.py annotate --beg 20191231 --end 20221231
  cygl1_coherency_filter.py report   --beg 20200101 --end 20201231 --lo 0.40 --hi 2.16
  cygl1_coherency_filter.py write    --beg 20191231 --end 20210101 --lo 0.40 --hi 2.16 \\
      --dst /discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_coh040_216
"""
import argparse
import os
from datetime import datetime, timedelta

import netCDF4 as nc
import numpy as np
import pandas as pd

SRC_ROOT = '/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/'
QC_ROOT = '/discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images/'
OUT_DIR = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/coherency_filter/'
EPOCH = datetime(2020, 1, 1)
DT = 10800


def src_path(day):
    return f'{SRC_ROOT}{day:Y%Y/M%m}/cygnss_l1_ddm3x5_crop_scalar_m36_{day:%Y%m%d}_all_cyg.nc4'


def days(beg, end):
    d, e = datetime.strptime(beg, '%Y%m%d'), datetime.strptime(end, '%Y%m%d')
    while d <= e:
        yield d
        d += timedelta(days=1)


_qc = {}
def qc_table(ymd, sc):
    k = (ymd, sc)
    if k not in _qc:
        fn = f'{QC_ROOT}cygnss_qc_m36_window_counts_{ymd}_cyg{sc:02d}/cygnss_l1_qc_pass_{ymd}_cyg{sc:02d}.csv'
        _qc[k] = (pd.read_csv(fn, usecols=['sample_id', 'ch_id', 'coherency_ratio'])
                  .drop_duplicates(['sample_id', 'ch_id']) if os.path.exists(fn) else None)
    return _qc[k]


def annotate(a):
    parts = []
    for day in days(a.beg, a.end):
        fn = src_path(day)
        if not os.path.exists(fn):
            print('missing', fn)
            continue
        with nc.Dataset(fn) as f:
            r = pd.DataFrame({k: np.asarray(f[k][:]) for k in
                              ['year', 'day', 'sc_num', 'sample_id', 'ch_id', 'sp_lon', 'sp_lat',
                               'ddm_timestamp_utc_sec', 'sp_nearest_tile_ig', 'sp_nearest_tile_jg',
                               'sp_nearest_tile_distance_km']})
        r['file_day'] = int(day.strftime('%Y%m%d'))
        r['obs_idx'] = np.arange(len(r))
        yr0 = r.year.map(lambda y: (datetime(int(y), 1, 1) - EPOCH).total_seconds())
        t = yr0 + (r.day - 1) * 86400 + r.ddm_timestamp_utc_sec
        r['T'] = np.ceil((t - DT / 2) / DT) * DT       # window center, half-open (T-DT/2, T+DT/2]
        r['lon'] = r.sp_lon.astype(np.float32)
        r['lat'] = r.sp_lat.astype(np.float32)
        out = []
        for (y, dy, sc), g in r.groupby(['year', 'day', 'sc_num']):
            ymd = (datetime(int(y), 1, 1) + timedelta(days=int(dy) - 1)).strftime('%Y%m%d')
            q = qc_table(ymd, int(sc))
            out.append(g.merge(q, on=['sample_id', 'ch_id'], how='left') if q is not None
                       else g.assign(coherency_ratio=np.nan))
        r = pd.concat(out).sort_values('obs_idx')
        parts.append(r.drop(columns=['sp_lon', 'sp_lat']))
        if day.day == 1:
            _qc.clear()
            print(f'{day:%Y-%m}: {len(r)} obs, no QC match {r.coherency_ratio.isna().sum()}', flush=True)
    d = pd.concat(parts, ignore_index=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    pq = f'{OUT_DIR}cygl1_obs_coherency_{a.beg}_{a.end}.parquet'
    d.to_parquet(pq)
    print('wrote', pq, len(d), 'rows; no QC match', int(d.coherency_ratio.isna().sum()))


def load_annot(a):
    pq = sorted(p for p in os.listdir(OUT_DIR) if p.startswith('cygl1_obs_coherency_'))
    d = pd.concat([pd.read_parquet(OUT_DIR + p) for p in pq]).drop_duplicates(['file_day', 'obs_idx'])
    d['keep'] = (d.coherency_ratio >= a.lo) & (d.coherency_ratio <= a.hi)
    return d


def report(a):
    d = load_annot(a)
    lo = (datetime.strptime(a.beg, '%Y%m%d') - EPOCH).total_seconds()
    hi = (datetime.strptime(a.end, '%Y%m%d') + timedelta(days=1) - EPOCH).total_seconds()
    # a row can sit in two daily files (previous-day tail); the reader reads both, keep one
    d = d.drop_duplicates(['T', 'lon', 'lat'])
    d = d[(d['T'] > lo) & (d['T'] <= hi)]
    key = ['T', 'sp_nearest_tile_ig', 'sp_nearest_tile_jg']
    sel = d.sort_values('sp_nearest_tile_distance_km').drop_duplicates(key)
    self_ = d[d.keep].sort_values('sp_nearest_tile_distance_km').drop_duplicates(key)
    m = sel.merge(self_[key + ['obs_idx', 'file_day']], on=key, how='left', suffixes=('', '_f'))
    same = (m.obs_idx == m.obs_idx_f) & (m.file_day == m.file_day_f)
    lost = m.obs_idx_f.isna()
    repl = ~same & ~lost
    n = len(m)
    print(f'{a.beg}-{a.end}, band [{a.lo}, {a.hi}]: raw rows {len(d)}, in band {d.keep.mean():.3f}, no QC {d.coherency_ratio.isna().mean():.4f}')
    print(f'  reader-selected (owner tile x window) obs: {n}')
    print(f'  unchanged {same.sum()} ({same.mean():.3f})  replaced by in-band obs {repl.sum()} ({repl.mean():.3f})'
          f'  lost {lost.sum()} ({lost.mean():.3f})  -> assimilated obs kept {(n - lost.sum()) / n:.3f}')
    m['month'] = pd.to_datetime(m['T'] - 1, unit='s', origin=EPOCH).dt.month
    print(m.assign(same=same, repl=repl, lost=lost).groupby('month')[['same', 'repl', 'lost']].mean().round(3).to_string())


def write(a):
    d = load_annot(a)
    os.makedirs(a.dst, exist_ok=True)
    tot_in = tot_out = 0
    for day in days(a.beg, a.end):
        fn = src_path(day)
        if not os.path.exists(fn):
            print('missing', fn)
            continue
        g = d[d.file_day == int(day.strftime('%Y%m%d'))]
        with nc.Dataset(fn) as src:
            n_src = len(src.dimensions['obs'])
            if len(g) != n_src:
                raise ValueError(f'{fn}: annotation has {len(g)} rows, file {n_src}')
            kept = np.sort(g.obs_idx[g.keep].to_numpy())
            tile_start = src['tile_start'][:]
            tile_count = src['tile_count'][:]
            new_start = np.zeros(len(kept), dtype=tile_start.dtype)
            take, cur = [], 0
            for k, i in enumerate(kept):
                new_start[k] = cur
                s0 = int(tile_start[i])
                take.append(np.arange(s0, s0 + int(tile_count[i])))
                cur += int(tile_count[i])
            take = np.concatenate(take) if take else np.array([], dtype=int)
            dst = f'{a.dst}/{day:Y%Y/M%m}/{os.path.basename(fn)}'
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with nc.Dataset(dst, 'w', format=src.file_format) as o:
                o.createDimension('obs', len(kept))
                o.createDimension('support', len(take))
                for name, var in src.variables.items():
                    if var.dimensions == ('obs',):
                        data = new_start if name == 'tile_start' else var[:][kept]
                    elif var.dimensions == ('support',):
                        data = var[:][take]
                    else:
                        raise ValueError(f'unexpected dims for {name}: {var.dimensions}')
                    v = o.createVariable(name, var.dtype, var.dimensions)
                    v.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
                    v[:] = data
                o.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
                o.setncattr('coherency_filter', f'keep {a.lo} <= coherency_ratio <= {a.hi} (NaN dropped)')
                o.setncattr('coherency_filter_script',
                            'geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/cygl1_coherency_filter.py')
                o.setncattr('coherency_filter_source_file', fn)
                o.setncattr('coherency_filter_n_obs_original', n_src)
                o.setncattr('coherency_filter_n_obs_kept', len(kept))
        tot_in += n_src
        tot_out += len(kept)
        if day.day == 1:
            print(f'{day:%Y%m%d}: {n_src} -> {len(kept)}', flush=True)
    print(f'total {tot_in} -> {tot_out} ({tot_out / tot_in:.3f} kept) into {a.dst}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmd', choices=['annotate', 'report', 'write'])
    p.add_argument('--beg', required=True, help='yyyymmdd (file day)')
    p.add_argument('--end', required=True, help='yyyymmdd inclusive')
    p.add_argument('--lo', type=float, default=0.40)
    p.add_argument('--hi', type=float, default=2.16)
    p.add_argument('--dst')
    a = p.parse_args()
    {'annotate': annotate, 'report': report, 'write': write}[a.cmd](a)


if __name__ == '__main__':
    main()
