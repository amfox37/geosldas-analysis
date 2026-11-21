#!/usr/bin/env python3
"""
Part A3: IMS check of VIIRS (and MODIS) snow cover, per VIIRS Basic_QA class (plan: docs/viirs_obs_study_plan.md).

IMS (NOAA/NIC Interactive Multisensor Snow and Ice, 4 km, Northern Hemisphere only) is a snow / no-snow
reference: an analyst product, not derived from the SCF files used here, but it draws on many sources including
MODIS and VIIRS imagery, so it is not fully independent.
Files: get_ims.py. File ims<YYYY><DDD>_00UTC is paired with obs day <YYYY><DDD> by name, as projects/IMS does
(its header says "Julian day of IMS data log: <DDD-1>", so there may be a one-day offset; checked with --shift 1).

IMS pixel -> 0.05 deg CMG cell: the nearest IMS pixel to the cell center (polar stereographic, true at 60N, lon0 = -80,
a = 6378137, 1/f = 291.505347349177; pixel centers at x = 4000*c - 12286000 m, y = 12286000 - 4000*r; checked against
the IMS lat/lon file to 4 m). Orientation: the ASCII rows start at the lower left, so they are flipped.
IMS categories: 1 sea, 2 land, 3 sea ice, 4 snow (0 outside). Only land and snow pixels are used (sea / sea ice cells
are dropped). Two IMS states per CMG cell:
  nn   : nearest pixel is snow
  dec  : 3x3 pixel window around it is all snow or all no snow among its land pixels (state = snow); others dropped
         (removes the edge of the snow line, where 4 km and 5 km pixels cannot be expected to agree)

Southern Hemisphere cells (and P6, SH winter) have no IMS reference: only cells north of the equator are used.

Counts written to <data>/viirs_obs_study/out/ims/:
  vs_ims.csv     per period x product x QA class x latitude band x threshold x IMS def: 2x2 counts of
                 (obs snow, IMS snow) for all kept cells with an IMS state. Products: VIIRS (QA class 0..3, the
                 reader drops class 2) and MODIS Aqua / Terra (class -1) as a benchmark.
                 thresholds: 'p' obs SCF > 0 (adds snow), 't' obs SCF >= SCF_ANA_BETA = 0.55 (below it the update removes snow)
  disagree.csv   per period x VIIRS product x QA class x band x threshold x IMS def x category (both snow / VIIRS only /
                 MODIS Aqua only / neither): cells where VIIRS and MODIS Aqua are both kept, with the number of
                 these cells where IMS has snow

usage: ims_check.py [--workers N] [--shift DAYS] [P1 P2 ...]     (default: all periods, 4 workers, shift 0)
"""

import argparse
import csv
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from viirs_snow_lib import scf_qc  # noqa: E402
from periods import DATA_DIR, MODIS, PERIODS, VIIRS, cmg_file, days, products  # noqa: E402
from coverage import BANDS, band_index  # noqa: E402
from get_ims import ims_file  # noqa: E402

OUT = DATA_DIR.parent / 'out' / 'ims'
IMS_N = 6144
IMS_EDGE = 12286000.0
IMS_DX = 4000.0
SCF_ANA_BETA = np.float32(0.55)
NH_ROWS = scf_qc.CMG_N_LAT // 2          # CMG rows are north first: rows 0..1799 are north of the equator
NB = 4                                   # NH latitude bands of BANDS: N60-90, N45-60, N30-45, N00-30
THR = ('p', 't')
DEFS = ('nn', 'dec')
CATS = ('both', 'v', 'm', 'none')        # (VIIRS snow, MODIS snow): both, VIIRS only, MODIS only, neither


def cell_pixels():
    """Row / column (top-first, as the lat/lon file) of the IMS pixel nearest to each NH CMG cell center."""
    cache = DATA_DIR.parent / 'ims' / 'cmg_nh_to_ims_pixel.npz'
    if cache.exists():
        d = np.load(cache)
        return d['r'], d['c']
    from pyproj import Transformer
    lon, lat = scf_qc.cmg_lon_lat()
    tr = Transformer.from_crs('EPSG:4326', '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +a=6378137 '
                              '+rf=291.505347349177 +x_0=0 +y_0=0 +units=m', always_xy=True)
    lo, la = np.meshgrid(lon.astype(np.float64), lat[:NH_ROWS].astype(np.float64))
    x, y = tr.transform(lo, la)
    c = np.rint((x + IMS_EDGE) / IMS_DX).astype(np.int32)
    r = np.rint((IMS_EDGE - y) / IMS_DX).astype(np.int32)
    ok = (c >= 0) & (c < IMS_N) & (r >= 0) & (r < IMS_N)
    r[~ok] = -1
    c[~ok] = -1
    np.savez_compressed(cache, r=r, c=c)
    return r, c


def read_ims(path):
    """IMS category grid (6144 x 6144 uint8, top row first, as the lat/lon file), or None if unreadable."""
    with gzip.open(path, 'rt', errors='replace') as f:
        rows = [ln for ln in f.read().split('\n') if len(ln) == IMS_N and ln.isdigit()]
    if len(rows) != IMS_N:
        return None
    a = np.frombuffer(''.join(rows).encode(), dtype=np.uint8).reshape(IMS_N, IMS_N) - 48
    return a[::-1]


def ims_state(path, r, c):
    """For NH CMG cells: dict nn / dec -> (valid, snow) boolean arrays (1800, 7200)."""
    from scipy.ndimage import uniform_filter
    a = read_ims(path)
    if a is None:
        return None
    land = (a == 2) | (a == 4)
    snow = a == 4
    ok = r >= 0
    rr, cc = np.where(ok, r, 0), np.where(ok, c, 0)
    out = {'nn': (land[rr, cc] & ok, snow[rr, cc] & ok)}
    nland = uniform_filter(land.astype(np.float32), size=3, mode='nearest')
    nsnow = uniform_filter(snow.astype(np.float32), size=3, mode='nearest')
    fs = np.where(nland > 0, nsnow / np.maximum(nland, 1e-6), np.nan)[rr, cc]
    eps = np.float32(1e-4)
    decisive = ok & land[rr, cc] & ((fs < eps) | (fs > 1 - eps))
    out['dec'] = (decisive, decisive & (fs > 0.5))
    return out


def _cont(idx, obs_snow, ims_snow, size):
    """2x2 counts per index: (obs snow, IMS snow) -> bins: both, obs only, IMS only, neither."""
    k = 2 * (~obs_snow).astype(np.int64) + (~ims_snow).astype(np.int64)
    return np.bincount(idx * 4 + k, minlength=size * 4).reshape(size, 4)


def process(args):
    period, day, shift = args
    from datetime import timedelta
    path = ims_file(day + timedelta(days=shift))
    if not path.exists():
        return period, day, None, None
    r, c = cell_pixels()
    st = ims_state(path, r, c)
    if st is None:
        return period, day, None, None
    lat = scf_qc.cmg_lon_lat()[1]
    band_row = band_index(lat[:NH_ROWS])                       # (1800,)
    band_all = np.broadcast_to(band_row[:, None], (NH_ROWS, scf_qc.CMG_N_LON))
    vs, dis = {}, {}

    def scf_of(prod, f):
        if prod in VIIRS:
            raw = scf_qc.read_viirs_cmg(f)
            keep, qa = scf_qc.viirs_qc_without_qa(raw)
            keep = keep & (qa >= 0) & (qa < 4)
            sc, ci = raw['Snow_Cover'].astype(np.float32), raw['Clear_Index'].astype(np.float32)
            scf = np.full(sc.shape, np.nan, np.float32)
            scf[keep] = np.minimum(sc[keep] / ci[keep], np.float32(scf_qc.QC_VIIRS['scf_cap']))
            return scf, keep, qa.astype(np.int64)
        scf, keep = scf_qc.modis_scf(scf_qc.read_modis_cmg(f))
        return scf, keep, np.full(keep.shape, -1, np.int64)

    data = {}
    for prod in products(period):
        f = cmg_file(prod, day)
        if f.exists():
            data[prod] = scf_of(prod, f)

    for prod, (scf, keep, qa) in data.items():
        scf, keep, qa = scf[:NH_ROWS], keep[:NH_ROWS], qa[:NH_ROWS]
        for d in DEFS:
            valid, isnow = st[d]
            m = keep & valid & (band_all >= 0) & (band_all < NB)
            b, q, x, s = band_all[m].astype(np.int64), qa[m] + 1, scf[m], isnow[m]     # qa -1..3 -> 0..4
            for t, thr in enumerate(THR):
                v = x > 0 if thr == 'p' else x >= SCF_ANA_BETA
                idx = b * 5 + q
                vs[(prod, d, thr)] = _cont(idx, v, s, NB * 5)

    # disagreement with MODIS Aqua (VIIRS products only)
    if 'MYD10C1' in data:
        mscf, mkeep, _ = data['MYD10C1']
        mscf, mkeep = mscf[:NH_ROWS], mkeep[:NH_ROWS]
        for prod in VIIRS:
            if prod not in data:
                continue
            scf, keep, qa = data[prod]
            scf, keep, qa = scf[:NH_ROWS], keep[:NH_ROWS], qa[:NH_ROWS]
            for d in DEFS:
                valid, isnow = st[d]
                m = keep & mkeep & valid & (band_all >= 0) & (band_all < NB)
                b, q, x, y, s = band_all[m].astype(np.int64), qa[m], scf[m], mscf[m], isnow[m]
                for t, thr in enumerate(THR):
                    v = x > 0 if thr == 'p' else x >= SCF_ANA_BETA
                    w = y > 0 if thr == 'p' else y >= SCF_ANA_BETA
                    cat = 2 * (~v).astype(np.int64) + (~w).astype(np.int64)        # 0 both, 1 v, 2 m, 3 none
                    idx = (b * 4 + q) * 4 + cat
                    size = NB * 4 * 4
                    n = np.bincount(idx, minlength=size)
                    ns = np.bincount(idx, weights=s.astype(np.float64), minlength=size)
                    dis[(prod, d, thr)] = np.stack([n, ns], axis=-1).reshape(NB, 4, 4, 2)
    return period, day, vs, dis


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('periods', nargs='*', default=list(PERIODS))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--shift', type=int, default=0, help='use the IMS file of day + SHIFT')
    a = ap.parse_args()
    out = OUT if a.shift == 0 else OUT.with_name(f'ims_shift{a.shift}')
    out.mkdir(parents=True, exist_ok=True)
    cell_pixels()                                           # build the cache before the workers start
    jobs = [(p, d, a.shift) for p in a.periods for d in days(p)]
    vs_tot, dis_tot, missing = {}, {}, []
    with ProcessPoolExecutor(a.workers) as ex:
        for period, day, vs, dis in ex.map(process, jobs):
            if vs is None:
                missing.append((period, day.isoformat()))
                continue
            for key, v in vs.items():
                vs_tot[(period,) + key] = vs_tot.get((period,) + key, 0) + v
            for key, v in dis.items():
                dis_tot[(period,) + key] = dis_tot.get((period,) + key, 0) + v
            print(period, day, 'ok', flush=True)

    with open(out / 'vs_ims.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['period', 'product', 'ims_def', 'thr', 'band', 'qa', 'both', 'obs_only', 'ims_only', 'neither'])
        for (period, prod, d, thr), arr in sorted(vs_tot.items()):
            for bi in range(NB):
                for qi in range(5):
                    n = arr[bi * 5 + qi]
                    if n.sum():
                        w.writerow([period, prod, d, thr, BANDS[bi][0], qi - 1, *map(int, n)])
    with open(out / 'disagree.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['period', 'product', 'ims_def', 'thr', 'band', 'qa', 'category', 'n', 'n_ims_snow'])
        for (period, prod, d, thr), arr in sorted(dis_tot.items()):
            for bi in range(NB):
                for qi in range(4):
                    for ci, cat in enumerate(CATS):
                        n, ns = arr[bi, qi, ci]
                        if n:
                            w.writerow([period, prod, d, thr, BANDS[bi][0], qi, cat, int(n), int(ns)])
    print('wrote', out, '; missing IMS days:', missing)


if __name__ == '__main__':
    main()
