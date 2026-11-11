#!/usr/bin/env python3
"""
Part A1: coverage of VIIRS Basic_QA class 1 (plan: docs/viirs_obs_study_plan.md).

For every VIIRS daily file of every period (periods.py), with all reader tests except
Basic_QA (scf_qc.viirs_qc_without_qa) and the QA class kept as a label:

cells.csv   cells per latitude band x SCF class x QA class (0, 1, 2, 3)
            (class 2 is counted for reference; the reader drops it)
tiles.csv   per tile space (M36, CF0360) and tile latitude band: tiles with obs from
            QA {0,3} and from {0,1,3}; tiles with obs only through class-1 cells; the
            same for snowy tiles (tile-mean SCF with {0,1,3} > 0); tiles with both
            class-1 and {0,3} cells (their tile mean changes if class 1 is dropped)
maps_<period>_<product>.npz
            0.5 deg maps summed over the period's days: kept cells {0,1,3}, class-1
            cells, snowy kept cells, snowy class-1 cells

A tile "has an obs" on a day if at least one kept CMG cell of that day's file maps to
it (lib/obs2tile.py, as GEOSldas). The reader picks the daily file per window by the
13:30 local overpass, so each tile gets its obs from one daily file per day; treating
the daily file as the tile's obs for that day is exact up to the tiles at the date line.

Output: <data>/viirs_obs_study/out/coverage/ (git-ignored).

usage: coverage.py [--workers N] [P1 P2 ...]     (default: all periods, 4 workers)
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from lib import obs2tile, scf_qc  # noqa: E402
from periods import DATA_DIR, PERIODS, VIIRS, cmg_file, days, products  # noqa: E402

TILECOORD = {
    'M36': Path('/Users/amfox/Desktop/LS_DAv8_M36_as_test3/output/SMAP_EASEv2_M36_GLOBAL/rc_out/'
                'LS_DAv8_M36_as_test3.ldas_tilecoord.bin'),
    # tile space of the M21C_test_CF0360_jpss_5d4d636_* runs (copied from Discover, run A2); the land
    # tiles are used (obs2tile.land_tiles). Not the older local M21C_test_CF0360 file: same tiles, other boxes.
    'CF0360': DATA_DIR.parent / 'validation_cf0360' / 'tiles' / 'M21C_test_CF0360_jpss_5d4d636_A2' / 'output' /
              'CF0360x6C_GLOBAL' / 'rc_out' / 'M21C_test_CF0360_jpss_5d4d636_A2.ldas_tilecoord.bin',
}
CACHE = DATA_DIR.parent / 'tilemap'
OUT = DATA_DIR.parent / 'out' / 'coverage'

# latitude bands (plan): label, lower, upper (deg); a cell/tile at lat goes to lower <= lat < upper
BANDS = [('N60-90', 60, 90), ('N45-60', 45, 60), ('N30-45', 30, 45), ('N00-30', 0, 30),
         ('S00-30', -30, 0), ('S30-45', -45, -30), ('S45-60', -60, -45), ('S60-90', -90, -60)]
SCF_CLASSES = ('0', '0-1', '1')      # SCF == 0, 0 < SCF < 1, SCF == 1 (after the cap)
MAP_RES = 0.5                        # deg, for the class-1 maps


def band_index(lat):
    idx = np.full(np.shape(lat), -1, dtype=np.int8)
    for k, (_, lo, hi) in enumerate(BANDS):
        idx[(lat >= lo) & (lat < hi)] = k
    idx[np.asarray(lat) >= 90] = 0
    return idx


def scf_class(scf):
    c = np.full(scf.shape, 1, dtype=np.int8)
    c[scf == 0] = 0
    c[scf >= 1] = 2
    return c


_STATE = {}


def _init():
    lat = scf_qc.cmg_lon_lat()[1]
    _STATE['cell_band'] = band_index(lat)                      # (3600,)
    for kind, f in TILECOORD.items():
        m, tc = obs2tile.cmg_tile_map(f, kind, CACHE)
        _STATE[kind] = (m, tc['N_tile'], band_index(tc['com_lat']))


def process(args):
    period, product, day = args
    path = cmg_file(product, day)
    if not path.exists():
        return None
    raw = scf_qc.read_viirs_cmg(path)
    keep, qa = scf_qc.viirs_qc_without_qa(raw)
    sc = raw['Snow_Cover'].astype(np.float32)
    ci = raw['Clear_Index'].astype(np.float32)
    scf = np.full(sc.shape, np.nan, np.float32)
    scf[keep] = np.minimum(sc[keep] / ci[keep], np.float32(scf_qc.QC_VIIRS['scf_cap']))
    tag = dict(period=period, product=product, day=day.isoformat())

    # cells: band x scf class x qa class
    rows = np.broadcast_to(_STATE['cell_band'][:, None], keep.shape)
    b, c, q = rows[keep], scf_class(scf[keep]), qa[keep]
    ok = (b >= 0) & (q >= 0) & (q <= 3)
    counts = np.zeros((len(BANDS), 3, 4), np.int64)
    np.add.at(counts, (b[ok], c[ok], q[ok]), 1)
    cells = [dict(tag, band=BANDS[i][0], scf_class=SCF_CLASSES[j], qa=k, n=int(counts[i, j, k]))
             for i in range(len(BANDS)) for j in range(3) for k in range(4)]

    # tiles
    k03 = keep & ((qa == 0) | (qa == 3))
    k1 = keep & (qa == 1)
    k013 = k03 | k1
    tiles = []
    for kind in TILECOORD:
        m, nt, tband = _STATE[kind]
        n03, _ = obs2tile.tile_mean(m, scf, k03, nt)
        n1, _ = obs2tile.tile_mean(m, scf, k1, nt)
        n013, mean013 = obs2tile.tile_mean(m, scf, k013, nt)
        snowy = (n013 > 0) & (mean013 > 0)
        only1 = (n03 == 0) & (n1 > 0)
        mixed = (n03 > 0) & (n1 > 0)
        for i, (bl, _, _) in enumerate(BANDS):
            inb = tband == i
            tiles.append(dict(tag, grid=kind, band=bl,
                              obs_03=int((inb & (n03 > 0)).sum()),
                              obs_013=int((inb & (n013 > 0)).sum()),
                              only_class1=int((inb & only1).sum()),
                              snowy_013=int((inb & snowy).sum()),
                              snowy_only_class1=int((inb & snowy & only1).sum()),
                              mixed=int((inb & mixed).sum()),
                              snowy_mixed=int((inb & snowy & mixed).sum())))

    # 0.5 deg maps
    f = int(round(MAP_RES / scf_qc.CMG_DLAT))
    def coarse(mask):
        return mask.reshape(scf_qc.CMG_N_LAT // f, f, scf_qc.CMG_N_LON // f, f).sum(axis=(1, 3)).astype(np.int32)
    snow = keep & (scf > 0)
    maps = dict(kept=coarse(k013), class1=coarse(k1), snowy_kept=coarse(k013 & snow), snowy_class1=coarse(k1 & snow))
    return cells, tiles, (period, product, maps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('periods', nargs='*', default=list(PERIODS))
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()

    jobs = [(p, prod, d) for p in a.periods for prod in products(p) if prod in VIIRS for d in days(p)]
    _init()      # build/load the tile maps once before forking
    OUT.mkdir(parents=True, exist_ok=True)
    all_cells, all_tiles, maps = [], [], {}
    missing = []
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        for job, res in zip(jobs, ex.map(process, jobs)):
            if res is None:
                missing.append(job)
                continue
            cells, tiles, (p, prod, m) = res
            all_cells += cells
            all_tiles += tiles
            acc = maps.setdefault((p, prod), {k: np.zeros_like(v) for k, v in m.items()} | {'n_days': 0})
            for k, v in m.items():
                acc[k] += v
            acc['n_days'] += 1
            print(f'{job[0]} {job[1]} {job[2]}  kept {sum(c["n"] for c in cells if c["qa"] != 2):9d}', flush=True)

    for name, rows in (('cells.csv', all_cells), ('tiles.csv', all_tiles)):
        with open(OUT / name, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    for (p, prod), m in maps.items():
        np.savez_compressed(OUT / f'maps_{p}_{prod}.npz', res_deg=MAP_RES, **m)
    print(f'wrote {OUT}; missing files: {[(j[0], j[1], j[2].isoformat()) for j in missing]}')


if __name__ == '__main__':
    main()
