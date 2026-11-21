#!/usr/bin/env python3
"""
Part B: JPSS-1 (VJ110C1) vs JPSS-2 (VJ210C1), obs only (plan: docs/viirs_obs_study_plan.md), periods P3-P7.

QA class 1 is dropped (amfox, 2026-11-16): VIIRS cells are kept with Basic_QA in {0, 3}
(scf_qc.viirs_qc_without_qa + qa in (0, 3)); MODIS Aqua with scf_qc.modis_scf.

B1  cell level, per period x latitude band x J1 QA class x subset (agreement.compare, as A2):
      j1_j2     : J1 (x) vs J2 (y) on cells where both are kept
      j1_aqua3, j2_aqua3 : J1 vs Aqua and J2 vs Aqua on the cells where J1, J2 and Aqua are all kept
                  (same cells, so the two are directly comparable)
B2  tile level, per tile space (M36, CF0360; viirs_snow_lib/obs2tile.py) x tile latitude band, per day: tiles with an
      obs from J1 only, J2 only, both (a tile has an obs if >= 1 kept cell maps to it); and per period the
      number of days (0..7) a tile has an obs, for J1 alone and J1 or J2, for all tiles and for tiles that
      are snowy on >= 1 day. "snowy" = tile obs > 0, where the tile obs of J1+J2 is what GEOSldas gives the
      update: the average of the species' tile means ("average in case there are multiple asnow obs").
B3  update side, tile-days where J1 and J2 both have an obs: how often the J1+J2 tile obs is on the other side
      of a threshold than the J1 tile obs, and how often J1 and J2 are on different sides (thresholds: SCF > 0
      adds snow; SCF >= SCF_ANA_BETA = 0.55 keeps it, below it the update removes snow), plus mean |J1+J2 - J1|.

Output: <data>/viirs_obs_study/out/jpss/{b1_stats,b2_tiledays,b2_daydist,b3_update}.csv (git-ignored)

usage: jpss_compare.py [--workers N] [P3 P4 ...]     (default: P3..P7, one worker per period)
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
from viirs_snow_lib import obs2tile, scf_qc  # noqa: E402
from periods import DATA_DIR, PERIODS, cmg_file, days  # noqa: E402
from coverage import BANDS, CACHE, TILECOORD, band_index  # noqa: E402
from agreement import SUBSETS, compare  # noqa: E402

OUT = DATA_DIR.parent / 'out' / 'jpss'
J1, J2, AQUA = 'VJ110C1', 'VJ210C1', 'MYD10C1'
PERIODS_B = ['P3', 'P4', 'P5', 'P6', 'P7']
KEEP_QA = (0, 3)
BETA = np.float32(0.55)
STAT_COLS = ['n', 'sx', 'sy', 'sxx', 'syy', 'sxy', 'p_both', 'p_v', 'p_m', 'p_none', 't_both', 't_v', 't_m', 't_none']


def load_viirs(prod, day):
    raw = scf_qc.read_viirs_cmg(cmg_file(prod, day))
    keep, qa = scf_qc.viirs_qc_without_qa(raw)
    keep = keep & np.isin(qa, KEEP_QA)
    sc, ci = raw['Snow_Cover'].astype(np.float32), raw['Clear_Index'].astype(np.float32)
    scf = np.full(sc.shape, np.nan, np.float32)
    scf[keep] = np.minimum(sc[keep] / ci[keep], np.float32(scf_qc.QC_VIIRS['scf_cap']))
    return scf, keep, qa.astype(np.int64)


def run_period(period):
    lat = scf_qc.cmg_lon_lat()[1]
    cell_band = band_index(lat)
    nb = len(BANDS)
    maps = {}
    for kind, f in TILECOORD.items():
        m, tc = obs2tile.cmg_tile_map(f, kind, CACHE)
        maps[kind] = (m, tc['N_tile'], band_index(tc['com_lat']))

    b1 = {}
    b2_days = []                      # per day rows
    b3 = []
    ndays = {k: {'j1': np.zeros(maps[k][1], np.int8), 'any': np.zeros(maps[k][1], np.int8),
                 'snowy': np.zeros(maps[k][1], bool)} for k in maps}

    for day in days(period):
        if not (cmg_file(J1, day).exists() and cmg_file(J2, day).exists() and cmg_file(AQUA, day).exists()):
            print(period, day, 'missing file, skipped', flush=True)
            continue
        s1, k1, q1 = load_viirs(J1, day)
        s2, k2, q2 = load_viirs(J2, day)
        sa, ka = scf_qc.modis_scf(scf_qc.read_modis_cmg(cmg_file(AQUA, day)))
        band_cells = np.broadcast_to(cell_band[:, None], k1.shape)

        # B1
        def comp(name, m, x, y, q):
            rows = band_cells[m]
            ok = rows >= 0
            a = compare(x[m][ok], y[m][ok], q[m][ok], rows[ok].astype(np.int64))
            b1[name] = b1.get(name, 0) + a
        comp('j1_j2', k1 & k2, s1, s2, q1)
        m3 = k1 & k2 & ka
        comp('j1_aqua3', m3, s1, sa, q1)
        comp('j2_aqua3', m3, s2, sa, q2)

        # B2, B3
        for kind, (tmap, nt, tband) in maps.items():
            n1, mean1 = obs2tile.tile_mean(tmap, s1, k1, nt)
            n2, mean2 = obs2tile.tile_mean(tmap, s2, k2, nt)
            h1, h2 = n1 > 0, n2 > 0
            both = h1 & h2
            comb = np.where(both, 0.5 * (mean1 + mean2), np.where(h1, mean1, mean2))
            hany = h1 | h2
            snowy = hany & (comb > 0)
            ndays[kind]['j1'] += h1
            ndays[kind]['any'] += hany
            ndays[kind]['snowy'] |= snowy
            snow1 = h1 & (mean1 > 0)
            snow2 = h2 & (mean2 > 0)
            sn_any = (snow1 | snow2)
            for i, (bl, _, _) in enumerate(BANDS):
                inb = tband == i
                b2_days.append(dict(period=period, day=day.isoformat(), grid=kind, band=bl,
                                    j1_only=int((inb & h1 & ~h2).sum()), j2_only=int((inb & h2 & ~h1).sum()),
                                    both=int((inb & both).sum()),
                                    snowy_j1_only=int((inb & h1 & ~h2 & (mean1 > 0)).sum()),
                                    snowy_j2_only=int((inb & h2 & ~h1 & (mean2 > 0)).sum()),
                                    snowy_both=int((inb & both & (comb > 0)).sum())))
                for subset, sel in (('both', inb & both), ('both_snowy', inb & both & sn_any)):
                    n = int(sel.sum())
                    if n == 0:
                        b3.append(dict(period=period, day=day.isoformat(), grid=kind, band=bl, subset=subset, n=0,
                                       change_p=0, change_t=0, differ_p=0, differ_t=0, sum_absdiff=0.0))
                        continue
                    a1, a2, ac = mean1[sel], mean2[sel], comb[sel]
                    b3.append(dict(period=period, day=day.isoformat(), grid=kind, band=bl, subset=subset, n=n,
                                   change_p=int(((a1 > 0) != (ac > 0)).sum()),
                                   change_t=int(((a1 >= BETA) != (ac >= BETA)).sum()),
                                   differ_p=int(((a1 > 0) != (a2 > 0)).sum()),
                                   differ_t=int(((a1 >= BETA) != (a2 >= BETA)).sum()),
                                   sum_absdiff=float(np.abs(ac - a1).sum())))
        print(period, day, 'ok', flush=True)

    dist = []
    for kind, (tmap, nt, tband) in maps.items():
        for i, (bl, _, _) in enumerate(BANDS):
            inb = tband == i
            for variant in ('j1', 'any'):
                for sub, selt in (('all', inb), ('snowy', inb & ndays[kind]['snowy'])):
                    h = np.bincount(ndays[kind][variant][selt], minlength=8)
                    for k in range(8):
                        dist.append(dict(period=period, grid=kind, band=bl, variant=variant, tiles=sub, days=k, n=int(h[k])))
    return period, b1, b2_days, dist, b3


def write(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('periods', nargs='*', default=PERIODS_B)
    ap.add_argument('--workers', type=int, default=5)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    b1rows, b2, dist, b3 = [], [], [], []
    with ProcessPoolExecutor(a.workers) as ex:
        for period, b1, d2, dd, d3 in ex.map(run_period, a.periods):
            for name, arr in b1.items():
                for bi, (bl, _, _) in enumerate(BANDS):
                    for qi in range(arr.shape[1]):
                        for si, sub in enumerate(SUBSETS):
                            v = arr[bi, qi, si]
                            if v[0]:
                                b1rows.append(dict(period=period, comp=name, band=bl, qa_j1=qi, subset=sub,
                                                   **{c: float(x) for c, x in zip(STAT_COLS, v)}))
            b2 += d2
            dist += dd
            b3 += d3
    write(OUT / 'b1_stats.csv', b1rows)
    write(OUT / 'b2_tiledays.csv', b2)
    write(OUT / 'b2_daydist.csv', dist)
    write(OUT / 'b3_update.csv', b3)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
