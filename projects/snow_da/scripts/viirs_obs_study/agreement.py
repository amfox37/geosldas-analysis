#!/usr/bin/env python3
"""
Part A2: agreement of VIIRS SCF with MODIS SCF per VIIRS Basic_QA class (plan: docs/viirs_obs_study_plan.md).

For every day of every period (periods.py), every VIIRS product and both MODIS references
(MYD10C1 Aqua, main; MOD10C1 Terra, second): the 0.05 deg CMG cells where both pass their
reader QC (VIIRS: all reader tests except Basic_QA, scf_qc.viirs_qc_without_qa, with the QA
class kept as a label; MODIS: scf_qc.modis_scf), same day, same cell. Class 2 (dropped by
the reader) is included for reference.

stats.csv   sums per period x product x ref x latitude band x QA class x subset, where
            subset is the MODIS SCF class ('ref0': SCF = 0, 'ref0-1': 0 < SCF < 1,
            'ref1': SCF = 1; together all matched cells) or 'snowy' (VIIRS or MODIS
            SCF > 0; overlaps the others). Columns:
              n, sx, sy, sxx, syy, sxy     x = VIIRS SCF, y = MODIS SCF (for bias, RMSD, r)
              p_*  snow / no snow at SCF > 0 (with no model snow, any obs SCF > 0 adds snow):
                   p_both, p_v (VIIRS only), p_m (MODIS only), p_none
              t_*  the same at SCF >= SCF_ANA_BETA = 0.55 (below it the update removes snow)
hist_<ref>.npz
            2-D histograms of (VIIRS, MODIS) SCF, 0.05 bins, per period x product x QA class
maps_<period>.npz
            0.5 deg maps, JPSS-1 vs MODIS Aqua, summed over the period's days: matched cells
            and cells that disagree at SCF > 0 and at 0.55, for QA classes 0, 1 and 3 (for A3)

Bias = mean(VIIRS - MODIS). The thresholds are SCF_ANA_ALPHA = 0.60 and SCF_ANA_BETA = 0.55
(clsm_ensupd_glob_param.F90): snow is added if the forecast asnow < ALPHA * obs SCF and
removed if obs SCF < BETA.

Output: <data>/viirs_obs_study/out/agreement/ (git-ignored).

usage: agreement.py [--workers N] [P1 P2 ...]     (default: all periods, 4 workers)
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
from lib import scf_qc  # noqa: E402
from periods import DATA_DIR, MODIS, PERIODS, VIIRS, cmg_file, days, products  # noqa: E402
from coverage import BANDS, band_index  # noqa: E402

OUT = DATA_DIR.parent / 'out' / 'agreement'
SCF_ANA_BETA = np.float32(0.55)
SUBSETS = ('ref0', 'ref0-1', 'ref1', 'snowy')
SUMS = ('n', 'sx', 'sy', 'sxx', 'syy', 'sxy')
CATS = ('both', 'v', 'm', 'none')
NQA = 4
HBINS = np.linspace(0, 1, 21)
MAP_RES = 0.5
MAP_QA = (0, 1, 3)


def _sums(idx, x, y, size):
    """Per index: n, sum x, sum y, sum x^2, sum y^2, sum xy (float64)."""
    x, y = x.astype(np.float64), y.astype(np.float64)
    return np.stack([np.bincount(idx, w, minlength=size) for w in
                     (np.ones_like(x), x, y, x * x, y * y, x * y)], axis=-1)


def _cats(idx, a, b, size):
    """Per index: 2x2 counts (both, a only, b only, none) of two boolean arrays."""
    k = 2 * (~a).astype(np.int64) + (~b).astype(np.int64)      # 0 both, 1 a only, 2 b only, 3 none
    return np.bincount(idx * 4 + k, minlength=size * 4).reshape(size, 4)


def compare(x, y, qa, band):
    """Stats array (band, qa, subset, 6 sums + 8 counts) for matched cells x (VIIRS), y (MODIS)."""
    nb = len(BANDS)
    out = np.zeros((nb, NQA, len(SUBSETS), len(SUMS) + 8))
    refcls = np.where(y == 0, 0, np.where(y >= 1, 2, 1))
    snowy = (x > 0) | (y > 0)
    for sub, sel in ((None, np.ones(x.shape, bool)), (3, snowy)):
        b, q, xs, ys = band[sel], qa[sel], x[sel], y[sel]
        if sub is None:
            idx = (b * NQA + q) * len(SUBSETS) + refcls[sel]
        else:
            idx = (b * NQA + q) * len(SUBSETS) + sub
        size = nb * NQA * len(SUBSETS)
        s = _sums(idx, xs, ys, size)
        c0 = _cats(idx, xs > 0, ys > 0, size)
        c1 = _cats(idx, xs >= SCF_ANA_BETA, ys >= SCF_ANA_BETA, size)
        out.reshape(size, -1)[:] += np.concatenate([s, c0, c1], axis=1)
    return out


def process(args):
    period, day = args
    lat = scf_qc.cmg_lon_lat()[1]
    band_row = band_index(lat)
    refs = {}
    for ref in MODIS:
        f = cmg_file(ref, day)
        if f.exists():
            refs[ref] = scf_qc.modis_scf(scf_qc.read_modis_cmg(f))
    res = {}
    hists = {}
    maps = {}
    f_map = int(round(MAP_RES / scf_qc.CMG_DLAT))

    def coarse(mask):
        return mask.reshape(scf_qc.CMG_N_LAT // f_map, f_map, scf_qc.CMG_N_LON // f_map, f_map).sum(axis=(1, 3)).astype(np.int32)

    for prod in products(period):
        if prod not in VIIRS:
            continue
        f = cmg_file(prod, day)
        if not f.exists():
            continue
        raw = scf_qc.read_viirs_cmg(f)
        keep, qa = scf_qc.viirs_qc_without_qa(raw)
        keep &= (qa >= 0) & (qa < NQA)
        sc = raw['Snow_Cover'].astype(np.float32)
        ci = raw['Clear_Index'].astype(np.float32)
        scf = np.full(sc.shape, np.nan, np.float32)
        scf[keep] = np.minimum(sc[keep] / ci[keep], np.float32(scf_qc.QC_VIIRS['scf_cap']))
        for ref, (mscf, mkeep) in refs.items():
            m = keep & mkeep
            rows = np.broadcast_to(band_row[:, None], m.shape)[m]
            ok = rows >= 0
            x, y, q = scf[m][ok], mscf[m][ok], qa[m][ok].astype(np.int64)
            res[(prod, ref)] = compare(x, y, q, rows[ok].astype(np.int64))
            for k in range(NQA):
                s = q == k
                hists[(prod, ref, k)] = np.histogram2d(x[s], y[s], bins=[HBINS, HBINS])[0].astype(np.int64)
            if prod == 'VJ110C1' and ref == 'MYD10C1':
                v0, m0 = scf > 0, mscf > 0
                v5, m5 = scf >= SCF_ANA_BETA, mscf >= SCF_ANA_BETA
                for k in MAP_QA:
                    mk = m & (qa == k)
                    maps[f'matched_qa{k}'] = coarse(mk)
                    maps[f'disagree0_qa{k}'] = coarse(mk & (v0 != m0))
                    maps[f'disagree55_qa{k}'] = coarse(mk & (v5 != m5))
    return period, day, res, hists, maps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('periods', nargs='*', default=list(PERIODS))
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, d) for p in a.periods for d in days(p)]
    acc, hacc, macc, ndays = {}, {}, {}, {}
    with ProcessPoolExecutor(a.workers) as ex:
        for period, day, res, hists, maps in ex.map(process, jobs):
            for (prod, ref), arr in res.items():
                key = (period, prod, ref)
                acc[key] = acc.get(key, 0) + arr
                ndays[key] = ndays.get(key, 0) + 1
            for (prod, ref, k), h in hists.items():
                hk = (ref, period, prod, k)
                hacc[hk] = hacc.get(hk, 0) + h
            if maps:
                m = macc.setdefault(period, {k: np.zeros_like(v) for k, v in maps.items()} | {'n_days': 0})
                for k, v in maps.items():
                    m[k] += v
                m['n_days'] += 1
            print(f'{period} {day}  pairs {len(res)}', flush=True)

    cols = list(SUMS) + [f'p_{c}' for c in CATS] + [f't_{c}' for c in CATS]
    with open(OUT / 'stats.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['period', 'product', 'ref', 'days', 'band', 'qa', 'subset'] + cols)
        for (period, prod, ref), arr in sorted(acc.items()):
            for i, (bl, _, _) in enumerate(BANDS):
                for k in range(NQA):
                    for j, sub in enumerate(SUBSETS):
                        v = arr[i, k, j]
                        w.writerow([period, prod, ref, ndays[(period, prod, ref)], bl, k, sub] +
                                   [int(round(t)) if c in ('n',) or c[:2] in ('p_', 't_') else repr(float(t))
                                    for c, t in zip(cols, v)])
    for ref in MODIS:
        hs = {f'{p}_{prod}_qa{k}': h for (r, p, prod, k), h in hacc.items() if r == ref}
        if hs:
            np.savez_compressed(OUT / f'hist_{ref}.npz', edges=HBINS, **hs)
    for period, m in macc.items():
        np.savez_compressed(OUT / f'maps_{period}.npz', res_deg=MAP_RES, **m)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
