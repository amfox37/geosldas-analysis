#!/usr/bin/env python3
"""
Task tile_whole_period, sensitivity checks on the key tile associations (reads tile_table.csv written by
tile_whole_period_link.py):
  1. partial Spearman with an extended confounder set: sfmc_mean, n_l1_obs, lat, lon (base)
     + alt_med + slope_med (topography; rs_smterm is ~0.63 rank-correlated with both)
     + OL O-F stdv of the target species group (room-for-improvement / % normalisation)
  2. Spearman within sfmc_mean terciles (aridity strata)
Block bootstrap as in the main script.
Usage: tile_whole_period_sensitivity.py --nboot 500
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tile_whole_period_link import P, ARMS, GROUPS, CONF, Boot, spear, pspear  # noqa: E402

KEY_X = ['rs_smterm', 'r_SMAPh', 'frac_out15', 'R_impl', 'rs_geom', 'r_L3']
KEY_Y = [f'{arm}:{g}' for arm in ['L1', 'L1err39', 'L3'] for g in ['SMAP_Tb', 'SMOS_Tb', 'L3', 'ASCAT', 'L1']] + \
        [f'L1-L3:{g}' for g in ['SMAP_Tb', 'SMOS_Tb', 'ASCAT']]


def ol_std(stats_dir, tag, period, nmin):
    ol = xr.open_dataset(stats_dir / f'temporal_stats_OL_paired_monitor_xmask_{tag}_{period}.nc4')
    o, n = ol.OmF_stdv.values, ol.N_data.values
    o = np.where(n >= nmin, o, np.nan)
    with np.errstate(all='ignore'):
        return {g: np.nanmean(o[:, idx], axis=1) for g, idx in GROUPS.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stats', type=Path, default=P / 'output/postproc_paired_density/stats_output')
    ap.add_argument('--period', default='20200101_20221231')
    ap.add_argument('--nmin', type=int, default=50)
    ap.add_argument('--nboot', type=int, default=500)
    ap.add_argument('--out', type=Path, default=P / 'output/operator_skill/da_link/tile_whole_period')
    a = ap.parse_args()
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
    t = pd.read_csv(a.out / 'tile_table.csv')
    boot = Boot(t.block.values, a.nboot, seed=2)
    rows = []
    for y in KEY_Y:
        arm, g = y.split(':')
        tag = ARMS['L1' if arm.startswith('L1-') else arm]
        t['ol_std'] = ol_std(a.stats, tag, a.period, a.nmin)[g]
        for x in KEY_X:
            for cname, conf in [('base', CONF), ('base+topo', CONF + ['alt_med', 'slope_med']),
                                ('base+topo+OLstd', CONF + ['alt_med', 'slope_med', 'ol_std'])]:
                cols = [x, y] + conf
                m = t.l1area.values & np.isfinite(t[cols].values).all(1)
                xv, yv, Z = t[x].values, t[y].values, t[conf].values
                pr = pspear(xv[m], yv[m], Z[m])
                lo, hi = boot.run(lambda i: [pspear(xv[i], yv[i], Z[i])], m)
                rows.append(dict(target=y, predictor=x, test='partial|' + cname, N=int(m.sum()),
                                 r=pr, lo=lo[0], hi=hi[0]))
            q = pd.qcut(t.sfmc_mean.where(t.l1area), 3, labels=False).values
            for k in range(3):
                m = t.l1area.values & (q == k) & np.isfinite(t[[x, y]].values).all(1)
                xv, yv = t[x].values, t[y].values
                lo, hi = boot.run(lambda i: [spear(xv[i], yv[i])], m)
                rows.append(dict(target=y, predictor=x, test=f'sfmc_tercile{k + 1}', N=int(m.sum()),
                                 r=spear(xv[m], yv[m]), lo=lo[0], hi=hi[0]))
    r = pd.DataFrame(rows)
    r['sig'] = (r.lo > 0) | (r.hi < 0)
    r.to_csv(a.out / 'sensitivity_confounders_strata.csv', index=False)
    r['s'] = r.apply(lambda z: f'{z.r:+.2f}{"*" if z.sig else " "}', axis=1)
    print(r.pivot_table(index=['predictor', 'target'], columns='test', values='s', aggfunc='first').to_string())


if __name__ == '__main__':
    main()
