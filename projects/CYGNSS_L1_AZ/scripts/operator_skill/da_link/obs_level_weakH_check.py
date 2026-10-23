#!/usr/bin/env python3
"""
obs_level (task D), step 3: is the monitor harm after weak-H L1 outliers caused by those obs, or by
the tile/season they occur in?  And how much of the pooled monitor dMSE do they and the
largest-increment L1 obs account for?

  (1) context-matched comparison: restrict to (tile, calendar month) cells that contain >= 1
      monitor obs attributed to a weak-H dominant L1 obs; compare monitor noise/gain for
      weak-H-dominant vs normal-dominant monitor obs inside those same cells.
  (2) contribution: dMSE contribution of each attribution class to the POOLED monitor dMSE
      (in % of the pooled OL MSE), and the pooled dMSE with that class's dF set to 0.
Usage:
  obs_level_weakH_check.py --tag 202001_202212 --inputs 202001_202012,202101_202112,202201_202212
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obs_level_analyze as A  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tag', required=True)
    p.add_argument('--inputs', required=True)
    p.add_argument('--data-dir', default=A.OUT0 + 'data/')
    p.add_argument('--boot', type=int, default=500)
    a = p.parse_args()
    out = f'{A.OUT0}{a.tag}/'
    rng = np.random.default_rng(2)
    L, M = A.load(a)
    d = A.l1_matched(L)
    d, dom, rec, cnt, nweak = A.attribute(M, d)
    M['has_L1'] = dom >= 0
    M['dom_weakH'] = pd.Series(d.weakH.values[np.maximum(dom, 0)]).where(dom >= 0).values
    q5 = d.inc.abs().quantile(0.8)
    M['dom_bigInc'] = pd.Series(d.inc.abs().values[np.maximum(dom, 0)] > q5).where(dom >= 0).values
    M['ym'] = M.year * 100 + M.month
    rows = []
    for grp, g in M[M.grp2.isin(['SMAP', 'L3', 'SMOS'])].groupby('grp2'):
        cells = g.loc[g.dom_weakH == 'weak-H outlier', ['tile', 'ym']].drop_duplicates()
        gc = g.merge(cells, on=['tile', 'ym'])
        for lab, sub in [('weak-H cells: weak-H dominant', gc[gc.dom_weakH == 'weak-H outlier']),
                         ('weak-H cells: normal dominant', gc[gc.dom_weakH == 'normal']),
                         ('weak-H cells: no L1 in 72h', gc[~gc.has_L1])]:
            rows.append(dict(grp=grp, test='context-matched', cls=lab, **A.boot_stats(sub, A.MON_SUMS, a.boot, rng)))
        tot_b2 = g.b2.sum()
        tot = 100 * (g.dd.sum() - 2 * g.di.sum()) / tot_b2
        for lab, msk in [('weak-H dominant', g.dom_weakH == 'weak-H outlier'),
                         (f'largest-|A-F| quintile dominant (|A-F|>{q5:.2f} dB)', g.dom_bigInc == True),  # noqa: E712
                         ('no L1 in 72h', ~g.has_L1)]:
            s = g[msk]
            contrib = 100 * (s.dd.sum() - 2 * s.di.sum()) / tot_b2
            rows.append(dict(grp=grp, test='contribution', cls=lab, N=int(msk.sum()), pooled_dmse=tot,
                             contrib_dmse=contrib, pooled_dmse_without_class_dF=tot - contrib))
    r = pd.DataFrame(rows)
    r.to_csv(f'{out}q3_weakH_context_and_contribution.csv', index=False)
    pd.set_option('display.width', 250)
    print(r.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))


if __name__ == '__main__':
    main()
