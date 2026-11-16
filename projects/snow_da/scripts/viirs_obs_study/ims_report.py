#!/usr/bin/env python3
"""
Part A3 tables (stdout, markdown) from ims_check.py output (<data>/viirs_obs_study/out/ims[_shift1]/).

usage: ims_report.py [--dir ims|ims_shift1]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from periods import DATA_DIR  # noqa: E402

NAME = {'VJ110C1': 'JPSS-1', 'VJ210C1': 'JPSS-2', 'VNP10C1': 'S-NPP', 'MYD10C1': 'MODIS Aqua', 'MOD10C1': 'MODIS Terra'}


def md(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(f'{x:.0f}' if isinstance(x, float) and x.is_integer() and abs(x) > 1000 else str(x) for x in r) + ' |')
    return '\n'.join(out)


def score(g):
    n = g[['both', 'obs_only', 'ims_only', 'neither']].sum()
    N = n.sum()
    return pd.Series({'cells': int(N), 'accuracy %': 100 * (n.both + n.neither) / N,
                      'false snow % of obs snow': 100 * n.obs_only / max(n.both + n.obs_only, 1),
                      'missed snow % of IMS snow': 100 * n.ims_only / max(n.both + n.ims_only, 1)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='ims')
    a = ap.parse_args()
    d = DATA_DIR.parent / 'out' / a.dir
    v = pd.read_csv(d / 'vs_ims.csv')
    dis = pd.read_csv(d / 'disagree.csv')
    v['label'] = v.apply(lambda r: NAME[r['product']] if r['qa'] < 0 else f'{NAME[r["product"]]} class {r["qa"]}', axis=1)

    for thr, tname in (('t', 'SCF >= 0.55 (snow removal threshold)'), ('p', 'SCF > 0 (snow addition)')):
        for ims, iname in (('nn', 'nearest IMS pixel'), ('dec', 'decisive 3x3 IMS window')):
            sel = v[(v.thr == thr) & (v.ims_def == ims)]
            t = sel.groupby(['label']).apply(score, include_groups=False).reset_index().round(2)
            print(f'### Obs vs IMS, {tname}; IMS = {iname}; all periods, NH, kept cells\n')
            print(md(t), '\n')

    sel = v[(v.thr == 't') & (v.ims_def == 'nn') & v['product'].isin(['VJ110C1', 'VNP10C1', 'VJ210C1']) & v.qa.isin([0, 1, 3])]
    t = sel.groupby(['band', 'qa']).apply(score, include_groups=False).reset_index().round(2)
    print('### VIIRS classes vs IMS by latitude band (0.55, nearest pixel, all VIIRS products)\n')
    print(md(t), '\n')

    sel = v[(v.thr == 't') & (v.ims_def == 'nn') & v['product'].isin(['VJ110C1', 'VNP10C1', 'VJ210C1']) & v.qa.isin([0, 1, 3])]
    t = sel.groupby(['period', 'qa']).apply(score, include_groups=False).reset_index()
    t = t.pivot(index='period', columns='qa', values=['cells', 'false snow % of obs snow']).round(1)
    t.columns = [f'{a} (class {b})' for a, b in t.columns]
    print('### VIIRS class 1 / 0 / 3 vs IMS by period (0.55, nearest pixel)\n')
    print(md(t.reset_index()), '\n')

    for thr, tname in (('t', '0.55'), ('p', '> 0')):
        for ims in ('nn', 'dec'):
            sel = dis[(dis.thr == thr) & (dis.ims_def == ims) & dis.qa.isin([0, 1, 3])]
            g = sel.groupby(['qa', 'category'])[['n', 'n_ims_snow']].sum()
            g['IMS snow %'] = 100 * g.n_ims_snow / g.n
            g = g.reset_index()
            g['category'] = g.category.map({'both': 'both snow', 'v': 'VIIRS snow only', 'm': 'MODIS Aqua snow only', 'none': 'neither'})
            print(f'### Cells with VIIRS and MODIS Aqua kept: share of cells with IMS snow, threshold {tname}, IMS def {ims}; all VIIRS products and periods\n')
            print(md(g.round({'IMS snow %': 1})), '\n')

    sel = dis[(dis.thr == 't') & (dis.ims_def == 'nn') & (dis.qa == 1) & dis.category.isin(['v', 'm'])]
    g = sel.groupby(['period', 'category'])[['n', 'n_ims_snow']].sum()
    g['IMS snow %'] = 100 * g.n_ims_snow / g.n
    print('### Class 1 disagreements with MODIS Aqua by period (0.55, nearest pixel)\n')
    print(md(g.reset_index().round(1)), '\n')


if __name__ == '__main__':
    main()
