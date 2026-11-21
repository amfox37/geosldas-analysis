#!/usr/bin/env python3
"""Part B tables (stdout, markdown) from jpss_compare.py output (<data>/viirs_obs_study/out/jpss/)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from periods import DATA_DIR  # noqa: E402

D = DATA_DIR.parent / 'out' / 'jpss'
NAMES = {'j1_j2': 'J1 vs J2', 'j1_aqua3': 'J1 vs Aqua (J1, J2, Aqua kept)', 'j2_aqua3': 'J2 vs Aqua (J1, J2, Aqua kept)'}


def md(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(f'{x:,.0f}' if isinstance(x, (int, np.integer)) or (isinstance(x, float) and x.is_integer() and abs(x) >= 1000) else str(x) for x in r) + ' |')
    return '\n'.join(out)


def stats(g):
    s = g[['n', 'sx', 'sy', 'sxx', 'syy', 'sxy']].sum()
    n = s.n
    vx, vy = s.sxx / n - (s.sx / n) ** 2, s.syy / n - (s.sy / n) ** 2
    return pd.Series({'cells': n, 'bias': (s.sx - s.sy) / n, 'RMSD': np.sqrt(max((s.sxx + s.syy - 2 * s.sxy) / n, 0)),
                      'r': (s.sxy / n - s.sx * s.sy / n / n) / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan})


def main():
    b1 = pd.read_csv(D / 'b1_stats.csv')
    allsub = b1[b1.subset.isin(['ref0', 'ref0-1', 'ref1'])]
    snowy = b1[b1.subset == 'snowy']
    print('### B1. Agreement, all matched cells and snowy cells (either > 0); periods P3-P7 pooled\n')
    rows = []
    for comp, nm in NAMES.items():
        a = stats(allsub[allsub.comp == comp]); s = stats(snowy[snowy.comp == comp])
        rows.append({'comparison': nm, 'cells': a.cells, 'bias': a.bias, 'RMSD': a.RMSD, 'r': a.r,
                     'snowy cells': s.cells, 'snowy bias': s.bias, 'snowy RMSD': s.RMSD, 'snowy r': s['r']})
    print(md(pd.DataFrame(rows).round(4)), '\n')

    print('### B1. J1 vs J2 snow / no-snow agreement of the cells (SCF > 0 and >= 0.55), pooled\n')
    j = b1[(b1.comp == 'j1_j2') & b1.subset.isin(['ref0', 'ref0-1', 'ref1'])]
    rows = []
    for t, nm in (('p', 'SCF > 0'), ('t', 'SCF >= 0.55')):
        s = j[[f'{t}_both', f'{t}_v', f'{t}_m', f'{t}_none']].sum()
        n = s.sum(); snow = s.iloc[0] + s.iloc[1] + s.iloc[2]
        rows.append({'threshold': nm, 'both snow': s.iloc[0], 'J1 only': s.iloc[1], 'J2 only': s.iloc[2], 'neither': s.iloc[3],
                     'disagree % of cells': 100 * (s.iloc[1] + s.iloc[2]) / n,
                     'disagree % of snow cells': 100 * (s.iloc[1] + s.iloc[2]) / snow})
    print(md(pd.DataFrame(rows).round(3)), '\n')

    print('### B1. By latitude band: J1 vs J2 and each vs Aqua (RMSD; snowy cells in brackets)\n')
    rows = []
    for band in b1.band.unique():
        r = {'band': band}
        for comp, short in (('j1_j2', 'J1-J2'), ('j1_aqua3', 'J1-Aqua'), ('j2_aqua3', 'J2-Aqua')):
            a = stats(allsub[(allsub.comp == comp) & (allsub.band == band)]); s = stats(snowy[(snowy.comp == comp) & (snowy.band == band)])
            r[f'{short} cells (M)'] = round(a.cells / 1e6, 2) if comp == 'j1_j2' else None
            r[f'{short} RMSD'] = f'{a.RMSD:.3f} ({s.RMSD:.3f})'
        r.pop('J1-Aqua cells (M)'); r.pop('J2-Aqua cells (M)')
        rows.append(r)
    print(md(pd.DataFrame(rows)), '\n')

    print('### B1. By period: bias and RMSD J1 vs J2, and each vs Aqua\n')
    rows = []
    for p in sorted(b1.period.unique()):
        r = {'period': p}
        for comp, short in (('j1_j2', 'J1-J2'), ('j1_aqua3', 'J1-Aqua'), ('j2_aqua3', 'J2-Aqua')):
            a = stats(allsub[(allsub.comp == comp) & (allsub.period == p)])
            r[f'{short} bias'] = round(a.bias, 4); r[f'{short} RMSD'] = round(a.RMSD, 4)
        rows.append(r)
    print(md(pd.DataFrame(rows)), '\n')

    t = pd.read_csv(D / 'b2_tiledays.csv')
    print('### B2. Added coverage: tile-days with an obs from J1 only / J2 only / both (P3-P7 pooled)\n')
    rows = []
    for grid, g in t.groupby('grid'):
        s = g[['j1_only', 'j2_only', 'both', 'snowy_j1_only', 'snowy_j2_only', 'snowy_both']].sum()
        rows.append({'grid': grid, 'J1 only': s.j1_only, 'J2 only': s.j2_only, 'both': s.both,
                     'J2 adds, % of J1 tile-days': 100 * s.j2_only / (s.j1_only + s.both),
                     'J1 adds over J2, %': 100 * s.j1_only / (s.j2_only + s.both),
                     'snowy J1 only': s.snowy_j1_only, 'snowy J2 only': s.snowy_j2_only, 'snowy both': s.snowy_both,
                     'snowy: J2 adds, % of J1': 100 * s.snowy_j2_only / (s.snowy_j1_only + s.snowy_both)})
    print(md(pd.DataFrame(rows).round(2)), '\n')
    for grid in ('CF0360', 'M36'):
        print(f'### B2. {grid}: J2 adds, % of J1 tile-days, by period and band (all tile-days / snowy tile-days)\n')
        g = t[t.grid == grid].groupby(['period', 'band'])[['j1_only', 'j2_only', 'both', 'snowy_j1_only', 'snowy_j2_only', 'snowy_both']].sum()
        g['all'] = 100 * g.j2_only / (g.j1_only + g.both)
        g['snowy'] = 100 * g.snowy_j2_only / (g.snowy_j1_only + g.snowy_both)
        tab = g[['all', 'snowy']].round(1)
        tab = tab['all'].astype(str) + ' / ' + tab['snowy'].astype(str)
        print(md(tab.unstack('band').reset_index()), '\n')

    dd = pd.read_csv(D / 'b2_daydist.csv')
    print('### B2. Days with an obs per tile over the 7-day period: J1 alone vs J1 or J2 (mean days; share of tiles with 7 days / with 0 days)\n')
    rows = []
    for (grid, sub), g in dd.groupby(['grid', 'tiles']):
        for variant, vn in (('j1', 'J1'), ('any', 'J1 or J2')):
            h = g[g.variant == variant].groupby('days').n.sum()
            tot = h.sum()
            rows.append({'grid': grid, 'tiles': sub, 'obs': vn, 'tile-periods': tot, 'mean days': (h.index * h).sum() / tot,
                         '% with 7 days': 100 * h.get(7, 0) / tot, '% with 0 days': 100 * h.get(0, 0) / tot})
    print(md(pd.DataFrame(rows).round(2)), '\n')

    u = pd.read_csv(D / 'b3_update.csv')
    print('### B3. Tile-days with J1 and J2 both: update side of J1 alone vs J1+J2 (averaged), and J1 vs J2\n')
    rows = []
    for (grid, sub), g in u.groupby(['grid', 'subset']):
        s = g[['n', 'change_p', 'change_t', 'differ_p', 'differ_t', 'sum_absdiff']].sum()
        rows.append({'grid': grid, 'tile-days': s.n, 'changes side at >0, %': 100 * s.change_p / s.n, 'changes side at 0.55, %': 100 * s.change_t / s.n,
                     'J1, J2 differ at >0, %': 100 * s.differ_p / s.n, 'J1, J2 differ at 0.55, %': 100 * s.differ_t / s.n,
                     'mean abs(J1+J2 - J1)': s.sum_absdiff / s.n, 'subset': sub})
    print(md(pd.DataFrame(rows).round(3)), '\n')
    g = u[u.subset == 'both_snowy'].groupby(['grid', 'period'])[['n', 'change_p', 'change_t']].sum()
    g['change >0 %'] = 100 * g.change_p / g.n; g['change 0.55 %'] = 100 * g.change_t / g.n
    print('### B3. Snowy tile-days with both, by period\n')
    print(md(g.reset_index()[['grid', 'period', 'n', 'change >0 %', 'change 0.55 %']].round(2)), '\n')


if __name__ == '__main__':
    main()
