#!/usr/bin/env python3
"""
Part A1 tables and figures from coverage.py output (for docs/viirs_obs_study_report.md).

Prints markdown tables to stdout and writes figures to
docs/viirs_obs_study_report_figures/:
  figA1_1_class1_only_tiles.png   snowy tile-days with obs only through class 1, per period
                                  and product, CF0360 and M36, with the 2% line of the
                                  decision rule (dot = period, bar = daily range)
  figA1_2_class1_maps.png         where class 1 occurs: class-1 cell-days per 0.5 deg cell,
                                  JPSS-1 (VJ110C1), one panel per period

usage: coverage_report.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from coverage import OUT  # noqa: E402
from periods import PERIODS  # noqa: E402

FIG = HERE.parents[1] / 'docs' / 'viirs_obs_study_report_figures'

# reference palette (dataviz skill references/palette.md), as in the reader-test report figures
SURFACE, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
CAT = {'VJ110C1': '#2a78d6', 'VJ210C1': '#eb6834', 'VNP10C1': '#1baf7a'}
NAME = {'VJ110C1': 'JPSS-1', 'VJ210C1': 'JPSS-2', 'VNP10C1': 'S-NPP'}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'text.color': INK, 'axes.labelcolor': INK2,
                     'xtick.color': MUTED, 'ytick.color': INK2, 'axes.edgecolor': AXIS,
                     'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE})


def md(df, fmt='{:.2f}'):
    cols = [str(c) for c in df.columns]
    lines = ['| ' + ' | '.join([df.index.name or ''] + cols) + ' |', '|' + '---|' * (len(cols) + 1)]
    for i, r in df.iterrows():
        lines.append('| ' + ' | '.join([str(i)] + ['' if pd.isna(v) else fmt.format(v) for v in r]) + ' |')
    return '\n'.join(lines)


def tables(c, t):
    k = c[c.qa != 2]
    g = k.groupby(['period', 'product', 'qa']).n.sum().unstack()
    sn = k[k.scf_class != '0'].groupby(['period', 'product', 'qa']).n.sum().unstack()
    ndays = k.groupby(['period', 'product']).day.nunique()
    cells = pd.DataFrame({'days': ndays, 'kept cells / day (M)': g.sum(1) / ndays / 1e6,
                          'class 3 % of kept': 100 * g[3] / g.sum(1),
                          'class 1 % of kept': 100 * g[1] / g.sum(1),
                          'class 1 % of snowy': 100 * sn[1] / sn.sum(1)})
    cells.index = [f'{p} {NAME[q]}' for p, q in cells.index]
    cells.index.name = 'period, satellite'
    print('### Cells\n')
    print(md(cells))

    gb = k[k['product'] == 'VJ110C1'].groupby(['band', 'period', 'qa']).n.sum().unstack()
    snb = k[(k['product'] == 'VJ110C1') & (k.scf_class != '0')].groupby(['band', 'period', 'qa']).n.sum().unstack()
    for title, x in (('class 1 % of kept cells', 100 * gb[1] / gb.sum(1)),
                     ('class 1 % of snowy cells', 100 * snb[1] / snb.sum(1)),
                     ('snowy kept cells per day (thousands)', snb.sum(1) / 7e3)):
        x = x.unstack('period')
        x.index.name = 'band'
        print(f'\n### JPSS-1, {title}, by latitude band\n')
        print(md(x, '{:.0f}' if 'thousands' in title else '{:.1f}'))

    s = t.groupby(['grid', 'period', 'product'])[['obs_013', 'only_class1', 'snowy_013',
                                                    'snowy_only_class1', 'snowy_mixed']].sum()
    d = t.groupby(['grid', 'period', 'product', 'day'])[['snowy_013', 'snowy_only_class1']].sum()
    dmax = (100 * d.snowy_only_class1 / d.snowy_013).groupby(['grid', 'period', 'product']).max()
    tiles = pd.DataFrame({'tile-days with obs only via class 1 (% of tile-days with obs)': 100 * s.only_class1 / s.obs_013,
                          'snowy: only via class 1 (%)': 100 * s.snowy_only_class1 / s.snowy_013,
                          'snowy: max day (%)': dmax,
                          'snowy: class 1 mixed with 0/3 (%)': 100 * s.snowy_mixed / s.snowy_013})
    for grid in ('CF0360', 'M36'):
        x = tiles.loc[grid].copy()
        x.index = [f'{p} {NAME[q]}' for p, q in x.index]
        x.index.name = 'period, satellite'
        print(f'\n### Tiles, {grid}\n')
        print(md(x))
    return tiles


def fig_tiles(tiles, t):
    d = t.groupby(['grid', 'period', 'product', 'day'])[['snowy_013', 'snowy_only_class1']].sum()
    p = 100 * d.snowy_only_class1 / d.snowy_013
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    periods = list(PERIODS)
    off = {'VJ110C1': -0.2, 'VJ210C1': 0., 'VNP10C1': 0.2}
    for ax, grid in zip(axs, ('CF0360', 'M36')):
        ax.axhline(2, color=INK2, lw=1, ls=(0, (4, 3)))
        ax.text(len(periods) - 0.5, 2.05, 'decision-rule line (2%)', ha='right', va='bottom', color=INK2, fontsize=8)
        for prod in CAT:
            for i, per in enumerate(periods):
                key = (grid, per, prod)
                if key not in tiles.index:
                    continue
                daily = p.loc[key]
                x = i + off[prod]
                ax.plot([x, x], [daily.min(), daily.max()], color=CAT[prod], lw=2, solid_capstyle='round', alpha=0.45)
                ax.plot(x, tiles.loc[key, 'snowy: only via class 1 (%)'], 'o', ms=6, color=CAT[prod],
                        mec=SURFACE, mew=1.5, label=NAME[prod] if i == 2 else None)
        ax.set_xticks(range(len(periods)), periods)
        ax.set_title(f'{grid} tiles', loc='left', fontsize=10, color=INK)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.set_ylim(0, 2.8)
    axs[0].set_ylabel('snowy tile-days with obs\nonly via class 1 (%)')
    axs[1].legend(frameon=False, loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / 'figA1_1_class1_only_tiles.png', dpi=150)
    plt.close(fig)


def fig_maps():
    import cartopy.crs as ccrs
    periods = list(PERIODS)
    bounds = [1, 3, 10, 30, 60, 100, 300, 700]
    cmap = ListedColormap(SEQ)
    norm = BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(figsize=(9, 9.2))
    for i, per in enumerate(periods):
        m = np.load(OUT / f'maps_{per}_VJ110C1.npz')
        c1 = np.ma.masked_less(m['class1'].astype(float), 1)
        ax = fig.add_subplot(4, 2, i + 1, projection=ccrs.Robinson())
        ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
        ax.coastlines(lw=0.3, color=MUTED)
        res = float(m['res_deg'])
        lon = np.arange(-180, 180 + res / 2, res)
        lat = np.arange(90, -90 - res / 2, -res)
        im = ax.pcolormesh(lon, lat, c1, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), rasterized=True)
        share = 100 * m['class1'].sum() / m['kept'].sum()
        first, last, _ = PERIODS[per]
        ax.set_title(f'{per}  {first} .. {last[5:]}   class 1 = {share:.2f}% of kept cells',
                     loc='left', fontsize=8.5, color=INK)
    cax = fig.add_axes([0.58, 0.1, 0.33, 0.015])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label('class-1 cell-days per 0.5° cell (7 days, max 700)', color=INK2, fontsize=8)
    cb.outline.set_visible(False)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.04, hspace=0.25, wspace=0.05)
    fig.savefig(FIG / 'figA1_2_class1_maps.png', dpi=150)
    plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(OUT / 'cells.csv')
    t = pd.read_csv(OUT / 'tiles.csv')
    tiles = tables(c, t)
    fig_tiles(tiles, t)
    fig_maps()
    print(f'\nfigures in {FIG}')


if __name__ == '__main__':
    main()
