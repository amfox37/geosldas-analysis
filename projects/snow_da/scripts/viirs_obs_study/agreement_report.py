#!/usr/bin/env python3
"""
Part A2 tables and figures from agreement.py output (for docs/viirs_obs_study_report.md).

Prints markdown tables to stdout and writes figures to
docs/viirs_obs_study_report_figures/:
  figA2_1_rmsd_by_qa.png        RMSD vs MODIS Aqua per period and QA class (0, 1, 3), all
                                matched cells and snowy ones (dot = JPSS-1, bar = range
                                over the satellites of the period)
  figA2_2_scatter_density.png   VIIRS vs MODIS Aqua SCF, 2-D histograms per QA class,
                                JPSS-1, all periods pooled
  figA2_3_class1_disagree_maps.png class-1 cells where VIIRS JPSS-1 and MODIS Aqua fall on
                                different sides of SCF_ANA_BETA = 0.55, per 0.5 deg cell,
                                one panel per period

Needs coverage.py output too (cells.csv: kept cells per QA class, for the matched share).

usage: agreement_report.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from agreement import OUT  # noqa: E402
from coverage import OUT as COV_OUT  # noqa: E402
from coverage_report import AXIS, FIG, GRID, INK, INK2, MUTED, NAME, SEQ, SURFACE, md  # noqa: E402
from periods import PERIODS  # noqa: E402

SUMS = ['n', 'sx', 'sy', 'sxx', 'syy', 'sxy']
CNT = ['p_both', 'p_v', 'p_m', 'p_none', 't_both', 't_v', 't_m', 't_none']
# QA classes: categorical slots 1-3 of the reference palette (fixed order)
QA_COL = {0: '#2a78d6', 1: '#eb6834', 3: '#1baf7a'}
QA_NAME = {0: 'class 0 (good)', 1: 'class 1 (poor)', 2: 'class 2 (bad)', 3: 'class 3 (other, low sun)'}


def stats(g):
    s = g[SUMS + CNT].sum()
    n = s.n
    mx, my = s.sx / n, s.sy / n
    vx, vy = s.sxx / n - mx ** 2, s.syy / n - my ** 2
    return pd.Series({'n': n, 'mean_v': mx, 'mean_m': my, 'bias': mx - my,
                      'rmsd': np.sqrt((s.sxx - 2 * s.sxy + s.syy) / n),
                      'r': (s.sxy / n - mx * my) / np.sqrt(vx * vy),
                      'agree0': (s.p_both + s.p_none) / n, 'v_only0': s.p_v / n, 'm_only0': s.p_m / n,
                      'agree55': (s.t_both + s.t_none) / n, 'v_only55': s.t_v / n, 'm_only55': s.t_m / n})


def per(d, keys):
    with np.errstate(invalid='ignore', divide='ignore'):
        return d.groupby(keys).apply(stats, include_groups=False)


def tables(d, cov):
    allm = d[d.subset != 'snowy']
    snowy = d[d.subset == 'snowy']
    out = {}
    for ref in ('MYD10C1', 'MOD10C1'):
        a = per(allm[allm.ref == ref], ['period', 'product', 'qa'])
        s = per(snowy[snowy.ref == ref], ['period', 'product', 'qa'])
        out[ref] = (a, s)

    # matched share: matched cells / kept cells, per QA class (Aqua)
    kept = cov.groupby(['period', 'product', 'qa']).n.sum()
    a, s = out['MYD10C1']
    t = pd.DataFrame({
        'class 1 matched (k)': a.n.xs(1, level='qa') / 1e3,
        'matched % of kept, class 0': 100 * a.n.xs(0, level='qa') / kept.xs(0, level='qa'),
        'class 1': 100 * a.n.xs(1, level='qa') / kept.xs(1, level='qa'),
        'class 3': 100 * a.n.xs(3, level='qa') / kept.xs(3, level='qa'),
        'RMSD 0': a.rmsd.xs(0, level='qa'), 'RMSD 1': a.rmsd.xs(1, level='qa'), 'RMSD 3': a.rmsd.xs(3, level='qa'),
        '1 / 3': a.rmsd.xs(1, level='qa') / a.rmsd.xs(3, level='qa'),
        '1 / 3 snowy': s.rmsd.xs(1, level='qa') / s.rmsd.xs(3, level='qa'),
        '1 / 0': a.rmsd.xs(1, level='qa') / a.rmsd.xs(0, level='qa'),
        'bias 1': a.bias.xs(1, level='qa'),
    })
    t.index = [f'{p} {NAME[q]}' for p, q in t.index]
    t.index.name = 'period, satellite'
    print('### Vs MODIS Aqua, per period (all matched cells unless "snowy")\n')
    print(md(t, '{:.2f}'))
    for ref in ('MYD10C1', 'MOD10C1'):
        a, s = out[ref]
        r_all = a.rmsd.xs(1, level='qa') / a.rmsd.xs(3, level='qa')
        r_sn = s.rmsd.xs(1, level='qa') / s.rmsd.xs(3, level='qa')
        print(f'\n{ref}: class-1/class-3 RMSD ratio > 1.5 in {(r_all > 1.5).sum()}/{r_all.size} period-satellites '
              f'(all cells; range {r_all.min():.2f}-{r_all.max():.2f}), {(r_sn > 1.5).sum()}/{r_sn.size} (snowy; '
              f'range {r_sn.min():.2f}-{r_sn.max():.2f})')

    a, _ = out['MOD10C1']
    t = pd.DataFrame({'RMSD 0': a.rmsd.xs(0, level='qa'), 'RMSD 1': a.rmsd.xs(1, level='qa'),
                      'RMSD 3': a.rmsd.xs(3, level='qa'),
                      '1 / 3': a.rmsd.xs(1, level='qa') / a.rmsd.xs(3, level='qa'),
                      'bias 1': a.bias.xs(1, level='qa')})
    t.index = [f'{p} {NAME[q]}' for p, q in t.index]
    t.index.name = 'period, satellite'
    print('\n### Vs MODIS Terra, per period (all matched cells)\n')
    print(md(t, '{:.2f}'))

    j = d[(d.ref == 'MYD10C1') & (d['product'] == 'VJ110C1') & d.qa.isin([0, 1, 3])]
    t = per(j, ['subset', 'qa'])[['n', 'mean_v', 'mean_m', 'bias', 'rmsd']]
    t['n'] = t.n / 1e3
    t.index = [f'{sub}, {QA_NAME[q]}' for sub, q in t.index]
    t.index.name = 'MODIS SCF class, VIIRS QA'
    t.columns = ['n (k)', 'mean VIIRS', 'mean MODIS', 'bias', 'RMSD']
    print('\n### JPSS-1 vs MODIS Aqua by MODIS SCF class, all periods\n')
    print(md(t, '{:.3f}'))

    t = per(allm[(allm.ref == 'MYD10C1') & (allm['product'] == 'VJ110C1')], ['period', 'qa'])
    x = pd.concat({'agree SCF>0 / =0 (%)': 100 * t.agree0.unstack(),
                   'VIIRS snow, MODIS none (%)': 100 * t.v_only0.unstack(),
                   'agree at 0.55 (%)': 100 * t.agree55.unstack()}, axis=1)
    x = x.loc[:, (slice(None), [0, 1, 3])]
    x.columns = [f'{a}, class {q}' for a, q in x.columns]
    x.index.name = 'period'
    print('\n### JPSS-1 vs MODIS Aqua, snow / no-snow agreement (all matched cells)\n')
    print(md(x, '{:.1f}'))

    t = per(allm[(allm.ref == 'MYD10C1') & (allm['product'] == 'VJ110C1') & allm.qa.isin([0, 1, 3])], ['band', 'qa'])
    x = pd.concat({'n (k)': t.n.unstack() / 1e3, 'RMSD': t.rmsd.unstack(), 'bias': t.bias.unstack()}, axis=1)
    x.columns = [f'{a}, class {q}' for a, q in x.columns]
    x.index.name = 'band'
    print('\n### JPSS-1 vs MODIS Aqua by latitude band, all periods (all matched cells)\n')
    print(md(x, '{:.2f}'))
    return out


def style(ax):
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def fig_rmsd(out):
    a, s = out['MYD10C1']
    periods = list(PERIODS)
    off = {0: -0.18, 1: 0., 3: 0.18}
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for ax, df, title in ((axs[0], a, 'all matched cells'), (axs[1], s, 'snowy (VIIRS or MODIS SCF > 0)')):
        for q in (0, 1, 3):
            for i, p in enumerate(periods):
                v = df.rmsd.loc[p].xs(q, level='qa')
                x = i + off[q]
                ax.plot([x, x], [v.min(), v.max()], color=QA_COL[q], lw=2, alpha=0.45, solid_capstyle='round')
                ax.plot(x, v.loc['VJ110C1'], 'o', ms=6, color=QA_COL[q], mec=SURFACE, mew=1.5,
                        label=QA_NAME[q] if i == 0 else None)
        ax.set_xticks(range(len(periods)), periods)
        ax.set_title(title, loc='left', fontsize=10, color=INK)
        style(ax)
    axs[0].set_ylabel('RMSD vs MODIS Aqua (SCF)')
    axs[0].set_ylim(0, 0.6)
    axs[0].legend(frameon=False, loc='upper left', fontsize=8, ncol=1)
    fig.tight_layout()
    fig.savefig(FIG / 'figA2_1_rmsd_by_qa.png', dpi=150)
    plt.close(fig)


def fig_density():
    h = np.load(OUT / 'hist_MYD10C1.npz')
    e = h['edges']
    fig, axs = plt.subplots(1, 3, figsize=(9, 3.3), sharey=True)
    cmap = ListedColormap(SEQ)
    for ax, q in zip(axs, (0, 1, 3)):
        tot = sum(h[k] for k in h.files if k.endswith(f'VJ110C1_qa{q}'))
        frac = np.ma.masked_equal(tot / tot.sum(), 0)
        im = ax.pcolormesh(e, e, frac.T, cmap=cmap, norm=LogNorm(1e-5, 1), rasterized=True)
        ax.plot([0, 1], [0, 1], color=INK2, lw=0.8, ls=(0, (4, 3)))
        ax.set_aspect('equal')
        ax.set_title(f'{QA_NAME[q]}\nn = {tot.sum() / 1e6:.2f} M', loc='left', fontsize=9, color=INK)
        ax.set_xlabel('VIIRS JPSS-1 SCF')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    axs[0].set_ylabel('MODIS Aqua SCF')
    cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.02)
    cb.set_label('share of the class\'s matched cells', color=INK2, fontsize=8)
    cb.outline.set_visible(False)
    fig.savefig(FIG / 'figA2_2_scatter_density.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def fig_disagree():
    import cartopy.crs as ccrs
    periods = list(PERIODS)
    bounds = [1, 2, 3, 5, 10, 20, 40]
    cmap = ListedColormap(SEQ)
    norm = BoundaryNorm(bounds, cmap.N, extend='max')
    fig = plt.figure(figsize=(9, 9.2))
    for i, p in enumerate(periods):
        m = np.load(OUT / f'maps_{p}.npz')
        res = float(m['res_deg'])
        c = np.ma.masked_less(m['disagree55_qa1'].astype(float), 1)
        ax = fig.add_subplot(4, 2, i + 1, projection=ccrs.Robinson())
        ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
        ax.coastlines(lw=0.3, color=MUTED)
        lon = np.arange(-180, 180 + res / 2, res)
        lat = np.arange(90, -90 - res / 2, -res)
        im = ax.pcolormesh(lon, lat, c, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), rasterized=True)
        nm = max(m['matched_qa1'].sum(), 1)
        first, last, _ = PERIODS[p]
        ax.set_title(f'{p}  {first} .. {last[5:]}   differ at 0.55: {100 * m["disagree55_qa1"].sum() / nm:.0f}%, '
                     f'at 0: {100 * m["disagree0_qa1"].sum() / nm:.0f}%', loc='left', fontsize=8.5, color=INK)
    cax = fig.add_axes([0.58, 0.1, 0.33, 0.015])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label('class-1 cell-days on the other side of SCF 0.55 than MODIS Aqua\n'
                 '(per 0.5\u00b0 cell, 7 days; % = share of matched class-1 cells)', color=INK2, fontsize=8)
    cb.outline.set_visible(False)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.04, hspace=0.25, wspace=0.05)
    fig.savefig(FIG / 'figA2_3_class1_disagree_maps.png', dpi=150)
    plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(OUT / 'stats.csv')
    cov = pd.read_csv(COV_OUT / 'cells.csv')
    out = tables(d, cov)
    fig_rmsd(out)
    fig_density()
    fig_disagree()
    print(f'\nfigures in {FIG}')


if __name__ == '__main__':
    main()
