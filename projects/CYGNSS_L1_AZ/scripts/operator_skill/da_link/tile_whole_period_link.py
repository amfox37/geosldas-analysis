#!/usr/bin/env python3
"""
Task tile_whole_period: does the per-tile DA benefit depend on per-tile CYGNSS L1 operator skill?

Per tile and arm: % change of O-F stdv vs the paired OL (same tag, same period), (DA-OL)/OL*100,
N_data>=--nmin per species, group value = nanmean over the group's species (as paper Fig. 7 tile_pct()).
Negative = DA better.  Operator metrics: output/operator_skill/analysis_202001_202212/tile_metrics.csv
(owner tile ig,jg; 2020-22 OL obs).  Primary sample = L1 area (n_l1_obs >= --l1min) with metrics.

Uncertainty: spatial block bootstrap, blocks = --block x --block tiles in (ig, jg); 95% percentile CIs.

Outputs (in --out):
  tile_table.csv                     per-tile metrics + all arm/group % changes (whole period + years)
  spearman_whole_period.csv          Spearman(target, predictor), raw and partial | confounders, block-boot CI
  multivar_rank_regression.csv       standardized rank-regression betas (operator set + confounders)
  landcover_groups.csv               DA benefit by MODIS land-cover class
  bins_<metric>_<q>.csv              tercile/quintile tables (+ top-bottom difference) with block-boot CIs
  yearly_spearman.csv                key associations per calendar year
  maps_*.png, bins_*.png, SUMMARY numbers are printed to stdout (log)
Usage:
  tile_whole_period_link.py --nboot 1000
"""
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

P = Path('/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ')
ARMS = {'L1': 'full_xc015_coh040216_fixedop', 'L1err39': 'full_xc015_coh040216_err39_fixedop',
        'L3': 'l3_fixedop', 'SMAP': 'smap_fixedop'}
ARM_COL = {'L1': '#D55E00', 'L1err39': '#CC79A7', 'L3': '#0072B2', 'SMAP': '#009E73'}
GROUPS = {'SMOS_Tb': [0, 1, 2, 3], 'SMAP_Tb': [4, 5, 6, 7], 'ASCAT': [8, 9, 10], 'L3': [11], 'L1': [12]}
OWN = {('L1', 'L1'), ('L1err39', 'L1'), ('L3', 'L3'), ('SMAP', 'SMAP_Tb')}
OPER = ['rs_H', 'rs_smterm', 'rs_geom', 'pr_geom_inc', 'R_impl', 'inn_std', 'snr', 'dyn_ratio', 'frac_out15',
        'r_SMAPh', 'r_L3']
ATTR = ['sfmc_mean', 'sfmc_std', 'slope_med', 'alt_med', 'opac_med', 'coh_tail_frac', 'n_l1_obs', 'lat', 'lon']
CONF = ['sfmc_mean', 'n_l1_obs', 'lat', 'lon']
MV_SET = ['rs_smterm', 'rs_geom', 'frac_out15', 'R_impl', 'r_L3']
KEY_X = ['rs_smterm', 'frac_out15', 'R_impl', 'rs_geom', 'r_L3', 'sfmc_mean']
YEARS = ['20200101_20201231', '20210101_20211231', '20220101_20221231']


def tile_pct(stats_dir, tag, period, nmin):
    da = xr.open_dataset(stats_dir / f'temporal_stats_DA_paired_{tag}_{period}.nc4')
    ol = xr.open_dataset(stats_dir / f'temporal_stats_OL_paired_monitor_xmask_{tag}_{period}.nc4')
    d, o, n = da.OmF_stdv.values, ol.OmF_stdv.values, da.N_data.values
    out = {}
    with np.errstate(all='ignore'):
        pct = np.where((n >= nmin) & (o > 0), 100.0 * (d - o) / o, np.nan)
        for g, idx in GROUPS.items():
            out[g] = np.nanmean(pct[:, idx], axis=1)
    return pd.DataFrame(out)


def ranks(v):
    return stats.rankdata(v, axis=0)


def spear(x, y):
    return np.corrcoef(ranks(x), ranks(y))[0, 1]


def resid(r, Z):
    A = np.column_stack([np.ones(len(r)), Z])
    return r - A @ np.linalg.lstsq(A, r, rcond=None)[0]


def pspear(x, y, Z):
    rz = ranks(Z)
    return np.corrcoef(resid(ranks(x), rz), resid(ranks(y), rz))[0, 1]


def mv_betas(y, X):
    """standardized OLS betas of rank(y) on rank(X columns)"""
    ry, rX = ranks(y), ranks(X)
    ry = (ry - ry.mean()) / ry.std()
    rX = (rX - rX.mean(0)) / rX.std(0)
    A = np.column_stack([np.ones(len(ry)), rX])
    return np.linalg.lstsq(A, ry, rcond=None)[0][1:]


class Boot:
    """spatial block bootstrap over a tile subset"""
    def __init__(self, blocks, nboot, seed=1):
        self.blocks, self.nboot = np.asarray(blocks), nboot
        self.rng = np.random.default_rng(seed)

    def run(self, fn, mask):
        idx = np.flatnonzero(mask)
        b = self.blocks[idx]
        ub = np.unique(b)
        members = [idx[b == k] for k in ub]
        out = []
        for _ in range(self.nboot):
            pick = self.rng.integers(0, len(ub), len(ub))
            out.append(fn(np.concatenate([members[k] for k in pick])))
        out = np.asarray(out, float)
        return np.nanpercentile(out, 2.5, axis=0), np.nanpercentile(out, 97.5, axis=0)


def build_table(a):
    tiles = pd.read_csv(P / 'output/operator_skill/da_link/tile_axis_ig_jg.csv')
    tm = pd.read_csv(P / 'output/operator_skill/analysis_202001_202212/tile_metrics.csv').drop(columns=['lon', 'lat'])
    t = tiles.merge(tm, on=['ig', 'jg'], how='left').sort_values('tile_axis_index').reset_index(drop=True)
    assert (t.tile_axis_index.values == np.arange(len(t))).all()
    t = t.rename(columns={'com_lon': 'lon', 'com_lat': 'lat'})
    for per in [a.period] + YEARS:
        sfx = '' if per == a.period else '@' + per[:4]
        for arm, tag in ARMS.items():
            p = tile_pct(a.stats, tag, per, a.nmin)
            for g in GROUPS:
                t[f'{arm}:{g}{sfx}'] = p[g].values
        for arm in ['L1', 'L1err39']:
            for g in GROUPS:
                t[f'{arm}-L3:{g}{sfx}'] = t[f'{arm}:{g}{sfx}'] - t[f'L3:{g}{sfx}']
    t['block'] = (t.ig // a.block) * 10000 + (t.jg // a.block)
    t['l1area'] = (t.n_l1_obs >= a.l1min) & t.N.notna()
    return t


def spearman_table(t, boot, ycols, xcols, conf, label=''):
    rows = []
    sub = t.l1area.values
    Z_all = t[conf].values
    for y in ycols:
        for x in xcols:
            m = sub & np.isfinite(t[y].values) & np.isfinite(t[x].values) & np.isfinite(Z_all).all(1)
            if m.sum() < 30:
                continue
            xv, yv = t[x].values, t[y].values
            zc = [c for c in conf if c != x]
            Z = t[zc].values
            r = spear(xv[m], yv[m])
            pr = pspear(xv[m], yv[m], Z[m])
            lo, hi = boot.run(lambda i: [spear(xv[i], yv[i]), pspear(xv[i], yv[i], Z[i])], m)
            rows.append(dict(period=label, target=y, predictor=x, N=int(m.sum()), rs=r, rs_lo=lo[0], rs_hi=hi[0],
                             prs=pr, prs_lo=lo[1], prs_hi=hi[1],
                             rs_sig=(lo[0] > 0) | (hi[0] < 0), prs_sig=(lo[1] > 0) | (hi[1] < 0)))
    return pd.DataFrame(rows)


def bin_table(t, boot, x, q, ycols):
    m0 = t.l1area.values & np.isfinite(t[x].values)
    lab = np.full(len(t), -1)
    lab[m0] = pd.qcut(t.loc[m0, x].rank(method='first'), q, labels=False)
    edges = [t.loc[m0 & (lab == k), x].agg(['min', 'max']).values for k in range(q)]
    rows = []
    for y in ycols:
        m = m0 & np.isfinite(t[y].values)
        yv = t[y].values

        def fn(i):
            li, yi = lab[i], yv[i]
            med = [np.median(yi[li == k]) if (li == k).any() else np.nan for k in range(q)]
            return med + [med[-1] - med[0]]
        lo, hi = boot.run(fn, m)
        full = fn(np.flatnonzero(m))
        for k in range(q + 1):
            sel = m & (lab == k) if k < q else m
            rows.append(dict(target=y, bin=(k + 1 if k < q else f'Q{q}-Q1'),
                             x_min=edges[k][0] if k < q else np.nan, x_max=edges[k][1] if k < q else np.nan,
                             N=int(sel.sum()), median=full[k], lo=lo[k], hi=hi[k],
                             mean=np.mean(yv[sel]) if k < q else np.nan))
    return pd.DataFrame(rows)


def mv_table(t, boot, ycols):
    xs = MV_SET + CONF
    rows = []
    for y in ycols:
        m = t.l1area.values & np.isfinite(t[[y] + xs].values).all(1)
        X, yv = t[xs].values, t[y].values
        b = mv_betas(yv[m], X[m])
        lo, hi = boot.run(lambda i: mv_betas(yv[i], X[i]), m)
        for k, x in enumerate(xs):
            rows.append(dict(target=y, predictor=x, N=int(m.sum()), beta=b[k], lo=lo[k], hi=hi[k],
                             sig=(lo[k] > 0) | (hi[k] < 0)))
    return pd.DataFrame(rows)


def landcover_table(t, boot, ycols):
    rows = []
    for y in ycols:
        for lc in [7, 10, 16]:
            m = t.l1area.values & (t.landcover.values == lc) & np.isfinite(t[y].values)
            yv = t[y].values
            lo, hi = boot.run(lambda i: [np.median(yv[i])], m)
            rows.append(dict(target=y, landcover=lc, N=int(m.sum()), median=np.median(yv[m]), lo=lo[0], hi=hi[0]))
    return pd.DataFrame(rows)


def map_axes(fig, pos):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    ax = fig.add_subplot(*pos, projection=ccrs.PlateCarree())
    ax.set_extent((-118, -106, 29, 40), crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor='0.35')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.6, edgecolor='0.25')
    return ax


def edges(c):
    c = np.unique(np.round(np.asarray(c, float), 4))
    mid = 0.5 * (c[:-1] + c[1:])
    return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def squares(ax, t, v, cmap, vmin, vmax):
    import cartopy.crs as ccrs
    xc, xe = edges(t.lon)
    yc, ye = edges(t.lat)
    ix, iy = np.searchsorted(xc, np.round(t.lon, 4)), np.searchsorted(yc, np.round(t.lat, 4))
    verts = [[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])] for i, j in zip(ix, iy)]
    pc = PolyCollection(verts, transform=ccrs.PlateCarree(), linewidths=0, cmap=cmap)
    pc.set_array(np.ma.masked_invalid(np.asarray(v, float)))
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def make_maps(t, out):
    tt = t.copy()
    for c in [c for c in tt.columns if ':' in c] + OPER + ['sfmc_mean']:
        tt.loc[~tt.l1area, c] = np.nan     # L1 area only, so all panels share one tile set
    panels = [
        [('rs_smterm', 'rs_smterm (SM-term skill)', 'Purples', 0, 0.6),
         ('frac_out15', 'frac_out15 (|inn|>15 dB)', 'Purples', 0, 0.015),
         ('L1:SMAP_Tb', 'L1 DA: SMAP Tb %', 'RdBu_r', -8, 8),
         ('L1:L3', 'L1 DA: CYGNSS L3 %', 'RdBu_r', -8, 8),
         ('L1:ASCAT', 'L1 DA: ASCAT %', 'RdBu_r', -8, 8)],
        [('sfmc_mean', 'sfmc_mean (OL SM, m3/m3)', 'Purples', 0.03, 0.15),
         ('R_impl', 'R_impl (dB^2)', 'Purples', 0, 20),
         ('L3:SMAP_Tb', 'L3 DA: SMAP Tb %', 'RdBu_r', -8, 8),
         ('L3:L3', 'L3 DA: CYGNSS L3 % (assimilated)', 'RdBu_r', -8, 8),
         ('L3:ASCAT', 'L3 DA: ASCAT %', 'RdBu_r', -8, 8)],
        [('L1:L1', 'L1 DA: CYGNSS L1 % (assimilated)', 'RdBu_r', -3, 3),
         ('L1:SMOS_Tb', 'L1 DA: SMOS Tb %', 'RdBu_r', -8, 8),
         ('L1-L3:SMAP_Tb', 'L1 minus L3 arm: SMAP Tb (pp)', 'RdBu_r', -8, 8),
         ('L1-L3:SMOS_Tb', 'L1 minus L3 arm: SMOS Tb (pp)', 'RdBu_r', -8, 8),
         ('L1-L3:ASCAT', 'L1 minus L3 arm: ASCAT (pp)', 'RdBu_r', -8, 8)]]
    fig = plt.figure(figsize=(19, 11))
    for i, row in enumerate(panels):
        for j, (c, title, cmap, vmin, vmax) in enumerate(row):
            ax = map_axes(fig, (3, 5, i * 5 + j + 1))
            pc = squares(ax, tt, tt[c], cmap, vmin, vmax)
            ax.set_title(title, fontsize=9)
            ax.text(0.02, 0.03, f'median {np.nanmedian(tt[c]):+.3g}, N={np.isfinite(tt[c]).sum()}',
                    transform=ax.transAxes, fontsize=7, bbox=dict(fc='white', ec='none', alpha=0.8))
            plt.colorbar(pc, ax=ax, shrink=0.75, extend='both')
    fig.suptitle('L1-area tiles (n_l1_obs>=100), 2020-2022: operator skill vs % change in O-F stdv vs OL '
                 '(negative = DA better)')
    fig.tight_layout()
    fig.savefig(out / 'maps_operator_vs_da_benefit.png', dpi=120, facecolor='white')
    plt.close(fig)


def plot_bins(bt, x, q, out):
    gs = ['SMAP_Tb', 'SMOS_Tb', 'L3', 'ASCAT', 'L1']
    fig, axs = plt.subplots(1, len(gs), figsize=(17, 3.6), sharex=True)
    for ax, g in zip(axs, gs):
        for k, arm in enumerate(['L1', 'L1err39', 'L3']):   # SMAP arm (own-fit -12%, L3/ASCAT +) only in CSV
            s = bt[(bt.target == f'{arm}:{g}') & (bt.bin.astype(str).str.isdigit())]
            if s.empty:
                continue
            xb = s.bin.astype(int).values + (k - 1) * 0.08
            own = ' (assim.)' if (arm, g) in OWN else ''
            ax.errorbar(xb, s['median'], yerr=[s['median'] - s.lo, s.hi - s['median']], fmt='o-', ms=4, lw=1.5,
                        capsize=2, color=ARM_COL[arm], label=arm + own)
        ax.axhline(0, color='0.6', lw=0.7)
        ax.set_title(g.replace('_', ' '), fontsize=10)
        ax.set_xlabel(f'{x} {"tercile" if q == 3 else "quintile"} (low -> high)')
        ax.legend(fontsize=7, frameon=False)
        ax.grid(axis='y', color='0.92')
    axs[0].set_ylabel('tile median % change O-F stdv vs OL\n(95% block-bootstrap CI)')
    fig.tight_layout()
    fig.savefig(out / f'bins_{x}_q{q}.png', dpi=120, facecolor='white')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stats', type=Path, default=P / 'output/postproc_paired_density/stats_output')
    ap.add_argument('--period', default='20200101_20221231')
    ap.add_argument('--nmin', type=int, default=50)
    ap.add_argument('--l1min', type=int, default=100)
    ap.add_argument('--block', type=int, default=3)
    ap.add_argument('--nboot', type=int, default=1000)
    ap.add_argument('--out', type=Path, default=P / 'output/operator_skill/da_link/tile_whole_period')
    ap.add_argument('--replot-bins', action='store_true', help='only redraw bins_*.png from existing bins_*.csv')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.replot_bins:
        for x in ['rs_smterm', 'frac_out15', 'R_impl', 'sfmc_mean']:
            for q in [3, 5]:
                plot_bins(pd.read_csv(a.out / f'bins_{x}_q{q}.csv'), x, q, a.out)
        return
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30); pd.set_option('display.max_rows', 500)

    t = build_table(a)
    t.to_csv(a.out / 'tile_table.csv', index=False)
    boot = Boot(t.block.values, a.nboot)
    print(f'L1-area tiles with operator metrics: {t.l1area.sum()}; blocks: {t.loc[t.l1area, "block"].nunique()} '
          f'({a.block}x{a.block}); nboot {a.nboot}')

    ycols = [c for c in t.columns if ':' in c and '@' not in c]
    print('\n== tile medians of % change (L1 area) ==')
    print(t.loc[t.l1area, ycols].agg(['median', 'mean', 'count']).T.round(2).to_string())

    sp = spearman_table(t, boot, ycols, OPER + ATTR, CONF, a.period)
    sp.to_csv(a.out / 'spearman_whole_period.csv', index=False)
    print('\n== Spearman (raw | partial on sfmc_mean,n_l1_obs,lat,lon) for key predictors ==')
    key = sp[sp.predictor.isin(KEY_X + ['rs_H', 'snr', 'inn_std', 'dyn_ratio', 'r_SMAPh', 'pr_geom_inc'])]
    print(key[['target', 'predictor', 'N', 'rs', 'rs_lo', 'rs_hi', 'prs', 'prs_lo', 'prs_hi', 'prs_sig']]
          .round(3).to_string(index=False))

    mv = mv_table(t, boot, ycols)
    mv.to_csv(a.out / 'multivar_rank_regression.csv', index=False)
    print('\n== multivariate rank regression (standardized betas) ==')
    print(mv.pivot_table(index='target', columns='predictor', values='beta')[MV_SET + CONF].round(2).to_string())
    print(mv[mv.sig & mv.predictor.isin(MV_SET)].round(3).to_string(index=False))

    lc = landcover_table(t, boot, ycols)
    lc.to_csv(a.out / 'landcover_groups.csv', index=False)
    print('\n== by landcover (7 open shrub, 10 grass, 16 barren) ==')
    print(lc.round(2).to_string(index=False))

    for x in ['rs_smterm', 'frac_out15', 'R_impl', 'sfmc_mean']:
        for q in [3, 5]:
            bt = bin_table(t, boot, x, q, ycols)
            bt.to_csv(a.out / f'bins_{x}_q{q}.csv', index=False)
            plot_bins(bt, x, q, a.out)
            if q == 5 or x in ['rs_smterm', 'frac_out15']:
                print(f'\n== {x} q{q} (tile median % change, 95% block CI) ==')
                bt['s'] = bt.apply(lambda r: f'{r["median"]:+.2f} [{r.lo:+.2f},{r.hi:+.2f}]', axis=1)
                print(bt.pivot_table(index='target', columns='bin', values='s', aggfunc='first').to_string())
                if q == 5:
                    print(bt.drop_duplicates('bin')[['bin', 'x_min', 'x_max', 'N']].to_string(index=False))

    yr = []
    for per in YEARS:
        y4 = '@' + per[:4]
        ycy = [f'{arm}:{g}{y4}' for arm in ['L1', 'L1err39', 'L3', 'SMAP'] for g in GROUPS] + \
              [f'L1-L3:{g}{y4}' for g in GROUPS]
        yr.append(spearman_table(t, boot, ycy, KEY_X, CONF, per))
    yr = pd.concat(yr)
    yr['target'] = yr.target.str.split('@').str[0]
    yr.to_csv(a.out / 'yearly_spearman.csv', index=False)
    whole = sp[sp.predictor.isin(KEY_X)].assign(period='2020-22')
    both = pd.concat([whole, yr])
    both['s'] = both.apply(lambda r: f'{r.rs:+.2f}{"*" if r.rs_sig else " "}/{r.prs:+.2f}{"*" if r.prs_sig else " "}',
                           axis=1)
    print('\n== yearly robustness: raw/partial Spearman (* = 95% block CI excludes 0) ==')
    print(both.pivot_table(index=['predictor', 'target'], columns='period', values='s', aggfunc='first').to_string())

    make_maps(t, a.out)
    print('\ndone ->', a.out)


if __name__ == '__main__':
    main()
