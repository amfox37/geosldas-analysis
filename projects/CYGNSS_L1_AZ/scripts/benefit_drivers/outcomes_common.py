"""Shared helpers for the outcome-structure analysis (EXPLORATORY). Run with python -I from this directory's parent or by path."""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
BD = ROOT / 'output' / 'benefit_drivers'
OUT = BD / 'analysis_outcomes'
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))

SMAP_FAM = ['out_iv_smapL3', 'out_omf_smap', 'out_omf_smos']
ASCAT_FAM = ['out_tc_dR2', 'out_iv_ascatH121', 'out_omf_ascat']
INDEP = SMAP_FAM + ASCAT_FAM
EXTRA = ['out_tc_dSig2pct', 'out_omf_l3obs', 'out_omf_l1obs']
LABEL = {'out_tc_dR2': 'TC dR2 (model)', 'out_tc_dSig2pct': 'TC dSigma2 %', 'out_iv_smapL3': 'IV dR SMAP L3',
         'out_iv_ascatH121': 'IV dR ASCAT H121', 'out_omf_smap': 'O-F SMAP', 'out_omf_smos': 'O-F SMOS',
         'out_omf_ascat': 'O-F ASCAT', 'out_omf_l3obs': 'O-F L3 obs', 'out_omf_l1obs': 'O-F L1 obs (own)',
         'out_omf_indep': 'O-F indep (SMAP+SMOS)', 'out_combined': 'combined (z)'}

def load():
    a = pd.read_csv(BD / 'analysis_table.csv')
    t = pd.read_csv(BD / 'discover' / 'tc_iv_per_tile.csv')
    keep = ['tile_axis_index', 'i_indg', 'j_indg'] + [c for c in t.columns if c.startswith('tc_')]
    a = a.merge(t[keep], on='tile_axis_index', how='left')
    a['block'] = (np.floor(a.com_lon / 2).astype(int).astype(str) + '_' + np.floor(a.com_lat / 2).astype(int).astype(str))
    return a

def m_only(a):
    return a[a.in_l1_mask.astype(bool)].copy()

def spearman_pairs(df, cols):
    return df[cols].corr(method='spearman')

def block_boot_idx(df, rng):
    blocks = df.block.values
    ub = np.unique(blocks)
    groups = {b: np.where(blocks == b)[0] for b in ub}
    pick = rng.choice(ub, len(ub), replace=True)
    return np.concatenate([groups[b] for b in pick])

def boot_corr(df, cols, nboot=1000, seed=1):
    rng = np.random.default_rng(seed)
    ub = np.unique(df.block.values)
    groups = {b: np.where(df.block.values == b)[0] for b in ub}
    X = df[cols].reset_index(drop=True)
    res = np.full((nboot, len(cols), len(cols)), np.nan)
    for k in range(nboot):
        pick = rng.choice(ub, len(ub), replace=True)
        idx = np.concatenate([groups[b] for b in pick])
        res[k] = X.iloc[idx].corr(method='spearman').values
    return res

def boot_stat(df, fn, nboot=1000, seed=2):
    """block bootstrap of scalar fn(df_resampled); returns (est, lo, hi)."""
    rng = np.random.default_rng(seed)
    ub = np.unique(df.block.values)
    groups = {b: np.where(df.block.values == b)[0] for b in ub}
    d = df.reset_index(drop=True)
    vals = []
    for k in range(nboot):
        pick = rng.choice(ub, len(ub), replace=True)
        idx = np.concatenate([groups[b] for b in pick])
        vals.append(fn(d.iloc[idx]))
    vals = np.array(vals, float)
    return fn(d), np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)

def tile_map(ax, a, col, title, mask_col='in_l1_mask', cmap='RdBu', sym=True, vmin=None, vmax=None, s=22, cat=None):
    import cartopy.crs as ccrs, cartopy.feature as cfeature
    inm = a[mask_col].astype(bool)
    bad = a[~(inm & a[col].notna())] if cat is None else a[~inm]
    ax.scatter(bad.com_lon, bad.com_lat, s=s, marker='s', c='0.8', lw=0, transform=ccrs.PlateCarree())
    g = a[inm & a[col].notna()]
    if cat is None:
        if sym and vmin is None:
            vmax = np.nanpercentile(np.abs(g[col]), 95); vmin = -vmax
        sc = ax.scatter(g.com_lon, g.com_lat, s=s, marker='s', c=g[col], cmap=cmap, vmin=vmin, vmax=vmax, lw=0, transform=ccrs.PlateCarree())
        plt.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
    else:
        for k, (c, v) in cat.items():
            gg = g[g[col] == k]
            ax.scatter(gg.com_lon, gg.com_lat, s=s, marker='s', c=c, lw=0, label=f'{k} ({len(gg)})', transform=ccrs.PlateCarree())
        ax.legend(fontsize=6, loc='lower left')
    ax.add_feature(cfeature.STATES, lw=0.4, edgecolor='0.3')
    ax.set_title(title, fontsize=8)
    ax.set_extent([a.com_lon.min() - .3, a.com_lon.max() + .3, a.com_lat.min() - .3, a.com_lat.max() + .3])

def mapfig(n, ncols=3, size=(4.2, 3.6)):
    import cartopy.crs as ccrs
    nr = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nr, ncols, figsize=(size[0] * ncols, size[1] * nr), subplot_kw=dict(projection=ccrs.PlateCarree()))
    return fig, np.atleast_1d(axs).ravel()

def scale_only(x, mask):
    return x / x[mask].std()
