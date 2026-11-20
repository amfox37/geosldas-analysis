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

_X0, _Y0, _CS, _NX = -17367530.45, 7314540.83, 36032.22, 964   # EASEv2 M36 corner and cell size (EPSG:6933)

def cell_mesh(a):
    """Cell-edge lon/lat meshes and (j, i) indices for the tiles in `a` (one tile per M36 cell, idx0 = j*964 + i)."""
    from pyproj import Transformer
    i, j = (a.idx0 % _NX).astype(int).values, (a.idx0 // _NX).astype(int).values
    i0, i1, j0, j1 = i.min(), i.max(), j.min(), j.max()
    X, Y = np.meshgrid(_X0 + np.arange(i0, i1 + 2) * _CS, _Y0 - np.arange(j0, j1 + 2) * _CS)
    lon, lat = Transformer.from_crs(6933, 4326, always_xy=True).transform(X, Y)
    return lon, lat, j - j0, i - i0

def tile_map(ax, a, col, title, mask_col='in_l1_mask', cmap='RdBu', sym=True, vmin=None, vmax=None, s=None, cat=None):
    """Draw tiles as filled M36 cells (pcolormesh, no gaps). Grey = outside the mask or no valid value."""
    import cartopy.crs as ccrs, cartopy.feature as cfeature
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    lon, lat, jj, ii = cell_mesh(a)
    shp = (lon.shape[0] - 1, lon.shape[1] - 1)
    inm = a[mask_col].astype(bool).values
    ok = (inm & a[col].notna().values) if cat is None else (inm & a[col].isin(list(cat)).values)
    kw = dict(transform=ccrs.PlateCarree(), shading='flat', edgecolors='face', linewidth=0.2, rasterized=True)
    grey = np.full(shp, np.nan); grey[jj[~ok], ii[~ok]] = 1
    ax.pcolormesh(lon, lat, np.ma.masked_invalid(grey), cmap=ListedColormap(['0.82']), **kw)
    v = np.full(shp, np.nan)
    if cat is None:
        x = a[col].values
        v[jj[ok], ii[ok]] = x[ok]
        if sym and vmin is None:
            vmax = np.nanpercentile(np.abs(x[ok]), 95); vmin = -vmax
        pc = ax.pcolormesh(lon, lat, np.ma.masked_invalid(v), cmap=cmap, vmin=vmin, vmax=vmax, **kw)
        plt.colorbar(pc, ax=ax, shrink=0.75, pad=0.02)
    else:
        keys = list(cat)
        code = {k: n for n, k in enumerate(keys)}
        v[jj[ok], ii[ok]] = [code[k] for k in a[col].values[ok]]
        ax.pcolormesh(lon, lat, np.ma.masked_invalid(v), cmap=ListedColormap([cat[k][0] for k in keys]), vmin=-0.5, vmax=len(keys) - 0.5, **kw)
        ax.legend(handles=[Patch(color=cat[k][0], label=f'{k} ({int((a[col].values[ok] == k).sum())})') for k in keys], fontsize=6, loc='lower left')
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
