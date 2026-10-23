#!/usr/bin/env python3
"""
TASK C (four_corners_ascat): does the CYGNSS L1 operator explain the ASCAT O-F penalty in the Four Corners area?

Inputs (read-only):
  output/postproc_paired_density/<exptag>/Y*/M*/<exptag>.ens_avg.ldas_ObsFcstAna_sums.<yyyymm>.nc4
      monthly per-tile raw sums (N, sum O, F, O^2, F^2, OF) for DA_paired_<tag> and OL_paired_monitor_xmask_<tag>;
      pooled over any set of months -> per-tile O-F stdv exactly as the toolkit's calc_temporal_stats_from_sums
      (checked against the whole-period temporal_stats nc4 at start-up).
  output/operator_skill/da_link/tile_axis_ig_jg.csv          tile axis -> ig, jg, lon, lat, n_l1_obs
  output/operator_skill/analysis_202001_202212/tile_metrics.csv   OL operator metrics per owner tile
  output/operator_skill/obs_table/ol_l1_obs_<yyyymm>.parquet      per-obs OL table (innovations, state, QC)
  output/ismn_fixedop_20200101_20221231/ismn_skill_stations.csv   ISMN station skill (for overlay/adjudication)

Region definition (--region hotspot, default): per-tile ASCAT % change of the primary L1 arm (2020-2022, N>=50 per
  species, group = mean over MetOp-A/B/C), smoothed with the mean over the (2r+1)x(2r+1) ig/jg neighbourhood of
  L1-area tiles (>= --smooth-min-valid valid), thresholded at the --hot-quantile of the smoothed field over the
  L1 area; the region is the 8-connected component of above-threshold L1-area tiles nearest to --seed (Four
  Corners, -109.05/37.0).  --region box uses --box lon0 lon1 lat0 lat1 instead.
Controls: rest of the L1 area (primary); rest minus a 1-tile buffer; same-latitude band outside the buffer;
  attribute-matched tiles (nearest neighbour on standardized altitude, slope, opacity, mean and std SFMC).
Statistics: stratified spatial block bootstrap (blocks of --block x --block tiles in ig/jg) for every
  region-vs-control difference.

Outputs in --out-dir: CSV tables, four_corners_ascat_map.png, SUMMARY numbers printed to the log.
Usage:
  four_corners_ascat.py [--region hotspot|box] [--hot-quantile 0.75] [--box -111.5 -107 34.5 37.5] [--nboot 2000]
"""
import argparse
import glob
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats, ndimage
from netCDF4 import Dataset

PROJ = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/'
PP = PROJ + 'output/postproc_paired_density/'
TILES = PROJ + 'output/operator_skill/da_link/tile_axis_ig_jg.csv'
TMET = PROJ + 'output/operator_skill/analysis_202001_202212/tile_metrics.csv'
OBST = PROJ + 'output/operator_skill/obs_table/ol_l1_obs_{}.parquet'
ISMN = PROJ + 'output/ismn_fixedop_20200101_20221231/ismn_skill_stations.csv'
ARMS = {'full_xc015_coh040216_fixedop': 'L1', 'full_xc015_coh040216_err39_fixedop': 'L1 err39',
        'l3_fixedop': 'L3', 'smap_fixedop': 'SMAP'}
PRIMARY = 'full_xc015_coh040216_fixedop'
GROUPS = {'SMOS Tb': [0, 1, 2, 3], 'SMAP Tb': [4, 5, 6, 7], 'ASCAT': [8, 9, 10], 'CYGNSS L3': [11], 'CYGNSS L1': [12]}
MONTHS = pd.period_range('2020-01', '2022-12', freq='M')
SEASON = {12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA', 7: 'JJA', 8: 'JJA',
          9: 'SON', 10: 'SON', 11: 'SON'}
COLD = {11, 12, 1, 2, 3}
DOMAIN = (-118.0, -106.0, 29.0, 40.0)
ORBIT_LAT = 37.9          # CYGNSS specular-point northern limit (project README)
SEL_BUFFER_KM = 200.      # preprocessing region = Arizona + 200 km
L1_MIN = 100
FREEZE = 273.15


# ------------------------------------------------------------------ O-F sums
def load_sums(exptag):
    keys = ['N_data', 'obs_obs_sum', 'obs_fcst_sum', 'obs_obs2_sum', 'obs_fcst2_sum', 'obsxfcst_sum']
    out = {k: [] for k in keys}
    for m in MONTHS:
        fn = f'{PP}{exptag}/Y{m.year}/M{m.month:02d}/{exptag}.ens_avg.ldas_ObsFcstAna_sums.{m.year}{m.month:02d}.nc4'
        with Dataset(fn) as nc:
            for k in keys:
                out[k].append(np.ma.filled(nc[k][:].astype(float), np.nan))
    return {k: np.nan_to_num(np.stack(v)) for k, v in out.items()}      # (month, tile, species)


def pool(s, msel):
    """Per-tile stats pooled over the selected months (same formulas as the toolkit, population moments)."""
    g = {k: v[msel].sum(axis=0) for k, v in s.items()}
    N = g['N_data']
    with np.errstate(all='ignore'):
        mO, mF = g['obs_obs_sum'] / N, g['obs_fcst_sum'] / N
        vO = g['obs_obs2_sum'] / N - mO ** 2
        vF = g['obs_fcst2_sum'] / N - mF ** 2
        cov = g['obsxfcst_sum'] / N - mO * mF
        omf = np.sqrt(np.maximum(vO + vF - 2 * cov, 0))
        r = cov / np.sqrt(vO * vF)
    bad = N < 1
    for a in (mO, mF, vO, vF, omf, r):
        a[bad] = np.nan
    return dict(N=N, O_mean=mO, F_mean=mF, F_std=np.sqrt(vF), O_std=np.sqrt(vO), OmF_std=omf, OmF_var=omf ** 2,
                OmF_mean=mO - mF, r=r)


def pct(da, ol, idx, nmin):
    d, o, n = da['OmF_std'][:, idx], ol['OmF_std'][:, idx], da['N'][:, idx]
    p = np.where((n >= nmin) & (o > 0), 100. * (d - o) / o, np.nan)
    return np.nanmean(p, axis=1)


def grp_mean(x, idx, n, nmin):
    v = np.where(n[:, idx] >= nmin, x[:, idx], np.nan)
    return np.nanmean(v, axis=1)


def pooled_pct(da, ol, idx, tiles_mask, nmin=1):
    """Region-pooled within-tile O-F stdv change: sqrt(sum N var_DA / sum N var_OL) - 1, mean over species."""
    out = []
    for i in idx:
        n = da['N'][:, i]
        ok = tiles_mask & (n >= nmin) & np.isfinite(da['OmF_var'][:, i]) & np.isfinite(ol['OmF_var'][:, i])
        if ok.sum() == 0:
            continue
        vd = (n[ok] * da['OmF_var'][ok, i]).sum()
        vo = (ol['N'][ok, i] * ol['OmF_var'][ok, i]).sum() * n[ok].sum() / ol['N'][ok, i].sum()
        out.append(100. * (np.sqrt(vd / vo) - 1))
    return np.mean(out) if out else np.nan


def decompose(da, ol, idx, nmin):
    """Split the O-F stdv change into amplitude and correlation parts, per tile (mean over species):
    OmF^2 = sO^2 + sF^2 - 2 r sO sF with the OL obs stdv sO.
      amp_only : DA forecast stdv, OL correlation   corr_only: OL forecast stdv, DA correlation"""
    sO, sF0, r0 = ol['O_std'][:, idx], ol['F_std'][:, idx], ol['r'][:, idx]
    sF1, r1 = da['F_std'][:, idx], da['r'][:, idx]
    base = ol['OmF_std'][:, idx]
    ok = (da['N'][:, idx] >= nmin) & (base > 0)
    f = lambda sF, r: np.sqrt(np.maximum(sO ** 2 + sF ** 2 - 2 * r * sO * sF, 0))
    out = {}
    for k, v in [('amp_only', f(sF1, r0)), ('corr_only', f(sF0, r1)), ('both_recon', f(sF1, r1))]:
        out[k] = np.nanmean(np.where(ok, 100 * (v / base - 1), np.nan), axis=1)
    out['dFstd_pct'] = np.nanmean(np.where(ok, 100 * (sF1 / sF0 - 1), np.nan), axis=1)
    out['dr'] = np.nanmean(np.where(ok, r1 - r0, np.nan), axis=1)
    return out


# ------------------------------------------------------------------ snow / frozen from the OL daily output
OLCAT = '/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_AZ_fixedop/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/'
SNOW_MIN = 1.0   # kg m-2 snow mass for a "snow day"


def snow_frozen(cache, t):
    """Per tile and month: fraction of days with OL snow mass > SNOW_MIN, daily TSOIL1 < 0 C, TSURF < 0 C."""
    if os.path.exists(cache):
        z = np.load(cache)
        return z['snow'], z['fsoil'], z['fsurf']
    snow, fsoil, fsurf = (np.zeros((len(MONTHS), len(t))) for _ in range(3))
    for k, m in enumerate(MONTHS):
        fs = sorted(glob.glob(f'{OLCAT}Y{m.year}/M{m.month:02d}/*.tavg24_1d_lnd_Nt.*.nc4'))
        for fn in fs:
            with Dataset(fn) as nc:
                if k == 0 and fn == fs[0]:
                    assert np.abs(nc['lon'][:] - t.com_lon.values).max() < 1e-3
                snow[k] += nc['SNOMASLAND'][0, :] > SNOW_MIN
                fsoil[k] += nc['TSOIL1'][0, :] < FREEZE
                fsurf[k] += nc['TSURFLAND'][0, :] < FREEZE
        snow[k] /= len(fs); fsoil[k] /= len(fs); fsurf[k] /= len(fs)
    np.savez(cache, snow=snow, fsoil=fsoil, fsurf=fsurf)
    return snow, fsoil, fsurf


# ------------------------------------------------------------------ geometry
def haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371. * 2 * np.arcsin(np.sqrt(a))


def arizona_polygon():
    import cartopy.io.shapereader as shp
    for rec in shp.Reader(shp.natural_earth('50m', 'cultural', 'admin_1_states_provinces_lakes')).records():
        if rec.attributes.get('name') == 'Arizona' and rec.attributes.get('admin') == 'United States of America':
            return rec.geometry
    raise SystemExit('Arizona polygon not found')


def distances(t):
    """Distance [km] of each tile centre to the preprocessing boundary (AZ + 200 km), the orbit limit, the
    domain box, and the empirical L1-area edge (nearest grid cell that is not an L1-area tile)."""
    import pyproj
    from shapely.ops import transform
    from shapely.geometry import Point
    az = arizona_polygon()
    aeqd = pyproj.Proj(proj='aeqd', lat_0=34.2, lon_0=-111.7, datum='WGS84')
    fwd = pyproj.Transformer.from_proj(pyproj.Proj('epsg:4326'), aeqd, always_xy=True).transform
    inv = pyproj.Transformer.from_proj(aeqd, pyproj.Proj('epsg:4326'), always_xy=True).transform
    sel = transform(fwd, az).buffer(SEL_BUFFER_KM * 1e3, resolution=64)
    xy = [fwd(lo, la) for lo, la in zip(t.com_lon, t.com_lat)]
    inside = np.array([sel.contains(Point(p)) for p in xy])
    dsel = np.array([sel.exterior.distance(Point(p)) / 1e3 for p in xy])
    t['d_sel_km'] = np.where(inside, dsel, -dsel)                      # >0 inside the preprocessing region
    t['d_orbit_km'] = (ORBIT_LAT - t.com_lat) * 111.2
    t['d_cov_km'] = np.minimum(t.d_sel_km, t.d_orbit_km)               # coverage boundary (preproc + orbit)
    lon0, lon1, lat0, lat1 = DOMAIN
    t['d_box_km'] = np.minimum.reduce([(t.com_lon - lon0) * 111.2 * np.cos(np.radians(t.com_lat)),
                                       (lon1 - t.com_lon) * 111.2 * np.cos(np.radians(t.com_lat)),
                                       (t.com_lat - lat0) * 111.2, (lat1 - t.com_lat) * 111.2])
    # empirical L1 edge on the ig/jg grid, padded by one cell (outside-domain and water cells count as edge)
    lon_of = t.groupby('ig').com_lon.mean()
    lat_of = t.groupby('jg').com_lat.mean()
    igs = np.arange(t.ig.min() - 1, t.ig.max() + 2)
    jgs = np.arange(t.jg.min() - 1, t.jg.max() + 2)
    glon = np.interp(igs, lon_of.index, lon_of.values, left=np.nan, right=np.nan)
    dlon = np.diff(lon_of.values).mean()
    glon[0], glon[-1] = lon_of.values[0] - dlon, lon_of.values[-1] + dlon
    glat = np.interp(jgs, lat_of.index, lat_of.values)
    dlat0, dlat1 = lat_of.values[1] - lat_of.values[0], lat_of.values[-1] - lat_of.values[-2]
    glat[0], glat[-1] = lat_of.values[0] - dlat0, lat_of.values[-1] + dlat1
    inl1 = set(zip(t.ig[t.l1], t.jg[t.l1]))
    gi, gj = np.meshgrid(igs, jgs, indexing='ij')
    edge = np.array([(i, j) not in inl1 for i, j in zip(gi.ravel(), gj.ravel())])
    elon = np.array([glon[i - igs[0]] for i in gi.ravel()[edge]])
    elat = np.array([glat[j - jgs[0]] for j in gj.ravel()[edge]])
    d_edge = np.array([haversine(lo, la, elon, elat).min() for lo, la in zip(t.com_lon, t.com_lat)])
    l1lon, l1lat = t.com_lon[t.l1].values, t.com_lat[t.l1].values
    d_l1 = np.array([haversine(lo, la, l1lon, l1lat).min() for lo, la in zip(t.com_lon, t.com_lat)])
    t['d_l1edge_km'] = np.where(t.l1, d_edge, -d_l1)                   # >0 inside the L1 area
    # buffer polygon in lon/lat for the map
    return transform(inv, sel)


# ------------------------------------------------------------------ region
def grid_of(t, v):
    G = np.full((t.ig.max() + 1, t.jg.max() + 1), np.nan)
    G[t.ig.values, t.jg.values] = v
    return G


def define_region(t, a):
    if a.region == 'box':
        lon0, lon1, lat0, lat1 = a.box
        reg = t.l1 & t.com_lon.between(lon0, lon1) & t.com_lat.between(lat0, lat1)
        return reg.values, f'box lon {lon0}..{lon1}, lat {lat0}..{lat1}, L1-area tiles', np.nan
    v = np.where(t.l1, t.asc_L1, np.nan)
    G = grid_of(t, v)
    r = a.smooth
    sm = np.full(len(t), np.nan)
    for k, (i, j) in enumerate(zip(t.ig.values, t.jg.values)):
        if not t.l1.iloc[k]:
            continue
        w = G[max(i - r, 0):i + r + 1, max(j - r, 0):j + r + 1]
        if np.isfinite(w).sum() >= a.smooth_min_valid:
            sm[k] = np.nanmean(w)
    t['asc_L1_smooth'] = sm
    thr = np.nanquantile(sm[t.l1.values], a.hot_quantile)
    hot = t.l1.values & (sm >= thr)
    H = np.zeros((t.ig.max() + 2, t.jg.max() + 2), bool)
    H[t.ig.values[hot], t.jg.values[hot]] = True
    lab, nlab = ndimage.label(H, structure=np.ones((3, 3)))
    tl = lab[t.ig.values, t.jg.values]
    dseed = haversine(t.com_lon.values, t.com_lat.values, a.seed[0], a.seed[1])
    best = min(range(1, nlab + 1), key=lambda c: dseed[tl == c].min())
    reg = tl == best
    sizes = sorted([(tl == c).sum() for c in range(1, nlab + 1)], reverse=True)
    desc = (f'hotspot: {(2 * r + 1)}x{(2 * r + 1)}-tile smoothed L1-arm ASCAT % change >= q{a.hot_quantile:.2f} '
            f'({thr:+.2f}%) of L1-area tiles; 8-connected component nearest ({a.seed[0]}, {a.seed[1]}); '
            f'{nlab} components above threshold, sizes {sizes[:6]}')
    return reg, desc, thr


def matched_control(t, reg, pool_mask, cols, rng):
    z = t[cols].astype(float)
    z = (z - z[pool_mask | reg].mean()) / z[pool_mask | reg].std()
    cand = np.where(pool_mask & z.notna().all(axis=1).values)[0]
    out = []
    for k in np.where(reg & z.notna().all(axis=1).values)[0]:
        d = ((z.values[cand] - z.values[k]) ** 2).sum(axis=1)
        out.append(cand[np.argmin(d)])
    return np.array(out, int)


# ------------------------------------------------------------------ bootstrap
class Boot:
    """Stratified spatial block bootstrap of group means (and their difference)."""
    def __init__(self, blocks, nboot, seed):
        self.blocks, self.nboot, self.seed = blocks, nboot, seed

    def _prep(self, vals, idx):
        v, b = vals[idx], self.blocks[idx]
        ok = np.isfinite(v)
        v, b = v[ok], b[ok]
        ub, inv = np.unique(b, return_inverse=True)
        return np.bincount(inv, weights=v, minlength=len(ub)), np.bincount(inv, minlength=len(ub)).astype(float), len(v)

    def draws(self, vals, idx, rng):
        s, c, n = self._prep(vals, idx)
        if n == 0:
            return np.full(self.nboot, np.nan), np.nan, 0
        k = rng.integers(0, len(s), size=(self.nboot, len(s)))
        return s[k].sum(1) / c[k].sum(1), s.sum() / c.sum(), n

    def diff(self, vals, idx_a, idx_b):
        rng = np.random.default_rng(self.seed)
        da, ma, na = self.draws(vals, idx_a, rng)
        db, mb, nb = self.draws(vals, idx_b, rng)
        q = lambda x: np.nanpercentile(x, [2.5, 97.5])
        return dict(mean_reg=ma, reg_lo=q(da)[0], reg_hi=q(da)[1], n_reg=na,
                    mean_ctl=mb, ctl_lo=q(db)[0], ctl_hi=q(db)[1], n_ctl=nb,
                    diff=ma - mb, diff_lo=q(da - db)[0], diff_hi=q(da - db)[1])


def boot_corr(x, y, blocks, idx, nboot, seed):
    """Pooled correlation of (tile-demeaned) x, y over obs idx, block bootstrap over spatial blocks."""
    ok = idx & np.isfinite(x) & np.isfinite(y)
    x, y, b = x[ok], y[ok], blocks[ok]
    if len(x) < 30:
        return np.nan, np.nan, np.nan, len(x)
    ub, inv = np.unique(b, return_inverse=True)
    S = np.stack([np.bincount(inv, weights=w, minlength=len(ub)) for w in (np.ones_like(x), x, y, x * x, y * y, x * y)])
    def r_of(T):
        n, sx, sy, sxx, syy, sxy = T
        return (sxy / n - sx * sy / n ** 2) / np.sqrt((sxx / n - (sx / n) ** 2) * (syy / n - (sy / n) ** 2))
    rng = np.random.default_rng(seed)
    k = rng.integers(0, len(ub), size=(nboot, len(ub)))
    rb = r_of(S[:, k].sum(axis=2))
    return r_of(S.sum(axis=1)), np.nanpercentile(rb, 2.5), np.nanpercentile(rb, 97.5), len(x)


# ------------------------------------------------------------------ obs table
def load_obs():
    cols = ['sp_nearest_tile_ig', 'sp_nearest_tile_jg', 'time', 'obs', 'fcst', 'inn_full', 'inn_coh', 'inn_ASCAT',
            'inn_SMAPh', 'inn_SMAPv', 'inn_L3', 'TSOIL1_FCST', 'TSURF_FCST', 'modis_land_cover', 'coherency_ratio',
            'w_in_domain']
    import pyarrow.parquet as pq
    have = set(pq.read_schema(OBST.format('202001')).names)
    parts = [pd.read_parquet(OBST.format(f'{m.year}{m.month:02d}'), columns=[c for c in cols if c in have])
             for m in MONTHS]
    d = pd.concat(parts, ignore_index=True)
    d = d.rename(columns={'sp_nearest_tile_ig': 'ig', 'sp_nearest_tile_jg': 'jg'})
    d['month'] = pd.to_datetime(d.time).dt.month
    d['ym'] = pd.to_datetime(d.time).dt.to_period('M')
    return d


def tile_obs_metrics(d, weak_h_q):
    qh = d.fcst.quantile(weak_h_q)
    d = d.assign(weakH=d.fcst < qh, frozen_soil=d.TSOIL1_FCST < FREEZE, frozen_surf=d.TSURF_FCST < FREEZE,
                 lc_shrub=d.modis_land_cover.isin([6, 7]), lc_grass=d.modis_land_cover == 10,
                 lc_woody=d.modis_land_cover.isin([1, 2, 3, 4, 5, 8, 9]), lc_barren=d.modis_land_cover == 16,
                 lc_crop=d.modis_land_cover.isin([12, 14]), coh_lo=d.coherency_ratio < 0.40,
                 coh_hi=d.coherency_ratio > 2.16, edge_fp=d.w_in_domain < 0.999)

    def c(g, x, y, nmin=30):
        ok = g[x].notna() & g[y].notna()
        return np.corrcoef(g.loc[ok, x], g.loc[ok, y])[0, 1] if ok.sum() >= nmin else np.nan
    pairs = [('inn_full', 'inn_ASCAT'), ('inn_coh', 'inn_ASCAT'), ('inn_full', 'inn_SMAPh'), ('inn_full', 'inn_L3'),
             ('inn_ASCAT', 'inn_L3'), ('inn_ASCAT', 'inn_SMAPh'), ('inn_L3', 'inn_SMAPh')]
    rows = []
    for (i, j), g in d.groupby(['ig', 'jg']):
        r = dict(ig=i, jg=j, n_obs_tab=len(g), frac_weakH=g.weakH.mean(), frac_frozen_soil=g.frozen_soil.mean(),
                 frac_frozen_surf=g.frozen_surf.mean(), lc_shrub=g.lc_shrub.mean(), lc_grass=g.lc_grass.mean(),
                 lc_woody=g.lc_woody.mean(), lc_barren=g.lc_barren.mean(), lc_crop=g.lc_crop.mean(),
                 coh_lo_frac=g.coh_lo.mean(), coh_hi_frac=g.coh_hi.mean(), frac_edge_fp=g.edge_fp.mean(),
                 n_ascat_pairs=int((g.inn_full.notna() & g.inn_ASCAT.notna()).sum()))
        for x, y in pairs:
            r[f'r_{x[4:]}_{y[4:]}'] = c(g, x, y)
        rows.append(r)
    return pd.DataFrame(rows), qh


# ------------------------------------------------------------------ map
def map_panel(fig, pos, t, v, title, vmin, vmax, cmap, reg, ctl, sel_poly, st=None, cb_label=''):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.collections import PolyCollection, LineCollection
    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent(DOMAIN, crs=ccrs.PlateCarree())
    xc, xe = edges(t.com_lon)
    yc, ye = edges(t.com_lat)
    ix = np.searchsorted(xc, np.round(t.com_lon, 4))
    iy = np.searchsorted(yc, np.round(t.com_lat, 4))
    verts = [[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])] for i, j in zip(ix, iy)]
    nl1 = ~t.l1.values
    bg = PolyCollection([verts[k] for k in np.where(nl1)[0]], facecolor='0.93', linewidths=0, transform=ccrs.PlateCarree())
    ax.add_collection(bg)
    pc = PolyCollection(verts, array=np.ma.masked_invalid(np.where(nl1, np.nan, v)), cmap=cmap, linewidths=0,
                        transform=ccrs.PlateCarree())
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.4, edgecolor='0.45')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='0.3')
    # region outline: tile edges whose neighbour is not in the region
    cell = {(i, j): k for k, (i, j) in enumerate(zip(ix, iy))}
    segs = []
    for k in np.where(reg)[0]:
        i, j = ix[k], iy[k]
        for di, dj, s in [(-1, 0, [(xe[i], ye[j]), (xe[i], ye[j + 1])]), (1, 0, [(xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1])]),
                          (0, -1, [(xe[i], ye[j]), (xe[i + 1], ye[j])]), (0, 1, [(xe[i], ye[j + 1]), (xe[i + 1], ye[j + 1])])]:
            nb = cell.get((i + di, j + dj))
            if nb is None or not reg[nb]:
                segs.append(s)
    ax.add_collection(LineCollection(segs, colors='k', linewidths=1.6, transform=ccrs.PlateCarree()))
    if ctl is not None and len(ctl):
        ax.scatter(t.com_lon.values[ctl], t.com_lat.values[ctl], s=7, marker='s', facecolor='none', edgecolor='k',
                   linewidths=0.6, transform=ccrs.PlateCarree(), zorder=6)
    if sel_poly is not None:
        x, y = sel_poly.exterior.xy
        ax.plot(x, y, color='#7B3294', lw=1.0, ls='--', transform=ccrs.PlateCarree())
        ax.plot([DOMAIN[0], DOMAIN[1]], [ORBIT_LAT, ORBIT_LAT], color='#7B3294', lw=1.0, ls=':', transform=ccrs.PlateCarree())
    if st is not None:
        ax.scatter(st.lon, st.lat, marker='^', s=14, facecolor='white', edgecolor='k', linewidths=0.6,
                   transform=ccrs.PlateCarree(), zorder=7)
    ax.set_title(title, fontsize=9)
    cb = fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.8, pad=0.03, extend='both')
    cb.set_label(cb_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return ax


def edges(centers):
    c = np.unique(np.round(np.asarray(centers, dtype=float), 4))
    mid = 0.5 * (c[:-1] + c[1:])
    return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


# ------------------------------------------------------------------ main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--region', choices=['hotspot', 'box'], default='hotspot')
    p.add_argument('--hot-quantile', type=float, default=0.75)
    p.add_argument('--smooth', type=int, default=1, help='smoothing half-width in tiles (1 -> 3x3)')
    p.add_argument('--smooth-min-valid', type=int, default=4)
    p.add_argument('--seed', type=float, nargs=2, default=[-109.05, 37.0], metavar=('LON', 'LAT'))
    p.add_argument('--box', type=float, nargs=4, default=[-111.5, -107.0, 34.5, 37.5])
    p.add_argument('--nmin', type=int, default=50, help='min N per tile/species (as Fig. 7)')
    p.add_argument('--nmin-season', type=int, default=30, help='min N per tile/species for seasonal splits')
    p.add_argument('--block', type=int, default=3, help='bootstrap block size in tiles')
    p.add_argument('--nboot', type=int, default=2000)
    p.add_argument('--weak-h-quantile', type=float, default=0.125)
    p.add_argument('--out-dir', default=PROJ + 'output/operator_skill/da_link/four_corners_ascat/')
    p.add_argument('--tag', default='', help='suffix for output file names (sensitivity runs)')
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sfx = f'_{a.tag}' if a.tag else ''
    out = lambda n: f'{a.out_dir}{n}{sfx}'
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)

    t = pd.read_csv(TILES)
    t['l1'] = t.n_l1_obs >= L1_MIN
    t['block'] = (t.ig // a.block) * 1000 + t.jg // a.block
    blocks = t.block.values
    print(f'{len(t)} tiles, {t.l1.sum()} in the L1 area (n_l1_obs >= {L1_MIN})')

    # ---- sums, verification
    sums = {}
    for tag in ARMS:
        for kind in ['DA_paired', 'OL_paired_monitor_xmask']:
            sums[(kind, tag)] = load_sums(f'{kind}_{tag}')
    allm = np.ones(len(MONTHS), bool)
    chk = pool(sums[('DA_paired', PRIMARY)], allm)
    ref = xr.open_dataset(f'{PP}stats_output/temporal_stats_DA_paired_{PRIMARY}_20200101_20221231.nc4')
    dn = np.nanmax(np.abs(chk['N'] - ref.N_data.values))
    ok = (ref.N_data.values >= 50)
    rel = np.nanmax(np.abs(chk['OmF_std'][ok] / ref.OmF_stdv.values[ok] - 1))
    print(f'verification vs temporal_stats nc4 (2020-2022, DA {PRIMARY}): max |dN| {dn:.0f}, '
          f'max rel. diff OmF_stdv (N>=50) {rel:.2e}')
    m21 = np.array([m.year >= 2021 for m in MONTHS])
    ref21 = xr.open_dataset(f'{PP}stats_output/temporal_stats_OL_paired_monitor_xmask_l3_fixedop_20210101_20221231.nc4')
    c21 = pool(sums[('OL_paired_monitor_xmask', 'l3_fixedop')], m21)
    ok = ref21.N_data.values >= 50
    print(f'  2021-2022 OL l3: max rel. diff {np.nanmax(np.abs(c21["OmF_std"][ok] / ref21.OmF_stdv.values[ok] - 1)):.2e}')

    periods = {'2020-2022': allm, '2021-2022': m21}
    P = {(per, tag): (pool(sums[('DA_paired', tag)], ms), pool(sums[('OL_paired_monitor_xmask', tag)], ms))
         for per, ms in periods.items() for tag in ARMS}
    for tag, lab in ARMS.items():
        da, ol = P[('2020-2022', tag)]
        for g, idx in GROUPS.items():
            t[f'pct_{lab}_{g}'] = pct(da, ol, idx, a.nmin)
    t['asc_L1'] = t['pct_L1_ASCAT']
    da, ol = P[('2020-2022', PRIMARY)]
    nasc = ol['N'][:, GROUPS['ASCAT']]
    t['N_ascat'] = nasc.sum(1)
    t['olr_ascat'] = grp_mean(ol['r'], GROUPS['ASCAT'], ol['N'], a.nmin)            # OL corr(O, F), ASCAT
    t['olr_smap'] = grp_mean(ol['r'], GROUPS['SMAP Tb'], ol['N'], a.nmin)
    t['olr_l3'] = grp_mean(ol['r'], GROUPS['CYGNSS L3'], ol['N'], a.nmin)
    with np.errstate(all='ignore'):
        fo = ol['F_std'] / ol['O_std']
    t['olFO_ascat'] = grp_mean(fo, GROUPS['ASCAT'], ol['N'], a.nmin)      # OL forecast/obs stdv ratio (scaled obs)
    t['olFO_smap'] = grp_mean(fo, GROUPS['SMAP Tb'], ol['N'], a.nmin)
    t['olFO_l3'] = grp_mean(fo, GROUPS['CYGNSS L3'], ol['N'], a.nmin)
    t['olomf_ascat'] = grp_mean(ol['OmF_std'], GROUPS['ASCAT'], ol['N'], a.nmin)
    t['olO_std_ascat'] = grp_mean(ol['O_std'], GROUPS['ASCAT'], ol['N'], a.nmin)
    t['olF_std_ascat'] = grp_mean(ol['F_std'], GROUPS['ASCAT'], ol['N'], a.nmin)
    for tag, lab in ARMS.items():
        da, ol = P[('2020-2022', tag)]
        t[f'dr_{lab}_ascat'] = grp_mean(da['r'] - ol['r'], GROUPS['ASCAT'], ol['N'], a.nmin)
        t[f'dFstd_{lab}_ascat'] = grp_mean(100 * (da['F_std'] / ol['F_std'] - 1), GROUPS['ASCAT'], ol['N'], a.nmin)
        t[f'dOmFmean_{lab}_ascat'] = grp_mean(da['OmF_mean'] - ol['OmF_mean'], GROUPS['ASCAT'], ol['N'], a.nmin)

    # ---- operator + obs-table metrics
    tm = pd.read_csv(TMET).drop(columns=['lon', 'lat'])
    t = t.merge(tm, on=['ig', 'jg'], how='left')
    d = load_obs()
    to, qh = tile_obs_metrics(d, a.weak_h_quantile)
    t = t.merge(to, on=['ig', 'jg'], how='left')
    print(f'obs table: {len(d)} obs; weak-H threshold (q{a.weak_h_quantile}) H < {qh:.1f} dB')
    print(f'  max fraction of L1 obs with TSOIL1_FCST < 0 C per tile: {to.frac_frozen_soil.max():.3f} '
          f'(L1 obs are model-QC screened for frozen/snow, so the obs table cannot measure frozen conditions)')
    snow, fsoil, fsurf = snow_frozen(a.out_dir + 'ol_snow_frozen_monthly.npz', t)
    t['snow_day_frac'], t['frozen_soil_day_frac'], t['frozen_surf_day_frac'] = snow.mean(0), fsoil.mean(0), fsurf.mean(0)
    cold_m = np.isin([m.month for m in MONTHS], list(COLD))
    t['snow_day_frac_cold'], t['frozen_soil_day_frac_cold'] = snow[cold_m].mean(0), fsoil[cold_m].mean(0)

    # ---- distances, region, controls
    sel_poly = distances(t)
    reg, desc, thr = define_region(t, a)
    t['region'] = reg
    L1 = t.l1.values
    rest = L1 & ~reg
    G = grid_of(t, reg.astype(float))
    near = np.zeros(len(t), bool)
    for k, (i, j) in enumerate(zip(t.ig.values, t.jg.values)):
        near[k] = np.nansum(G[max(i - 1, 0):i + 2, max(j - 1, 0):j + 2]) > 0
    rest_buf = rest & ~near
    lat_band = rest_buf & t.com_lat.between(t.com_lat[reg].min(), t.com_lat[reg].max()).values
    mcols = ['alt_med', 'slope_med', 'opac_med', 'sfmc_mean', 'sfmc_std']
    reg_asc = reg & np.isfinite(t.asc_L1.values)
    pool_asc = rest_buf & np.isfinite(t.asc_L1.values)
    matched = matched_control(t, reg_asc, pool_asc, mcols, None)
    print(f'\nREGION: {desc}')
    mt = t.iloc[np.unique(matched)]
    print(f'  matched-control tiles: lon {mt.com_lon.min():.1f}..{mt.com_lon.max():.1f}, lat {mt.com_lat.min():.1f}..'
          f'{mt.com_lat.max():.1f}; median alt {mt.alt_med.median():.0f} m (region {t.alt_med[reg].median():.0f}), '
          f'lc_grass {mt.lc_grass.mean():.2f} (region {t.lc_grass[reg].mean():.2f})')
    print('  matched tiles (lon, lat, times used):', [(round(x, 2), round(y, 2), int(c)) for x, y, c in
          zip(mt.com_lon, mt.com_lat, np.bincount(matched)[np.unique(matched)])])
    print(f'  {reg.sum()} tiles ({reg_asc.sum()} with ASCAT N>=50), lon {t.com_lon[reg].min():.2f}..{t.com_lon[reg].max():.2f}, '
          f'lat {t.com_lat[reg].min():.2f}..{t.com_lat[reg].max():.2f}; rest of L1 area {rest.sum()} tiles '
          f'({(rest & np.isfinite(t.asc_L1.values)).sum()} with ASCAT); rest minus 1-tile buffer {rest_buf.sum()}; '
          f'same-lat band {lat_band.sum()}; matched control {len(matched)} draws, {len(np.unique(matched))} unique tiles')
    t['matched_count'] = np.bincount(matched, minlength=len(t))
    t['ctl_latband'] = lat_band
    t.to_csv(out('tile_table') + '.csv', index=False)

    B = Boot(blocks, a.nboot, 1)
    ctl_sets = {'rest': np.where(rest)[0], 'rest_buffer': np.where(rest_buf)[0], 'lat_band': np.where(lat_band)[0],
                'matched': matched}
    ridx = np.where(reg)[0]

    # ---- Q1: % change by arm and species group, region vs controls
    rows = []
    for per in periods:
        for tag, lab in ARMS.items():
            da, ol = P[(per, tag)]
            for g, idx in GROUPS.items():
                v = pct(da, ol, idx, a.nmin)
                for cn, cidx in ctl_sets.items():
                    r = B.diff(v, ridx, cidx)
                    r.update(period=per, arm=lab, group=g, control=cn)
                    if cn == 'rest':
                        r['pooled_reg'] = pooled_pct(da, ol, idx, reg)
                        r['pooled_rest'] = pooled_pct(da, ol, idx, rest)
                        r['mean_L1area'] = np.nanmean(v[L1])
                    rows.append(r)
    q1 = pd.DataFrame(rows)
    q1.to_csv(out('q1_pct_change_region_vs_control') + '.csv', index=False)
    print('\nQ1: O-F stdv % change (DA-OL)/OL, mean of per-tile values, region vs rest of L1 area (95% block-bootstrap CI)')
    show = q1[q1.control == 'rest'][['period', 'arm', 'group', 'n_reg', 'mean_reg', 'reg_lo', 'reg_hi', 'n_ctl', 'mean_ctl',
                                     'diff', 'diff_lo', 'diff_hi', 'pooled_reg', 'pooled_rest']]
    print(show.to_string(index=False, float_format=lambda x: f'{x:7.2f}'))
    print('\nQ1b: ASCAT and SMAP Tb, region vs the other controls (2020-2022)')
    sh = q1[(q1.period == '2020-2022') & (q1.group.isin(['ASCAT', 'SMAP Tb', 'CYGNSS L3'])) & (q1.control != 'rest')]
    print(sh[['arm', 'group', 'control', 'n_reg', 'mean_reg', 'n_ctl', 'mean_ctl', 'diff', 'diff_lo', 'diff_hi']]
          .to_string(index=False, float_format=lambda x: f'{x:7.2f}'))

    rows = []
    for tag, lab in ARMS.items():
        da, ol = P[('2020-2022', tag)]
        for g in ['ASCAT', 'SMAP Tb', 'CYGNSS L3']:
            dc = decompose(da, ol, GROUPS[g], a.nmin)
            for k, v in dc.items():
                if g == 'ASCAT':
                    t[f'{k}_{lab}_ascat'] = v
                r = B.diff(v, ridx, ctl_sets['rest'])
                r.update(arm=lab, group=g, part=k)
                rows.append(r)
    q1d = pd.DataFrame(rows)
    q1d.to_csv(out('q1_decomposition_amplitude_vs_correlation') + '.csv', index=False)
    print('\nQ1d: O-F stdv % change split into forecast-amplitude and correlation parts (2020-2022; sO, sF, r per tile/species):')
    print('     amp_only = DA forecast stdv with OL corr(O,F); corr_only = OL forecast stdv with DA corr; dFstd_pct = % change of forecast stdv')
    print(q1d[['arm', 'group', 'part', 'n_reg', 'mean_reg', 'n_ctl', 'mean_ctl', 'diff', 'diff_lo', 'diff_hi']]
          .to_string(index=False, float_format=lambda x: f'{x:7.3f}'))

    # ---- Q2: operator / attribute diagnostics
    diag = ['N', 'n_l1_obs', 'rs_H', 'r_sc', 'rs_smterm', 'rs_geom', 'pr_geom_inc', 'frac_weakH', 'frac_out15',
            'bias_raw', 'dyn_ratio', 'R_impl', 'snr', 'inn_std', 'inc_slope', 'inc_med', 'slope_med', 'alt_med',
            'opac_med', 'sfmc_mean', 'sfmc_std', 'coh_tail_frac', 'coh_lo_frac', 'coh_hi_frac', 'snr_db_med',
            'n_eff_med', 'w_owner_med', 'frac_edge_fp', 'frac_frozen_soil', 'snow_day_frac', 'snow_day_frac_cold',
            'frozen_soil_day_frac', 'frozen_soil_day_frac_cold', 'frozen_surf_day_frac', 'lc_shrub', 'lc_grass',
            'lc_woody', 'lc_barren', 'lc_crop', 'r_SMAPh', 'r_L3',
            'r_full_ASCAT', 'r_coh_ASCAT', 'r_full_SMAPh', 'r_full_L3', 'r_ASCAT_L3', 'r_ASCAT_SMAPh', 'r_L3_SMAPh',
            'N_ascat', 'olr_ascat', 'olFO_ascat', 'olFO_smap', 'olFO_l3', 'sd_smterm', 'sd_geom', 'olomf_ascat', 'olO_std_ascat', 'olF_std_ascat', 'olr_smap', 'olr_l3',
            'dr_L1_ascat', 'dFstd_L1_ascat', 'dOmFmean_L1_ascat', 'dr_L3_ascat', 'dFstd_L3_ascat',
            'dr_SMAP_ascat', 'dFstd_SMAP_ascat', 'd_cov_km', 'd_sel_km', 'd_orbit_km', 'd_l1edge_km', 'd_box_km']
    rows = []
    for c in diag:
        v = t[c].values.astype(float)
        for cn in ['rest', 'matched']:
            r = B.diff(v, ridx, ctl_sets[cn])
            r.update(metric=c, control=cn,
                     median_reg=np.nanmedian(v[ridx]), median_ctl=np.nanmedian(v[ctl_sets[cn]]))
            rows.append(r)
    q2 = pd.DataFrame(rows)
    q2.to_csv(out('q2_operator_diagnostics_region_vs_control') + '.csv', index=False)
    print('\nQ2: tile diagnostics, region vs rest of L1 area (mean, 95% block-bootstrap CI of difference)')
    print(q2[q2.control == 'rest'][['metric', 'n_reg', 'mean_reg', 'median_reg', 'n_ctl', 'mean_ctl', 'median_ctl',
                                    'diff', 'diff_lo', 'diff_hi']].to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    lc = t[L1].groupby('region').landcover.value_counts(normalize=True).unstack(0).round(3)
    print('\nland-cover mode (MODIS IGBP) fraction of tiles, region (True) vs rest (False):')
    print(lc.to_string())

    attrs = ['alt_med', 'slope_med', 'lc_grass', 'lc_shrub', 'sfmc_mean', 'sfmc_std', 'opac_med', 'snow_day_frac',
             'frozen_soil_day_frac', 'N_ascat', 'olr_ascat', 'olFO_smap', 'sd_smterm', 'olr_smap', 'olr_l3', 'n_l1_obs', 'rs_H', 'rs_smterm',
             'rs_geom', 'dyn_ratio', 'bias_raw', 'R_impl', 'snr', 'frac_weakH', 'coh_tail_frac', 'r_SMAPh', 'r_L3',
             'r_full_ASCAT', 'r_ASCAT_L3', 'r_L3_SMAPh', 'dFstd_L1_ascat', 'd_cov_km', 'd_sel_km', 'com_lat', 'com_lon']
    rows = []
    rng = np.random.default_rng(4)
    for y in ['pct_L1_ASCAT', 'pct_SMAP_ASCAT', 'pct_L3_ASCAT']:
        for at in attrs:
            for nm, m in [('L1 area', L1), ('rest', rest), ('region', reg)]:
                ok = m & np.isfinite(t[y].values) & np.isfinite(t[at].values.astype(float))
                yy, xx, bb = t[y].values[ok], t[at].values[ok].astype(float), blocks[ok]
                rs = stats.spearmanr(yy, xx)[0] if ok.sum() > 10 else np.nan
                lo = hi = np.nan
                if nm != 'region' and ok.sum() > 10:
                    ub, inv = np.unique(bb, return_inverse=True)
                    mem = [np.where(inv == b)[0] for b in range(len(ub))]
                    bs = []
                    for _ in range(min(a.nboot, 500)):
                        pk = np.concatenate([mem[k] for k in rng.integers(0, len(ub), len(ub))])
                        bs.append(stats.spearmanr(yy[pk], xx[pk])[0])
                    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
                rows.append(dict(y=y, attr=at, area=nm, n=ok.sum(), spearman=rs, lo=lo, hi=hi))
    q2a = pd.DataFrame(rows)
    q2a.to_csv(out('q2_attribute_spearman') + '.csv', index=False)
    print('\nQ2a: Spearman of per-tile ASCAT % change with tile attributes (95% block-bootstrap CI, <=500 draws)')
    w = q2a[q2a.y == 'pct_L1_ASCAT'].pivot_table(index='attr', columns='area', values=['spearman', 'lo', 'hi'])
    print(w.round(2).to_string())
    w = q2a[q2a.y == 'pct_SMAP_ASCAT'].pivot_table(index='attr', columns='area', values='spearman')
    print('  same for the SMAP arm (spearman):')
    print(w.round(2).to_string())

    # pooled within-tile innovation correlations (obs level), tile-demeaned
    d = d.merge(t[['ig', 'jg', 'region', 'l1', 'block']], on=['ig', 'jg'], how='left')
    for c in ['inn_full', 'inn_coh', 'inn_ASCAT', 'inn_SMAPh', 'inn_L3']:
        d[c + '_w'] = d[c] - d.groupby(['ig', 'jg'])[c].transform('mean')
    rows = []
    inreg = (d.region == True).values
    inrest = ((d.region == False) & (d.l1 == True)).values
    for x, y in [('inn_full', 'inn_ASCAT'), ('inn_coh', 'inn_ASCAT'), ('inn_full', 'inn_SMAPh'), ('inn_full', 'inn_L3'),
                 ('inn_ASCAT', 'inn_L3'), ('inn_ASCAT', 'inn_SMAPh'), ('inn_L3', 'inn_SMAPh')]:
        for nm, m in [('region', inreg), ('rest', inrest)]:
            for cold, mm in [('all', np.ones(len(d), bool)), ('cold Nov-Mar', d.month.isin(COLD).values),
                             ('warm Apr-Oct', ~d.month.isin(COLD).values)]:
                r, lo, hi, n = boot_corr(d[x + '_w'].values, d[y + '_w'].values, d.block.values, m & mm, a.nboot, 2)
                rows.append(dict(x=x, y=y, area=nm, months=cold, r=r, lo=lo, hi=hi, n=n))
    q2c = pd.DataFrame(rows)
    q2c.to_csv(out('q2_innovation_correlations_pooled') + '.csv', index=False)
    print('\nQ2c: pooled within-tile correlation of innovations at L1 obs times (95% block-bootstrap CI)')
    print('     expected signs if both see the same SM anomaly: L1~ASCAT +, L1~L3 +, L1~SMAPh -, ASCAT~L3 +, ASCAT~SMAPh -')
    print(q2c.to_string(index=False, float_format=lambda x: f'{x:7.3f}'))

    # ---- Q3: season
    rows = []
    mon = np.array([m.month for m in MONTHS])
    seas_sets = {s: np.array([SEASON[m] == s for m in mon]) for s in ['DJF', 'MAM', 'JJA', 'SON']}
    seas_sets['cold Nov-Mar'] = np.isin(mon, list(COLD))
    seas_sets['warm Apr-Oct'] = ~np.isin(mon, list(COLD))
    for sn, ms in seas_sets.items():
        for tag, lab in ARMS.items():
            da, ol = pool(sums[('DA_paired', tag)], ms), pool(sums[('OL_paired_monitor_xmask', tag)], ms)
            for g in ['ASCAT', 'SMAP Tb', 'CYGNSS L3', 'CYGNSS L1']:
                v = pct(da, ol, GROUPS[g], a.nmin_season)
                r = B.diff(v, ridx, ctl_sets['rest'])
                r.update(season=sn, arm=lab, group=g, pooled_reg=pooled_pct(da, ol, GROUPS[g], reg),
                         pooled_rest=pooled_pct(da, ol, GROUPS[g], rest))
                if g == 'ASCAT' and tag == PRIMARY:
                    r['snow_reg'], r['snow_rest'] = snow[ms][:, reg].mean(), snow[ms][:, rest].mean()
                    r['fsoil_reg'], r['fsoil_rest'] = fsoil[ms][:, reg].mean(), fsoil[ms][:, rest].mean()
                    n = ol['N'][:, GROUPS['ASCAT']].sum(1)
                    r['ascat_N_per_tile_reg'] = n[reg].mean() / ms.sum()
                    r['ascat_N_per_tile_rest'] = n[rest].mean() / ms.sum()
                rows.append(r)
    q3 = pd.DataFrame(rows)
    q3.to_csv(out('q3_season_region_vs_rest') + '.csv', index=False)
    print(f'\nQ3: seasonal % change (per-tile N>={a.nmin_season}), region vs rest; pooled = region-pooled within-tile stdv')
    print(q3[q3.group.isin(['ASCAT', 'SMAP Tb'])][['season', 'arm', 'group', 'n_reg', 'mean_reg', 'n_ctl', 'mean_ctl', 'diff',
                                                   'diff_lo', 'diff_hi', 'pooled_reg', 'pooled_rest',
                                                   'ascat_N_per_tile_reg', 'ascat_N_per_tile_rest', 'snow_reg', 'snow_rest',
                                                   'fsoil_reg', 'fsoil_rest']]
          .to_string(index=False, float_format=lambda x: f'{x:7.2f}'))
    # monthly pooled series + frozen fraction + ASCAT N
    rows = []
    for k, m in enumerate(MONTHS):
        ms = np.zeros(len(MONTHS), bool)
        ms[k] = True
        r = dict(month=str(m))
        for tag, lab in ARMS.items():
            da, ol = pool(sums[('DA_paired', tag)], ms), pool(sums[('OL_paired_monitor_xmask', tag)], ms)
            for g in ['ASCAT', 'SMAP Tb', 'CYGNSS L3']:
                r[f'{lab}|{g}|reg'] = pooled_pct(da, ol, GROUPS[g], reg)
                r[f'{lab}|{g}|rest'] = pooled_pct(da, ol, GROUPS[g], rest)
            if tag == PRIMARY:
                n = ol['N'][:, GROUPS['ASCAT']].sum(1)
                r['ascat_N_reg'], r['ascat_N_rest'] = n[reg].mean(), n[rest].mean()
                r['ascat_tilefrac_reg'] = (n[reg] > 0).mean()
                r['ascat_tilefrac_rest'] = (n[rest] > 0).mean()
        r['snow_reg'], r['snow_rest'] = snow[k][reg].mean(), snow[k][rest].mean()
        r['fsoil_reg'], r['fsoil_rest'] = fsoil[k][reg].mean(), fsoil[k][rest].mean()
        dm = d[d.ym == m]
        r['frozen_soil_reg'] = (dm.TSOIL1_FCST[dm.region == True] < FREEZE).mean()
        r['frozen_soil_rest'] = (dm.TSOIL1_FCST[(dm.region == False) & (dm.l1 == True)] < FREEZE).mean()
        r['frozen_surf_reg'] = (dm.TSURF_FCST[dm.region == True] < FREEZE).mean()
        r['L1_obs_reg'] = int((dm.region == True).sum())
        rows.append(r)
    q3m = pd.DataFrame(rows)
    q3m.to_csv(out('q3_monthly_pooled_region_vs_rest') + '.csv', index=False)
    print('\nQ3m: monthly region-pooled within-tile ASCAT/SMAP Tb O-F stdv % change, ASCAT N, frozen fraction at L1 obs')
    cols = ['month', 'L1|ASCAT|reg', 'L1|ASCAT|rest', 'L3|ASCAT|reg', 'L3|ASCAT|rest', 'L1|SMAP Tb|reg', 'L1|SMAP Tb|rest',
            'ascat_N_reg', 'ascat_N_rest', 'ascat_tilefrac_reg', 'snow_reg', 'snow_rest', 'fsoil_reg', 'fsoil_rest', 'L1_obs_reg']
    print(q3m[cols].to_string(index=False, float_format=lambda x: f'{x:7.2f}'))
    cm = q3m.assign(cal=[int(s[-2:]) for s in q3m.month]).groupby('cal')[cols[1:]].mean()
    print('\ncalendar-month means of the monthly values:')
    print(cm.to_string(float_format=lambda x: f'{x:7.2f}'))
    # correlation of monthly region penalty with frozen fraction / ASCAT N
    x = q3m['L1|ASCAT|reg'] - q3m['L1|ASCAT|rest']
    for c in ['snow_reg', 'fsoil_reg', 'ascat_N_reg', 'ascat_tilefrac_reg', 'L1_obs_reg']:
        ok = x.notna() & q3m[c].notna()
        rs = stats.spearmanr(x[ok], q3m[c][ok])
        print(f'  Spearman(monthly region-minus-rest L1 ASCAT %, {c}) = {rs[0]:+.2f} (p={rs[1]:.3f}, N={ok.sum()} months)')

    # ---- Q4: distance to the coverage / preprocessing boundary
    rows = []
    v = t.asc_L1.values
    for dist in ['d_cov_km', 'd_sel_km', 'd_orbit_km', 'd_l1edge_km', 'd_box_km']:
        dd = t[dist].values
        for nm, m in [('L1 area', L1), ('region', reg), ('rest', rest)]:
            ok = m & np.isfinite(v) & np.isfinite(dd)
            rs = stats.spearmanr(v[ok], dd[ok])
            rows.append(dict(dist=dist, area=nm, n=ok.sum(), spearman=rs[0], p_naive=rs[1]))
    q4s = pd.DataFrame(rows)
    # regression ASCAT% ~ region + d (per 100 km), block bootstrap over all L1-area blocks
    rng = np.random.default_rng(3)
    ok = L1 & np.isfinite(v)
    X = lambda dist: np.column_stack([np.ones(ok.sum()), reg[ok], t[dist].values[ok] / 100.])
    ub, inv = np.unique(blocks[ok], return_inverse=True)
    members = [np.where(inv == b)[0] for b in range(len(ub))]
    reg_rows = []
    for dist in ['d_cov_km', 'd_l1edge_km', 'd_sel_km', 'd_orbit_km']:
        Xm, y = X(dist), v[ok]
        beta = np.linalg.lstsq(Xm, y, rcond=None)[0]
        bb = []
        for _ in range(a.nboot):
            pick = np.concatenate([members[k] for k in rng.integers(0, len(ub), len(ub))])
            bb.append(np.linalg.lstsq(Xm[pick], y[pick], rcond=None)[0])
        bb = np.array(bb)
        # same within the rest of the L1 area only (does distance matter away from Four Corners?)
        okr = ~reg[ok]
        b_rest = np.linalg.lstsq(Xm[okr][:, [0, 2]], y[okr], rcond=None)[0]
        bbr = []
        for _ in range(a.nboot):
            pick = np.concatenate([members[k] for k in rng.integers(0, len(ub), len(ub))])
            pick = pick[okr[pick]]
            bbr.append(np.linalg.lstsq(Xm[pick][:, [0, 2]], y[pick], rcond=None)[0][1])
        reg_rows.append(dict(dist=dist, b_region=beta[1], b_region_lo=np.percentile(bb[:, 1], 2.5),
                             b_region_hi=np.percentile(bb[:, 1], 97.5), b_dist_per100km=beta[2],
                             b_dist_lo=np.percentile(bb[:, 2], 2.5), b_dist_hi=np.percentile(bb[:, 2], 97.5),
                             b_dist_rest_only=b_rest[1], b_dist_rest_lo=np.percentile(bbr, 2.5),
                             b_dist_rest_hi=np.percentile(bbr, 97.5), n=ok.sum()))
    q4r = pd.DataFrame(reg_rows)
    # binned
    bins = [-1e9, 0, 75, 150, 250, 1e9]
    labs = ['outside', '0-75', '75-150', '150-250', '>250']
    t['d_cov_bin'] = pd.cut(t.d_cov_km, bins, labels=labs)
    t['d_l1edge_bin'] = pd.cut(t.d_l1edge_km, [-1e9, -100, 0, 75, 150, 250, 1e9],
                               labels=['outside>100', 'outside0-100', '0-75', '75-150', '150-250', '>250'])
    t['area'] = np.where(reg, 'region', np.where(L1, 'rest', 'non-L1'))
    rows = []
    for bcol in ['d_cov_bin', 'd_l1edge_bin']:
        for (ar, bn), g in t.groupby(['area', bcol], observed=True):
            gv = g[['pct_L1_ASCAT', 'pct_L1 err39_ASCAT', 'pct_L3_ASCAT', 'pct_SMAP_ASCAT', 'pct_L1_SMAP Tb', 'n_l1_obs']]
            rows.append(dict(binning=bcol, area=ar, bin=bn, n_tiles=len(g), n_ascat=int(g.pct_L1_ASCAT.notna().sum()),
                             **gv.mean().round(2).to_dict()))
    q4b = pd.DataFrame(rows)
    # tiles just north of the region outside the L1 area (spill-over check)
    north = (~L1) & t.com_lon.between(t.com_lon[reg].min() - 0.4, t.com_lon[reg].max() + 0.4).values & \
        (t.com_lat.values > t.com_lat[reg].max())
    q4s.to_csv(out('q4_distance_spearman') + '.csv', index=False)
    q4r.to_csv(out('q4_distance_regression') + '.csv', index=False)
    q4b.to_csv(out('q4_distance_bins') + '.csv', index=False)
    print('\nQ4: distance to boundary vs per-tile L1-arm ASCAT % change')
    print(q4s.to_string(index=False, float_format=lambda x: f'{x:7.3f}'))
    print('\nregression ASCAT% = a + b_region*region + b_dist*d/100km (L1 area), and b_dist within the rest only; block-bootstrap CIs')
    print(q4r.to_string(index=False, float_format=lambda x: f'{x:7.2f}'))
    print('\nbinned means (km):')
    print(q4b.to_string(index=False))
    nb = t[north]
    print(f'\nnon-L1 tiles directly north of the region (lat > {t.com_lat[reg].max():.2f}): {north.sum()} tiles, '
          f'{nb.pct_L1_ASCAT.notna().sum()} with ASCAT; L1-arm ASCAT % mean {nb.pct_L1_ASCAT.mean():+.2f}, '
          f'L3 {nb.pct_L3_ASCAT.mean():+.2f}, SMAP {nb.pct_SMAP_ASCAT.mean():+.2f}; L1 SMAP Tb {nb["pct_L1_SMAP Tb"].mean():+.2f}')
    nb2 = t[north & (t.com_lat.values < t.com_lat[reg].max() + 0.8)]
    print(f'  first 2 rows north ({len(nb2)} tiles): L1-arm ASCAT {nb2.pct_L1_ASCAT.mean():+.2f}, L3 {nb2.pct_L3_ASCAT.mean():+.2f}')

    # ---- ISMN
    st = pd.read_csv(ISMN)
    st = st.merge(t[['tile_axis_index', 'region', 'l1', 'n_l1_obs']], left_on='tile_index', right_on='tile_axis_index')
    sloc = st.drop_duplicates(['station_key', 'domain'])
    print('\nISMN stations (unique, by domain): in region / in L1 area / all')
    for dom in ['surface', 'rz']:
        s = sloc[sloc.domain == dom]
        print(f'  {dom}: {int(s.region.sum())} stations on {s[s.region].tile_index.nunique()} tiles / '
              f'{int(s.l1.sum())} / {len(s)}')
    wide = st[st.region].pivot_table(index=['station_key', 'domain', 'network'], columns='run', values='R')
    if len(wide):
        wide = wide.assign(**{f'dR_{c}': wide[c] - wide['OL'] for c in wide.columns if c != 'OL'})
        print(wide.round(3).to_string())
        wide.to_csv(out('ismn_region_stations') + '.csv')

    # ---- one summary table: % change region vs rest (2020-2022 and 2021-2022) + amplitude/correlation split
    sm = q1[q1.control == 'rest'][['period', 'arm', 'group', 'n_reg', 'mean_reg', 'reg_lo', 'reg_hi', 'n_ctl', 'mean_ctl',
                                   'ctl_lo', 'ctl_hi', 'diff', 'diff_lo', 'diff_hi', 'pooled_reg', 'pooled_rest']]
    mc = q1[q1.control == 'matched'][['period', 'arm', 'group', 'mean_ctl', 'diff', 'diff_lo', 'diff_hi']]
    mc.columns = ['period', 'arm', 'group', 'mean_matched', 'diff_matched', 'diff_matched_lo', 'diff_matched_hi']
    sm = sm.merge(mc, on=['period', 'arm', 'group'])
    dc = q1d.pivot_table(index=['arm', 'group'], columns='part', values=['mean_reg', 'mean_ctl']).reset_index()
    dc.columns = ['arm', 'group'] + [f'{p}_{"reg" if v == "mean_reg" else "rest"}' for v, p in dc.columns[2:]]
    dc['period'] = '2020-2022'
    sm = sm.merge(dc, on=['period', 'arm', 'group'], how='left')
    sm.to_csv(out('summary_table') + '.csv', index=False, float_format='%.3f')

    # ---- figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    stn = sloc[(sloc.domain == 'surface') & sloc.l1]
    fig = plt.figure(figsize=(13.5, 9.4))
    gs = fig.add_gridspec(2, 3, left=0.03, right=0.98, top=0.95, bottom=0.10, wspace=0.14, hspace=0.22)
    ctlm = np.unique(matched)
    map_panel(fig, gs[0, 0], t, t.pct_L1_ASCAT.values, '(a) L1 DA: ASCAT O-F stdv change', -10, 10, 'RdBu_r', reg, ctlm,
              sel_poly, stn, '% (DA-OL)/OL; negative = better')
    map_panel(fig, gs[0, 1], t, t.pct_L3_ASCAT.values, '(b) L3 DA: ASCAT O-F stdv change', -10, 10, 'RdBu_r', reg, None,
              sel_poly, None, '% (DA-OL)/OL')
    map_panel(fig, gs[0, 2], t, t.pct_SMAP_ASCAT.values, '(c) SMAP DA: ASCAT O-F stdv change', -40, 40, 'RdBu_r', reg,
              None, sel_poly, None, '% (DA-OL)/OL')
    map_panel(fig, gs[1, 0], t, t['pct_L1_SMAP Tb'].values, '(d) L1 DA: SMAP Tb O-F stdv change', -5, 5, 'RdBu_r', reg,
              None, sel_poly, None, '% (DA-OL)/OL')
    map_panel(fig, gs[1, 1], t, t.dFstd_L1_ascat.values, '(e) L1 DA: change of forecast SM stdv at ASCAT obs', -10, 10,
              'PuOr_r', reg, None, sel_poly, None, '% (DA-OL)/OL of F stdv')
    sg = gs[1, 2].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax = fig.add_subplot(sg[0])
    ax2 = fig.add_subplot(sg[1], sharex=ax)
    tm_ = pd.PeriodIndex(q3m.month, freq='M').to_timestamp() + pd.Timedelta(days=14)
    cols_ = {'L1': '#D55E00', 'L3': '#0072B2'}
    for lab, c in cols_.items():
        ax.plot(tm_, q3m[f'{lab}|ASCAT|reg'], '-o', ms=3, lw=1.5, color=c, label=f'{lab} DA, region')
        ax.plot(tm_, q3m[f'{lab}|ASCAT|rest'], '--', lw=1.2, color=c, label=f'{lab} DA, rest of L1 area')
    ax.axhline(0, color='0.5', lw=0.7)
    for k, mt in enumerate(tm_):
        if MONTHS[k].month in COLD:
            ax.axvspan(mt - pd.Timedelta(days=15), mt + pd.Timedelta(days=15), color='0.92', lw=0, zorder=0)
    ax.set_ylabel('ASCAT O-F stdv change (%)')
    ax.set_title('(f) Monthly ASCAT change, region-pooled (grey = Nov-Mar)', fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc='upper left')
    plt.setp(ax.get_xticklabels(), visible=False)
    ax2.bar(tm_, q3m.ascat_N_reg, width=25, color='0.55')
    ax2.set_ylabel('ASCAT obs\nper tile, region', fontsize=7)
    ax2.tick_params(labelsize=7)
    import matplotlib.dates as mdates
    ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.setp(ax2.get_xticklabels(), rotation=30, ha='right')
    fig.text(0.01, 0.005, f'Black outline: penalty region ({reg.sum()} tiles; {desc.split(";")[0]}). Open squares (a): '
             f'attribute-matched control tiles. Triangles (a): ISMN surface stations on L1-area tiles. Purple dashed: Arizona + 200 km '
             f'preprocessing boundary; dotted: {ORBIT_LAT} N orbit limit. Grey: tiles outside the L1 area.',
             fontsize=7, wrap=True)
    fig.savefig(out('four_corners_ascat_map') + '.png', dpi=150, facecolor='white')
    plt.close(fig)
    print(f'\nwrote {a.out_dir}')


if __name__ == '__main__':
    main()
