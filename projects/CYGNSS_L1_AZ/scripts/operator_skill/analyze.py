#!/usr/bin/env python3
"""
Where does the CYGNSS L1 operator do well, and where not?  (steps 2-3 of the plan)

Reads the per-obs OL table from extract_obs_table.py and writes:
  tile_metrics.csv     one row per owner tile (N >= --nmin):
      r_sc        Pearson(scaled O, H(x))            -> does H(x) track O over time
      rs_H        Spearman(raw O, H(x))              \
      rs_sm_own   Spearman(raw O, owner-tile SFMC)    > operator vs operator-free baselines
      rs_sm_fp    Spearman(raw O, footprint SFMC)    /   (same obs, rank-based, so the dB/SFMC
                                                          nonlinearity does not matter)
      rs_geom     Spearman(raw O, geometry term 10log10 sum C_t)        H = geometry term
      rs_smterm   Spearman(raw O, soil-moisture term H - geometry term)     + soil-moisture term
      rs_inc, pr_geom_inc  Spearman(O, incidence) and rank corr(O, geometry term | incidence):
                  does the terrain/geometry part add information beyond incidence angle
      sd_geom, sd_smterm  within-tile std of each term [dB]
      bias_raw    mean(O - H(x)) unscaled [dB]
      dyn_ratio   median(m_std / o_std) from the clim: <1 means O varies more than H(x)
      inn_std     std of scaled O - H(x)
      R_impl      var(scaled O - H(x)) - mean(fcstvar) : implied obs+operator error variance
      snr         mean(fcstvar) / R_impl
      r_SMAPh, r_L3  corr(L1 innovation, SMAP Tbh / L3 innovation): soil-moisture signal in O
      geometry/terrain/vegetation medians, land-cover mode, coherency-tail fraction
      inc_slope   within-tile OLS slope of scaled innovation on incidence angle [dB/deg] (+ t)
  binned_<var>.csv     per-obs stats in quantile bins of each stratifying variable
  tile_attr_spearman.csv   rank correlation between tile skill metrics and tile attributes
  maps_*.png, binned_*.png
Usage:
  analyze.py --start 202001 --end 202212
"""
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

TILE_FILE = ('/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_AZ_fixedop/output/SMAP_EASEv2_M36_GLOBAL/'
             'cat/ens_avg/Y2020/M06/OLv8_M36_AZ_fixedop.inst3_1d_lndfcstana_Nt.20200601_0000z.nc4')
HERE = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/operator_skill/'
COH_LO, COH_HI = 0.40, 2.16
SEASON = {12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
          6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'}
BIN_VARS = ['sp_inc_angle', 'geom_db', 'n_eff', 'w_max', 'w_owner', 'sp_nearest_tile_distance_km', 'srtm_slope',
            'coefficient_weighted_opacity', 'SFMC_FCST', 'fcst', 'coherency_ratio', 'ddm_snr',
            'sp_rx_gain', 'sc_num', 'month']
TILE_ATTRS = ['inc_med', 'n_eff_med', 'w_max_med', 'w_owner_med', 'slope_med', 'alt_med', 'opac_med',
              'sfmc_mean', 'sfmc_std', 'coh_tail_frac', 'snr_db_med', 'lat', 'lon']
SKILL = ['r_sc', 'rs_H', 'rs_geom', 'pr_geom_inc', 'rs_smterm', 'sd_geom', 'sd_smterm', 'rs_H_minus_sm_fp', 'bias_raw', 'dyn_ratio', 'inn_std', 'R_impl', 'snr',
         'r_SMAPh', 'r_L3', 'inc_slope']


def load(a):
    months = pd.period_range(pd.Period(a.start, 'M'), pd.Period(a.end, 'M'), freq='M')
    fs = [f'{HERE}obs_table/ol_l1_obs_{m.strftime("%Y%m")}.parquet' for m in months]
    missing = [f for f in fs if not os.path.exists(f)]
    if missing:
        raise SystemExit(f'missing {len(missing)} month files, e.g. {missing[0]}')
    d = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    d['month'] = pd.to_datetime(d.time).dt.month
    d['season'] = d.month.map(SEASON)
    d['dyn'] = d.m_std_full / d.o_std_full
    d['coh_tail'] = ~d.coherency_ratio.between(COH_LO, COH_HI)
    d['tile'] = list(zip(d.sp_nearest_tile_ig, d.sp_nearest_tile_jg))
    return d


def corr(x, y, f=np.corrcoef, nmin=30):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < nmin:
        return np.nan
    return f(x[ok], y[ok])[0, 1] if f is np.corrcoef else f(x[ok], y[ok])[0]


def partial_rank(y, x, z, nmin=50):
    """Rank correlation of y and x after removing a linear rank dependence on z."""
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(z)
    if ok.sum() < nmin:
        return np.nan
    ry, rx, rz = (stats.rankdata(v[ok]) for v in (y, x, z))
    res = lambda v: v - np.polyval(np.polyfit(rz, v, 1), rz)
    return np.corrcoef(res(ry), res(rx))[0, 1]


def tile_metrics(g):
    O, H = g.obs.values, g.fcst.values
    inn = g.inn_full.values
    ok = np.isfinite(inn)
    R = np.nanvar(inn, ddof=1) - g.fcstvar[ok].mean()
    rs_H = corr(O, H, stats.spearmanr)
    rs_fp = corr(O, g.sfmc_fp.values, stats.spearmanr)
    slope, tstat = np.nan, np.nan
    x = g.sp_inc_angle.values
    if ok.sum() >= 50 and np.nanstd(x[ok]) > 0.5:
        lr = stats.linregress(x[ok], inn[ok])
        slope, tstat = lr.slope, lr.slope / lr.stderr
    return pd.Series(dict(
        N=len(g), r_sc=corr(g.obs_sc_full.values, H), rs_H=rs_H,
        rs_sm_own=corr(O, g.SFMC_FCST.values, stats.spearmanr), rs_sm_fp=rs_fp,
        rs_H_minus_sm_fp=rs_H - rs_fp,
        rs_geom=corr(O, g.geom_db.values, stats.spearmanr),
        rs_inc=corr(O, g.sp_inc_angle.values, stats.spearmanr),
        pr_geom_inc=partial_rank(O, g.geom_db.values, g.sp_inc_angle.values),
        rs_smterm=corr(O, H - g.geom_db.values, stats.spearmanr),
        sd_geom=np.nanstd(g.geom_db.values), sd_smterm=np.nanstd(H - g.geom_db.values),
        bias_raw=np.nanmean(O - H), dyn_ratio=np.nanmedian(g.dyn), inn_std=np.nanstd(inn, ddof=1),
        R_impl=R, snr=g.fcstvar[ok].mean() / R if R > 0 else np.nan,
        r_SMAPh=corr(inn, g.inn_SMAPh.values), r_L3=corr(inn, g.inn_L3.values),
        inc_slope=slope, inc_slope_t=tstat,
        inc_med=np.nanmedian(x), n_eff_med=g.n_eff.median(), w_max_med=g.w_max.median(),
        w_owner_med=g.w_owner.median(), slope_med=g.srtm_slope.median(), alt_med=g.srtm_dem_alt.median(),
        opac_med=g.coefficient_weighted_opacity.median(), sfmc_mean=g.SFMC_FCST.mean(),
        sfmc_std=g.SFMC_FCST.std(), coh_tail_frac=g.coh_tail.mean(), snr_db_med=g.ddm_snr_db.median(),
        landcover=g.modis_land_cover.mode().iloc[0] if g.modis_land_cover.notna().any() else np.nan,
        frac_out15=(np.abs(inn) > 15).mean()))


def binned(d, var, nb=8):
    x = d[var]
    if x.nunique() <= 12:
        b = x
    else:
        b = pd.qcut(x, nb, duplicates='drop')
    rows = []
    for k, g in d.groupby(b, observed=True):
        inn = g.inn_full
        rows.append(dict(bin=str(k), N=len(g), x_med=g[var].median(), inn_mean=inn.mean(),
                         inn_std=inn.std(), R_impl=inn.var() - g.fcstvar[inn.notna()].mean(),
                         fcstvar=g.fcstvar.mean(), bias_raw=(g.obs - g.fcst).mean(),
                         r_SMAPh=corr(inn.values, g.inn_SMAPh.values),
                         r_L3=corr(inn.values, g.inn_L3.values),
                         # within-tile: demean innovations by tile so between-tile offsets do not leak in
                         inn_mean_within=g.inn_within.mean()))
    return pd.DataFrame(rows)


def tile_centres():
    from netCDF4 import Dataset
    with Dataset(TILE_FILE) as nc:
        return pd.DataFrame({'ig': np.asarray(nc['IG'][:]), 'jg': np.asarray(nc['JG'][:]),
                             'lon': np.asarray(nc['lon'][:]), 'lat': np.asarray(nc['lat'][:])})


def tile_polys(t):
    # cell edges from the full domain grid, so any subset of tiles gets the same cells
    allt = tile_centres()
    def edges(c):
        c = np.unique(np.round(c, 4))
        mid = 0.5 * (c[:-1] + c[1:])
        return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])
    xc, xe = edges(allt.lon.values)
    yc, ye = edges(allt.lat.values)
    ix, iy = np.searchsorted(xc, np.round(t.lon, 4)), np.searchsorted(yc, np.round(t.lat, 4))
    return [[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])] for i, j in zip(ix, iy)]


def state_borders():
    try:
        import cartopy.io.shapereader as shp
        geoms = shp.Reader(shp.natural_earth('50m', 'cultural', 'admin_1_states_provinces_lakes')).geometries()
        out = []
        for g in geoms:
            for poly in getattr(g, 'geoms', [g]):
                xy = np.asarray(poly.exterior.coords)
                if ((xy[:, 0] > -120) & (xy[:, 0] < -104) & (xy[:, 1] > 27) & (xy[:, 1] < 42)).any():
                    out.append(xy)
        return out
    except Exception as e:      # no Natural Earth data available: maps without borders
        print('no state borders:', e)
        return []


BORDERS = []


def maps(t, cols, out, title):
    n = len(cols)
    fig, axs = plt.subplots(int(np.ceil(n / 3)), 3, figsize=(13, 3.9 * np.ceil(n / 3)), squeeze=False)
    verts = tile_polys(t)
    for ax, (c, (vmin, vmax, cmap)) in zip(axs.flat, cols.items()):
        pc = PolyCollection(verts, array=np.ma.masked_invalid(t[c].values.astype(float)), cmap=cmap, linewidths=0)
        pc.set_clim(vmin, vmax)
        ax.add_collection(pc)
        for b in BORDERS:
            ax.plot(b[:, 0], b[:, 1], color='0.3', lw=0.5)
        ax.set_xlim(-118, -106); ax.set_ylim(29, 40); ax.set_aspect(1.15)
        ax.set_title(c); plt.colorbar(pc, ax=ax, shrink=0.8)
    for ax in axs.flat[n:]:
        ax.axis('off')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--start', required=True)
    p.add_argument('--end', required=True)
    p.add_argument('--nmin', type=int, default=100, help='min obs per tile for tile metrics')
    a = p.parse_args()
    global BORDERS
    BORDERS = state_borders()
    out = f'{HERE}analysis_{a.start}_{a.end}/'
    os.makedirs(out, exist_ok=True)
    d = load(a)
    d['inn_within'] = d.inn_full - d.groupby('tile').inn_full.transform('mean')
    print(f'{len(d)} obs, {d.tile.nunique()} owner tiles, {a.start}-{a.end}')
    o = d[d.inn_full.abs() > 15]
    q = d.fcst.quantile(0.125)
    print(f'|scaled innovation| > 15 dB: {len(o)} obs on {o.tile.nunique()} tiles; {(o.fcst < q).mean():.0%} in the '
          f'lowest-H(x) octile; median H(x) {o.fcst.median():.1f} dB vs O {o.obs.median():.1f} dB '
          f'(all obs: H {d.fcst.median():.1f}, O {d.obs.median():.1f}); median coherency {o.coherency_ratio.median():.2f}')
    print('  by month:', o.groupby('month').size().to_dict())

    # tile metrics, all seasons and by season
    t = d.groupby('tile').filter(lambda g: len(g) >= a.nmin).groupby('tile').apply(tile_metrics).reset_index()
    t[['ig', 'jg']] = pd.DataFrame(t.tile.tolist(), index=t.index)
    t = t.merge(tile_centres(), on=['ig', 'jg'], how='left')
    t.drop(columns='tile').to_csv(f'{out}tile_metrics.csv', index=False)
    seas = []
    for s, ds in d.groupby('season'):
        ts = ds.groupby('tile').filter(lambda g: len(g) >= a.nmin // 4).groupby('tile').apply(tile_metrics)
        ts['season'] = s
        seas.append(ts.reset_index())
    seas = pd.concat(seas)
    seas[['ig', 'jg']] = pd.DataFrame(seas.tile.tolist(), index=seas.index)
    seas.drop(columns='tile').to_csv(f'{out}tile_metrics_by_season.csv', index=False)
    print(f'\n{len(t)} tiles with N >= {a.nmin}; medians:')
    print(t[SKILL + ['rs_sm_own', 'rs_sm_fp']].median().round(3).to_string())
    print('\nby season (tile medians):')
    print(seas.groupby('season')[['N', 'r_sc', 'rs_H', 'rs_geom', 'rs_smterm', 'sd_geom', 'sd_smterm', 'rs_sm_fp', 'R_impl', 'snr', 'r_SMAPh', 'r_L3', 'dyn_ratio']]
          .median().round(3).to_string())

    # which tile attributes go with skill
    sp = pd.DataFrame({s: {at: stats.spearmanr(t[s], t[at], nan_policy='omit')[0] for at in TILE_ATTRS}
                       for s in SKILL}).round(2)
    sp.to_csv(f'{out}tile_attr_spearman.csv')
    print('\nSpearman, tile skill (columns) vs tile attributes (rows):')
    print(sp.to_string())

    # per-obs binned stats
    pd.set_option('display.width', 250)
    for v in BIN_VARS:
        b = binned(d, v)
        b.to_csv(f'{out}binned_{v}.csv', index=False)
        print(f'\n--- {v} ---')
        print(b.drop(columns='bin').to_string(index=False, float_format=lambda x: f'{x:8.3f}'))

    # figures
    maps(t, {'r_sc': (0, 0.8, 'viridis'), 'rs_geom': (0, 0.8, 'viridis'), 'rs_smterm': (-0.2, 0.6, 'RdBu_r'),
             'sd_smterm': (0, 1.5, 'viridis'),
             'dyn_ratio': (0.3, 1.2, 'viridis'), 'bias_raw': (-10, 2, 'RdBu_r'),
             'R_impl': (0, 12, 'magma_r'), 'snr': (0, 0.6, 'viridis'),
             'r_SMAPh': (-0.5, 0.2, 'RdBu'), 'r_L3': (-0.2, 0.5, 'RdBu_r'), 'inc_slope': (-0.1, 0.1, 'RdBu_r')},
         f'{out}maps_tile_skill.png', f'L1 operator skill per owner tile, OL {a.start}-{a.end} (N >= {a.nmin})')
    maps(t, {'n_eff_med': (1, 3, 'viridis'), 'w_owner_med': (0, 1, 'viridis'), 'slope_med': (0, 10, 'viridis'),
             'opac_med': (0, 0.3, 'viridis'), 'sfmc_mean': (0.03, 0.2, 'viridis'), 'coh_tail_frac': (0, 0.5, 'viridis')},
         f'{out}maps_tile_attrs.png', 'Tile attributes')
    fig, axs = plt.subplots(3, 5, figsize=(17, 9.5))
    for ax, v in zip(axs.flat, BIN_VARS):
        b = pd.read_csv(f'{out}binned_{v}.csv')
        x = np.arange(len(b))
        ax.plot(x, b.inn_mean_within, 'o-', label='mean inn (within tile)')
        ax.plot(x, b.R_impl / 10, 's-', label='R_impl / 10')
        ax.plot(x, b.r_L3, '^-', label='r(L1, L3 inn)')
        ax.plot(x, -b.r_SMAPh, 'v-', label='-r(L1, SMAPh inn)')
        ax.set_xticks(x); ax.set_xticklabels([f'{v:.3g}' for v in b.x_med], rotation=60, fontsize=7)
        ax.axhline(0, color='0.6', lw=0.6); ax.set_title(v, fontsize=9)
    axs.flat[0].legend(fontsize=7)
    for ax in axs.flat[len(BIN_VARS):]:
        ax.axis('off')
    fig.tight_layout(); fig.savefig(f'{out}binned_summary.png', dpi=120); plt.close(fig)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
