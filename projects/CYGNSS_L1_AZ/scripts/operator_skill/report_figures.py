#!/usr/bin/env python3
"""
Figures for docs/cygl1_operator_skill_report.md (where the CYGNSS L1 operator does well/badly, and how that links to DA benefit).

Reads only finished outputs:
  output/operator_skill/obs_table/ol_l1_obs_<yyyymm>.parquet        per-obs OL table (extract_obs_table.py)
  output/operator_skill/analysis_<start>_<end>/tile_metrics.csv      per-tile operator skill (analyze.py)
  output/operator_skill/da_link/<task>/                              step-4 DA-link outputs (scripts/operator_skill/da_link/)

Usage:
  report_figures.py --figs all
  report_figures.py --figs op_maps weakh
"""
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature

PROJ = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/'
OPS = PROJ + 'output/operator_skill/'
LINK = OPS + 'da_link/'
OUT = PROJ + 'docs/cygl1_operator_skill_report_figures/'
DOMAIN = (-118.0, -106.0, 29.0, 40.0)
L1_AREA_MIN_OBS = 100
COH_LO, COH_HI = 0.40, 2.16
# Okabe-Ito, as in paper/make_paper_figures.py; the monitors take the colour of the arm that assimilates them.
# SMOS (4th slot) is below the CVD target next to SMAP, so it always also gets its own marker.
MON = {'SMAP': ('#009E73', 'o'), 'L3': ('#0072B2', 's'), 'SMOS': ('#CC79A7', '^')}
DPI = 150

plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8, 'figure.titlesize': 11, 'axes.spines.top': False, 'axes.spines.right': False,
})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(OUT + name, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('wrote', OUT + name)


# ---------------------------------------------------------------------------------------------
# maps
# ---------------------------------------------------------------------------------------------
def map_axes(fig, pos):
    ax = fig.add_subplot(*pos, projection=ccrs.PlateCarree())
    ax.set_extent(DOMAIN, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor='0.35')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.6, edgecolor='0.25')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, edgecolor='0.25')
    return ax


def _edges(c):
    c = np.unique(np.round(np.asarray(c, float), 4))
    mid = 0.5 * (c[:-1] + c[1:])
    return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def _cells(lon, lat):
    allt = pd.read_csv(LINK + 'tile_axis_ig_jg.csv')
    xc, xe = _edges(allt.com_lon)
    yc, ye = _edges(allt.com_lat)
    ix = np.searchsorted(xc, np.round(np.asarray(lon, float), 4))
    iy = np.searchsorted(yc, np.round(np.asarray(lat, float), 4))
    return ix, iy, xe, ye


def tile_squares(ax, lon, lat, values, **kw):
    """Each M36 tile as its cell rectangle; edges from the full 909-tile grid so any subset gets the same cells."""
    ix, iy, xe, ye = _cells(lon, lat)
    verts = [[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])] for i, j in zip(ix, iy)]
    vmin, vmax = kw.pop('vmin', None), kw.pop('vmax', None)
    pc = PolyCollection(verts, transform=ccrs.PlateCarree(), linewidths=0, **kw)
    pc.set_array(np.ma.masked_invalid(np.asarray(values, float)))
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def outline(ax, lon, lat, mask, **kw):
    """Draw the outer boundary of the tiles where mask is True."""
    ix, iy, xe, ye = _cells(lon, lat)
    cells = set(zip(ix[mask], iy[mask]))
    for i, j in cells:
        for di, dj, seg in [(-1, 0, [(xe[i], ye[j]), (xe[i], ye[j + 1])]), (1, 0, [(xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1])]),
                            (0, -1, [(xe[i], ye[j]), (xe[i + 1], ye[j])]), (0, 1, [(xe[i], ye[j + 1]), (xe[i + 1], ye[j + 1])])]:
            if (i + di, j + dj) not in cells:
                ax.plot(*zip(*seg), transform=ccrs.PlateCarree(), **{'color': 'k', 'lw': 1.4, **kw})


def map_grid(t, panels, ncol, name, suptitle, lon='lon', lat='lat', width=4.6, region=None, extra=None):
    """panels: list of (column, title, vmin, vmax, cmap, colorbar label). region: bool mask to outline on every map.
    extra: callable(ax) drawing one more, non-map panel after the maps."""
    nrow = int(np.ceil((len(panels) + (extra is not None)) / ncol))
    fig = plt.figure(figsize=(width * ncol, 3.9 * nrow))
    for k, (c, title, vmin, vmax, cmap, lab) in enumerate(panels):
        ax = map_axes(fig, (nrow, ncol, k + 1))
        v = t[c].values.astype(float)
        pc = tile_squares(ax, t[lon], t[lat], v, cmap=cmap, vmin=vmin, vmax=vmax)
        if region is not None:
            outline(ax, t[lon], t[lat], np.asarray(region))
        ax.set_title(f'({chr(97 + k)}) {title}', fontsize=9, loc='left')
        ax.text(0.02, 0.03, f'tile median {np.nanmedian(v):.2f}  (N={np.isfinite(v).sum()})', transform=ax.transAxes,
                fontsize=7, bbox=dict(fc='white', ec='none', alpha=0.8))
        cb = fig.colorbar(pc, ax=ax, shrink=0.75, extend='both', pad=0.02)
        cb.set_label(lab, fontsize=8)
    if extra is not None:
        ax = fig.add_subplot(nrow, ncol, len(panels) + 1)
        extra(ax)
        title = ax.get_title()
        ax.set_title('')
        ax.set_title(f'({chr(97 + len(panels))}) {title}', fontsize=9, loc='left')
    fig.suptitle(suptitle)
    fig.tight_layout()
    save(fig, name)


def load_tiles(a):
    t = pd.read_csv(f'{OPS}analysis_{a.start}_{a.end}/tile_metrics.csv')
    m = pd.read_csv(LINK + 'tile_axis_ig_jg.csv')
    return t.merge(m[['ig', 'jg', 'n_l1_obs']], on=['ig', 'jg'], how='left')


def fig_op_maps(a):
    t = load_tiles(a)
    t['frac_out15_pct'] = 100 * t.frac_out15
    panels = [
        ('rs_H', 'Spearman(O, H(x))', 0, 0.8, 'Blues', 'rank correlation'),
        ('pr_geom_inc', 'geometry term, incidence removed', 0, 0.8, 'Blues', 'partial rank correlation'),
        ('rs_smterm', 'soil-moisture term: Spearman(O, H - geom)', -0.6, 0.6, 'RdBu_r', 'rank correlation'),
        ('R_impl', 'implied obs-error variance', 0, 20, 'Purples', 'dB$^2$'),
        ('frac_out15_pct', 'weak-H outliers (|innovation| > 15 dB)', 0, 2, 'Oranges', '% of tile obs'),
        ('r_L3', 'corr(L1 innovation, L3 innovation)', -0.5, 0.5, 'RdBu_r', 'Pearson r'),
    ]
    map_grid(t, panels, 3, 'fig01_operator_skill_maps.png',
             f'CYGNSS L1 operator skill per owner tile, open loop {a.start}-{a.end} (tiles with >= 100 obs)')


# ---------------------------------------------------------------------------------------------
# per-obs table helpers
# ---------------------------------------------------------------------------------------------
def load_obs(a, cols):
    months = pd.period_range(pd.Period(a.start, 'M'), pd.Period(a.end, 'M'), freq='M')
    d = pd.concat([pd.read_parquet(f'{OPS}obs_table/ol_l1_obs_{m.strftime("%Y%m")}.parquet', columns=cols)
                   for m in months], ignore_index=True)
    d['tile'] = d.sp_nearest_tile_ig.astype(int) * 10000 + d.sp_nearest_tile_jg.astype(int)
    return d


def cluster_boot_mean(y, x_bin, tile, nboot=300, seed=1):
    """Mean of y per bin with a tile-cluster bootstrap 95% CI (tiles resampled with replacement)."""
    ok = np.isfinite(y)
    df = pd.DataFrame({'y': y[ok], 'b': x_bin[ok], 't': tile[ok]})
    S = df.pivot_table(index='t', columns='b', values='y', aggfunc='sum', fill_value=0, observed=True)
    N = df.pivot_table(index='t', columns='b', values='y', aggfunc='count', fill_value=0, observed=True)
    rng = np.random.default_rng(seed)
    w = rng.multinomial(len(S), np.full(len(S), 1 / len(S)), size=nboot)
    boot = (w @ S.values) / np.maximum(w @ N.values, 1)
    est = S.values.sum(0) / N.values.sum(0)
    return pd.DataFrame({'bin': S.columns, 'mean': est, 'lo': np.percentile(boot, 2.5, 0),
                         'hi': np.percentile(boot, 97.5, 0), 'N': N.values.sum(0)})


def fig_geometry(a):
    cols = ['sp_nearest_tile_ig', 'sp_nearest_tile_jg', 'inn_full', 'sp_inc_angle', 'ddm_snr_db',
            'coefficient_weighted_opacity', 'sc_num', 'srtm_slope']
    d = load_obs(a, cols)
    d['inn_within'] = d.inn_full - d.groupby('tile').inn_full.transform('mean')
    specs = [('sp_inc_angle', 'incidence angle (deg)'), ('ddm_snr_db', 'DDM SNR (dB)'),
             ('coefficient_weighted_opacity', 'vegetation opacity'), ('srtm_slope', 'terrain slope (deg)'),
             ('sc_num', 'CYGNSS spacecraft')]
    fig, axs = plt.subplots(1, len(specs), figsize=(3.1 * len(specs), 3.0), sharey=True)
    rows = []
    for ax, (v, lab) in zip(axs, specs):
        if v == 'sc_num':
            b = d[v].values
        else:
            b = pd.qcut(d[v], 10, labels=False, duplicates='drop').values
        r = cluster_boot_mean(d.inn_within.values, b, d.tile.values)
        r['x'] = d.groupby(b)[v].median().reindex(r.bin).values
        r['var'] = v
        rows.append(r)
        ax.axhline(0, color='0.6', lw=0.8)
        if v == 'sc_num':   # categories: no line, no band
            ax.errorbar(r.x, r['mean'], yerr=[r['mean'] - r.lo, r.hi - r['mean']], fmt='o', color='#D55E00', ms=5,
                        capsize=2)
            ax.set_xticks(r.x)
        else:
            ax.fill_between(r.x, r.lo, r.hi, color='#D55E00', alpha=0.2, lw=0)
            ax.plot(r.x, r['mean'], '-o', color='#D55E00', lw=2, ms=4)
        ax.set_xlabel(lab)
    axs[0].set_ylabel('within-tile mean scaled innovation (dB)')
    fig.suptitle(f'Residual dependence of the scaled L1 innovation O - H(x) after removing each tile\'s mean, '
                 f'{a.start}-{a.end} (deciles; band = 95% tile-cluster bootstrap)', fontsize=10)
    fig.tight_layout()
    save(fig, 'fig02_residual_geometry.png')
    pd.concat(rows).to_csv(OUT + 'fig02_residual_geometry.csv', index=False)


def fig_weakh(a):
    cols = ['sp_nearest_tile_ig', 'sp_nearest_tile_jg', 'obs', 'fcst', 'inn_full', 'coherency_ratio', 'time']
    d = load_obs(a, cols)
    d['weak'] = d.inn_full.abs() > 15
    d['coh_ok'] = d.coherency_ratio.between(COH_LO, COH_HI)
    fig = plt.figure(figsize=(14, 3.9))
    # (a) O vs H(x)
    ax = fig.add_subplot(1, 3, 1)
    ok = np.isfinite(d.obs) & np.isfinite(d.fcst)
    hb = ax.hexbin(d.fcst[ok], d.obs[ok], gridsize=70, extent=(-45, 0, -30, 0), bins='log', cmap='Greys', mincnt=1)
    w = d[d.weak]
    ax.scatter(w.fcst, w.obs, s=3, color='#D55E00', alpha=0.5, lw=0, label=f'|innovation| > 15 dB (N={len(w)})')
    ax.plot([-30, 0], [-30, 0], color='0.4', lw=0.8, ls='--')
    ax.set_xlim(-45, 0); ax.set_ylim(-30, 0)
    ax.set_xlabel('H(x), open loop (dB)'); ax.set_ylabel('observed O (dB)')
    ax.set_title('(a) observed vs simulated L1 (raw dB)', loc='left', fontsize=9)
    ax.legend(loc='lower right', markerscale=4)
    fig.colorbar(hb, ax=ax, shrink=0.8, label='obs per cell')
    # (b) seasonality and the coherency filter
    ax = fig.add_subplot(1, 3, 2)
    w = w.assign(month=pd.to_datetime(w.time).dt.month)
    c = w.groupby(['month', 'coh_ok']).size().unstack(fill_value=0).reindex(range(1, 13), fill_value=0)
    ax.bar(c.index, c.get(False, 0), color='0.7', width=0.8, label='removed by coherency filter')
    ax.bar(c.index, c.get(True, 0), bottom=c.get(False, 0), color='#D55E00', width=0.8,
           label='passes filter (can reach the DA)')
    ax.set_xticks(range(1, 13)); ax.set_xticklabels('JFMAMJJASOND')
    ax.set_xlabel('month'); ax.set_ylabel('weak-H outliers, 2020-22')
    ax.set_title(f'(b) when they occur ({100 * (~w.coh_ok).mean():.0f}% removed by coherency 0.40-2.16)',
                 loc='left', fontsize=9)
    ax.legend()
    # (c) monitor harm after weak-H vs normal L1 analyses in the same tile x calendar-month cells (task D)
    ax = fig.add_subplot(1, 3, 3)
    q = pd.read_csv(LINK + 'obs_level/202001_202212/q3_weakH_context_and_contribution.csv')
    q = q[q.test == 'context-matched']
    cls = ['weak-H cells: weak-H dominant', 'weak-H cells: normal dominant', 'weak-H cells: no L1 in 72h']
    labs = ['after weak-H\nL1 analysis', 'after normal\nL1 analysis', 'no L1 analysis\nin 72 h']
    x = np.arange(len(cls))
    for k, (g, (col, mk)) in enumerate(MON.items()):
        s = q[q.grp == g].set_index('cls').reindex(cls)
        ax.errorbar(x + (k - 1) * 0.22, s.dmse, yerr=[s.dmse - s.dmse_lo, s.dmse_hi - s.dmse], fmt=mk, color=col,
                    ms=6, capsize=2, lw=1.2, label=f'{g} Tb' if g != 'L3' else 'CYGNSS L3')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labs)
    ax.set_ylabel('monitor MSE change vs OL (%)')
    ax.set_title('(c) monitor O-F in the same tile x month cells', loc='left', fontsize=9)
    ax.legend(loc='upper right')
    fig.suptitle('The weak-H failure mode: H(x) predicts almost no coherent return while the observation is normal')
    fig.tight_layout()
    save(fig, 'fig03_weak_h_failure_mode.png')


# ---------------------------------------------------------------------------------------------
# step 4, task D: obs-level DA effect by operator class
# ---------------------------------------------------------------------------------------------
def fig_obs_level(a):
    q = pd.read_csv(LINK + 'obs_level/202001_202212/q3_monitor_by_L1_class.csv')
    q = q[q.rule == 'dominant72']
    specs = [('fcst_ol_tab', 'H(x) octile of the L1 obs', None),
             ('smterm', 'SM-term octile of the L1 obs', None),
             ('tile_smskill', 'tile SM-term skill tercile', ['T1 low', 'T2 mid', 'T3 high']),
             ('period', 'period', ['spring window', 'other', 'Oct-Jan'])]
    fig, axs = plt.subplots(2, len(specs), figsize=(3.4 * len(specs), 5.6), sharey='row')
    for j, (v, lab, order) in enumerate(specs):
        s = q[q['var'] == v]
        for k, (g, (col, mk)) in enumerate(MON.items()):
            sg = s[s.grp == g]
            sg = sg.set_index('bin').reindex(order).reset_index() if order else sg.assign(
                b=pd.to_numeric(sg.bin)).sort_values('b')
            x = np.arange(len(sg)) + (k - 1) * 0.12
            name = f'{g} Tb' if g != 'L3' else 'CYGNSS L3'
            for i, (m, lo, hi) in enumerate([('dmse', 'dmse_lo', 'dmse_hi'), ('alpha', 'alpha_lo', 'alpha_hi')]):
                ax = axs[i, j]
                ax.errorbar(x, sg[m], yerr=[sg[m] - sg[lo], sg[hi] - sg[m]], fmt='-' + mk, color=col, ms=5,
                            capsize=2, lw=1.2, label=name)
        for i in range(2):
            ax = axs[i, j]
            ax.axhline(0 if i == 0 else 1, color='0.6', lw=0.8, ls='-' if i == 0 else '--')
            n = len(order) if order else s.bin.nunique()
            ax.set_xticks(range(n))
            if order:
                ax.set_xticklabels([o.replace(' window', '') for o in order])
            else:
                xm = s[s.grp == 'SMAP'].assign(b=pd.to_numeric(s.bin)).sort_values('b').x_med
                ax.set_xticklabels([f'{x:.1f}' for x in xm], fontsize=7, rotation=45)
        axs[1, j].set_xlabel(lab + ('\n(bin median, dB)' if not order else ''))
        axs[0, j].set_title(f'({chr(97 + j)})', loc='left', fontsize=9)
    axs[0, 0].set_ylabel('monitor MSE change vs OL (%)\nnegative = DA better')
    axs[1, 0].set_ylabel('alpha_opt (1 = increments\nright size, 0 = no signal)')
    axs[1, 0].set_ylim(-0.5, 2.5)
    axs[0, 0].legend(loc='lower left')
    fig.suptitle('Monitor O-F after an L1 analysis, by operator class of the dominant L1 obs in the preceding 72 h '
                 '(primary L1 arm, 2020-22, 95% tile-cluster bootstrap)', fontsize=10)
    fig.tight_layout()
    save(fig, 'fig04_obs_level_monitor_by_operator_class.png')


# ---------------------------------------------------------------------------------------------
# step 4, task C: Four Corners ASCAT penalty
# ---------------------------------------------------------------------------------------------
def fig_four_corners(a):
    d = LINK + 'four_corners_ascat/'
    t = pd.read_csv(d + 'tile_table.csv')
    t = t[t.l1.astype(bool)]           # L1 area
    reg = t.region.astype(bool).values
    q = pd.read_csv(d + 'q2_innovation_correlations_pooled.csv')
    q = q[(q.months == 'all') & (q.x != 'inn_coh')]
    # Tb innovations anti-correlate with SM innovations (wetter -> colder), so Tb pairs are sign-flipped: + = consistent.
    pairs = [('inn_full', 'inn_ASCAT', 'L1 ~ ASCAT', 1), ('inn_full', 'inn_SMAPh', 'L1 ~ SMAP Tb', -1),
             ('inn_full', 'inn_L3', 'L1 ~ L3', 1), ('inn_ASCAT', 'inn_L3', 'ASCAT ~ L3', 1),
             ('inn_ASCAT', 'inn_SMAPh', 'ASCAT ~ SMAP Tb', -1), ('inn_L3', 'inn_SMAPh', 'L3 ~ SMAP Tb', -1)]

    def bars(ax):
        y = np.arange(len(pairs))
        for k, (area, col) in enumerate([('region', '#D55E00'), ('rest', '0.55')]):
            r = [q[(q.x == x) & (q.y == yv) & (q.area == area)].iloc[0] for x, yv, _, _ in pairs]
            v = np.array([s * row.r for row, (_, _, _, s) in zip(r, pairs)])
            lo = np.array([min(s * row.lo, s * row.hi) for row, (_, _, _, s) in zip(r, pairs)])
            hi = np.array([max(s * row.lo, s * row.hi) for row, (_, _, _, s) in zip(r, pairs)])
            ax.barh(y + (0.2 if k == 0 else -0.2), v, height=0.38, color=col, xerr=[v - lo, hi - v], capsize=2,
                    label='penalty region' if area == 'region' else 'rest of L1 area')
        ax.set_yticks(y); ax.set_yticklabels([p[2] for p in pairs]); ax.invert_yaxis()
        ax.axvline(0, color='0.6', lw=0.8)
        ax.set_xlabel('within-tile innovation correlation\n(Tb pairs x -1, so + = consistent)')
        ax.set_title('who agrees with whom (95% CI)')
        ax.legend(loc='lower right')

    panels = [('pct_L1_ASCAT', 'L1 DA: ASCAT O-F stdv change', -10, 10, 'RdBu_r', '% (DA-OL)/OL'),
              ('pct_L3_ASCAT', 'L3 DA: ASCAT O-F stdv change', -10, 10, 'RdBu_r', '% (DA-OL)/OL'),
              ('pct_SMAP_ASCAT', 'SMAP DA: ASCAT O-F stdv change', -40, 40, 'RdBu_r', '% (DA-OL)/OL (wider scale)'),
              ('rs_smterm', 'L1 operator: SM-term skill', -0.6, 0.6, 'RdBu_r', 'rank correlation'),
              ('dFstd_L1_ascat', 'L1 DA: forecast SM stdv at ASCAT obs', -10, 10, 'RdBu_r', '% change vs OL')]
    map_grid(t, panels, 3, 'fig06_four_corners_ascat.png',
             f'Four Corners ASCAT penalty, 2020-22 (outline: penalty region, {reg.sum()} tiles; L1-area tiles only)',
             lon='com_lon', lat='com_lat', region=reg, extra=bars)


# ---------------------------------------------------------------------------------------------
# step 4, task B: spring windows
# ---------------------------------------------------------------------------------------------
SPRING = [('2020-05', '2020-06'), ('2021-04', '2021-06'), ('2022-02', '2022-05')]   # as paper/make_paper_figures.py


def shade_spring(ax):
    for a0, b0 in SPRING:
        ax.axvspan(pd.Timestamp(a0 + '-01'), pd.Timestamp(b0 + '-01') + pd.offsets.MonthEnd(1), color='0.9', zorder=0, lw=0)


def fig_spring_monthly(a):
    d = LINK + 'spring_windows/'
    m = pd.read_csv(d + 'monthly_series_merged.csv')
    c = pd.read_csv(d + 'monthly_correlations.csv')
    m['t'] = pd.to_datetime(m.ym.astype(str), format='%Y%m') + pd.Timedelta(days=14)

    def rlab(x):
        r = c[(c.y == 'SMAP_prim_minus_ctrl') & (c.x == x) & (c.lag == 0)].iloc[0]
        return f'r = {r.r:+.2f} with L1 - L3 (calendar-month anomalies {r.r_anom:+.2f})'

    fig, axs = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
    ax = axs[0]
    ax.plot(m.t, m.dom_SMAP_prim, '-o', color='#D55E00', ms=3.5, lw=2, label='L1 DA')
    ax.plot(m.t, m.dom_SMAP_ctrl, '-s', color='#0072B2', ms=3.5, lw=2, label='L3 DA (control)')
    ax.bar(m.t, m.SMAP_prim_minus_ctrl, width=20, color='0.55', label='L1 - L3', zorder=1)
    ax.set_ylabel('SMAP Tb O-F stdv\nchange vs OL (%)')
    ax.set_title('(a) SMAP Tb O-F change (domain); negative = DA better', loc='left', fontsize=9)
    ax.legend(ncol=3, loc='lower left')
    series = [('rs_smterm_w', '(b) L1 operator SM-term skill (within-tile Spearman, assimilated obs)', 'rank corr', '#D55E00'),
              ('r_SMAPh_w', '(c) corr(L1 innovation, SMAP Tbh innovation); more negative = more consistent', 'Pearson r',
               '#009E73'),
              ('fcstvar', '(d) ensemble spread of H(x): sets the gain (high in spring 2021 and Apr-May 2022)', 'dB$^2$', '0.3')]
    for ax, (v, title, lab, col) in zip(axs[1:], series):
        ax.plot(m.t, m[v], '-o', color=col, ms=3.5, lw=2)
        ax.set_ylabel(lab)
        ax.set_title(title + ('' if v == 'fcstvar' else f'\n{rlab(v)}'), loc='left', fontsize=9)
    for ax in axs:
        shade_spring(ax)
        ax.axhline(0, color='0.6', lw=0.8) if ax is axs[0] else None
    fig.suptitle('Spring windows (grey): the L1-specific SMAP penalty tracks the collapse of the operator\'s SM sensitivity',
                 fontsize=10)
    fig.tight_layout()
    save(fig, 'fig07_spring_monthly.png')


def fig_spring_maps(a):
    d = LINK + 'spring_windows/'
    w = pd.read_csv(d + 'window_tile_da_pct.csv')
    w = w[(w.year == 'all') & (w.cls == 'spring') & (w.n_l1_obs >= L1_AREA_MIN_OBS)]
    w.loc[w.N_SMAP < 20, 'SMAP'] = np.nan
    t = w[w.arm == 'full_xc015_coh040216_fixedop'][['ig', 'jg', 'com_lon', 'com_lat', 'SMAP']].merge(
        w[w.arm == 'l3_fixedop'][['ig', 'jg', 'SMAP']], on=['ig', 'jg'], suffixes=('_L1', '_L3'))
    k = pd.read_csv(d + 'tile_link_table.csv')
    k = k[(k.year == 'all') & (k.variant == 'da')][['ig', 'jg', 'r_SMAPh_spr', 'rs_smterm_spr', 'rs_smterm_dgood']]
    t = t.merge(k, on=['ig', 'jg'], how='left')
    t['dSMAP'] = t.SMAP_L1 - t.SMAP_L3
    panels = [('SMAP_L1', 'L1 DA: SMAP Tb change, 3 springs', -10, 10, 'RdBu_r', '% (DA-OL)/OL'),
              ('SMAP_L3', 'L3 DA (control): SMAP Tb change, 3 springs', -10, 10, 'RdBu_r', '% (DA-OL)/OL'),
              ('rs_smterm_dgood', 'SM-term skill, spring minus good months', -0.6, 0.6, 'RdBu_r', 'rank corr change'),
              ('r_SMAPh_spr', 'spring corr(L1 inn, SMAP Tbh inn)', -0.5, 0.5, 'RdBu_r', 'Pearson r')]

    def scatter(ax):
        ok = t[['r_SMAPh_spr', 'dSMAP']].notna().all(axis=1)
        ax.scatter(t.r_SMAPh_spr[ok], t.dSMAP[ok], s=8, color='#D55E00', alpha=0.6, lw=0)
        ax.axhline(0, color='0.6', lw=0.8)
        rho = pd.Series(t.r_SMAPh_spr[ok]).rank().corr(pd.Series(t.dSMAP[ok]).rank())
        ax.text(0.03, 0.95, f'Spearman {rho:+.2f} (N={ok.sum()})', transform=ax.transAxes, va='top', fontsize=8)
        ax.set_xlabel('spring corr(L1 inn, SMAP Tbh inn)')
        ax.set_ylabel('SMAP Tb change, L1 DA - L3 DA (points)')
        ax.set_title('L1-specific spring penalty per tile')

    map_grid(t, panels, 3, 'fig08_spring_tile_maps.png',
             'Spring windows pooled (May-Jun 2020, Apr-Jun 2021, Feb-May 2022), L1-area tiles with N >= 20 SMAP obs',
             lon='com_lon', lat='com_lat', extra=scatter)


FIGS = {'op_maps': fig_op_maps, 'geometry': fig_geometry, 'weakh': fig_weakh, 'obs_level': fig_obs_level, 'four_corners': fig_four_corners,
        'spring_monthly': fig_spring_monthly, 'spring_maps': fig_spring_maps}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--figs', nargs='+', default=['all'], help=f'any of {list(FIGS)} or all')
    p.add_argument('--start', default='202001')
    p.add_argument('--end', default='202212')
    a = p.parse_args()
    for f in (FIGS if 'all' in a.figs else a.figs):
        FIGS[f](a)


if __name__ == '__main__':
    main()
