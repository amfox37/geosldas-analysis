#!/usr/bin/env python3
"""
spring_windows: do the CYGNSS L1 operator-skill diagnostics explain the yearly spring degradation of L1 DA?

Steps (all outputs in --out-dir):
  1. Window operator metrics. Every month of each year is classed spring (per-year window, --spring),
     good (--good-months, default Oct-Jan = Jan + Oct-Dec of the same calendar year) or rest.
     From the per-obs OL table (scripts/operator_skill/extract_obs_table.py) compute
       - per-tile metrics via analyze.tile_metrics (+ extra monitor correlations), tile medians over
         L1-area tiles, and paired spring-minus-good tile differences with spatial-block-bootstrap CIs
       - pooled and pooled-within-tile metrics (pooled_metrics below)
     Two obs populations ("variants"):
       da  : obs the primary arm assimilated (obs_da present: coherency 0.40-2.16), scaled with the
             coherency-screened clim (= what DA_L1_full_xc015_coh040216_fixedop used). PRIMARY.
       all : all OL L1 obs, full-stream clim (inn_full), as in analyze.py.
  2. Per-tile DA change. Monthly ObsFcstAna_sums pooled exactly over each window (same algebra as
     postproc_ObsFcstAna.calc_temporal_stats_from_sums), verified against the monthly / yearly /
     whole-period temporal_stats nc4. Per-tile % (DA-OL)/OL OmF_stdv per species, group = nanmean over
     the group's species. Then tile-level Spearman between DA change and operator metrics (block bootstrap).
  3. Monthly (36-month) series of operator metrics vs monthly SMAP Tb % change (domain table and
     L1-area tile mean) and noise-gain alpha_opt / corr; contemporaneous and lagged correlations, raw
     and as calendar-month anomalies, with autocorrelation-adjusted effective N; plus the Feb-Jun
     cross-year test (the months where the window moves).
  4. Figure.

Usage:
  spring_windows_analysis.py                       # defaults = paper SPRING list, Oct-Jan good
  spring_windows_analysis.py --spring 2020-05:2020-06,2021-04:2021-06,2022-02:2022-05 --good-months 1,10,11,12
"""
import argparse
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from netCDF4 import Dataset

P = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/'
sys.path.insert(0, P + 'scripts/operator_skill')
import analyze as AZ  # noqa: E402  (reuse load / tile_metrics / corr / tile_polys)

SUMS = P + 'output/postproc_paired_density/'
STATS = SUMS + 'stats_output/'
TILE_AXIS = P + 'output/operator_skill/da_link/tile_axis_ig_jg.csv'
NOISE_GAIN = P + 'output/full_period_2020_2022/noise_gain_{tag}_202001_202212.csv'
MONTH_TABLES = [P + 'output/month_arm_table_2020/table.csv',
                P + 'output/month_arm_table_2021_2022/month_arm_omf_stdv_pct_vs_OL_2021_2022.csv']
GROUPS = {'SMOS': [0, 1, 2, 3], 'SMAP': [4, 5, 6, 7], 'ASCAT': [8, 9, 10], 'L3': [11], 'L1': [12]}
MONS = ['SMAPh', 'SMAPv', 'SMOSh', 'SMOSv', 'ASCAT', 'L3']
CLASSES = ['spring', 'good', 'rest']


# ----------------------------------------------------------------------------- helpers
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--start', default='202001')
    p.add_argument('--end', default='202212')
    p.add_argument('--spring', default='2020-05:2020-06,2021-04:2021-06,2022-02:2022-05',
                   help='comma list of yyyy-mm:yyyy-mm spring windows, one per year')
    p.add_argument('--good-months', default='1,10,11,12', help='calendar months of the "good" window')
    p.add_argument('--arms', default='full_xc015_coh040216_fixedop,l3_fixedop', help='first = primary, second = control')
    p.add_argument('--tile-nmin-window', type=int, default=30, help='min L1 obs per tile per window for tile metrics')
    p.add_argument('--tile-nmin-month', type=int, default=10, help='min L1 obs per tile per month (monthly tile medians)')
    p.add_argument('--da-nmin', type=int, default=20, help='min N_data per tile and species for window O-F stats')
    p.add_argument('--da-nmin-month', type=int, default=5, help='min N_data per tile and species for single-month O-F stats')
    p.add_argument('--l1-area-nmin', type=int, default=100, help='L1 area = tiles with n_l1_obs >= this')
    p.add_argument('--nboot', type=int, default=1000)
    p.add_argument('--nboot-link', type=int, default=500, help='bootstrap reps for tile-link correlations')
    p.add_argument('--block', type=int, default=3, help='spatial bootstrap block size in ig/jg tiles')
    p.add_argument('--out-dir', default=P + 'output/operator_skill/da_link/spring_windows/')
    p.add_argument('--skip-obs', action='store_true', help='reuse cached step-1 CSVs')
    return p.parse_args()


def month_class(a):
    """dict yyyymm -> (class, year) for every month in [start, end]."""
    spring = set()
    for w in a.spring.split(','):
        m0, m1 = (pd.Period(x, 'M') for x in w.split(':'))
        spring |= {m.strftime('%Y%m') for m in pd.period_range(m0, m1, freq='M')}
    good = {int(x) for x in a.good_months.split(',')}
    out = {}
    for m in pd.period_range(pd.Period(a.start, 'M'), pd.Period(a.end, 'M'), freq='M'):
        k = m.strftime('%Y%m')
        out[k] = 'spring' if k in spring else ('good' if m.month in good else 'rest')
    return out


def block_id(ig, jg, b):
    return (np.asarray(ig) // b) * 10000 + (np.asarray(jg) // b)


def block_boot(df, fn, nboot, b, seed=0):
    """Spatial block bootstrap of statistic fn(df); df has ig, jg. Returns (est, lo, hi) 95%."""
    est = fn(df)
    if not np.isfinite(est):
        return est, np.nan, np.nan
    blk = block_id(df.ig, df.jg, b)
    ub = np.unique(blk)
    idx = {u: np.flatnonzero(blk == u) for u in ub}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(nboot):
        pick = rng.choice(ub, len(ub), replace=True)
        rows = np.concatenate([idx[u] for u in pick])
        vals.append(fn(df.iloc[rows]))
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return est, lo, hi


def spear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return stats.spearmanr(x[ok], y[ok])[0] if ok.sum() >= 10 else np.nan


# ----------------------------------------------------------------------------- step 1: obs-side metrics
def variant(d, v):
    """Population + clim: 'da' = obs assimilated by the primary arm with the coherency clim."""
    if v == 'all':
        return d
    g = d[d.obs_da.notna()].copy()
    g['inn_full'] = g.inn_coh
    g['obs_sc_full'] = g.obs_sc_coh
    g['dyn'] = g.m_std_coh / g.o_std_coh
    return g


def tile_metrics_ext(g):
    s = AZ.tile_metrics(g)
    inn = g.inn_full.values
    for m in ['SMAPv', 'SMOSh', 'SMOSv', 'ASCAT']:
        s[f'r_{m}'] = AZ.corr(inn, g[f'inn_{m}'].values)
    s['inn_mean'] = np.nanmean(inn)
    s['fcstvar_mean'] = g.fcstvar.mean()
    s['o_std_clim'] = g.o_std_full.median() if 'o_std_full' in g else np.nan
    return s


def demean(x, key):
    return x - x.groupby(key).transform('mean')


def pooled_metrics(g):
    """Pooled (all obs) and pooled-within-tile (tile means removed) metrics for one set of obs."""
    key = g.tile
    O, H, G = g.obs, g.fcst, g.geom_db
    S = H - G
    inn = g.inn_full
    ok = inn.notna()
    out = dict(N=len(g), n_tiles=g.tile.nunique(),
               rs_H=spear(O.values, H.values), rs_smterm=spear(O.values, S.values), rs_geom=spear(O.values, G.values),
               rs_H_w=spear(demean(O, key).values, demean(H, key).values),
               rs_smterm_w=spear(demean(O, key).values, demean(S, key).values),
               rs_geom_w=spear(demean(O, key).values, demean(G, key).values),
               r_sc=AZ.corr(g.obs_sc_full.values, H.values),
               r_sc_w=AZ.corr(demean(g.obs_sc_full, key).values, demean(H, key).values),
               frac_out15=(inn.abs() > 15)[ok].mean(), inn_mean=inn.mean(), inn_std=inn.std(),
               R_impl=inn.var() - g.fcstvar[ok].mean(), fcstvar=g.fcstvar[ok].mean(),
               sd_geom_w=demean(G, key).std(), sd_smterm_w=demean(S, key).std(),
               sfmc_mean=g.SFMC_FCST.mean(), sfmc_sd_w=demean(g.SFMC_FCST, key).std(),
               sfmc_sd=g.SFMC_FCST.std(), rzmc_mean=g.RZMC_FCST.mean(),
               opac_med=g.coefficient_weighted_opacity.median(), coh_med=g.coherency_ratio.median(),
               dyn_med=g.dyn.median(), o_std_clim_med=g.o_std_full.median())
    out['snr'] = out['fcstvar'] / out['R_impl'] if out['R_impl'] > 0 else np.nan
    for m in MONS:
        x = g[f'inn_{m}']
        out[f'r_{m}'] = AZ.corr(inn.values, x.values)
        out[f'r_{m}_w'] = AZ.corr(demean(inn, key).values, demean(x, key).values)
        out[f'mean_inn_{m}'] = x.mean()
    return out


TILE_COLS = ['N', 'rs_H', 'rs_smterm', 'rs_geom', 'r_sc', 'sd_geom', 'sd_smterm', 'inn_std', 'R_impl', 'snr',
             'r_SMAPh', 'r_SMAPv', 'r_SMOSh', 'r_SMOSv', 'r_L3', 'r_ASCAT', 'frac_out15', 'inn_mean',
             'sfmc_mean', 'sfmc_std', 'opac_med', 'dyn_ratio', 'fcstvar_mean']


def per_tile(d, nmin):
    t = d.groupby('tile').filter(lambda g: len(g) >= nmin)
    if len(t) == 0:
        return pd.DataFrame()
    t = t.groupby('tile').apply(tile_metrics_ext).reset_index()
    t[['ig', 'jg']] = pd.DataFrame(t.tile.tolist(), index=t.index)
    return t.drop(columns='tile')


def step1(a, d, l1area, out):
    win_rows, tile_rows, mon_rows, mtile_rows = [], [], [], []
    for v in ['da', 'all']:
        dv = variant(d, v)
        dv = dv[dv.tile.isin(l1area)]
        # windows: per year x class, and all years pooled
        for (yr, cl), g in list(dv.groupby(['year_', 'cls'])) + [(('all', c), dv[dv.cls == c]) for c in CLASSES]:
            win_rows.append(dict(variant=v, year=yr, cls=cl, **pooled_metrics(g)))
            t = per_tile(g, a.tile_nmin_window)
            t['variant'], t['year'], t['cls'] = v, yr, cl
            tile_rows.append(t)
            print(f'  {v} {yr} {cl}: {len(g)} obs, {len(t)} tiles', flush=True)
        for ym, g in dv.groupby('ym'):
            mon_rows.append(dict(variant=v, ym=ym, **pooled_metrics(g)))
            t = per_tile(g, a.tile_nmin_month)
            r = t[TILE_COLS].median().add_prefix('tmed_').to_dict()
            r.update(variant=v, ym=ym, tmed_ntiles=len(t))
            mtile_rows.append(r)
        print(f'  {v}: monthly done', flush=True)
    win = pd.DataFrame(win_rows)
    tiles = pd.concat(tile_rows, ignore_index=True)
    mon = pd.DataFrame(mon_rows).merge(pd.DataFrame(mtile_rows), on=['variant', 'ym'])
    win.to_csv(out + 'window_pooled_metrics.csv', index=False)
    tiles.to_csv(out + 'window_tile_metrics.csv', index=False)
    mon.to_csv(out + 'monthly_operator_metrics.csv', index=False)
    return win, tiles, mon


def window_tile_summary(a, tiles, out):
    """Tile medians per window and paired spring-minus-good / spring-minus-rest differences with CIs."""
    rows = []
    for (v, yr), g in tiles.groupby(['variant', 'year']):
        by = {c: g[g.cls == c].set_index(['ig', 'jg']) for c in CLASSES}
        for col in TILE_COLS:
            r = dict(variant=v, year=yr, metric=col)
            for c in CLASSES:
                r[f'med_{c}'] = by[c][col].median() if len(by[c]) else np.nan
                r[f'ntile_{c}'] = by[c][col].notna().sum() if len(by[c]) else 0
            for ref in ['good', 'rest']:
                j = by['spring'][[col]].join(by[ref][[col]], rsuffix='_ref', how='inner').dropna().reset_index()
                j['dd'] = j[col] - j[f'{col}_ref']
                est, lo, hi = block_boot(j, lambda x: x.dd.median(), a.nboot, a.block) if len(j) >= 20 else (np.nan,) * 3
                r.update({f'd_{ref}_med': est, f'd_{ref}_lo': lo, f'd_{ref}_hi': hi, f'd_{ref}_n': len(j)})
            rows.append(r)
    s = pd.DataFrame(rows)
    s.to_csv(out + 'window_tile_median_summary.csv', index=False)
    return s


# ----------------------------------------------------------------------------- step 2: DA per-tile
SUMVARS = ['N_data', 'obs_obs_sum', 'obs_obs2_sum', 'obs_fcst_sum', 'obs_fcst2_sum', 'obsxfcst_sum']


def read_sums(exptag, ym):
    f = f'{SUMS}{exptag}/Y{ym[:4]}/M{ym[4:]}/{exptag}.ens_avg.ldas_ObsFcstAna_sums.{ym}.nc4'
    with Dataset(f) as nc:
        return {k: np.ma.filled(nc[k][:], 0).astype(np.float64) for k in SUMVARS}


def stats_from_sums(s):
    """Same algebra as postproc_ObsFcstAna.calc_temporal_stats_from_sums (population moments)."""
    N = s['N_data'].copy()
    with np.errstate(invalid='ignore', divide='ignore'):
        Nn = np.where(N > 0, N, np.nan)
        Om, Fm = s['obs_obs_sum'] / Nn, s['obs_fcst_sum'] / Nn
        Ov, Fv = s['obs_obs2_sum'] / Nn - Om ** 2, s['obs_fcst2_sum'] / Nn - Fm ** 2
        oxf = s['obsxfcst_sum'] / Nn
        sd = np.sqrt(Ov + Fv - 2 * (oxf - Om * Fm))
    return N, sd


def pool(sums, months):
    out = {k: 0. for k in SUMVARS}
    for m in months:
        for k in SUMVARS:
            out[k] = out[k] + sums[m][k]
    return out


def read_nc4_stdv(exptag, d0, d1):
    f = f'{STATS}temporal_stats_{exptag}_{d0}_{d1}.nc4'
    if not os.path.exists(f):
        return None
    with Dataset(f) as nc:
        return np.ma.filled(nc['OmF_stdv'][:].astype(float), np.nan), np.ma.filled(nc['N_data'][:].astype(float), np.nan)


def verify(exptag, sums, months, out_rows, nmin):
    """Compare pooled-sums OmF_stdv with the toolkit nc4 for each month, each year and the whole period."""
    checks = [(m, [m]) for m in months]
    for y in sorted({m[:4] for m in months}):
        checks.append((y, [m for m in months if m[:4] == y]))
    checks.append(('all', months))
    for lab, ms in checks:
        p0, p1 = pd.Period(ms[0], 'M'), pd.Period(ms[-1], 'M')
        r = read_nc4_stdv(exptag, p0.start_time.strftime('%Y%m%d'), p1.end_time.strftime('%Y%m%d'))
        if r is None:
            continue
        N, sd = stats_from_sums(pool(sums, ms))
        ok = (N > 0) & np.isfinite(r[0]) & np.isfinite(sd)
        diff = np.abs(sd - r[0])[ok]
        rel = (diff / np.abs(r[0][ok]))
        okn = ok & (N >= nmin)
        out_rows.append(dict(exptag=exptag, period=lab, n_cells=int(ok.sum()), max_abs_diff=diff.max(),
                             max_rel_diff=rel.max(), N_match=bool(np.array_equal(N[N > 0], r[1][N > 0])),
                             nan_mismatch=int(((N > 0) & (np.isfinite(sd) != np.isfinite(r[0]))).sum()),
                             n_cells_Nmin=int(okn.sum()), max_abs_diff_Nmin=np.abs(sd - r[0])[okn].max(),
                             max_rel_diff_Nmin=(np.abs(sd - r[0])[okn] / np.abs(r[0][okn])).max(),
                             nan_mismatch_Nmin=int(((N >= nmin) & (np.isfinite(sd) != np.isfinite(r[0]))).sum())))


def within_pooled(sums, blocks):
    """Pool several windows with each window's own O-F mean removed: var = sum N_w var_w / sum N_w."""
    Nt, acc = 0., 0.
    for ms in blocks:
        N, sd = stats_from_sums(pool(sums, ms))
        acc = acc + np.where(N > 0, N * np.nan_to_num(sd) ** 2, 0.)
        Nt = Nt + N
    with np.errstate(invalid='ignore', divide='ignore'):
        return Nt, np.where(Nt > 0, np.sqrt(acc / Nt), np.nan)


def pct_change(sums_da, sums_ol, months, nmin, blocks=None):
    if blocks is None:
        Nd, sd_da = stats_from_sums(pool(sums_da, months))
        No, sd_ol = stats_from_sums(pool(sums_ol, months))
    else:
        Nd, sd_da = within_pooled(sums_da, blocks)
        No, sd_ol = within_pooled(sums_ol, blocks)
    pct = 100. * (sd_da - sd_ol) / sd_ol
    pct[(Nd < nmin) | (No < nmin)] = np.nan
    grp = {g: np.nanmean(pct[:, i], axis=1) for g, i in GROUPS.items()}
    nobs = {f'N_{g}': Nd[:, i].sum(axis=1) for g, i in GROUPS.items()}
    return pct, grp, nobs


def step2(a, mc, tax, out):
    months = sorted(mc)
    arms = a.arms.split(',')
    vrows, rows, mrows = [], [], []
    for arm in arms:
        sda = {m: read_sums(f'DA_paired_{arm}', m) for m in months}
        sol = {m: read_sums(f'OL_paired_monitor_xmask_{arm}', m) for m in months}
        verify(f'DA_paired_{arm}', sda, months, vrows, a.da_nmin)
        verify(f'OL_paired_monitor_xmask_{arm}', sol, months, vrows, a.da_nmin)
        wins = {}
        for y in sorted({m[:4] for m in months}):
            for c in CLASSES:
                wins[(y, c)] = [m for m in months if m[:4] == y and mc[m] == c]
        for (y, c), ms in list(wins.items()) + [(('all', c), None) for c in CLASSES]:
            if ms is None:   # all years: per-year windows pooled with each window's own mean removed
                blocks = [wins[(yy, c)] for yy in sorted({m[:4] for m in months}) if wins[(yy, c)]]
                ms = sum(blocks, [])
            else:
                blocks = None
            if not ms:
                continue
            _, grp, nobs = pct_change(sda, sol, ms, a.da_nmin, blocks)
            df = pd.DataFrame({**grp, **nobs})
            df['tile_axis_index'] = np.arange(len(df))
            df['arm'], df['year'], df['cls'], df['months'] = arm, y, c, ' '.join(ms)
            rows.append(df)
        for m in months:   # monthly per-tile series (lower N threshold: ~10 SMAP obs per species per tile-month)
            _, grp, _ = pct_change(sda, sol, [m], a.da_nmin_month)
            df = pd.DataFrame(grp)
            df['tile_axis_index'] = np.arange(len(df))
            df['arm'], df['ym'] = arm, m
            mrows.append(df)
    ver = pd.DataFrame(vrows)
    ver.to_csv(out + 'verify_sums_vs_nc4.csv', index=False)
    print(ver.groupby('exptag')[['max_abs_diff', 'max_rel_diff', 'nan_mismatch', 'max_abs_diff_Nmin', 'max_rel_diff_Nmin', 'nan_mismatch_Nmin']].max().to_string())
    print('N match all:', ver.N_match.all(), ' periods checked per exptag:', ver.groupby('exptag').size().to_dict())
    tile = pd.concat(rows, ignore_index=True).merge(tax, on='tile_axis_index')
    tile.to_csv(out + 'window_tile_da_pct.csv', index=False)
    mt = pd.concat(mrows, ignore_index=True).merge(tax, on='tile_axis_index')
    return tile, mt, ver


def da_window_summary(a, tile, out):
    """L1-area tile mean of % change per window, block-bootstrap CI."""
    rows = []
    l1 = tile[tile.n_l1_obs >= a.l1_area_nmin]
    for (arm, y, c), g in l1.groupby(['arm', 'year', 'cls']):
        for grp in GROUPS:
            gg = g[['ig', 'jg', grp]].dropna()
            est, lo, hi = block_boot(gg, lambda x: x[grp].mean(), a.nboot, a.block) if len(gg) >= 20 else (np.nan,) * 3
            rows.append(dict(arm=arm, year=y, cls=c, months=g.months.iloc[0], group=grp, ntiles=len(gg),
                             mean_pct=est, lo=lo, hi=hi, frac_tiles_worse=(gg[grp] > 0).mean()))
    s = pd.DataFrame(rows)
    s.to_csv(out + 'window_da_pct_l1area_summary.csv', index=False)
    return s


LINK_BOOT_Y = ['SMAP_spr_prim', 'SMAP_spr_ctrl', 'SMAP_spr_prim_minus_ctrl', 'SMOS_spr_prim', 'L3_spr_prim',
               'SMAP_spr_minus_good_prim', 'SMAP_spr_minus_good_ctrl', 'SMAP_spr_minus_good_prim_minus_ctrl']


def tile_link(a, tiles, tile_da, out):
    """Are the tiles that degrade in spring those whose operator metrics collapse in spring?"""
    arms = a.arms.split(',')
    prim, ctrl = arms[0], arms[1]
    da = tile_da[tile_da.n_l1_obs >= a.l1_area_nmin]
    metrics = ['rs_smterm', 'rs_H', 'r_SMAPh', 'r_SMAPv', 'r_SMOSh', 'r_L3', 'R_impl', 'inn_std', 'snr',
               'frac_out15', 'inn_mean', 'sd_smterm', 'sfmc_mean', 'sfmc_std']
    rows, keep = [], []
    for y in ['all'] + sorted(y for y in tiles.year.astype(str).unique() if y != 'all'):
        if y == 'all':
            ty = tiles[tiles.year.astype(str) == 'all']
        else:
            ty = tiles[tiles.year.astype(str) == y]
        for v in ['da', 'all']:
            tv = ty[ty.variant == v]
            sp = tv[tv.cls == 'spring'].set_index(['ig', 'jg'])
            gd = tv[tv.cls == 'good'].set_index(['ig', 'jg'])
            rs = tv[tv.cls == 'rest'].set_index(['ig', 'jg'])
            ops = sp[metrics].add_suffix('_spr').join(
                (sp[metrics] - gd[metrics]).add_suffix('_dgood')).join((sp[metrics] - rs[metrics]).add_suffix('_drest'))
            d_sp = {arm: da[(da.arm == arm) & (da.year.astype(str) == y) & (da.cls == 'spring')].set_index(['ig', 'jg'])
                    for arm in arms}
            d_gd = {arm: da[(da.arm == arm) & (da.year.astype(str) == y) & (da.cls == 'good')].set_index(['ig', 'jg'])
                    for arm in arms}
            ys = {}
            for grp in ['SMAP', 'SMOS', 'L3', 'L1']:
                ys[f'{grp}_spr_prim'] = d_sp[prim][grp]
                ys[f'{grp}_spr_ctrl'] = d_sp[ctrl][grp]
                ys[f'{grp}_spr_prim_minus_ctrl'] = d_sp[prim][grp] - d_sp[ctrl][grp]
                ys[f'{grp}_spr_minus_good_prim'] = d_sp[prim][grp] - d_gd[prim][grp]
                ys[f'{grp}_spr_minus_good_ctrl'] = d_sp[ctrl][grp] - d_gd[ctrl][grp]
                ys[f'{grp}_spr_minus_good_prim_minus_ctrl'] = ys[f'{grp}_spr_minus_good_prim'] - ys[f'{grp}_spr_minus_good_ctrl']
            J = ops.join(pd.DataFrame(ys), how='inner').reset_index()
            J['year'], J['variant'] = y, v
            keep.append(J)
            for yc in ys:
                for xc in [c for c in ops.columns]:
                    jj = J[['ig', 'jg', xc, yc]].dropna()
                    if len(jj) < 30:
                        continue
                    if v == 'da' and yc in LINK_BOOT_Y and not xc.endswith('_drest'):
                        est, lo, hi = block_boot(jj, lambda q: spear(q[xc].values, q[yc].values), a.nboot_link, a.block)
                    else:
                        est, lo, hi = spear(jj[xc].values, jj[yc].values), np.nan, np.nan
                    rows.append(dict(year=y, variant=v, y=yc, x=xc, n=len(jj), rho=est, lo=lo, hi=hi))
    L = pd.DataFrame(rows)
    L.to_csv(out + 'tile_link_spearman.csv', index=False)
    J = pd.concat(keep, ignore_index=True)
    J.to_csv(out + 'tile_link_table.csv', index=False)
    return L, J


# ----------------------------------------------------------------------------- step 3: monthly series
def neff(x, y):
    def r1(z):
        z = z - np.nanmean(z)
        return np.nansum(z[1:] * z[:-1]) / np.nansum(z * z)
    ok = np.isfinite(x) & np.isfinite(y)
    n = ok.sum()
    rr = r1(np.where(ok, x, np.nan)) * r1(np.where(ok, y, np.nan))
    return n * (1 - rr) / (1 + rr) if rr > -1 else n


def corr_p(x, y, dof_scale=1.0):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 8:
        return np.nan, np.nan, np.nan, 0
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    ne = max(min(neff(x, y), ok.sum()) * dof_scale, 3.01)
    t = r * np.sqrt((ne - 2) / max(1 - r * r, 1e-12))
    return r, 2 * stats.t.sf(abs(t), ne - 2), ne, int(ok.sum())


def step3(a, mc, mon, mt, out):
    months = sorted(mc)
    arms = a.arms.split(',')
    # DA side: domain tables, L1-area tile-mean, noise-gain
    tabs = pd.concat([pd.read_csv(f) for f in MONTH_TABLES])
    tabs = tabs[tabs.period.astype(str).str.len() == 6]
    tabs['ym'] = tabs.period.astype(str)
    y = pd.DataFrame({'ym': months})
    y['cls'] = y.ym.map(mc)
    for arm, lab in zip(arms, ['prim', 'ctrl']):
        tt = tabs[tabs.arm == arm].drop_duplicates('ym').set_index('ym')
        for g in ['SMAP', 'SMOS', 'L3', 'L1']:
            y[f'dom_{g}_{lab}'] = y.ym.map(tt[g])
        m1 = mt[(mt.arm == arm) & (mt.n_l1_obs >= a.l1_area_nmin)].groupby('ym')[['SMAP', 'SMOS', 'L3', 'L1']].mean()
        for g in ['SMAP', 'SMOS', 'L3', 'L1']:
            y[f'l1a_{g}_{lab}'] = y.ym.map(m1[g])
        ng = pd.read_csv(NOISE_GAIN.format(tag=arm))
        ng = ng[(ng.grp == 'SMAP') & (ng.month != 'all')]
        ng['ym'] = ng.month.str.replace('-', '')
        ng = ng.set_index('ym')
        for c in ['alpha_opt', 'corr', 'noise', 'gain', 'dMSE']:
            y[f'ng_{c}_{lab}'] = y.ym.map(ng[c])
    y['SMAP_prim_minus_ctrl'] = y.dom_SMAP_prim - y.dom_SMAP_ctrl
    y.to_csv(out + 'monthly_da_series.csv', index=False)

    md = mon[mon.variant == 'da'].copy()
    md['ym'] = md.ym.astype(str)
    S = y.merge(md, on='ym', how='left')
    S['cal'] = S.ym.str[4:].astype(int)
    S.to_csv(out + 'monthly_series_merged.csv', index=False)
    xcols = ['rs_H', 'rs_H_w', 'rs_smterm', 'rs_smterm_w', 'rs_geom_w', 'r_sc', 'r_sc_w', 'frac_out15', 'inn_mean',
             'inn_std', 'R_impl', 'snr', 'fcstvar', 'sd_geom_w', 'sd_smterm_w', 'sfmc_mean', 'sfmc_sd_w', 'rzmc_mean',
             'opac_med', 'coh_med', 'dyn_med', 'o_std_clim_med', 'N',
             'r_SMAPh', 'r_SMAPv', 'r_SMOSh', 'r_SMOSv', 'r_L3', 'r_ASCAT',
             'r_SMAPh_w', 'r_SMAPv_w', 'r_L3_w', 'mean_inn_SMAPh', 'mean_inn_SMAPv', 'mean_inn_L3',
             'tmed_rs_H', 'tmed_rs_smterm', 'tmed_r_SMAPh', 'tmed_r_L3', 'tmed_R_impl', 'tmed_sd_smterm', 'tmed_inn_mean']
    ycols = ['dom_SMAP_prim', 'l1a_SMAP_prim', 'SMAP_prim_minus_ctrl', 'dom_L3_prim', 'ng_alpha_opt_prim', 'ng_corr_prim',
             'ng_noise_prim', 'dom_SMAP_ctrl']
    anom = lambda s: s - s.groupby(S.cal).transform('mean')
    feb_jun = S.cal.between(2, 6).values
    rows = []
    for yc in ycols:
        for xc in xcols:
            for lag in [0, 1, -1]:     # lag 1: operator metric leads DA change by one month
                xr = S[xc].shift(lag).values
                yr = S[yc].values
                r, p, ne, n = corr_p(xr, yr)
                xa, ya = anom(S[xc]).shift(lag).values, anom(S[yc]).values
                ra, pa, nea, _ = corr_p(xa, ya, dof_scale=24 / 36)
                row = dict(y=yc, x=xc, lag=lag, n=n, r=r, p=p, neff=ne, r_anom=ra, p_anom=pa, neff_anom=nea)
                if lag == 0:
                    ok = feb_jun & np.isfinite(xr) & np.isfinite(yr)
                    row['r_febjun'] = np.corrcoef(xr[ok], yr[ok])[0, 1] if ok.sum() > 5 else np.nan
                    row['rs_febjun'] = stats.spearmanr(xr[ok], yr[ok])[0] if ok.sum() > 5 else np.nan
                    row['p_febjun'] = stats.spearmanr(xr[ok], yr[ok])[1] if ok.sum() > 5 else np.nan
                    ra2, pa2, _, _ = corr_p(xa[feb_jun], ya[feb_jun], dof_scale=10 / 15)
                    row['r_febjun_anom'], row['p_febjun_anom'] = ra2, pa2
                    sp = S.cls.values == 'spring'
                    row['spring_mean'] = np.nanmean(S[xc].values[sp])
                    row['good_mean'] = np.nanmean(S[xc].values[S.cls.values == 'good'])
                    row['rest_mean'] = np.nanmean(S[xc].values[S.cls.values == 'rest'])
                    row['spring_anom_mean'] = np.nanmean(xa[sp])
                    row['spring_anom_t'] = stats.ttest_1samp(xa[sp & np.isfinite(xa)], 0)[0]
                rows.append(row)
    C = pd.DataFrame(rows)
    C.to_csv(out + 'monthly_correlations.csv', index=False)
    return S, C


# ----------------------------------------------------------------------------- step 4: figure
def figure(a, S, J, tile_da, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    BLUE, ORANGE, AQUA, GREY = '#2a78d6', '#eb6834', '#1baf7a', '0.45'
    t = pd.to_datetime(S.ym, format='%Y%m') + pd.Timedelta(days=14)
    spr = []
    for w in a.spring.split(','):
        m0, m1 = (pd.Period(x, 'M') for x in w.split(':'))
        spr.append((m0.start_time, m1.end_time))
    panels = [
        ('SMAP Tb O-F stdv change vs OL [%]  (negative = DA better)',
         [('dom_SMAP_prim', 'L1 arm', BLUE, '-'), ('dom_SMAP_ctrl', 'L3-only control', ORANGE, '-'),
          ('SMAP_prim_minus_ctrl', 'L1 arm minus control', GREY, '--')], True),
        ('SMAP noise-gain, L1 arm: corr(dF, inn)', [('ng_corr_prim', 'corr(dF, inn)', BLUE, '-')], True),
        ('SMAP noise-gain, L1 arm: noise = sum dF^2 [% of OL MSE]', [('ng_noise_prim', 'noise', BLUE, '-')], False),
        ('Operator SM-term skill: Spearman(O, H - geometry term)',
         [('rs_smterm_w', 'pooled within-tile', BLUE, '-'), ('tmed_rs_smterm', 'tile median', AQUA, ':')], True),
        ('corr(L1 innovation, SMAP Tbh innovation), within-tile  (more negative = more shared SM signal)',
         [('r_SMAPh_w', 'SMAP Tbh', BLUE, '-')], True),
        ('OL surface SM at L1 obs: within-tile std [m3/m3]', [('sfmc_sd_w', 'within-tile std', BLUE, '-')], False),
        ('Ensemble variance of H(x) [dB^2]', [('fcstvar', 'fcstvar', BLUE, '-')], False),
    ]
    fig = plt.figure(figsize=(14, 17))
    gs = fig.add_gridspec(len(panels) + 2, 3, height_ratios=[1] * len(panels) + [0.15, 2.2])
    for i, (title, series, zero) in enumerate(panels):
        ax = fig.add_subplot(gs[i, :])
        for a0, a1 in spr:
            ax.axvspan(a0, a1, color='0.9', lw=0, zorder=0)
        for col, lab, c, ls in series:
            ax.plot(t, S[col], ls, color=c, lw=2, marker='o', ms=3, label=lab)
        if zero:
            ax.axhline(0, color='0.6', lw=0.7)
        ax.set_title(title, loc='left', fontsize=10)
        if len(series) > 1:
            ax.legend(fontsize=8, ncol=3, frameon=False, loc='upper left', bbox_to_anchor=(1.0, 1.0))
        ax.spines[['top', 'right']].set_visible(False)
        if i < len(panels) - 1:
            ax.set_xticklabels([])
    # bottom row: per-tile pooled-spring map + 2 scatters
    jj = J[(J.year == 'all') & (J.variant == 'da')]
    ax = fig.add_subplot(gs[-1, 0])
    m = tile_da[(tile_da.arm == a.arms.split(',')[0]) & (tile_da.year.astype(str) == 'all') & (tile_da.cls == 'spring')
                & (tile_da.n_l1_obs >= a.l1_area_nmin)].dropna(subset=['SMAP'])
    tt = m.merge(AZ.tile_centres(), on=['ig', 'jg'], how='inner')   # grid-cell centres for tile_polys
    pc = PolyCollection(AZ.tile_polys(tt), array=tt.SMAP.values, cmap='RdBu_r', linewidths=0)
    pc.set_clim(-10, 10)
    ax.add_collection(pc)
    for b in AZ.BORDERS:
        ax.plot(b[:, 0], b[:, 1], color='0.3', lw=0.5)
    ax.set_xlim(-118, -106); ax.set_ylim(29, 40); ax.set_aspect(1.15)
    ax.set_title(f'L1 arm SMAP Tb change, all spring windows [%]\n(L1-area tiles, N>={a.da_nmin}; red = DA worse)', fontsize=9)
    plt.colorbar(pc, ax=ax, shrink=0.7)
    for k, (xc, yp, yc, lab, ylab) in enumerate([
            ('r_SMAPh_spr', 'SMAP_spr_prim', 'SMAP_spr_ctrl', 'corr(L1, SMAPh innovation) in spring', 'spring SMAP Tb change [%]'),
            ('sfmc_mean_dgood', 'SMAP_spr_minus_good_prim', 'SMAP_spr_minus_good_ctrl', 'OL SFMC, spring - good [m3/m3]',
             'SMAP Tb change, spring - good [% points]')]):
        ax = fig.add_subplot(gs[-1, k + 1])
        q = jj[[xc, yp, yc]].dropna()
        ax.scatter(q[xc], q[yp], s=9, color=BLUE, alpha=0.6, label='L1 arm', edgecolors='none')
        ax.scatter(q[xc], q[yc], s=9, color=ORANGE, alpha=0.4, label='L3-only control', edgecolors='none')
        ax.axhline(0, color='0.6', lw=0.7)
        ax.set_xlabel(lab); ax.set_ylabel(ylab)
        ax.set_title(f'per tile, 3 springs, n={len(q)}\nSpearman: L1 arm {spear(q[xc].values, q[yp].values):+.2f}, '
                     f'control {spear(q[xc].values, q[yc].values):+.2f}', fontsize=9)
        ax.set_ylim(-25, 25)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('CYGNSS L1 operator metrics (OL, assimilated obs, L1-area tiles) vs DA SMAP Tb change; '
                 'grey = spring windows', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out + 'spring_windows_summary.png', dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    a = parse_args()
    out = a.out_dir.rstrip('/') + '/'
    os.makedirs(out, exist_ok=True)
    AZ.BORDERS = AZ.state_borders()
    mc = month_class(a)
    tax = pd.read_csv(TILE_AXIS)
    l1area = set(zip(tax.ig[tax.n_l1_obs >= a.l1_area_nmin], tax.jg[tax.n_l1_obs >= a.l1_area_nmin]))
    print(f'L1 area: {len(l1area)} tiles; month classes:', {c: sum(v == c for v in mc.values()) for c in CLASSES})

    if a.skip_obs and os.path.exists(out + 'window_tile_metrics.csv'):
        win = pd.read_csv(out + 'window_pooled_metrics.csv')
        tiles = pd.read_csv(out + 'window_tile_metrics.csv')
        mon = pd.read_csv(out + 'monthly_operator_metrics.csv', dtype={'ym': str})
    else:
        d = AZ.load(a)
        tt = pd.to_datetime(d.time)
        d['ym'] = tt.dt.strftime('%Y%m')
        d['year_'] = tt.dt.strftime('%Y')
        d['cls'] = d.ym.map(mc)
        print(f'{len(d)} obs; assimilated by primary arm: {d.obs_da.notna().sum()}', flush=True)
        win, tiles, mon = step1(a, d, l1area, out)
    tiles['year'] = tiles.year.astype(str)
    window_tile_summary(a, tiles, out)

    tile_da, mt, ver = step2(a, mc, tax, out)
    da_window_summary(a, tile_da, out)
    L, J = tile_link(a, tiles, tile_da, out)
    S, C = step3(a, mc, mon, mt, out)
    figure(a, S, J, tile_da, out)
    print('wrote', out)


if __name__ == '__main__':
    main()
