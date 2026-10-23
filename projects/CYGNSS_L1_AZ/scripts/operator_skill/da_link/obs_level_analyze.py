#!/usr/bin/env python3
"""
obs_level (task D), step 2: link the L1 operator-skill bins to what the DA did with each L1 obs,
and to the independent monitors that follow it.

Input: obs_level_l1_*.parquet / obs_level_mon_*.parquet from obs_level_extract.py.

Definitions (DA arm = DA_L1_full_xc015_coh040216_fixedop, OL = OLv8_M36_AZ_fixedop; O = the DA
arm's scaled obs for both, matched on T/species/tile/lon/lat):
  omf_da = O - F_DA, omf_ol = O - F_OL, oma = O - A_DA, inc = A_DA - F_DA,
  ninn = (O - F_DA) / sqrt(obsvar + fcstvar_DA)
  dF = F_DA - F_OL, inn = O - F_OL                     (cygl1_noise_gain_split.py)
  noise = 100 sum dF^2 / sum inn^2, gain = -200 sum dF inn / sum inn^2, dMSE = noise + gain,
  alpha_opt = sum dF inn / sum dF^2
  dstd = 100 (std(omf_da) - std(omf_ol)) / std(omf_ol)   (pooled over the obs in the bin)
  K_eff = sum inc*omf_da / sum omf_da^2 (regression of the increment on the innovation),
  rms_ratio = rms(O-A)/rms(O-F_DA)
All %-changes are (DA-OL)/OL, negative = DA better.  CIs: 2.5-97.5 % of a tile-cluster bootstrap
(owner tile resampled with replacement, B replicates, per bin).

Operator classes of an L1 obs (from the OL operator table, extract_obs_table.py):
  weakH        |inn_full| > 15 dB (scaled with the full-stream clim, as in analyze.py's headline)
  coh_tail     coherency_ratio outside [0.40, 2.16]
  H octile     octiles of the OL H(x) over ALL OL-table obs (2020-22)
  smterm       H(x) - geom_db (soil-moisture/reflectivity part of H)
  tile_smskill tercile over the 556 tiles of tile_metrics.csv rs_smterm (owner tile)
  continuous bins: octiles over all OL-table obs

Monitor attribution (question 3): a monitor obs (SMAP Tb, L3, SMOS, ASCAT) in tile k at time Tm
is attributed to the assimilated L1 obs with the LARGEST |A-F| among those in the same tile with
Tm - 72 h <= T_L1 < Tm (strictly earlier window: the OFA fcst at Tm is pre-analysis).  The monitor
obs takes that L1 obs's operator class.  Monitor obs in L1-area tiles (n_l1_obs >= 100) with no L1
analysis in the preceding 72 h form the 'no L1 in 72h' control.  Sensitivity: the most recent L1
obs within 24 h.

Usage:
  obs_level_analyze.py --tag 202006_202006 --inputs 202006_202006
  obs_level_analyze.py --tag 202001_202212 --inputs 202001_202012,202101_202112,202201_202212
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

HERE = '/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/operator_skill/'
OUT0 = HERE + 'da_link/obs_level/'
COH_LO, COH_HI = 0.40, 2.16
SEASON = {12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
          6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'}
SPRING = {2020: [5, 6], 2021: [4, 5, 6], 2022: [2, 3, 4, 5]}
CONT = ['fcst_ol_tab', 'smterm', 'geom_db', 'sp_inc_angle', 'coherency_ratio', 'ddm_snr_db', 'srtm_slope',
        'coefficient_weighted_opacity', 'n_eff', 'w_max']
CAT = ['weakH', 'tile_smskill', 'sc_num', 'season', 'period']
LABEL = {'fcst_ol_tab': 'H(x) octile [dB]', 'smterm': 'SM term H - geom [dB]', 'geom_db': 'geometry term [dB]',
         'sp_inc_angle': 'incidence [deg]', 'coherency_ratio': 'coherency ratio', 'ddm_snr_db': 'DDM SNR [dB]',
         'srtm_slope': 'SRTM slope [deg]', 'coefficient_weighted_opacity': 'opacity', 'n_eff': 'n_eff',
         'w_max': 'w_max', 'weakH': 'weak-H class', 'tile_smskill': 'tile SM-term skill tercile',
         'sc_num': 'spacecraft', 'season': 'season', 'period': 'period', 'inc_q': '|A-F| quintile (L1)'}
B_DEFAULT = 500
MON_GRP = {'SMAPh': 'SMAP', 'SMAPv': 'SMAP', 'SMOSh': 'SMOS', 'SMOSv': 'SMOS', 'L3': 'L3', 'ASCAT': 'ASCAT'}


# ------------------------------------------------------------------ statistics
L1_SUMS = ['n', 'a1', 'a2', 'b1', 'b2', 'dd', 'di', 'oa2', 'kf', 'ainc', 'nn2']
MON_SUMS = ['n', 'a1', 'a2', 'b1', 'b2', 'dd', 'di']


def per_obs_sums(omf_da, omf_ol, oma=None, inc=None, ninn=None):
    dF = omf_ol - omf_da          # F_DA - F_OL
    s = dict(n=np.ones_like(omf_da), a1=omf_da, a2=omf_da ** 2, b1=omf_ol, b2=omf_ol ** 2,
             dd=dF ** 2, di=dF * omf_ol)
    if oma is not None:
        s.update(oa2=oma ** 2, kf=inc * omf_da, ainc=np.abs(inc), nn2=ninn ** 2)
    return pd.DataFrame(s)


def stats_from(S):
    """S: dict of sums (scalars or arrays over bootstrap replicates)."""
    n = S['n']
    sd = lambda s1, s2: np.sqrt(np.maximum(s2 / n - (s1 / n) ** 2, 0) * n / np.maximum(n - 1, 1))
    o = dict(dstd=100 * (sd(S['a1'], S['a2']) - sd(S['b1'], S['b2'])) / sd(S['b1'], S['b2']),
             dmse=100 * (S['a2'] - S['b2']) / S['b2'], noise=100 * S['dd'] / S['b2'],
             gain=-200 * S['di'] / S['b2'], alpha=S['di'] / S['dd'])
    if 'oa2' in S:
        o.update(rms_ratio=np.sqrt(S['oa2'] / S['a2']), K_eff=S['kf'] / S['a2'], mean_abs_inc=S['ainc'] / n,
                 rms_ninn=np.sqrt(S['nn2'] / n))
    return o


def boot_stats(df, sums_cols, B, rng, cluster='tile'):
    """Point estimate + tile-cluster bootstrap 95% CI of stats_from for one group."""
    M = df.groupby(cluster)[sums_cols].sum()
    pt = stats_from({c: M[c].sum() for c in sums_cols})
    nt = len(M)
    out = {'N': int(M['n'].sum()), 'N_tiles': nt}
    if nt >= 5 and B > 0:
        W = rng.multinomial(nt, np.full(nt, 1. / nt), size=B).astype(float)
        SB = W @ M.values
        bs = stats_from({c: SB[:, i] for i, c in enumerate(sums_cols)})
    for k, v in pt.items():
        out[k] = v
        if nt >= 5 and B > 0:
            lo, hi = np.nanpercentile(bs[k], [2.5, 97.5])
            out[k + '_lo'], out[k + '_hi'] = lo, hi
    return out


# ------------------------------------------------------------------ data
def load(a):
    L = pd.concat([pd.read_parquet(f'{a.data_dir}obs_level_l1_{t}.parquet') for t in a.inputs.split(',')],
                  ignore_index=True)
    M = pd.concat([pd.read_parquet(f'{a.data_dir}obs_level_mon_{t}.parquet') for t in a.inputs.split(',')],
                  ignore_index=True)
    tm = pd.read_csv(f'{HERE}analysis_202001_202212/tile_metrics.csv')
    ta = pd.read_csv(f'{HERE}da_link/tile_axis_ig_jg.csv')
    ta['tilenum'] = ta.tile_axis_index + 1          # OFA tilenum = tile_axis_index + 1 (verified)
    q = tm.rs_smterm.quantile([1 / 3, 2 / 3]).values
    tm['tile_smskill'] = np.where(tm.rs_smterm < q[0], 'T1 low', np.where(tm.rs_smterm < q[1], 'T2 mid', 'T3 high'))
    tm = tm.merge(ta[['ig', 'jg', 'tilenum']], on=['ig', 'jg'])
    print(f'rs_smterm tercile edges {q.round(3)} over {len(tm)} tiles')

    # time / class columns
    t = pd.to_datetime(L['T'], unit='s', origin='2020-01-01')
    L['month'], L['year'] = t.dt.month, t.dt.year
    L['season'] = L.month.map(SEASON)
    L['period'] = np.where([m in SPRING.get(y, []) for y, m in zip(L.year, L.month)], 'spring window',
                           np.where(L.month.isin([10, 11, 12, 1]), 'Oct-Jan', 'other'))
    L['weakH'] = np.where(L.inn_full.abs() > 15, 'weak-H outlier', 'normal')
    L['coh_tail'] = ~L.coherency_ratio.between(COH_LO, COH_HI)
    L['smterm'] = L.fcst_ol_tab - L.geom_db
    tab = L[L.in_tab]
    for v in CONT:
        edges = np.unique(tab[v].quantile(np.linspace(0, 1, 9)).values)
        edges[0], edges[-1] = -np.inf, np.inf
        L[v + '_bin'] = pd.cut(L[v], edges, labels=False)
        mids = tab.groupby(pd.cut(tab[v], edges, labels=False))[v].median()
        L[v + '_med'] = L[v + '_bin'].map(mids)
    L = L.merge(tm[['tilenum', 'tile_smskill', 'rs_smterm']], on='tilenum', how='left')
    L['tile_smskill'] = L.tile_smskill.fillna('no tile metric')
    L['tile'] = L.tilenum
    L['l1_area'] = L.tilenum.isin(ta.tilenum[ta.n_l1_obs >= 100])

    M['grp2'] = M.grp.map(MON_GRP)
    M['tile'] = M.tilenum
    tm_ = pd.to_datetime(M['T'], unit='s', origin='2020-01-01')
    M['month'], M['year'] = tm_.dt.month, tm_.dt.year
    M['l1_area'] = M.tilenum.isin(ta.tilenum[ta.n_l1_obs >= 100])
    s = per_obs_sums((M.obs_da - M.fcst_da).values.astype(float), (M.obs_da - M.fcst_ol).values.astype(float))
    for c in s:
        M[c] = s[c].values
    return L, M


# ------------------------------------------------------------------ question 1
def q1(L, out):
    rows = []
    T = L[L.in_tab]
    def add(name, sub):
        d = sub[sub.in_da & sub.fcst_ol.notna()]
        r = dict(cls=name, N_OL_table=len(sub), N_in_DA=int(sub.in_da.sum()),
                 frac_in_DA=sub.in_da.mean() if len(sub) else np.nan,
                 N_assim_flag1=int((d.assim_flag == 1).sum()))
        if len(d):
            r.update(med_abs_inc=np.median(np.abs(d.ana_da - d.fcst_da)),
                     p90_abs_inc=np.percentile(np.abs(d.ana_da - d.fcst_da), 90),
                     med_abs_omf=np.median(np.abs(d.obs_da - d.fcst_da)),
                     med_abs_ninn=np.median(np.abs(d.obs_da - d.fcst_da) / np.sqrt(d.obsvar + d.fcstvar_da)),
                     p90_abs_ninn=np.percentile(np.abs(d.obs_da - d.fcst_da) / np.sqrt(d.obsvar + d.fcstvar_da), 90),
                     med_obsvar=d.obsvar.median(), med_fcstvar=d.fcstvar_da.median(),
                     med_coherency=d.coherency_ratio.median(), med_H=d.fcst_ol_tab.median(),
                     med_O_raw=d.obs_raw.median(), N_abs_innDA_gt15=int(((d.obs_da - d.fcst_ol).abs() > 15).sum()))
        rows.append(r)
    add('all', T)
    add('weak-H outlier (|inn_full|>15)', T[T.weakH == 'weak-H outlier'])
    add('weak-H outlier, coherency in [0.40,2.16]', T[(T.weakH == 'weak-H outlier') & ~T.coh_tail])
    add('normal', T[T.weakH == 'normal'])
    add('coherency tail (<0.40 or >2.16)', T[T.coh_tail])
    add('coherency < 0.40', T[T.coherency_ratio < COH_LO])
    add('coherency > 2.16', T[T.coherency_ratio > COH_HI])
    add('coherency NaN', T[T.coherency_ratio.isna()])
    add('coherency in range', T[~T.coh_tail & T.coherency_ratio.notna()])
    add('coherency in range, not in DA', T[~T.coh_tail & T.coherency_ratio.notna() & ~T.in_da])
    for b in range(8):
        add(f'H octile {b + 1}', T[T.fcst_ol_tab_bin == b])
    r = pd.DataFrame(rows)
    r.to_csv(f'{out}q1_assimilated_by_class.csv', index=False)
    print('\n=== Q1: which OL-table L1 obs reach the DA arm ===')
    print(r.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    print(f'DA-arm L1 obs not in OL table: {(~L.in_tab).sum()}; assim_flag values in DA L1:',
          L.assim_flag.value_counts().to_dict())
    return r


# ------------------------------------------------------------------ question 2
def l1_matched(L):
    d = L[L.in_da & L.in_tab & L.fcst_ol.notna()].copy()
    d['omf_da'] = d.obs_da - d.fcst_da
    d['omf_ol'] = d.obs_da - d.fcst_ol
    d['oma'] = d.obs_da - d.ana_da
    d['inc'] = d.ana_da - d.fcst_da
    d['ninn'] = d.omf_da / np.sqrt(d.obsvar + d.fcstvar_da)
    s = per_obs_sums(d.omf_da.values, d.omf_ol.values, d.oma.values, d.inc.values, d.ninn.values)
    for c in s:
        d['s_' + c] = s[c].values
    return d


def binned_l1(d, B, rng, sub_label='all'):
    rows = []
    cols = ['s_' + c for c in L1_SUMS]
    for v in CONT + CAT:
        key = v + '_bin' if v in CONT else v
        for k, g in d.groupby(key, observed=True):
            r = boot_stats(g.rename(columns=lambda c: c[2:] if c.startswith('s_') else c), L1_SUMS, B, rng)
            ratio = g.oma / g.omf_da
            r.update(var=v, bin=k if v not in CONT else int(k) + 1,
                     x_med=g[v + '_med'].iloc[0] if v in CONT else np.nan,
                     med_ratio_OmA_OmF=np.nanmedian(ratio[np.abs(g.omf_da) > 0.1]),
                     med_abs_ninn=np.median(np.abs(g.ninn)), subset=sub_label)
            rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ question 3
def attribute(M, d, win_h=72.):
    """Dominant (largest |A-F|) assimilated L1 obs in the same tile, Tm-win <= T_L1 < Tm."""
    d = d.sort_values(['tilenum', 'T']).reset_index(drop=True)
    kl = d.tilenum.values * 1e9 + d['T'].values
    km = M.tilenum.values * 1e9 + M['T'].values
    lo = np.searchsorted(kl, km - win_h * 3600, 'left')
    hi = np.searchsorted(kl, km, 'left')
    cnt = hi - lo
    mi = np.repeat(np.arange(len(M)), cnt)
    off = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    li = np.repeat(lo, cnt) + off
    pr = pd.DataFrame({'mi': mi, 'li': li, 'ainc': np.abs(d.inc.values[li])})
    pr = pr.sort_values(['mi', 'ainc'], ascending=[True, False]).drop_duplicates('mi')
    dom = np.full(len(M), -1)
    dom[pr.mi.values] = pr.li.values
    rec = np.where((cnt > 0) & (M['T'].values - d['T'].values[np.maximum(hi - 1, 0)] <= 24 * 3600), hi - 1, -1)
    nweak = np.bincount(mi, weights=(d.weakH.values[li] == 'weak-H outlier'), minlength=len(M))
    return d, dom, rec, cnt, nweak


def monitor_link(M, d, B, rng, out):
    d, dom, rec, cnt, nweak = attribute(M, d)
    M = M.copy()
    M['n_L1_72h'], M['n_weakH_72h'] = cnt, nweak
    M['has_L1'] = dom >= 0
    M['lag_h'] = np.where(dom >= 0, (M['T'].values - d['T'].values[np.maximum(dom, 0)]) / 3600., np.nan)
    d['inc_q'] = pd.qcut(d.inc.abs(), 5, labels=False)
    d['inc_q_med'] = d.inc_q.map(d.groupby('inc_q').inc.apply(lambda x: x.abs().median()))
    attrs = [v + '_bin' for v in CONT] + [v + '_med' for v in CONT] + CAT + ['inc_q', 'inc_q_med']
    for c in attrs:
        vals = d[c].values
        M['dom_' + c] = pd.Series(vals[np.maximum(dom, 0)]).where(dom >= 0).values
        M['rec_' + c] = pd.Series(vals[np.maximum(rec, 0)]).where(rec >= 0).values
    rows = []
    for grp, gm in M.groupby('grp2'):
        base = gm[gm.l1_area]
        for lab, sub in [('all obs (domain)', gm), ('L1 area', base),
                         ('L1 area, no L1 in 72h', base[~base.has_L1]), ('after an L1 analysis (72h)', gm[gm.has_L1]),
                         ('lag 0-24h', gm[gm.lag_h <= 24]), ('lag 24-72h', gm[gm.lag_h > 24]),
                         ('any weak-H in 72h', gm[gm.n_weakH_72h > 0])]:
            rows.append(dict(grp=grp, var='control', bin=lab, rule='dominant72', **boot_stats(sub, MON_SUMS, B, rng)))
        for rule, pre in [('dominant72', 'dom_'), ('recent24', 'rec_')]:
            for v in CONT + CAT + ['inc_q']:
                key = pre + (v + '_bin' if v in CONT else v)
                for k, g in gm.groupby(key, observed=True):
                    r = boot_stats(g, MON_SUMS, B, rng)
                    xm = g[pre + v + '_med'].iloc[0] if (v in CONT or v == 'inc_q') else np.nan
                    rows.append(dict(grp=grp, var=v, bin=int(k) + 1 if (v in CONT or v == 'inc_q') else k,
                                     x_med=xm, rule=rule, **r))
    r = pd.DataFrame(rows)
    r.to_csv(f'{out}q3_monitor_by_L1_class.csv', index=False)
    keep = ['grp', 'n_L1_72h', 'n_weakH_72h', 'has_L1', 'lag_h', 'T', 'tilenum', 'dom_weakH', 'dom_tile_smskill',
            'dom_fcst_ol_tab_bin', 'dom_inc_q', 'dom_period']
    M[keep].to_parquet(f'{out}q3_monitor_attribution.parquet')
    return r, M


# ------------------------------------------------------------------ figure
def figure(b2, r3, out, tag):
    blue, orange, aqua = '#2a78d6', '#eb6834', '#1baf7a'
    vars_ = ['fcst_ol_tab', 'smterm', 'geom_db', 'sp_inc_angle', 'coherency_ratio', 'ddm_snr_db', 'srtm_slope',
             'coefficient_weighted_opacity', 'n_eff', 'weakH', 'tile_smskill', 'season', 'period']
    cols = ['L1 own-fit change [% (DA-OL)/OL]', 'L1 alpha_opt (1 = optimal)', 'L1 increment: mean |A-F| [dB]',
            'Monitor dMSE after L1 analysis [% of OL MSE]']
    fig, axs = plt.subplots(len(vars_), 4, figsize=(17, 2.3 * len(vars_)), squeeze=False)
    b = b2[b2.subset == 'all']
    for i, v in enumerate(vars_):
        g = b[b['var'] == v].copy()
        m3 = r3[(r3['var'] == v) & (r3.rule == 'dominant72')]
        if v in CONT:
            xl = [f'{x:.3g}' for x in g.x_med]
        else:
            order = sorted(g.bin.astype(str).unique())
            g = g.set_index(g.bin.astype(str)).loc[order].reset_index(drop=True)
            xl = order
        x = np.arange(len(g))
        def ebar(ax, xx, y, lo, hi, c, lab, mk='o', dx=0.):
            ax.errorbar(xx + dx, y, yerr=[y - lo, hi - y], fmt=mk + '-', color=c, ms=5, lw=1.6, capsize=2, label=lab)
        ax = axs[i, 0]
        ebar(ax, x, g.dstd, g.dstd_lo, g.dstd_hi, blue, 'std(O-F)', dx=-0.08)
        ebar(ax, x, g.dmse, g.dmse_lo, g.dmse_hi, orange, 'MSE(O-F)', 's', dx=0.08)
        ax = axs[i, 1]
        ebar(ax, x, g.alpha, g.alpha_lo, g.alpha_hi, blue, 'alpha_opt')
        ax.axhline(1, color='0.5', lw=0.8, ls='--')
        ax = axs[i, 2]
        ebar(ax, x, g.mean_abs_inc, g.mean_abs_inc_lo, g.mean_abs_inc_hi, blue, 'mean |A-F|')
        ax2 = axs[i, 3]
        for grp, c, mk, dx in [('SMAP', blue, 'o', -0.1), ('L3', orange, 's', 0.0), ('SMOS', aqua, '^', 0.1)]:
            mm = m3[m3.grp == grp].copy()
            if not len(mm):
                continue
            if v in CONT:
                mm = mm.sort_values('bin'); xx = mm.bin.values.astype(int) - 1
            else:
                mm = mm.set_index(mm.bin.astype(str)).reindex(xl).reset_index(drop=True); xx = x
            ebar(ax2, xx, mm.dmse.values, mm.dmse_lo.values, mm.dmse_hi.values, c, grp, mk, dx)
        for j in range(4):
            a_ = axs[i, j]
            if j != 1:
                a_.axhline(0, color='0.6', lw=0.8)
            a_.set_xticks(x); a_.set_xticklabels(xl, rotation=30 if v in CONT else 15, fontsize=7)
            a_.tick_params(axis='y', labelsize=7)
            a_.grid(axis='y', color='0.92', lw=0.6)
            for s in ['top', 'right']:
                a_.spines[s].set_visible(False)
            if i == 0:
                a_.set_title(cols[j], fontsize=9)
        axs[i, 0].set_ylabel(LABEL[v], fontsize=8)
    axs[0, 0].legend(fontsize=7, frameon=False)
    axs[0, 3].legend(fontsize=7, frameon=False)
    fig.suptitle(f'CYGNSS L1 obs-level DA link, DA_L1_full_xc015_coh040216_fixedop vs OL, {tag} '
                 '(bars: tile-cluster bootstrap 95% CI; monitor = SMAP/L3/SMOS obs in the same tile within 72 h '
                 'after the L1 analysis, binned by the largest-|A-F| L1 obs)', fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(f'{out}obs_level_bins_{tag}.png', dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------ cross-checks
def crosscheck(d, M, out, tag):
    rows = []
    r = stats_from({c: d['s_' + c].sum() for c in L1_SUMS})
    rows.append(dict(grp='L1', N=len(d), **{k: r[k] for k in ['dstd', 'dmse', 'noise', 'gain', 'alpha']}))
    for g, gm in M.groupby('grp2'):
        r = stats_from({c: gm[c].sum() for c in MON_SUMS})
        rows.append(dict(grp=g, N=len(gm), **r))
    # per-tile OmF stdv, averaged over tiles (closer to the temporal_stats / report numbers)
    def tile_std(x):
        s = x.groupby('tile').agg(n=('n', 'sum'), a1=('a1', 'sum'), a2=('a2', 'sum'), b1=('b1', 'sum'), b2=('b2', 'sum'))
        s = s[s.n >= 10]
        sd = lambda s1, s2: np.sqrt(np.maximum(s2 / s.n - (s1 / s.n) ** 2, 0) * s.n / (s.n - 1))
        return 100 * (sd(s.a1, s.a2).mean() - sd(s.b1, s.b2).mean()) / sd(s.b1, s.b2).mean()
    dd = d.rename(columns=lambda c: c[2:] if c.startswith('s_') else c)
    rows[0]['dstd_tilemean'] = tile_std(dd)
    for i, (g, gm) in enumerate(M.groupby('grp2')):
        rows[i + 1]['dstd_tilemean'] = tile_std(gm)
    r = pd.DataFrame(rows)
    r.to_csv(f'{out}crosscheck_pooled_{tag}.csv', index=False)
    print('\n=== pooled cross-check (% (DA-OL)/OL; noise/gain/dMSE % of OL MSE) ===')
    print(r.to_string(index=False, float_format=lambda x: f'{x:9.3f}'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tag', required=True)
    p.add_argument('--inputs', required=True, help='comma list of <start>_<end> tags of extract outputs')
    p.add_argument('--data-dir', default=OUT0 + 'data/')
    p.add_argument('--boot', type=int, default=B_DEFAULT)
    p.add_argument('--seed', type=int, default=1)
    a = p.parse_args()
    out = f'{OUT0}{a.tag}/'
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    pd.set_option('display.width', 260); pd.set_option('display.max_columns', 40)
    L, M = load(a)
    print(f'{len(L)} L1 rows ({L.in_tab.sum()} OL table, {L.in_da.sum()} DA), {len(M)} monitor rows')
    q1(L, out)
    d = l1_matched(L)
    crosscheck(d, M, out, a.tag)

    b2 = [binned_l1(d, a.boot, rng)]
    for per in ['spring window', 'Oct-Jan']:
        b2.append(binned_l1(d[d.period == per], a.boot // 2, rng, per))
    b2 = pd.concat(b2, ignore_index=True)
    b2.to_csv(f'{out}q2_l1_by_bin.csv', index=False)
    show = ['subset', 'var', 'bin', 'x_med', 'N', 'N_tiles', 'dstd', 'dstd_lo', 'dstd_hi', 'dmse', 'dmse_lo', 'dmse_hi',
            'noise', 'gain', 'alpha', 'alpha_lo', 'alpha_hi', 'mean_abs_inc', 'K_eff', 'rms_ratio', 'med_ratio_OmA_OmF',
            'rms_ninn']
    print('\n=== Q2: L1 own fit / analysis draw / noise-gain by operator bin ===')
    print(b2[b2.subset == 'all'][show].to_string(index=False, float_format=lambda x: f'{x:7.3f}'))

    r3, M2 = monitor_link(M, d, a.boot, rng, out)
    show3 = ['grp', 'rule', 'var', 'bin', 'x_med', 'N', 'N_tiles', 'dstd', 'dstd_lo', 'dstd_hi', 'dmse', 'dmse_lo',
             'dmse_hi', 'noise', 'gain', 'alpha', 'alpha_lo', 'alpha_hi']
    print('\n=== Q3: monitor noise/gain by operator class of the preceding L1 analysis ===')
    sel = r3[r3.grp.isin(['SMAP', 'L3', 'SMOS']) & ((r3.rule == 'dominant72') |
                                                     r3['var'].isin(['weakH', 'tile_smskill', 'fcst_ol_tab']))]
    print(sel[show3].to_string(index=False, float_format=lambda x: f'{x:7.3f}'))
    figure(b2, r3, out, a.tag)
    print('wrote', out)


if __name__ == '__main__':
    main()
