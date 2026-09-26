#!/usr/bin/env python3
"""Analyze CYGNSS L1 close-pair innovation correlations (output of build_l1_pairs.py).

Innovation d = scaled O - F from the unscaled fixed-operator OL (OLv8_M36_AZ_fixedop), with the obs
scaled offline by the new L1 z-score clim. Pairs = two L1 obs in the same analysis cycle within 0.6 deg.
Correlations are Pearson over symmetrized pairs; CIs are day-block bootstrap (pairs within a day are not independent).
"""
import argparse
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--pairs', default='/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/'
               'output/obs_error_correlation/l1_pairs_3yr.parquet')
p.add_argument('--nboot', type=int, default=200)
args = p.parse_args()

df = pd.read_parquet(args.pairs)
df['day'] = df.t.str[:8]
df['month'] = df.t.str[4:6].astype(int)
df['acq'] = np.where(df.same_track, 'same-track', np.where(df.same_sc, 'same-sc other', 'different-sc'))
rng = np.random.default_rng(0)


def corr(x):
    if len(x) < 30:
        return np.nan
    a = np.r_[x.d1.values, x.d2.values]
    b = np.r_[x.d2.values, x.d1.values]
    return np.corrcoef(a, b)[0, 1]


def boot(x, nboot=args.nboot):
    """day-block bootstrap 90% CI of corr"""
    if len(x) < 30:
        return np.nan, np.nan
    days = x.day.unique()
    g = {d: v for d, v in x.groupby('day')}
    out = []
    for _ in range(nboot):
        s = pd.concat([g[d] for d in rng.choice(days, len(days))])
        out.append(corr(s))
    return np.nanpercentile(out, 5), np.nanpercentile(out, 95)


def row(x, ci=True):
    c = corr(x)
    lo, hi = boot(x) if ci else (np.nan, np.nan)
    return f'{len(x):8d}  {c:6.3f}  [{lo:5.3f}, {hi:5.3f}]  ovl {x.overlap.median():5.3f}'


dist_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6]
df['dbin'] = pd.cut(df.dist, dist_bins)
ivar = np.var(np.r_[df.d1, df.d2])

print('=== 0. dataset')
print(f'pairs {len(df):,}; days {df.day.nunique()}; innovation var {ivar:.2f} dB^2 (std {np.sqrt(ivar):.2f}); mean {np.mean(np.r_[df.d1, df.d2]):+.3f}')
print(df.acq.value_counts().to_string())

print('\n=== 1. corr vs distance, all pairs (N, corr, 90% CI, median overlap)')
for k, g in df.groupby('dbin', observed=True):
    print(f'{str(k):13s} {row(g)}  same-track {100*g.same_track.mean():5.1f}%')

print('\n=== 2. corr vs distance by acquisition type')
for acq in ['same-track', 'same-sc other', 'different-sc']:
    print(f'-- {acq}')
    for k, g in df[df.acq == acq].groupby('dbin', observed=True):
        print(f'   {str(k):13s} {row(g)}')

print('\n=== 3. same-track: corr vs time separation |dt| (s)')
st = df[df.same_track]
for k, g in st.groupby(pd.cut(st.dt_s, [0, 1.5, 2.5, 3.5, 5.5, 8.5, 12.5, 30]), observed=True):
    print(f'   dt {str(k):12s} {row(g)}  median dist {g.dist.median():.3f} deg')

print('\n=== 4. footprint overlap at fixed distance, split by acquisition (corr; N)')
ob = [-1e-9, 0.1, 0.3, 0.6, 0.8, 1.01]
for acq in ['same-track', 'different-sc']:
    x = df[(df.acq == acq) & (df.dist < 0.3)]
    t = x.groupby(['dbin', pd.cut(x.overlap, ob)], observed=True).apply(lambda g: f'{corr(g):.3f} ({len(g)})').unstack()
    print(f'-- {acq}\n{t.to_string()}')

print('\n=== 5. different-sc pairs within 0.15 deg: incidence-angle difference')
x = df[(df.acq == 'different-sc') & (df.dist < 0.15)]
for k, g in x.groupby(pd.cut(x.dinc, [0, 3, 6, 10, 15, 25, 90]), observed=True):
    print(f'   d_inc {str(k):10s} {row(g)}')
print('   (same, restricted to overlap > 0.6)')
x2 = x[x.overlap > 0.6]
for k, g in x2.groupby(pd.cut(x2.dinc, [0, 5, 15, 90]), observed=True):
    print(f'   d_inc {str(k):10s} {row(g, ci=False)}')

print('\n=== 6. season: corr at d<0.1 and 0.3-0.6 by season')
seas = {12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'}
df['season'] = df.month.map(seas)
for s in ['DJF', 'MAM', 'JJA', 'SON']:
    a = df[(df.season == s) & (df.dist < 0.1)]
    b = df[(df.season == s) & (df.dist >= 0.3)]
    a_st = a[a.same_track]
    a_ds = a[a.acq == 'different-sc']
    print(f'   {s}: d<0.1 all {corr(a):.3f} (N {len(a)}), same-track {corr(a_st):.3f}, different-sc {corr(a_ds):.3f};  d 0.3-0.6 {corr(b):.3f}')

print('\n=== 7. Hollingsworth-Loennberg style split (different-sc pairs only)')
x = df[df.acq == 'different-sc']
bins = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
xs, ys = [], []
for k, g in x.groupby(pd.cut(x.dist, bins), observed=True):
    xs.append(g.dist.median())
    ys.append(corr(g))
xs, ys = np.array(xs), np.array(ys)
m = xs >= 0.15
# fit c(d) = b * exp(-d/L) on d>=0.15 (assumed free of obs-error correlation); b = forecast-error share at d=0
A = np.polyfit(xs[m], np.log(ys[m]), 1)
L = -1 / A[0]
b = np.exp(A[1])
print('   bins (median d, corr):', ' '.join(f'({a:.3f},{c:.3f})' for a, c in zip(xs, ys)))
print(f'   fit c(d)=b*exp(-d/L) on d>=0.15: b={b:.3f} (forecast-error share of innovation var at d=0), L={L:.2f} deg')
print(f'   => obs-error share of innovation variance ~ {1-b:.3f}; ensemble-based estimate (Pf~1 dB^2 / {ivar:.1f}) obs share ~ {1-1.0/ivar:.3f}')
for acq in ['same-track', 'different-sc']:
    for lo, hi in [(0, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.3)]:
        g = df[(df.acq == acq) & (df.dist > lo) & (df.dist <= hi)]
        if len(g) < 30:
            continue
        d = g.dist.median()
        cf = b * np.exp(-d / L)
        rho_o = (corr(g) - cf) / (1 - b)
        print(f'   {acq:13s} d {lo:.2f}-{hi:.2f} (med {d:.3f}): innov corr {corr(g):.3f}, fcst part {cf:.3f} -> implied obs-err corr {rho_o:.3f}')

print('\n=== 8. R models vs observed innovation corr (all pairs), with obs share a=1-b')
a = 1 - b
print('   dist   observed | Gaussian 0.625  0.25   0.15   0.10 | 0.15+nugget0.5 | (model innov corr = a*rho_R + b*exp(-d/L))')
for k, g in df.groupby('dbin', observed=True):
    d = g.dist.median()
    cf = b * np.exp(-d / L)
    mods = [a * np.exp(-0.5 * (d / xc) ** 2) + cf for xc in (0.625, 0.25, 0.15, 0.10)]
    nug = a * 0.5 * np.exp(-0.5 * (d / 0.15) ** 2) + cf
    print(f'   {d:.3f}  {corr(g):6.3f}  | ' + '  '.join(f'{v:5.3f}' for v in mods) + f' | {nug:5.3f}')
