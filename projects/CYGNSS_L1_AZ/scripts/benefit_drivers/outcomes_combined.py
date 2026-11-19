"""Task 2: robustness of the combined score (EXPLORATORY)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from outcomes_common import *
from scipy import stats
a = load(); mask = a.in_l1_mask.astype(bool).values
def pr(x):  # centred percentile rank within mask
    r = pd.Series(np.where(mask, x, np.nan)).rank(); return ((r - .5) / r.notna().sum() - .5).values
def sc(x): return (x / np.nanstd(np.where(mask, x, np.nan))).values if hasattr(x, 'values') else x / np.nanstd(np.where(mask, x, np.nan))
for c in INDEP:
    a['r_' + c] = pr(a[c].values); a['s_' + c] = np.where(mask, a[c] / a.loc[mask, c].std(), np.nan)
def fam(prefix, cols, minn=1):
    X = a[[prefix + c for c in cols]]; v = X.mean(axis=1); v[X.notna().sum(axis=1) < minn] = np.nan; return v
a['Frank_smap'] = fam('r_', SMAP_FAM); a['Frank_ascat'] = fam('r_', ASCAT_FAM)
a['Fsc_smap'] = fam('s_', SMAP_FAM); a['Fsc_ascat'] = fam('s_', ASCAT_FAM)
a['Frank_ascat_noTC'] = fam('r_', ['out_iv_ascatH121', 'out_omf_ascat']); a['Fsc_ascat_noTC'] = fam('s_', ['out_iv_ascatH121', 'out_omf_ascat'])
a['V0_z_tc_indep'] = a.out_combined
a['V1_rank_tc_indep'] = (a.r_out_tc_dR2 + pr(a.out_omf_indep.values)) / 2
a['V2_famrank'] = (a.Frank_smap + a.Frank_ascat) / 2
a['V3_famscale'] = (a.Fsc_smap + a.Fsc_ascat) / 2
a['V2b_famrank_noTC'] = (a.Frank_smap + a.Frank_ascat_noTC) / 2
a['V3b_famscale_noTC'] = (a.Fsc_smap + a.Fsc_ascat_noTC) / 2
# conservative: all available measures positive (need >=1 per family and >=4 measures)
P = a[INDEP]; nav = P.notna().sum(axis=1); npos = (P > 0).sum(axis=1)
nav_s = P[SMAP_FAM].notna().sum(axis=1); nav_a = P[ASCAT_FAM].notna().sum(axis=1)
ok = mask & (nav >= 4) & (nav_s >= 1) & (nav_a >= 1)
a['V4_allpos'] = np.where(ok, np.where(npos == nav, 'all_pos', np.where(npos == 0, 'all_neg', 'mixed')), 'na')
a['n_meas_avail'] = nav; a['n_meas_pos'] = npos
# family categories from scale-only family scores (sign = better/worse)
def cat(s, t):
    out = np.full(len(a), 'na', dtype=object); v = s.notna() & t.notna()
    out[v & (s > 0) & (t > 0)] = 'both_better'; out[v & (s > 0) & (t <= 0)] = 'smap_only'; out[v & (s <= 0) & (t > 0)] = 'ascat_only'; out[v & (s <= 0) & (t <= 0)] = 'both_worse'
    return np.where(mask, out, 'na')
a['fam_cat'] = cat(a.Fsc_smap, a.Fsc_ascat); a['fam_cat_noTC'] = cat(a.Fsc_smap, a.Fsc_ascat_noTC)
m = a[mask].copy()
vers = ['V0_z_tc_indep', 'V1_rank_tc_indep', 'V2_famrank', 'V3_famscale', 'V2b_famrank_noTC', 'V3b_famscale_noTC', 'Frank_smap', 'Frank_ascat', 'Frank_ascat_noTC']
rows = [(v, m[v].notna().sum()) for v in vers]; print(rows)
# pairwise Spearman among versions with CIs (block bootstrap)
Rv = m[vers].corr(method='spearman'); Bv = boot_corr(m, vers, 500, seed=9)
lo, hi = np.nanpercentile(Bv, 2.5, 0), np.nanpercentile(Bv, 97.5, 0)
pd.DataFrame([(x, y, Rv.loc[x, y], lo[i, j], hi[i, j], m[[x, y]].dropna().shape[0]) for i, x in enumerate(vers) for j, y in enumerate(vers) if i < j],
             columns=['a', 'b', 'spearman', 'lo', 'hi', 'n']).round(3).to_csv(OUT / 'combined_versions_corr.csv', index=False)
print(Rv.round(2).to_string())
# TC-validity sensitivity
tcv = m.out_tc_dR2.notna(); fam_ok = m.Frank_smap.notna() & m.Frank_ascat_noTC.notna()
print('mask', len(m), 'TC valid', tcv.sum(), 'noTC-family score available', fam_ok.sum(), 'of which no TC', (fam_ok & ~tcv).sum())
both = m[tcv & m.V2_famrank.notna() & m.V2b_famrank_noTC.notna()]
print('on TC-valid tiles: V2 vs V2b noTC rho', stats.spearmanr(both.V2_famrank, both.V2b_famrank_noTC)[0], len(both))
ct = pd.crosstab(m.fam_cat, m.fam_cat_noTC); print(ct)
# agreement between categories: kappa
def kappa(x, y):
    t = pd.crosstab(x, y); n = t.values.sum(); po = np.trace(t.values) / n if t.shape[0] == t.shape[1] else np.nan
    if np.isnan(po): idx = t.index.intersection(t.columns); po = sum(t.loc[i, i] for i in idx) / n; pe = sum(t.loc[i].sum() * t[i].sum() for i in idx) / n ** 2
    else: pe = (t.sum(1).values * t.sum(0).values).sum() / n ** 2
    return (po - pe) / (1 - pe), po
v = m[(m.fam_cat != 'na') & (m.fam_cat_noTC != 'na') & m.out_tc_dR2.notna()]
print('kappa fam_cat vs fam_cat_noTC (TC-valid tiles)', kappa(v.fam_cat, v.fam_cat_noTC), len(v))
v2 = m[(m.fam_cat != 'na') & (m.out_agree != 'na')]
print('crosstab out_agree vs fam_cat'); print(pd.crosstab(v2.out_agree, v2.fam_cat))
pd.crosstab(m.out_agree, m.fam_cat).to_csv(OUT / 'crosstab_outagree_vs_famcat.csv')
pd.crosstab(m.fam_cat, m.V4_allpos).to_csv(OUT / 'crosstab_famcat_vs_allpositive.csv'); print(pd.crosstab(m.fam_cat, m.V4_allpos))
# measure-level sign disagreement counts
print('family cat counts', m.fam_cat.value_counts().to_dict(), '| noTC', m.fam_cat_noTC.value_counts().to_dict(), '| V4', m.V4_allpos.value_counts().to_dict())
# family scores correlation overall + CI
d = m[['Frank_smap', 'Frank_ascat']].dropna(); est = stats.spearmanr(d.Frank_smap, d.Frank_ascat)[0]
f = lambda df: stats.spearmanr(df.Frank_smap, df.Frank_ascat, nan_policy='omit')[0]
fr = boot_stat(m, f, 1000); print('family score rho', fr, len(d))
dn = m[['Frank_smap', 'Frank_ascat_noTC']].dropna(); print('no-TC family rho', stats.spearmanr(dn.iloc[:, 0], dn.iloc[:, 1])[0], len(dn))
# chance level of 4-class agreement given family sign: expected if independent
pS = (m.Fsc_smap.dropna() > 0).mean(); pA = (m.Fsc_ascat.dropna() > 0).mean()
print('P(SMAP fam>0)', pS, 'P(ASCAT fam>0)', pA, 'expected both_better if indep', pS * pA, 'observed', (m.fam_cat == 'both_better').sum() / (m.fam_cat != 'na').sum())
# --- split-half spatial reliability
rng = np.random.default_rng(11)
lon, lat = m.com_lon.values, m.com_lat.values; blocks = m.block.values; ub = np.unique(blocks)
D = np.hypot((lon[:, None] - lon[None]) * np.cos(np.radians(34)), lat[:, None] - lat[None])
comps = ['out_tc_dR2', 'out_iv_ascatH121', 'out_omf_ascat', 'out_iv_smapL3', 'out_omf_smap', 'out_omf_smos'] + ['V0_z_tc_indep', 'V1_rank_tc_indep', 'V2_famrank', 'V3_famscale', 'V2b_famrank_noTC', 'Frank_smap', 'Frank_ascat', 'Frank_ascat_noTC']
res = {c: [[], []] for c in comps}
for it in range(200):
    hb = rng.permutation(ub)[:len(ub) // 2]; inA = np.isin(blocks, hb)
    for c in comps:
        x = m[c].values
        for k, (dmin, dmax) in enumerate([(0.0, 0.9), (0.6, 2.2)]):  # near neighbours incl adjacent blocks / gapped
            nbr = []
            ia = np.where(inA & np.isfinite(x))[0]; ib = np.where(~inA & np.isfinite(x))[0]
            sel = (D[np.ix_(ia, ib)] > dmin) & (D[np.ix_(ia, ib)] <= dmax)
            xb = x[ib]; w = sel.astype(float); cnt = w.sum(1); ok_ = cnt >= 2
            pred = (w @ xb)[ok_] / cnt[ok_]
            if ok_.sum() > 15: res[c][k].append(stats.spearmanr(x[ia][ok_], pred)[0])
rel = pd.DataFrame([(c, np.nanmedian(res[c][0]), np.nanpercentile(res[c][0], 5), np.nanpercentile(res[c][0], 95), np.nanmedian(res[c][1]), np.nanpercentile(res[c][1], 5), np.nanpercentile(res[c][1], 95)) for c in comps],
                   columns=['score', 'near_median', 'near_p5', 'near_p95', 'gap_median', 'gap_p5', 'gap_p95']).round(3)
rel.to_csv(OUT / 'split_half_spatial_reliability.csv', index=False); print(rel.to_string())
# year-split reliability for O-F measures (same tile, different years)
yr = []
for g_ in ['smap', 'smos', 'ascat', 'l3obs']:
    for y1, y2 in [(2020, 2021), (2020, 2022), (2021, 2022)]:
        d = m[[f'out_omf_{g_}_{y1}', f'out_omf_{g_}_{y2}']].dropna(); yr.append((g_, y1, y2, stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0], len(d)))
yr = pd.DataFrame(yr, columns=['group', 'y1', 'y2', 'spearman', 'n']).round(3); yr.to_csv(OUT / 'year_split_reliability_omf.csv', index=False); print(yr.to_string())
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
x = np.arange(len(rel)); ax[0].bar(x - .2, rel.near_median, .4, label='near (<0.9 deg)'); ax[0].bar(x + .2, rel.gap_median, .4, label='gapped (0.6-2.2 deg)')
ax[0].set_xticks(x); ax[0].set_xticklabels(rel.score, rotation=70, ha='right', fontsize=7); ax[0].set_ylabel('Spearman tile vs mean of other-half neighbours (median of 200 splits)'); ax[0].legend(fontsize=7); ax[0].set_title('split-half spatial reliability (2x2-deg block halves)', fontsize=8)
ax[1].scatter(m.Frank_smap, m.Frank_ascat, s=8, c=np.where(m.out_tc_dR2.notna(), 'C0', 'C1')); ax[1].set_xlabel('SMAP-family rank score'); ax[1].set_ylabel('ASCAT-family rank score'); ax[1].axhline(0, c='.7'); ax[1].axvline(0, c='.7')
ax[1].set_title(f'family scores: rho={est:.2f} [{fr[1]:.2f},{fr[2]:.2f}] n={len(d)}', fontsize=8)
plt.tight_layout(); plt.savefig(OUT / 'fig06_combined_reliability.png', dpi=150); plt.close()
# maps of main combined variants
fig, axs = mapfig(6)
for ax, c in zip(axs, ['V0_z_tc_indep', 'V2_famrank', 'V2b_famrank_noTC', 'Frank_smap', 'Frank_ascat', 'V3_famscale']): tile_map(ax, a, c, c)
plt.tight_layout(); plt.savefig(OUT / 'fig07_maps_combined.png', dpi=130); plt.close()
a[['tile_axis_index', 'com_lon', 'com_lat', 'in_l1_mask', 'Frank_smap', 'Frank_ascat', 'Fsc_smap', 'Fsc_ascat', 'Frank_ascat_noTC', 'Fsc_ascat_noTC', 'V0_z_tc_indep', 'V1_rank_tc_indep', 'V2_famrank', 'V3_famscale', 'V2b_famrank_noTC', 'V4_allpos', 'n_meas_avail', 'n_meas_pos', 'fam_cat', 'fam_cat_noTC', 'out_agree']].round(4).to_csv(OUT / 'combined_scores_per_tile.csv', index=False)
