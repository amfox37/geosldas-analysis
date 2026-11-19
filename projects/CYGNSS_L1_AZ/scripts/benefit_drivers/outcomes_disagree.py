"""Task 3: characterise family-disagreement tiles (EXPLORATORY; descriptive, no multivariate modelling)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from outcomes_common import *
from scipy import stats
a = load(); s = pd.read_csv(OUT / 'combined_scores_per_tile.csv')[['tile_axis_index', 'fam_cat', 'fam_cat_noTC', 'V4_allpos', 'Fsc_smap', 'Fsc_ascat_noTC']]
a = a.merge(s, on='tile_axis_index'); fc = pd.read_csv(OUT / 'four_corners_flags.csv'); a = a.merge(fc[['tile_axis_index', 'fc_box', 'fc_hot']], on='tile_axis_index')
a['cls'] = a.fam_cat_noTC
# na split: missing ASCAT vs missing otherwise
m = a[a.in_l1_mask.astype(bool)].copy()
print(m.cls.value_counts().to_dict())
COVS = {'n_l1_obs': 'n_l1_obs', 'op_rs_smterm': 'op_rs_smterm', 'elev_m': 'st_tilecoord_elev_m', 'aridity_P_PET': 'cl_aridity_index_P_over_PET', 'sm_sfmc_std': 'sm_sfmc_std'}
IG = {1: 'ENF', 2: 'EBF', 3: 'DNF', 4: 'DBF', 5: 'MF', 6: 'CShrub', 7: 'OShrub', 8: 'WSav', 9: 'Sav', 10: 'Grass', 11: 'Wetl', 12: 'Crop', 13: 'Urban', 14: 'CropMos', 15: 'Snow', 16: 'Barren', 17: 'Water'}
order = ['both_better', 'smap_only', 'ascat_only', 'both_worse', 'na']
rows = []
for c in order:
    sub = m[m.cls == c]
    for k, col in COVS.items():
        est = boot_stat(sub, lambda d, col=col: d[col].median(), 1000, seed=3)
        # difference from rest (within the 552 mask)
        rest = m[m.cls != c]
        rows.append((c, k, len(sub), *est, rest[col].median()))
tab = pd.DataFrame(rows, columns=['class', 'covariate', 'n', 'median', 'ci_lo', 'ci_hi', 'median_rest_of_mask']).round(3)
tab.to_csv(OUT / 'disagree_covariate_medians.csv', index=False); print(tab.to_string())
# class-vs-rest median differences with block-bootstrap CI
def dd(df, c, col):
    return df.loc[df.cls == c, col].median() - df.loc[df.cls != c, col].median()
rows = []
for c in order[:4]:
    for k, col in COVS.items():
        rows.append((c, k, *boot_stat(m, lambda d, c=c, col=col: dd(d, c, col), 1000, seed=4)))
dif = pd.DataFrame(rows, columns=['class', 'covariate', 'median_diff_vs_rest', 'lo', 'hi']).round(3); dif['excludes0'] = (dif.lo > 0) | (dif.hi < 0)
dif.to_csv(OUT / 'disagree_covariate_median_diff_vs_rest.csv', index=False); print(dif.to_string())
# smap_only vs ascat_only directly
rows = []
for k, col in COVS.items():
    f = lambda d, col=col: d.loc[d.cls == 'smap_only', col].median() - d.loc[d.cls == 'ascat_only', col].median()
    rows.append((k, *boot_stat(m, f, 1000, seed=6)))
print(pd.DataFrame(rows, columns=['cov', 'smap_only_minus_ascat_only', 'lo', 'hi']).round(3).to_string())
# IGBP and Four Corners composition
m['igbp'] = m.ig_dominant_class.map(IG)
comp = pd.crosstab(m.cls, m.igbp).reindex(order); comp.to_csv(OUT / 'disagree_igbp_counts.csv'); print(comp.to_string())
fcr = m.groupby('cls')[['fc_box', 'fc_hot']].agg(['sum', 'mean']).reindex(order).round(3); fcr.to_csv(OUT / 'disagree_fourcorners_share.csv'); print(fcr.to_string())
# why na: which measures are missing
m['has_ascat'] = m[['out_iv_ascatH121', 'out_omf_ascat']].notna().any(axis=1); m['has_smap'] = m[SMAP_FAM].notna().any(axis=1)
print(pd.crosstab(m.cls, [m.has_smap, m.has_ascat]))
# measure-level families signs among classes: median of each outcome
print(m.groupby('cls')[INDEP].median().reindex(order).round(3).to_string()); m.groupby('cls')[INDEP].median().reindex(order).round(4).to_csv(OUT / 'disagree_outcome_medians_by_class.csv')
# L3 arm comparison: do the same tiles disagree for L3 DA?  (is the pattern arm-specific?)
for c in ['out_omf_ascat_L3', 'out_iv_ascatH121_L3', 'out_omf_smap_L3', 'out_omf_smos_L3', 'out_iv_smapL3_L3']:
    print(c, m.groupby('cls')[c].median().reindex(order[:4]).round(3).to_dict())
# maps
cat = {'both_better': ('tab:blue', ''), 'smap_only': ('tab:green', ''), 'ascat_only': ('tab:orange', ''), 'both_worse': ('tab:red', '')}
fig, axs = mapfig(2, ncols=2, size=(5.5, 4.8))
tile_map(axs[0], a, 'cls', 'family agreement (SMAP fam: IV L3, O-F SMAP, O-F SMOS; ASCAT fam: IV ASCAT, O-F ASCAT)', cat=cat)
a['ssc'] = a.Fsc_smap - a.Fsc_ascat_noTC
tile_map(axs[1], a, 'ssc', 'SMAP-family minus ASCAT-family scale-only score (green/blue = SMAP side)', cmap='PiYG')
plt.tight_layout(); plt.savefig(OUT / 'fig08_disagreement_map.png', dpi=140); plt.close()
# boxplots of covariates
fig, axs = plt.subplots(1, 5, figsize=(16, 3.8))
for ax, (k, col) in zip(axs, COVS.items()):
    ax.boxplot([m.loc[m.cls == c, col].dropna() for c in order], labels=[f'{c}\n{(m.cls == c).sum()}' for c in order], showfliers=False)
    ax.set_title(k, fontsize=9); ax.tick_params(axis='x', labelsize=6, rotation=40)
plt.tight_layout(); plt.savefig(OUT / 'fig09_disagree_covariates.png', dpi=140); plt.close()
