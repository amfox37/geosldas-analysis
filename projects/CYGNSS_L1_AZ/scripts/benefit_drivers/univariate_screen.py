#!/usr/bin/env python3
"""EXPLORATORY univariate screening of L1-DA benefit outcomes vs tile covariates (L1 mask only).

Spearman rho, 2x2 deg block-bootstrap 95% CI (2000 draws; ranks computed once on each pair's valid sample,
then weighted Pearson of ranks per draw), BH-FDR within outcome (block-bootstrap p primary, naive p secondary),
quintile-binned means with block-bootstrap CIs for top covariates per family, heatmap.
"""
import re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / 'output' / 'benefit_drivers'
OUT = ROOT / 'analysis_univariate'; OUT.mkdir(exist_ok=True)
B, MIN_N = 2000, 30
rng = np.random.default_rng(20261119)

d = pd.read_csv(ROOT / 'analysis_table.csv')
d = d[d.in_l1_mask.astype(bool)].reset_index(drop=True)
T = len(d)
blk = (np.floor(d.com_lon / 2).astype(int).astype(str) + '_' + np.floor(d.com_lat / 2).astype(int).astype(str))
bid = pd.factorize(blk)[0]; NB = bid.max() + 1
draw_blocks = rng.integers(0, NB, size=(B, NB))
BC = np.zeros((B, NB)); [np.add.at(BC[i], draw_blocks[i], 1) for i in range(B)]
W = BC[:, bid]                        # B x T tile weights per draw
print('tiles', T, 'blocks', NB)

# ---------- outcomes
def add_diff(name, a, b):
    d[name] = d[a] - d[b]
for m in ['omf_indep', 'omf_smap', 'omf_smos', 'omf_ascat', 'iv_smapL3', 'iv_ascatH121']:
    add_diff(f'out_{m}_L1mL3', f'out_{m}', f'out_{m}_L3')
d['out_tc_dR2_L1mL3'] = d['out_tc_dR2_L1_minus_L3']
UNITS = {'omf': 'pct O-F stddev reduction', 'iv': 'delta corr (R)', 'tc': 'delta R2', 'comb': 'z-score units'}
def unit(o): return UNITS['omf'] if 'omf' in o else UNITS['iv'] if '_iv_' in o else UNITS['tc'] if 'tc_' in o else UNITS['comb']
FAMILY = {'out_omf_indep': 'SMAP', 'out_omf_smap': 'SMAP', 'out_omf_smos': 'SMAP', 'out_iv_smapL3': 'SMAP',
          'out_omf_ascat': 'ASCAT', 'out_iv_ascatH121': 'ASCAT', 'out_tc_dR2': 'ASCAT', 'out_combined': 'combined'}
OUTCOMES = list(FAMILY) + [f'out_{m}_L1mL3' for m in ['omf_indep', 'omf_smap', 'omf_smos', 'omf_ascat', 'iv_smapL3', 'iv_ascatH121', 'tc_dR2']]
for o in OUTCOMES:
    FAMILY.setdefault(o, FAMILY[o.replace('_L1mL3', '')] + ' (L1-L3)')

# ---------- covariates
CATEG = {'cv_stat_modis_landcover_mode', 'cv_stat_mwrtm_vegcls', 'cv_stat_mwrtm_soilcls', 'op_landcover',
         'ig_dominant_class', 'st_mosaic_primary_veg_type', 'st_mosaic_secondary_veg_type'}
POST = [c for c in d if c.startswith('cv_sn_da_')]          # post-treatment (computed from the DA arm itself)
def group(c):
    n = c.split('_', 1)[1] if c[:3] in ('cv_', 'op_', 'st_', 'vg_', 'cl_', 'sm_', 'ig_') else c
    if c.startswith('ig_'): return 'landcover (IGBP)'
    if c.startswith('vg_') or c.startswith('cv_stat_opacity') or c.startswith('cv_stat_veg_geom') or c.startswith('cv_stat_attenuation') or c == 'op_opac_med' or 'canopy' in c: return 'vegetation/canopy'
    if c.startswith('cl_'): return 'climate (precip/aridity)'
    if c.startswith('sm_'): return 'OL soil-moisture climate'
    if re.search(r'elev|dem|slope_obs|slope_med|alt_med|cti|geo_area|tile_area|frac_atm|bcs_mean', c) or c in ('op_slope_med',): return 'terrain/topography'
    if 'z0' in c: return 'roughness (ASCAT z0)'
    if c.startswith('st_soil') or re.search(r'sand|clay|porosity|peat', c): return 'soil texture'
    if re.search(r'inc_angle|rs_inc|inc_slope|inc_med|pr_geom_inc', c): return 'incidence angle'
    if re.search(r'rs_geom|sd_geom|std_ratio_geom|brcs', c): return 'geometry term'
    if c.startswith('cv_sup_'): return 'obs supply/quality'
    if c.startswith('cv_sn_') or re.search(r'op_(snr|R_impl|inn_std|dyn|bias)|sfmc_(mean|std)@', c): return 'signal/noise'
    if c.startswith('op_') or c.startswith('cv_op_'): return 'operator skill'
    return 'other'
covs = []
for c in d.columns:
    if not re.match(r'(cv|op|st|vg|cl|sm|ig)_', c) or c in CATEG or c in POST: continue
    if not pd.api.types.is_numeric_dtype(d[c]): continue
    x = d[c]
    if x.notna().sum() < MIN_N or x.nunique() <= 4 or x.value_counts(normalize=True).iloc[0] > 0.9: continue
    covs.append(c)
GRP = {c: group(c) for c in covs}
pd.Series(GRP, name='group').rename_axis('covariate').to_csv(OUT / 'univariate_covariate_groups.csv')
print('covariates', len(covs), pd.Series(GRP).value_counts().to_dict())

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

def wrho(rx, ry, Wm):
    s = Wm.sum(1)
    mx = Wm @ rx / s; my = Wm @ ry / s
    cxy = Wm @ (rx * ry) / s - mx * my
    vx = Wm @ rx**2 / s - mx**2; vy = Wm @ ry**2 / s - my**2
    return cxy / np.sqrt(vx * vy)

rows = []
for o in OUTCOMES:
    y = d[o].values
    for c in covs:
        x = d[c].values; m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < MIN_N: continue
        rx, ry = stats.rankdata(x[m]), stats.rankdata(y[m])
        r = np.corrcoef(rx, ry)[0, 1]
        bs = wrho(rx, ry, W[:, m]); bs = bs[np.isfinite(bs)]
        pb = min(1, 2 * (min((bs <= 0).sum(), (bs >= 0).sum()) + 1) / (len(bs) + 1))
        pn = stats.spearmanr(x[m], y[m]).pvalue
        # effect size: outcome mean difference top vs bottom quintile of covariate
        q = pd.qcut(pd.Series(x[m]).rank(method='first'), 5, labels=False).values
        eff = y[m][q == 4].mean() - y[m][q == 0].mean()
        rows.append(dict(outcome=o, family=FAMILY[o], units=unit(o), covariate=c, group=GRP[c], n=int(m.sum()), rho=r,
                         ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5), p_block=pb, p_naive=pn,
                         effect_Q5_minus_Q1=eff, outcome_sd=np.std(y[m], ddof=1)))
R = pd.DataFrame(rows)
R['fdr_block'] = R.groupby('outcome').p_block.transform(lambda s: bh(s.values))
R['fdr_naive'] = R.groupby('outcome').p_naive.transform(lambda s: bh(s.values))
R['ci_excl0'] = (R.ci_lo > 0) | (R.ci_hi < 0)
R['abs_rho'] = R.rho.abs()
R = R.sort_values(['outcome', 'abs_rho'], ascending=[True, False])
R.drop(columns='abs_rho').to_csv(OUT / 'univariate_spearman_all.csv', index=False)
print(R.groupby('outcome').agg(ntests=('rho', 'size'), n_fdr05=('fdr_block', lambda s: (s < .05).sum()), n_naive05=('fdr_naive', lambda s: (s < .05).sum()),
                                max_abs=('abs_rho', 'max')).round(3).to_string())

# ---------- top covariates per family (greedy, drop collinear |rho|>0.9 with already chosen)
cm = d[covs].corr(method='spearman')
def top_for(prim, k=12):
    s = R[R.outcome.isin(prim)].groupby('covariate').agg(mr=('abs_rho', 'mean'), mn=('n', 'min'), nout=('rho', 'size'))
    s = s[(s.nout == len(prim)) & (s.mn >= 100)].sort_values('mr', ascending=False)
    pick = []
    for c in s.index:
        if all(abs(cm.loc[c, p]) < 0.9 for p in pick): pick.append(c)
        if len(pick) == k: break
    return pick
FAMS = {'SMAP': ['out_omf_indep', 'out_iv_smapL3', 'out_omf_smap', 'out_omf_smos'],
        'ASCAT': ['out_omf_ascat', 'out_iv_ascatH121', 'out_tc_dR2']}
HEAD = {'SMAP': ['out_omf_indep', 'out_iv_smapL3'], 'ASCAT': ['out_omf_ascat', 'out_iv_ascatH121', 'out_tc_dR2']}
TOP = {f: top_for(FAMS[f]) for f in FAMS}
pd.DataFrame([(f, i + 1, c, GRP[c]) for f, l in TOP.items() for i, c in enumerate(l)], columns=['family', 'rank', 'covariate', 'group']).to_csv(OUT / 'univariate_top_covariates_by_family.csv', index=False)

# ---------- quintile bins
def binned(x, y, nb=5):
    m = np.isfinite(x) & np.isfinite(y); xx, yy, Wm, bb = x[m], y[m], W[:, m], bid[m]
    q = pd.qcut(pd.Series(xx).rank(method='first'), nb, labels=False).values
    out = []
    for k in range(nb):
        s = q == k; w = Wm[:, s]; ok = w.sum(1) > 0
        bm = (w[ok] @ yy[s]) / w[ok].sum(1)
        out.append(dict(bin=k + 1, n=int(s.sum()), cov_mean=xx[s].mean(), cov_lo=xx[s].min(), cov_hi=xx[s].max(), mean=yy[s].mean(),
                        ci_lo=np.percentile(bm, 2.5), ci_hi=np.percentile(bm, 97.5)))
    return pd.DataFrame(out)
brows = []
for f, cl in TOP.items():
    for o in FAMS[f]:
        for c in cl:
            b = binned(d[c].values, d[o].values); b.insert(0, 'covariate', c); b.insert(0, 'outcome', o); b.insert(0, 'family', f); brows.append(b)
pd.concat(brows).to_csv(OUT / 'univariate_binned_means_quintile.csv', index=False)
B_ = pd.concat(brows)
for f in FAMS:
    for o in HEAD[f]:
        fig, axs = plt.subplots(3, 4, figsize=(16, 10)); axs = axs.ravel()
        for ax, c in zip(axs, TOP[f]):
            b = B_[(B_.outcome == o) & (B_.covariate == c)]
            ax.errorbar(b['bin'], b['mean'], yerr=[b['mean'] - b.ci_lo, b.ci_hi - b['mean']], fmt='o-', capsize=3, color='C0')
            ax.axhline(0, color='k', lw=.6)
            rr = R[(R.outcome == o) & (R.covariate == c)].iloc[0]
            ax.set_title(f'{c}\nrho={rr.rho:+.2f} [{rr.ci_lo:+.2f},{rr.ci_hi:+.2f}] n={rr.n}', fontsize=8)
            ax.set_xlabel('covariate quintile (1=low)'); ax.set_ylabel(unit(o), fontsize=7)
        fig.suptitle(f'EXPLORATORY: {o} by covariate quintile (block-bootstrap 95% CI); top covariates of {f} family', fontsize=10)
        fig.tight_layout(); fig.savefig(OUT / f'univariate_binned_{o}.png', dpi=110); plt.close(fig)

# ---------- heatmaps (median rho and strongest rho per group)
grp_order = sorted(set(GRP.values()))
outs = [o for o in OUTCOMES]
med = pd.DataFrame(index=outs, columns=grp_order, dtype=float); strong = med.copy(); sig = med.copy().astype(object)
for o in outs:
    ro = R[R.outcome == o]
    for g in grp_order:
        rg = ro[ro.group == g]
        if len(rg):
            med.loc[o, g] = rg.rho.median(); j = rg.rho.abs().idxmax(); strong.loc[o, g] = rg.loc[j, 'rho']
            sig.loc[o, g] = '*' if (rg.fdr_block < .05).any() else ''
med.to_csv(OUT / 'univariate_heatmap_median_rho_by_group.csv'); strong.to_csv(OUT / 'univariate_heatmap_strongest_rho_by_group.csv')
fig, axs = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
for ax, M, t in zip(axs, [med, strong], ['median rho within group', 'strongest |rho| covariate in group (* = FDR<0.05 block)']):
    im = ax.imshow(M.values.astype(float), cmap='RdBu_r', vmin=-.6, vmax=.6, aspect='auto')
    ax.set_xticks(range(len(grp_order))); ax.set_xticklabels(grp_order, rotation=50, ha='right', fontsize=8)
    ax.set_yticks(range(len(outs))); ax.set_yticklabels(outs, fontsize=8); ax.set_title(t, fontsize=10)
    for i in range(len(outs)):
        for j in range(len(grp_order)):
            v = M.values[i, j]
            if np.isfinite(v): ax.text(j, i, f'{v:+.2f}' + (sig.values[i, j] if M is strong else ''), ha='center', va='center', fontsize=6)
fig.colorbar(im, ax=axs, shrink=.7, label='Spearman rho')
fig.suptitle('EXPLORATORY: Spearman rho, outcomes x covariate group (groups internally collinear)', fontsize=11)
fig.savefig(OUT / 'univariate_heatmap_rho_outcome_by_group.png', dpi=110, bbox_inches='tight'); plt.close(fig)

# ---------- focus list (terrain / vegetation / incidence / roughness)
FOC = [c for c in covs if GRP[c] in ('terrain/topography', 'vegetation/canopy', 'incidence angle', 'roughness (ASCAT z0)', 'geometry term')]
Fo = R[R.covariate.isin(FOC) & R.outcome.isin(list(FAMILY)[:7] + ['out_combined'])]
Fo.pivot(index='covariate', columns='outcome', values='rho').assign(group=lambda x: x.index.map(GRP)).round(3).to_csv(OUT / 'univariate_focus_terrain_veg_incidence_rho.csv')

# ---------- family consistency: covariate rho in SMAP vs ASCAT family
a = R[R.outcome == 'out_omf_indep'].set_index('covariate').rho; b = R[R.outcome == 'out_omf_ascat'].set_index('covariate').rho
cc = pd.concat([a, b], axis=1, keys=['smap', 'ascat']).dropna()
print('rho-profile agreement indep vs ascat (Spearman across covariates):', stats.spearmanr(cc.smap, cc.ascat)[0].round(2), len(cc))
for f in FAMS:
    print('\nTOP', f); 
    for o in HEAD[f]:
        print(o); print(R[(R.outcome == o)].head(10)[['covariate', 'group', 'n', 'rho', 'ci_lo', 'ci_hi', 'fdr_block', 'fdr_naive', 'effect_Q5_minus_Q1']].round(3).to_string(index=False))
    print(TOP[f])
