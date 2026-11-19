"""EXPLORATORY multivariate attribution of L1 DA benefit (spatial block CV).
Run: python -I multivariate_attribution.py
Outputs -> output/benefit_drivers/analysis_multivariate/
"""
import re, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import GroupKFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import GradientBoostingRegressor
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[2]
BD = ROOT / 'output/benefit_drivers'
OUT = BD / 'analysis_multivariate'; OUT.mkdir(exist_ok=True, parents=True)
NREP, NIMP, NFOLD, MAXMISS, RHO = 10, 5, 5, 0.30, 0.90

d = pd.read_csv(BD / 'analysis_table.csv', na_values=['na', 'NA', ''])
tc = pd.read_csv(BD / 'discover/tc_iv_per_tile.csv')
base_cols = ['tc_R2mod_OL', 'iv_R_smapL3_OL', 'iv_R_ascatH121_OL']
d = d.merge(tc[['tile_id'] + base_cols], on='tile_id', how='left')
d = d[d.in_l1_mask.astype(bool)].reset_index(drop=True)
print('mask tiles', len(d))

# ---- outcomes
for k in ['out_omf_indep', 'out_omf_ascat', 'out_iv_smapL3', 'out_iv_ascatH121']:
    d[k + '_L1mL3'] = d[k] - d[k + '_L3']
OUTCOMES = ['out_omf_indep', 'out_omf_ascat', 'out_iv_smapL3', 'out_iv_ascatH121', 'out_tc_dR2', 'out_combined',
            'out_omf_indep_L1mL3', 'out_omf_ascat_L1mL3', 'out_iv_smapL3_L1mL3', 'out_iv_ascatH121_L1mL3',
            'out_tc_dR2_L1_minus_L3']

# ---- candidate covariates and groups
CODE_COLS = {'cv_stat_modis_landcover_mode', 'cv_stat_mwrtm_vegcls', 'cv_stat_mwrtm_soilcls', 'ig_dominant_class',
             'st_mosaic_primary_veg_type', 'st_mosaic_secondary_veg_type', 'op_landcover'}
pref = ('cv_', 'ig_', 'op_', 'st_', 'vg_', 'cl_', 'sm_')
cands = [c for c in d.columns if c.startswith(pref) and c not in CODE_COLS
         and not c.startswith('cv_sn_da_l1_')]           # post-DA (outcome-linked) diagnostics excluded
cands += base_cols

def grp(c):
    if c in base_cols: return 'baseline_skill'
    if re.search(r'dem_|slope_obs|^st_cti|canopy_height|roughness_z0|elev|inc_slope|op_inc_med|coh_tail|frac_coh_pass|coh_ratio|inc_angle|op_slope_med|op_alt_med', c): return 'terrain_roughness'
    if re.match(r'(sm_|op_sd_|op_dyn|op_inn|op_sfmc|op_std_ratio)', c) or c.startswith('cv_sn_'): return 'signal_noise'
    if re.match(r'(cv_op_|op_r|op_pr_|op_bias|op_snr|op_R_impl)', c): return 'operator_skill'
    if re.match(r'(cv_sup_|op_N|op_n_eff|op_w_|op_frac_out)', c): return 'obs_supply'
    if re.match(r'(vg_|ig_|op_opac)', c) or re.search(r'opacity|veg_geometry|attenuation|brcs|mosaic_.*frac', c): return 'vegetation'
    if c.startswith('cl_'): return 'climate'
    if re.search(r'soil|sand|clay|poros|peat', c): return 'soil'
    return 'other'
G = {c: grp(c) for c in cands}
pd.Series(G, name='group').rename_axis('covariate').to_csv(OUT / 'multivariate_covariate_groups.csv')
GROUPS = ['obs_supply', 'operator_skill', 'signal_noise', 'terrain_roughness', 'vegetation', 'climate', 'soil', 'baseline_skill', 'other']

blk = (np.floor(d.com_lon / 2).astype(int).astype(str) + '_' + np.floor(d.com_lat / 2).astype(int).astype(str))
blk_id = pd.factorize(blk)[0]
print('2x2 deg blocks in mask:', blk_id.max() + 1)

def preprocess(sub):
    """Unsupervised column screening (outcome not used). Returns kept cols + log."""
    log = []
    X = sub[cands]
    keep = []
    for c in cands:
        x = X[c]
        if x.isna().mean() > MAXMISS: log.append((c, 'missing>30%')); continue
        v = x.dropna()
        if v.std() < 1e-10 or v.nunique() < 2 or v.round(10).value_counts(normalize=True).iloc[0] > 0.95:
            log.append((c, 'near-constant')); continue
        keep.append(c)
    # |rho|>0.9 pruning (Spearman, pairwise); keep earlier in priority order: baseline, then lower missing
    order = sorted(keep, key=lambda c: (c not in base_cols, X[c].isna().mean()))
    R = X[order].corr(method='spearman').abs()
    kept = []
    for c in order:
        hit = [k for k in kept if R.loc[c, k] > RHO]
        if hit: log.append((c, f'|rho|>0.9 with {hit[0]} (rho={R.loc[c, hit[0]]:.2f})'))
        else: kept.append(c)
    return kept, log

def mk_gb(seed):
    return GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.04, subsample=0.7,
                                     max_features=0.33, min_samples_leaf=5, random_state=seed)

def fit_predict(model, Xtr, ytr, gtr, Xte, seed):
    """Fit on train (impute within fold) and return predictor closure."""
    imp = SimpleImputer(strategy='median').fit(Xtr)
    A = imp.transform(Xtr)
    if model == 'enet':
        sc = StandardScaler().fit(A); A = sc.transform(A)
        inner = list(GroupKFold(4, shuffle=True, random_state=seed).split(A, ytr, gtr))
        m = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], alphas=40, cv=inner, max_iter=5000, random_state=seed).fit(A, ytr)
        return lambda X: m.predict(sc.transform(imp.transform(X)))
    m = mk_gb(seed).fit(A, ytr)
    return lambda X: m.predict(imp.transform(X))

def cv_repeat(X, y, g, r, cols, model_list, do_imp, subsets):
    """One CV repeat; returns dict of OOF R2 + importance."""
    seed = 1000 + r
    rng = np.random.RandomState(seed)
    gkf = GroupKFold(NFOLD, shuffle=True, random_state=seed)
    folds = list(gkf.split(X, y, g))
    colidx = {c: i for i, c in enumerate(cols)}
    sst = ((y - y.mean()) ** 2).sum()
    res = {}
    def r2(p): return 1 - ((y - p) ** 2).sum() / sst
    yn = rng.permutation(y)
    oof = {k: np.zeros(len(y)) for k in
           [f'{m}' for m in model_list] + [f'null_{m}' for m in model_list] + ['latlon']}
    imp_rec = {m: [] for m in model_list}
    for tr, te in folds:
        for m in model_list:
            f = fit_predict(m, X[tr], y[tr], g[tr], X[te], seed)
            oof[m][te] = p0 = f(X[te])
            fn = fit_predict(m, X[tr], yn[tr], g[tr], X[te], seed)
            oof[f'null_{m}'][te] = fn(X[te])
            if do_imp:
                base = ((y[te] - p0) ** 2).sum()
                row = {}
                units = [[c] for c in cols] + [[c for c in cols if G[c] == gg] for gg in GROUPS]
                names = cols + ['GROUP:' + gg for gg in GROUPS]
                for nm, u in zip(names, units):
                    if not u: continue
                    ix = [colidx[c] for c in u]; acc = 0
                    for _ in range(3):
                        Xp = X[te].copy(); pi = rng.permutation(len(te))
                        Xp[:, ix] = X[te][pi][:, ix]
                        acc += ((y[te] - f(Xp)) ** 2).sum()
                    row[nm] = (acc / 3 - base)
                imp_rec[m].append(row)
        ll = fit_predict('gb', LL[tr], y[tr], g[tr], LL[te], seed)
        oof['latlon'][te] = ll(LL[te])
    for k, v in oof.items(): res['r2_' + k] = r2(v)
    res['oof_gb'] = oof['gb']
    if do_imp:
        for m in model_list:
            df = pd.DataFrame(imp_rec[m]).sum() / sst   # R2 loss, summed over folds
            res['imp_' + m] = df
    if subsets is not None:  # group ablation (GB and enet)
        ab = {}
        for m in model_list:
            for nm, cs in subsets.items():
                if not cs: continue
                ix = [colidx[c] for c in cs]; p = np.zeros(len(y))
                for tr, te in folds:
                    f = fit_predict(m, X[tr][:, ix], y[tr], g[tr], X[te][:, ix], seed); p[te] = f(X[te][:, ix])
                ab[(m, nm)] = r2(p)
        res['ablate'] = ab
    return res

LL = None
summary, group_imp_rows, cov_imp_rows, ablate_rows, drop_log, comp_rows = [], [], [], [], [], []
pdp_store = {}
for oc in OUTCOMES:
    sub = d[d[oc].notna()].reset_index(drop=True)
    kept, log = preprocess(sub)
    for c, why in log: drop_log.append(dict(outcome=oc, covariate=c, reason=why))
    y = sub[oc].values.astype(float); X = sub[kept].values.astype(float)
    g = pd.factorize(blk[d[oc].notna().values].values)[0]
    LL = sub[['com_lon', 'com_lat']].values
    cols = kept
    ngr = {gg: sum(G[c] == gg for c in kept) for gg in GROUPS}
    print(f'== {oc}: n={len(y)} blocks={g.max()+1} p={len(kept)}', ngr, flush=True)
    subsets = {}
    for gg in GROUPS:
        subsets['only_' + gg] = [c for c in kept if G[c] == gg]
        subsets['drop_' + gg] = [c for c in kept if G[c] != gg]
    subsets['SN'] = [c for c in kept if G[c] == 'signal_noise']
    subsets['SN+terrain'] = [c for c in kept if G[c] in ('signal_noise', 'terrain_roughness')]
    subsets['SN+baseline'] = [c for c in kept if G[c] in ('signal_noise', 'baseline_skill')]
    subsets['SN+baseline+terrain'] = [c for c in kept if G[c] in ('signal_noise', 'baseline_skill', 'terrain_roughness')]
    subsets['all'] = list(kept)
    subsets = {k: v for k, v in subsets.items() if v}
    models = ['enet', 'gb']
    out = Parallel(n_jobs=NREP, backend='loky')(delayed(cv_repeat)(X, y, g, r, cols, models, r < NIMP,
                                                  subsets if r < NIMP else None) for r in range(NREP))
    row = dict(outcome=oc, n=len(y), n_blocks=int(g.max() + 1), n_cov=len(kept), y_sd=y.std())
    for k in ['enet', 'gb', 'null_enet', 'null_gb', 'latlon']:
        v = np.array([o['r2_' + k] for o in out]); row[f'r2_{k}_mean'] = v.mean(); row[f'r2_{k}_sd'] = v.std(ddof=1)
    nul = np.array([o['r2_null_gb'] for o in out])
    row['gb_minus_latlon_mean'] = np.mean([o['r2_gb'] - o['r2_latlon'] for o in out])
    row['gb_minus_latlon_sd'] = np.std([o['r2_gb'] - o['r2_latlon'] for o in out], ddof=1)
    summary.append(row)
    for m in models:
        imps = [o['imp_' + m] for o in out if 'imp_' + m in o]
        I = pd.concat(imps, axis=1)
        for nm in I.index:
            rec = dict(outcome=oc, model=m, name=nm, imp_mean=I.loc[nm].mean(), imp_sd=I.loc[nm].std(ddof=1))
            (group_imp_rows if nm.startswith('GROUP:') else cov_imp_rows).append(
                {**rec, 'group': nm[6:] if nm.startswith('GROUP:') else G[nm]})
    abl = pd.DataFrame([{**dict(zip(['model', 'subset'], k)), 'r2': v, 'rep': i}
                        for i, o in enumerate(out) if 'ablate' in o for k, v in o['ablate'].items()])
    for (m, s), q in abl.groupby(['model', 'subset']):
        ablate_rows.append(dict(outcome=oc, model=m, subset=s, n_cov=len(subsets[s]), r2_mean=q.r2.mean(), r2_sd=q.r2.std(ddof=1)))
    for m in models:
        P = abl[abl.model == m].pivot(index='rep', columns='subset', values='r2')
        for a, b in [('SN', 'SN+terrain'), ('SN+baseline', 'SN+baseline+terrain')]:
            if a in P and b in P:
                dd = P[b] - P[a]
                comp_rows.append(dict(outcome=oc, model=m, base=a, plus_terrain=b, delta_r2_mean=dd.mean(),
                                      delta_r2_sd=dd.std(ddof=1), r2_base=P[a].mean(), r2_with_terrain=P[b].mean()))
    # PDP data: full-data GB (median-imputed), top-6 covariates by GB perm importance
    ci = pd.DataFrame([r for r in cov_imp_rows if r['outcome'] == oc and r['model'] == 'gb']).sort_values('imp_mean', ascending=False)
    top = ci.name.head(6).tolist()
    imp = SimpleImputer(strategy='median').fit(X); A = imp.transform(X)
    ms = [mk_gb(s).fit(A, y) for s in range(5)]
    pd_data = {}
    for c in top:
        j = cols.index(c); gridv = np.quantile(A[:, j], np.linspace(0.05, 0.95, 20)); pv = []
        for v in gridv:
            B = A.copy(); B[:, j] = v; pv.append(np.mean([m_.predict(B).mean() for m_ in ms]))
        pd_data[c] = (gridv, np.array(pv), A[:, j])
    pdp_store[oc] = pd_data

S = pd.DataFrame(summary); S.to_csv(OUT / 'multivariate_cv_r2_summary.csv', index=False)
pd.DataFrame(group_imp_rows).to_csv(OUT / 'multivariate_group_perm_importance.csv', index=False)
pd.DataFrame(cov_imp_rows).sort_values(['outcome', 'model', 'imp_mean'], ascending=[True, True, False]).to_csv(OUT / 'multivariate_covariate_perm_importance.csv', index=False)
pd.DataFrame(ablate_rows).to_csv(OUT / 'multivariate_group_ablation_r2.csv', index=False)
pd.DataFrame(comp_rows).to_csv(OUT / 'multivariate_terrain_vs_signalnoise.csv', index=False)
pd.DataFrame(drop_log).to_csv(OUT / 'multivariate_preprocessing_dropped.csv', index=False)

# ---- figures
fig, ax = plt.subplots(figsize=(12, 5)); w = 0.16; x = np.arange(len(S))
for i, (k, lab) in enumerate([('enet', 'elastic-net'), ('gb', 'grad. boosting'), ('latlon', 'lat/lon-only GB'),
                              ('null_gb', 'null GB (shuffled y)'), ('null_enet', 'null enet')]):
    ax.bar(x + (i - 2) * w, S[f'r2_{k}_mean'], w, yerr=S[f'r2_{k}_sd'], label=lab, capsize=2)
ax.axhline(0, color='k', lw=.8); ax.set_xticks(x); ax.set_xticklabels(S.outcome, rotation=40, ha='right')
ax.set_ylabel('spatial-block CV R2 (mean +/- sd over 10 repeats)'); ax.legend(fontsize=8); ax.set_ylim(-0.6, 0.5)
ax.set_title('EXPLORATORY: block-CV R2 vs null and lat/lon baseline'); fig.tight_layout(); fig.savefig(OUT / 'multivariate_cv_r2.png', dpi=130); plt.close()

for nm, df, fn in [('GB group permutation importance (R2 loss, held-out blocks)', pd.DataFrame(group_imp_rows), 'multivariate_group_perm_importance_gb.png')]:
    q = df[df.model == 'gb'].pivot(index='outcome', columns='group', values='imp_mean').reindex(columns=GROUPS).reindex(OUTCOMES)
    fig, ax = plt.subplots(figsize=(10, 5)); im = ax.imshow(q.values, cmap='viridis', aspect='auto')
    ax.set_xticks(range(len(q.columns))); ax.set_xticklabels(q.columns, rotation=40, ha='right')
    ax.set_yticks(range(len(q.index))); ax.set_yticklabels(q.index)
    plt.colorbar(im, label='R2 loss'); ax.set_title('EXPLORATORY: ' + nm); fig.tight_layout(); fig.savefig(OUT / fn, dpi=130); plt.close()
A_ = pd.DataFrame(ablate_rows)
q = A_[(A_.model == 'gb') & A_.subset.str.startswith('only_')].pivot(index='outcome', columns='subset', values='r2_mean').reindex(OUTCOMES)
q2 = A_[(A_.model == 'gb') & (A_.subset == 'all')].set_index('outcome').r2_mean
dq = pd.DataFrame({g_: q2 - A_[(A_.model == 'gb') & (A_.subset == 'drop_' + g_)].set_index('outcome').r2_mean for g_ in GROUPS}).reindex(OUTCOMES)
fig, ax = plt.subplots(figsize=(10, 5)); v = np.nanmax(np.abs(dq.values)); im = ax.imshow(dq.values, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto')
ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS, rotation=40, ha='right'); ax.set_yticks(range(len(OUTCOMES))); ax.set_yticklabels(OUTCOMES)
plt.colorbar(im, label='CV R2(all) - CV R2(drop group)'); ax.set_title('EXPLORATORY: drop-one-group, GB (positive = group helps)'); fig.tight_layout()
fig.savefig(OUT / 'multivariate_drop_one_group_gb.png', dpi=130); plt.close()
dq.to_csv(OUT / 'multivariate_drop_one_group_gb_delta.csv')

for oc, pdd in pdp_store.items():
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))
    for ax, (c, (gv, pv, xv)) in zip(axs.ravel(), pdd.items()):
        ax.plot(gv, pv, 'o-'); ax.set_title(f'{c}\n[{G[c]}]', fontsize=8); ax2 = ax.twinx(); ax2.hist(xv, bins=20, alpha=.15, color='gray'); ax2.set_yticks([])
    fig.suptitle(f'EXPLORATORY partial dependence (GB, full-data fit, 5 seeds) - {oc}', fontsize=10); fig.tight_layout()
    fig.savefig(OUT / f'multivariate_pdp_{oc}.png', dpi=110); plt.close()
print(S.round(3).to_string())
