"""Task 1: outcome structure (EXPLORATORY). Outputs to output/benefit_drivers/analysis_outcomes/."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from outcomes_common import *
from scipy import stats
a = load(); m = m_only(a)
cols = INDEP + EXTRA
# --- correlation matrix + block bootstrap CI
R = spearman_pairs(m, cols); B = boot_corr(m, cols, 1000)
lo, hi = np.nanpercentile(B, 2.5, axis=0), np.nanpercentile(B, 97.5, axis=0)
n = m[cols].notna().astype(int); Npair = n.T.dot(n)
rows = [(c1, c2, R.loc[c1, c2], lo[i, j], hi[i, j], Npair.loc[c1, c2]) for i, c1 in enumerate(cols) for j, c2 in enumerate(cols) if i < j]
pd.DataFrame(rows, columns=['a', 'b', 'spearman', 'ci_lo', 'ci_hi', 'n_pair']).round(3).to_csv(OUT / 'corr_outcomes_blockboot.csv', index=False)
fig, ax = plt.subplots(figsize=(7.5, 6.5))
im = ax.imshow(R.values, cmap='RdBu_r', vmin=-1, vmax=1)
ax.set_xticks(range(len(cols))); ax.set_yticks(range(len(cols)))
ax.set_xticklabels([LABEL[c] for c in cols], rotation=60, ha='right', fontsize=7); ax.set_yticklabels([LABEL[c] for c in cols], fontsize=7)
for i in range(len(cols)):
    for j in range(len(cols)):
        ax.text(j, i, f'{R.values[i,j]:.2f}' + ('' if i == j else f'\n[{lo[i,j]:.2f},{hi[i,j]:.2f}]'), ha='center', va='center', fontsize=5)
ax.axhline(2.5, c='k', lw=.8); ax.axvline(2.5, c='k', lw=.8)
plt.colorbar(im, shrink=.8); ax.set_title('Spearman, 552-tile mask, pairwise-complete; [95% 2x2-deg block-bootstrap CI]\nSMAP family (first 3) | ASCAT family (next 3) | extras', fontsize=8)
plt.tight_layout(); plt.savefig(OUT / 'fig01_corr_matrix.png', dpi=150); plt.close()
# family-level block means: within vs between
def fam_stat(df):
    c = df[INDEP].corr(method='spearman').values
    w_s = np.mean([c[0, 1], c[0, 2], c[1, 2]]); w_a = np.mean([c[3, 4], c[3, 5], c[4, 5]])
    bt = np.mean(c[:3, 3:]); return np.array([w_s, w_a, bt])
rng = np.random.default_rng(5); ub = m.block.unique(); g = {b: np.where(m.block.values == b)[0] for b in ub}
bs = np.array([fam_stat(m.iloc[np.concatenate([g[b] for b in rng.choice(ub, len(ub))])]) for _ in range(1000)])
est = fam_stat(m)
fam = pd.DataFrame({'quantity': ['mean within-SMAP-family r', 'mean within-ASCAT-family r', 'mean between-family r'], 'est': est,
                    'lo': np.percentile(bs, 2.5, 0), 'hi': np.percentile(bs, 97.5, 0)}).round(3)
fam.to_csv(OUT / 'family_corr_summary.csv', index=False); print(fam)
# --- PCA / factor on complete cases (ranks -> normal scores)
cc = m[INDEP].dropna(); print('PCA complete cases', len(cc))
Z = cc.rank().apply(lambda c: stats.norm.ppf((c - .5) / len(c)))
Zs = (Z - Z.mean()) / Z.std(); C = np.corrcoef(Zs.values.T)
ev, evec = np.linalg.eigh(C); o = np.argsort(ev)[::-1]; ev, evec = ev[o], evec[:, o]
def varimax(L, it=100):
    p, k = L.shape; Rm = np.eye(k); d = 0
    for _ in range(it):
        Lr = L @ Rm; u, s, vt = np.linalg.svd(L.T @ (Lr ** 3 - Lr @ np.diag((Lr ** 2).sum(0)) / p)); Rm = u @ vt
        if s.sum() < d * (1 + 1e-6): break
        d = s.sum()
    return L @ Rm
Lr = varimax(evec[:, :2] * np.sqrt(ev[:2]))
for k in range(2):  # sign: majority positive
    if Lr[:, k].sum() < 0: Lr[:, k] *= -1
pca = pd.DataFrame({'measure': INDEP, 'PC1': evec[:, 0] * np.sqrt(ev[0]), 'PC2': evec[:, 1] * np.sqrt(ev[1]),
                    'varimax_F1': Lr[:, 0], 'varimax_F2': Lr[:, 1]}).round(3)
pca.to_csv(OUT / 'pca_loadings.csv', index=False)
# parallel-analysis-ish: bootstrap eigenvalue CIs by block
def eigs(df):
    d = df[INDEP].dropna(); return np.sort(np.linalg.eigvalsh(d.corr(method='spearman').values))[::-1]
eb = np.array([eigs(m.iloc[np.concatenate([g[b] for b in rng.choice(ub, len(ub))])]) for _ in range(500)])
ev_tab = pd.DataFrame({'comp': range(1, 7), 'eigenvalue': ev, 'pct_var': ev / ev.sum() * 100, 'boot_lo': np.percentile(eb, 2.5, 0), 'boot_hi': np.percentile(eb, 97.5, 0)}).round(3)
ev_tab.to_csv(OUT / 'pca_eigenvalues.csv', index=False); print(ev_tab); print(pca)
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(range(1, 7), ev, 'o-'); ax[0].fill_between(range(1, 7), ev_tab.boot_lo, ev_tab.boot_hi, alpha=.3); ax[0].axhline(1, c='k', ls=':'); ax[0].set_xlabel('component'); ax[0].set_ylabel('eigenvalue'); ax[0].set_title(f'scree (n={len(cc)} complete tiles, 95% block-boot band)', fontsize=8)
for i, c in enumerate(INDEP):
    ax[1].scatter(Lr[i, 0], Lr[i, 1], c='C0' if c in SMAP_FAM else 'C3'); ax[1].annotate(LABEL[c], (Lr[i, 0], Lr[i, 1]), fontsize=7)
ax[1].set_xlabel('varimax F1'); ax[1].set_ylabel('varimax F2'); ax[1].axhline(0, c='.7'); ax[1].axvline(0, c='.7'); ax[1].set_title('rotated loadings (blue = SMAP family, red = ASCAT family)', fontsize=8)
plt.tight_layout(); plt.savefig(OUT / 'fig02_pca.png', dpi=150); plt.close()
# --- TC decomposition. R2s*R2a is model-independent -> exact: dlogR2m = .5 dlog r2(m,SMAP) + .5 dlog r2(m,ASCAT), r2(m,x)=R2m*R2x
t = m.copy()
for arm in ['OL', 'L1', 'L3']:
    t[f'r2ms_{arm}'] = t[f'tc_R2mod_{arm}'] * t[f'tc_R2_smapL3_{arm}']
    t[f'r2ma_{arm}'] = t[f'tc_R2mod_{arm}'] * t[f'tc_R2_ascatH121_{arm}']
    t[f'sa_{arm}'] = t[f'tc_R2_smapL3_{arm}'] * t[f'tc_R2_ascatH121_{arm}']
chk = (np.log(t.sa_L1) - np.log(t.sa_OL)).dropna()
print('check |dlog(R2s*R2a)| L1-OL  median/p95:', chk.abs().median(), chk.abs().quantile(.95), len(chk))
t['dlogR2'] = np.log(t.tc_R2mod_L1 / t.tc_R2mod_OL)
t['leg_smap'] = .5 * np.log(t.r2ms_L1 / t.r2ms_OL); t['leg_ascat'] = .5 * np.log(t.r2ma_L1 / t.r2ma_OL)
t['dr2ms'] = t.r2ms_L1 - t.r2ms_OL; t['dr2ma'] = t.r2ma_L1 - t.r2ma_OL
t['d_R2smapleg'] = t.tc_R2_smapL3_L1 - t.tc_R2_smapL3_OL; t['d_R2ascatleg'] = t.tc_R2_ascatH121_L1 - t.tc_R2_ascatH121_OL
t['d_sig2_smap'] = t.tc_sigma2_smapL3_L1 / t.tc_sigma2_smapL3_OL - 1; t['d_sig2_ascat'] = t.tc_sigma2_ascatH121_L1 / t.tc_sigma2_ascatH121_OL - 1
tt = t.dropna(subset=['dlogR2', 'leg_smap', 'leg_ascat', 'out_tc_dR2'])
print('n TC decomposable', len(tt), 'resid of identity', (tt.dlogR2 - tt.leg_smap - tt.leg_ascat).abs().median())
def shares(df):
    v = df.dlogR2.var(); return np.array([np.cov(df.dlogR2, df.leg_smap)[0, 1] / v, np.cov(df.dlogR2, df.leg_ascat)[0, 1] / v])
e = boot_stat(tt, lambda d: shares(d)[0]), boot_stat(tt, lambda d: shares(d)[1])
tt = tt.assign(resid=tt.dlogR2 - tt.leg_smap - tt.leg_ascat)
e3 = boot_stat(tt, lambda d: np.cov(d.dlogR2, d.resid)[0, 1] / d.dlogR2.var())
clean = tt[np.abs(np.log(tt.sa_L1 / tt.sa_OL)) < 0.01]
ec = [boot_stat(clean, lambda d, k=k: shares(d)[k]) for k in (0, 1)]
print('clean-identity tiles', len(clean), 'of', len(tt))
dec = [('variance share of dlogR2_mod from SMAP leg', *e[0]), ('variance share from ASCAT leg', *e[1]), ('variance share from identity residual (clamping/NaN-inconsistent tiles)', *e3),
       (f'clean-identity subset n={len(clean)}: SMAP leg share', *ec[0]), ('clean-identity subset: ASCAT leg share', *ec[1])]
rows = []
for x in ['dlogR2', 'leg_smap', 'leg_ascat', 'dr2ms', 'dr2ma', 'd_R2smapleg', 'd_R2ascatleg', 'd_sig2_smap', 'd_sig2_ascat']:
    for y in INDEP + ['out_tc_dSig2pct']:
        d = tt[[x, y]].dropna(); rows.append((x, y, stats.spearmanr(d[x], d[y])[0], len(d)))
pd.DataFrame(rows, columns=['tc_quantity', 'outcome', 'spearman', 'n']).round(3).to_csv(OUT / 'tc_decomposition_corr.csv', index=False)
pd.DataFrame(dec, columns=['quantity', 'est', 'lo', 'hi']).round(3).to_csv(OUT / 'tc_decomposition_shares.csv', index=False); print(dec)
summ = tt[['dlogR2', 'leg_smap', 'leg_ascat', 'out_tc_dR2']].describe().round(4); print(summ)
summ.to_csv(OUT / 'tc_decomposition_summary.csv')
print(tt[['out_tc_dR2', 'leg_smap', 'leg_ascat', 'out_iv_ascatH121', 'out_iv_smapL3', 'out_omf_indep', 'out_omf_ascat']].corr(method='spearman').round(2).to_string())
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for k, (c, ttl) in enumerate([('leg_smap', 'SMAP leg: 0.5*dlog r2(model,SMAP)'), ('leg_ascat', 'ASCAT leg: 0.5*dlog r2(model,ASCAT)')]):
    ax[k].scatter(tt[c], tt.dlogR2, s=8); ax[k].set_xlabel(ttl, fontsize=8); ax[k].set_ylabel('dlog R2_mod (L1-OL)'); r = stats.spearmanr(tt[c], tt.dlogR2)[0]; ax[k].set_title(f'rho={r:.2f}', fontsize=8)
    ax[k].axhline(0, c='.7'); ax[k].axvline(0, c='.7')
ax[2].scatter(tt.leg_smap, tt.leg_ascat, s=8, c=tt.dlogR2, cmap='RdBu', vmin=-.5, vmax=.5); ax[2].set_xlabel('SMAP leg'); ax[2].set_ylabel('ASCAT leg'); ax[2].set_title('legs coloured by dlogR2', fontsize=8)
plt.tight_layout(); plt.savefig(OUT / 'fig03_tc_legs.png', dpi=150); plt.close()
# --- regional: Four Corners box and hotspot
box = (a.com_lon.between(-111.5, -107)) & (a.com_lat.between(34.5, 37.5))
a['fc_box'] = box
# hotspot: 3x3 smoothing of omf_ascat penalty on grid; top quartile of L1-area, component around seed
g_ = a.set_index(['i_indg', 'j_indg'])
v = a.out_omf_ascat.where(a.in_l1_mask.astype(bool))
pos = {(i, j): k for k, (i, j) in enumerate(zip(a.i_indg, a.j_indg))}
sm = np.full(len(a), np.nan)
for k, (i, j) in enumerate(zip(a.i_indg, a.j_indg)):
    nb = [v.iloc[pos[(i + di, j + dj)]] for di in (-1, 0, 1) for dj in (-1, 0, 1) if (i + di, j + dj) in pos]
    sm[k] = np.nanmean(nb) if np.isfinite(nb).any() else np.nan
a['ascat_smooth'] = sm
thr = np.nanpercentile(sm[a.in_l1_mask.astype(bool)], 25)   # penalty = more negative out_omf_ascat (positive = better)
cand = set(np.where((sm <= thr) & a.in_l1_mask.astype(bool))[0])
cl = sorted(cand); seed = cl[int(np.argmin(((a.com_lon.values[cl] + 109.05) ** 2 + (a.com_lat.values[cl] - 37.0) ** 2)))]; comp = {seed}; fr = list(comp)
while fr:
    k = fr.pop()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            q = pos.get((a.i_indg[k] + di, a.j_indg[k] + dj))
            if q is not None and q in cand and q not in comp: comp.add(q); fr.append(q)
a['fc_hot'] = False; a.loc[list(comp), 'fc_hot'] = True
print('hotspot tiles', a.fc_hot.sum(), 'box tiles in mask', (a.fc_box & a.in_l1_mask.astype(bool)).sum())
a[['tile_axis_index', 'fc_box', 'fc_hot', 'ascat_smooth']].to_csv(OUT / 'four_corners_flags.csv', index=False)
m = m_only(a); rows = []
for rname in ['fc_box', 'fc_hot']:
    for nm, sub in [('inside', m[m[rname]]), ('outside', m[~m[rname]])]:
        for x, y in [('out_tc_dR2', 'out_iv_ascatH121'), ('out_tc_dR2', 'out_omf_ascat'), ('out_tc_dR2', 'out_omf_indep'), ('out_iv_ascatH121', 'out_omf_ascat'), ('out_iv_smapL3', 'out_omf_indep'), ('out_iv_ascatH121', 'out_iv_smapL3')]:
            d = sub[[x, y]].dropna(); rows.append((rname, nm, x, y, stats.spearmanr(d[x], d[y])[0] if len(d) > 5 else np.nan, len(d)))
rg = pd.DataFrame(rows, columns=['region', 'where', 'x', 'y', 'spearman', 'n']).round(3); rg.to_csv(OUT / 'family_corr_inside_outside_region.csv', index=False); print(rg.to_string())
# medians of outcomes inside/outside + share of summed negative contributions
rows = []
for rname in ['fc_box', 'fc_hot']:
    for c in INDEP:
        for nm, sub in [('inside', m[m[rname]]), ('outside', m[~m[rname]])]:
            x = sub[c].dropna(); rows.append((rname, c, nm, len(x), x.median(), x.mean(), (x < 0).mean()))
        x = m[c].dropna(); neg = x.clip(upper=0).sum(); ins = m.loc[x.index][rname]
        rows.append((rname, c, 'share_of_summed_penalty_inside', len(x), np.nan, x[ins].clip(upper=0).sum() / neg if neg < 0 else np.nan, np.nan))
rr = pd.DataFrame(rows, columns=['region', 'outcome', 'where', 'n', 'median', 'mean_or_share', 'frac_negative']).round(3); rr.to_csv(OUT / 'outcomes_by_region.csv', index=False); print(rr.to_string())
# leave-region-out: ascat-family correlations and PCA-ish
rows = []
for rname in ['fc_box', 'fc_hot']:
    sub = m[~m[rname]]
    for x, y in [('out_tc_dR2', 'out_iv_ascatH121'), ('out_tc_dR2', 'out_omf_ascat'), ('out_tc_dR2', 'out_omf_indep'), ('out_iv_smapL3', 'out_omf_indep')]:
        d = sub[[x, y]].dropna(); rows.append((rname, x, y, stats.spearmanr(d[x], d[y])[0], len(d), 'full', stats.spearmanr(*[m[[x, y]].dropna()[k] for k in (x, y)])[0]))
print(pd.DataFrame(rows, columns=['excl', 'x', 'y', 'rho_excl', 'n', '_', 'rho_full']).round(3).to_string())
# --- maps
fig, axs = mapfig(9)
for ax, c in zip(axs, INDEP + ['out_tc_dSig2pct', 'out_omf_l3obs', 'out_omf_l1obs']):
    tile_map(ax, a, c, LABEL[c] + ' (blue = L1 better)')
for ax in axs[9:]: ax.axis('off')
plt.tight_layout(); plt.savefig(OUT / 'fig04_maps_outcomes.png', dpi=130); plt.close()
fig, axs = mapfig(6)
tmp = a.copy()
tmp['leg_smap'] = np.nan; tmp['leg_ascat'] = np.nan; tmp['dlogR2'] = np.nan
tmp.loc[tt.index.map(lambda i: i) if False else a.index[a.tile_axis_index.isin(tt.tile_axis_index)], 'leg_smap'] = tt.set_index('tile_axis_index').leg_smap.reindex(a.tile_axis_index[a.tile_axis_index.isin(tt.tile_axis_index)]).values
tmp.loc[a.index[a.tile_axis_index.isin(tt.tile_axis_index)], 'leg_ascat'] = tt.set_index('tile_axis_index').leg_ascat.reindex(a.tile_axis_index[a.tile_axis_index.isin(tt.tile_axis_index)]).values
tmp.loc[a.index[a.tile_axis_index.isin(tt.tile_axis_index)], 'dlogR2'] = tt.set_index('tile_axis_index').dlogR2.reindex(a.tile_axis_index[a.tile_axis_index.isin(tt.tile_axis_index)]).values
tmp['hot'] = np.where(a.fc_hot, 1., np.where(a.in_l1_mask, 0., np.nan))
tmp['avail_ascat'] = a[ASCAT_FAM].notna().sum(axis=1).astype(float)
for ax, c, tl in zip(axs, ['dlogR2', 'leg_smap', 'leg_ascat', 'ascat_smooth', 'hot', 'avail_ascat'],
                     ['TC dlogR2_mod', 'SMAP-leg part', 'ASCAT-leg part', 'O-F ASCAT, 3x3 smoothed', 'Four Corners hotspot (1)', 'n ASCAT-family measures available']):
    tile_map(ax, tmp, c, tl, sym=(c not in ('hot', 'avail_ascat')), cmap='RdBu' if c not in ('hot', 'avail_ascat') else 'viridis', vmin=0 if c in ('hot', 'avail_ascat') else None, vmax=1 if c == 'hot' else (3 if c == 'avail_ascat' else None))
plt.tight_layout(); plt.savefig(OUT / 'fig05_maps_tc_legs_region.png', dpi=130); plt.close()
# ASCAT availability inside mask
print('ASCAT-fam availability in mask:', m[ASCAT_FAM].notna().sum().to_dict(), 'TC valid', m.out_tc_dR2.notna().sum())
print('ASCAT iv available fraction by hot', m.groupby('fc_hot').out_iv_ascatH121.apply(lambda s: s.notna().mean()).to_dict())
