"""EXPLORATORY: (A) time resolution of L1 benefit drivers, (B) roughness proxies vs benefit.

Association only. Block bootstrap = 2x2 deg blocks (resample blocks with replacement).
IGOT source is private and not used; the operator is described from repo docs only.
Outputs -> output/benefit_drivers/analysis_seasonal_roughness/
"""
import numpy as np, pandas as pd, warnings
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata
from scipy.spatial import cKDTree
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
BASE = "/Users/amfox/Desktop/geosldas-analysis/projects/CYGNSS_L1_AZ/output/benefit_drivers/"
OUT = BASE + "analysis_seasonal_roughness/"
B = 1000
SEAS = ["DJF", "MAM", "JJA", "SON"]

d = pd.read_csv(BASE + "analysis_table.csv")
d = d[d.in_l1_mask.astype(bool)].reset_index(drop=True)
d["block"] = (np.floor((d.com_lon + 180) / 2).astype(int) * 1000 + np.floor((d.com_lat + 90) / 2).astype(int))
blocks = d.block.values
ub = np.unique(blocks)
idx_by_block = {b: np.where(blocks == b)[0] for b in ub}
print("mask tiles", len(d), "blocks", len(ub))


def boot_idx():
    pick = rng.choice(ub, len(ub), replace=True)
    return np.concatenate([idx_by_block[b] for b in pick])


BOOT = [boot_idx() for _ in range(B)]  # shared resamples


def boot_stat(fn, *arrs):
    """fn(*arrays)->float; returns point, lo, hi (block bootstrap percentile 95%)."""
    arrs = [np.asarray(a, float) for a in arrs]
    pt = fn(*arrs)
    vals = []
    for bi in BOOT:
        try:
            v = fn(*[a[bi] for a in arrs])
        except Exception:
            v = np.nan
        vals.append(v)
    vals = np.array(vals)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 50:
        return pt, np.nan, np.nan
    return pt, *np.percentile(vals, [2.5, 97.5])


def spear(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 20:
        return np.nan
    return spearmanr(x[m], y[m])[0]


def resid_rank(y, Z):
    ry = rankdata(y)
    Zr = np.column_stack([np.ones(len(y))] + [rankdata(z) for z in Z.T])
    beta = np.linalg.lstsq(Zr, ry, rcond=None)[0]
    return ry - Zr @ beta


def pspear(x, y, Z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(1)
    if m.sum() < 30:
        return np.nan
    x, y, Z = x[m], y[m], Z[m]
    return np.corrcoef(resid_rank(x, Z), resid_rank(y, Z))[0, 1]


def ci(t):
    return f"{t[0]:+.2f} [{t[1]:+.2f},{t[2]:+.2f}]"


# ======================================================================
# PART A: TIME
# ======================================================================
A = []  # report lines
dom = pd.read_csv(BASE + "benefit_table_seasonal_domain.csv")
gain = {}
for g in ["SMAP", "SMOS", "ASCAT", "L3obs", "L1obs"]:
    l1 = dom[(dom.arm == "L1") & (dom.group == g)].set_index("season").omf_pct_better
    l3 = dom[(dom.arm == "L3") & (dom.group == g)].set_index("season").omf_pct_better
    gain[g] = l1.reindex(SEAS)
    gain[g + "_L1minusL3"] = (l1 - l3).reindex(SEAS)
gain = pd.DataFrame(gain)
gain.to_csv(OUT + "seasonal_domain_gain_used.csv")

# seasonal covariate panel (long, tile x season)
cov_def = {  # name: column template with {S} (as given) or {s} lower
    "rs_smterm": "op_rs_smterm@{S}", "rs_H": "op_rs_H@{S}", "rs_geom": "op_rs_geom@{S}",
    "sd_smterm": "op_sd_smterm@{S}", "sd_geom": "op_sd_geom@{S}", "R_impl": "op_R_impl@{S}",
    "snr": "op_snr@{S}", "dyn_ratio": "op_dyn_ratio@{S}", "frac_out15": "op_frac_out15@{S}",
    "op_sfmc_mean": "op_sfmc_mean@{S}", "N_obs": "op_N@{S}",
    "H_std_db": "cv_sn_H_std_db_{s}", "H_spread_db": "cv_sn_spread_H_db_{s}", "K_mean": "cv_sn_K_mean_{s}",
    "sfmc_std": "sm_sfmc_std_{S}", "sfmc_mean": "sm_sfmc_mean_{S}",
    "LAI": "vg_lai_{S}_mean", "NDVI": "vg_ndvi_{S}_mean", "veg_opacity": "vg_vegopacity_{S}_mean",
    "precip_mm": "cl_precip_{S}_mean_mm",
}
rows = []
for S in SEAS:
    t = pd.DataFrame({"tile": d.tile_axis_index, "block": d.block, "season": S})
    for k, tpl in cov_def.items():
        c = tpl.format(S=S, s=S.lower())
        t[k] = d[c].values if c in d else np.nan
    t["geom_var_frac"] = t.sd_geom**2 / (t.sd_geom**2 + t.sd_smterm**2)
    rows.append(t)
panel = pd.concat(rows, ignore_index=True)

# add roughness_diagnostics tile-season (tilenum -> tile_axis_index via nearest centre; offset verified)
rd = pd.read_csv(BASE + "roughness_diagnostics_by_tile_season.csv")
full = pd.read_csv(BASE + "analysis_table.csv")
kd = cKDTree(full[["com_lon", "com_lat"]].values)
cen = rd.groupby("tilenum")[["lon", "lat"]].median()
dist, ii = kd.query(cen.values)
mp = pd.Series(full.tile_axis_index.values[ii], index=cen.index)
mp[dist > 0.25] = np.nan
print("tilenum->tile_axis_index: fraction with tilenum-1 == matched axis index:", float(((mp.index - 1) == mp.values).mean()))
rd["tile"] = rd.tilenum.map(mp)
rd = rd.dropna(subset=["tile"]).astype({"tile": int})
panel = panel.merge(rd[["tile", "season", "n", "resid_std", "slope_inc_per10deg", "frac_diffuse_coh_lt05",
                        "resid_var_over_fcstvar"]].rename(columns={"n": "rd_n", "slope_inc_per10deg": "inc_slope_season"}),
                    on=["tile", "season"], how="left")
panel.to_csv(OUT + "seasonal_tile_season_panel.csv", index=False)

covs = [c for c in panel.columns if c not in ("tile", "block", "season", "rd_n")]
# A1: within-tile seasonal anomalies, domain mean per season, paired contrast with block CI
# anomaly = value - tile mean over the 4 seasons, restricted to tiles with all 4 seasons for that covariate
anom = {}
a1 = []
for c in covs:
    w = panel.pivot(index="tile", columns="season", values=c)[SEAS]
    ok = w.notna().all(1)
    w = w[ok]
    an = w.sub(w.mean(1), axis=0)
    anom[c] = an
    blk = d.set_index("tile_axis_index").loc[w.index, "block"].values
    # bootstrap by blocks, domain mean anomaly per season
    ubk = np.unique(blk)
    ib = {b: np.where(blk == b)[0] for b in ubk}
    anv = an.values
    bs = []
    for _ in range(500):
        pick = rng.choice(ubk, len(ubk), replace=True)
        ix = np.concatenate([ib[b] for b in pick])
        bs.append(anv[ix].mean(0))
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    for k, S in enumerate(SEAS):
        a1.append(dict(covariate=c, season=S, n_tiles=int(ok.sum()), domain_mean=w[S].mean(),
                       mean_anom=anv[:, k].mean(), anom_lo=lo[k], anom_hi=hi[k]))
a1 = pd.DataFrame(a1)
a1.to_csv(OUT + "A1_seasonal_anomalies_domain.csv", index=False)

# A1b: n=4 association between seasonal anomaly of covariate and domain gain (descriptive only)
a1b = []
for c in covs:
    sub = a1[a1.covariate == c].set_index("season").reindex(SEAS)
    for g in ["SMAP", "SMOS", "ASCAT", "SMAP_L1minusL3", "ASCAT_L1minusL3"]:
        r = spearmanr(sub.mean_anom.values, gain[g].values)[0]
        a1b.append(dict(covariate=c, gain=g, spearman_n4=r))
a1b = pd.DataFrame(a1b).pivot(index="covariate", columns="gain", values="spearman_n4")
a1b.to_csv(OUT + "A1b_domain_season_n4_rank_assoc.csv")

# A2: tile-level seasonal swing vs whole-period outcome (tile resolution; swing is spatial, not temporal)
outcomes = {"SMAP-family": ["out_omf_indep", "out_iv_smapL3"],
            "ASCAT-family": ["out_omf_ascat", "out_iv_ascatH121", "out_tc_dR2"]}
all_out = [o for v in outcomes.values() for o in v]
swing = pd.DataFrame({"tile": d.tile_axis_index})
for c in ["rs_smterm", "sd_smterm", "sd_geom", "snr", "sfmc_std", "LAI", "veg_opacity", "H_std_db", "R_impl", "frac_out15", "dyn_ratio"]:
    w = panel.pivot(index="tile", columns="season", values=c)[SEAS]
    ok = w.notna().sum(1) >= 3
    sw = (w.max(1) - w.min(1)).where(ok)
    swing[c + "_range"] = swing.tile.map(sw)
    swing[c + "_JJAminusDJF"] = swing.tile.map((w["JJA"] - w["DJF"]).where(ok))
    swing[c + "_MAMminusrest"] = swing.tile.map((w["MAM"] - w[["DJF", "JJA", "SON"]].mean(1)).where(ok))
nobs = d[["cv_sup_n_obs_djf", "cv_sup_n_obs_mam", "cv_sup_n_obs_jja", "cv_sup_n_obs_son"]].sum(1)
swing["obs_frac_MAMJJA"] = (d.cv_sup_n_obs_mam + d.cv_sup_n_obs_jja) / nobs
a2 = []
for sc in [c for c in swing.columns if c != "tile"]:
    for o in all_out:
        t = boot_stat(spear, swing[sc].values, d[o].values)
        a2.append(dict(swing=sc, outcome=o, n=int((swing[sc].notna() & d[o].notna()).sum()), rho=t[0], lo=t[1], hi=t[2]))
a2 = pd.DataFrame(a2)
a2["sig"] = (a2.lo > 0) | (a2.hi < 0)
a2.to_csv(OUT + "A2_tile_seasonal_swing_vs_outcome.csv", index=False)

# A3: per tile-year outcomes
yrs = [2020, 2021, 2022]
a3 = []
varrows = []
for fam in ["smos", "smap", "ascat", "l3obs", "l1obs"]:
    Y = d[[f"out_omf_{fam}_{y}" for y in yrs]].values
    N = d[[f"cv_sup_n_obs_{y}" for y in yrs]].values
    ok = np.isfinite(Y).all(1)
    Yo, No, bo = Y[ok], N[ok], blocks[ok]
    tm = Yo.mean(1, keepdims=True)
    an = Yo - tm
    # variance decomposition (balanced 3 years)
    vb = Yo.mean(1).var(ddof=1)
    vw = an.var(ddof=0) * 3 / 2  # within-tile variance with df correction
    icc = (vb - vw / 3) / (vb + 2 * vw / 3)
    # domain means by year with block CI
    dm = []
    for yi, y in enumerate(yrs):
        t = boot_stat(lambda a: np.nanmean(a), Yo[:, yi][:, None] if False else Yo[:, yi]) if False else None
    ym = []
    for yi, y in enumerate(yrs):
        vals = [Yo[np.isin(np.arange(len(Yo)), np.where(ok)[0][bi_]) , yi].mean() for bi_ in []]
    # year-mean CI via bootstrap on masked index
    okidx = np.where(ok)[0]
    pos = {t: i for i, t in enumerate(okidx)}
    bs = []
    for bi in BOOT[:500]:
        sel = [pos[t] for t in bi if t in pos]
        bs.append(Yo[sel].mean(0))
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    for yi, y in enumerate(yrs):
        varrows.append(dict(family=fam, year=y, n_tiles=int(ok.sum()), domain_mean_pct=Yo[:, yi].mean(), lo=lo[yi], hi=hi[yi],
                            icc_tile=icc, between_tile_sd=np.sqrt(vb), within_tile_year_sd=np.sqrt(vw)))
    # does year anomaly of outcome follow year anomaly of obs count?
    Nan = No - No.mean(1, keepdims=True)
    t = boot_stat(lambda a, b: spearmanr(a.ravel(), b.ravel())[0], an, Nan) if False else None
    # bootstrap needs tile-level arrays: use tile index arrays
    full_an = np.full((len(d), 3), np.nan); full_an[ok] = an
    full_N = np.full((len(d), 3), np.nan); full_N[ok] = Nan
    def f(a, b):
        m = np.isfinite(a).all(1)
        return spearmanr(a[m].ravel(), b[m].ravel())[0]
    pt = f(full_an, full_N)
    vals = [f(full_an[bi], full_N[bi]) for bi in BOOT[:500]]
    a3.append(dict(family=fam, test="year anomaly outcome vs year anomaly n_obs (tile-year)", rho=pt, lo=np.percentile(vals, 2.5), hi=np.percentile(vals, 97.5)))
    # year-to-year stability: between-year corr pairs
    for (i, j) in [(0, 1), (1, 2), (0, 2)]:
        a3.append(dict(family=fam, test=f"tile corr {yrs[i]} vs {yrs[j]}", rho=spearmanr(Yo[:, i], Yo[:, j])[0], lo=np.nan, hi=np.nan))
    # is year-instability (SD across years) related to tile covariates?
    sdy = pd.Series(np.nan, index=range(len(d))); sdy[ok] = Yo.std(1, ddof=1)
    for c in ["op_rs_smterm", "op_std_ratio_geom_to_smterm", "op_slope_med", "vg_vegopacity_clim_mean", "n_l1_obs", "st_canopy_height_m"]:
        t = boot_stat(spear, d[c].values, sdy.values)
        a3.append(dict(family=fam, test=f"SD across years vs {c}", rho=t[0], lo=t[1], hi=t[2]))
pd.DataFrame(varrows).to_csv(OUT + "A3_year_means_variance_decomp.csv", index=False)
pd.DataFrame(a3).to_csv(OUT + "A3_year_tests.csv", index=False)

# figure A: seasonal anomalies vs domain gain
show = ["rs_smterm", "sd_smterm", "snr", "sfmc_std", "LAI", "frac_out15"]
fig, ax = plt.subplots(2, 4, figsize=(15, 6.5))
for k, g in enumerate(["SMAP", "SMOS", "ASCAT", "SMAP_L1minusL3"]):
    ax[0, k].plot(SEAS, gain[g].values, "o-", color="C3")
    ax[0, k].axhline(0, color="k", lw=.5)
    ax[0, k].set_title(f"domain O-F gain, {g}")
    ax[0, k].set_ylabel("% better than OL")
for k, c in enumerate(show[:4]):
    s = a1[a1.covariate == c].set_index("season").reindex(SEAS)
    ax[1, k].errorbar(SEAS, s.mean_anom, yerr=[s.mean_anom - s.anom_lo, s.anom_hi - s.mean_anom], fmt="s-", capsize=3)
    ax[1, k].axhline(0, color="k", lw=.5)
    ax[1, k].set_title(f"within-tile seasonal anomaly: {c}")
fig.suptitle("EXPLORATORY: domain-season gain (4 points) vs within-tile covariate anomalies (block-bootstrap 95% CI)")
fig.tight_layout(); fig.savefig(OUT + "A_seasonal_gain_vs_anomalies.png", dpi=130); plt.close(fig)

# ======================================================================
# PART B: ROUGHNESS
# ======================================================================
d["abs_inc_slope"] = d.op_inc_slope.abs()
d["geom_var_frac"] = d.op_sd_geom**2 / (d.op_sd_geom**2 + d.op_sd_smterm**2)
d["log_nobs"] = np.log(d.n_l1_obs)
proxies = {
    # static terrain / vegetation
    "terrain: CTI std (st_)": "st_cti_std_logm",
    "terrain: DEM std over obs (cv_)": "cv_stat_dem_obs_std_m",
    "terrain: slope mean over obs deg (cv_)": "cv_stat_slope_obs_mean_deg",
    "terrain: slope p90 deg (cv_)": "cv_stat_slope_obs_p90_deg",
    "terrain: slope median (op_)": "op_slope_med",
    "terrain: altitude median (op_)": "op_alt_med",
    "veg: canopy height m (st_)": "st_canopy_height_m",
    "veg: opacity clim (vg_)": "vg_vegopacity_clim_mean",
    "veg: opacity median (op_)": "op_opac_med",
    "ASCAT z0 m (ASCAT's, not IGOT's)": "st_ascat_roughness_z0_m",
    # operator-residual diagnostics
    "resid: |inc slope of O-H| (op_)": "abs_inc_slope",
    "resid: inc slope signed": "op_inc_slope",
    "resid: coherency tail frac": "op_coh_tail_frac",
    "resid: frac |innov|>15 dB": "op_frac_out15",
    # SNR / variance partition
    "var: sd_geom (dB)": "op_sd_geom",
    "var: sd_smterm (dB)": "op_sd_smterm",
    "var: sd_geom/sd_smterm": "op_std_ratio_geom_to_smterm",
    "var: geom variance fraction": "geom_var_frac",
    "var: dyn_ratio": "op_dyn_ratio",
    "var: snr": "op_snr",
    "var: R_impl": "op_R_impl",
    # skill
    "skill: rs_smterm": "op_rs_smterm",
    "skill: rs_H": "op_rs_H",
    "skill: rs_geom": "op_rs_geom",
    "skill: partial geom|inc": "op_pr_geom_inc",
    # confounders / controls
    "ctl: log n_l1_obs": "log_nobs",
    "ctl: latitude": "com_lat",
    "ctl: coh ratio mean": "cv_sup_coh_ratio_mean",
}
b1 = []
for name, col in proxies.items():
    for fam, outs in outcomes.items():
        for o in outs:
            t = boot_stat(spear, d[col].values, d[o].values)
            b1.append(dict(family=fam, outcome=o, proxy=name, col=col, n=int((d[col].notna() & d[o].notna()).sum()),
                           rho=t[0], lo=t[1], hi=t[2]))
b1 = pd.DataFrame(b1)
b1["ci_excludes0"] = (b1.lo > 0) | (b1.hi < 0)
b1.to_csv(OUT + "B1_proxy_spearman.csv", index=False)

# B2: partial spearman controlling for rs_smterm + log n_obs (and a stricter set with lat, opacity)
ctlA = ["op_rs_smterm", "log_nobs"]
ctlB = ["op_rs_smterm", "log_nobs", "com_lat", "vg_vegopacity_clim_mean"]
b2 = []
test_prox = [k for k, c in proxies.items() if not k.startswith(("ctl", "skill: rs_smterm"))]
for name in test_prox:
    col = proxies[name]
    for fam, outs in outcomes.items():
        for o in outs:
            res = dict(family=fam, outcome=o, proxy=name)
            for lab, ctl in [("A", ctlA), ("B", ctlB)]:
                Z = d[ctl].values
                t = boot_stat(lambda x, y, Z: pspear(x, y, Z), d[col].values, d[o].values, Z)
                res[f"pr_{lab}"], res[f"lo_{lab}"], res[f"hi_{lab}"] = t
            b2.append(res)
b2 = pd.DataFrame(b2)
b2["A_excl0"] = (b2.lo_A > 0) | (b2.hi_A < 0)
b2["B_excl0"] = (b2.lo_B > 0) | (b2.hi_B < 0)
b2.to_csv(OUT + "B2_partial_spearman.csv", index=False)

# B3: tercile contrast (top minus bottom tercile mean outcome)
def terc(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    lo, hi = np.percentile(x, [33.3, 66.7])
    return y[x >= hi].mean() - y[x <= lo].mean()
b3 = []
for name in ["resid: |inc slope of O-H| (op_)", "var: sd_geom/sd_smterm", "var: geom variance fraction",
             "terrain: slope median (op_)", "veg: canopy height m (st_)", "resid: coherency tail frac"]:
    for fam, outs in outcomes.items():
        for o in outs:
            t = boot_stat(terc, d[proxies[name]].values, d[o].values)
            b3.append(dict(proxy=name, outcome=o, top_minus_bottom_tercile=t[0], lo=t[1], hi=t[2]))
pd.DataFrame(b3).to_csv(OUT + "B3_tercile_contrasts.csv", index=False)

# B4: incremental R2 of roughness block over base (rank OLS), in-sample and block-CV
static_rough = ["abs_inc_slope", "geom_var_frac", "op_slope_med", "st_cti_std_logm", "st_canopy_height_m", "st_ascat_roughness_z0_m"]
resid_diag = ["op_coh_tail_frac", "op_frac_out15", "op_dyn_ratio", "op_snr"]
base = ["op_rs_smterm", "log_nobs"]
conf = ["com_lat", "vg_vegopacity_clim_mean"]


def rank_mat(df, cols):
    return np.column_stack([rankdata(df[c].values) for c in cols])


def r2(X, y):
    X = np.column_stack([np.ones(len(y)), X])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    return 1 - e.var() / y.var()


def cv_r2(X, y, blk, K=5, rep=10):
    out = []
    ubk = np.unique(blk)
    for r in range(rep):
        f = {b: k for b, k in zip(rng.permutation(ubk), np.arange(len(ubk)) % K)}
        fold = np.array([f[b] for b in blk])
        pred = np.zeros(len(y))
        for k in range(K):
            tr, te = fold != k, fold == k
            Xt = np.column_stack([np.ones(tr.sum()), X[tr]])
            bb = np.linalg.lstsq(Xt, y[tr], rcond=None)[0]
            pred[te] = np.column_stack([np.ones(te.sum()), X[te]]) @ bb
        out.append(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return np.mean(out)


b4 = []
for fam, outs in outcomes.items():
    for o in outs:
        for lab, add, ctl in [("static roughness proxies | base", static_rough, base),
                              ("static roughness proxies | base+lat+opacity", static_rough, base + conf),
                              ("operator-residual diagnostics | base", resid_diag, base),
                              ("geometry share + |inc slope| only | base", ["abs_inc_slope", "geom_var_frac"], base)]:
            cols = list(dict.fromkeys(ctl + add + [o]))
            sub = d.dropna(subset=cols).reset_index(drop=True)
            if len(sub) < 60:
                continue
            y = rankdata(sub[o].values)
            Xc = rank_mat(sub, ctl); Xf = rank_mat(sub, ctl + add)
            dr2 = r2(Xf, y) - r2(Xc, y)
            dcv = cv_r2(Xf, y, sub.block.values) - cv_r2(Xc, y, sub.block.values)
            # block bootstrap of in-sample dR2
            bb_ = sub.block.values; ubk = np.unique(bb_); ib = {b: np.where(bb_ == b)[0] for b in ubk}
            vals = []
            for _ in range(500):
                ix = np.concatenate([ib[b] for b in rng.choice(ubk, len(ubk), replace=True)])
                vals.append(r2(Xf[ix], y[ix]) - r2(Xc[ix], y[ix]))
            # permutation null for dR2: shuffle added columns jointly (rows), keep base fixed
            null = []
            for _ in range(300):
                p = rng.permutation(len(y))
                null.append(r2(np.column_stack([Xc, Xf[p][:, Xc.shape[1]:]]), y) - r2(Xc, y))
            b4.append(dict(family=fam, outcome=o, added=lab, n=len(sub), R2_base=r2(Xc, y), R2_full=r2(Xf, y), dR2=dr2,
                           dR2_lo=np.percentile(vals, 2.5), dR2_hi=np.percentile(vals, 97.5),
                           perm_p=(np.sum(np.array(null) >= dr2) + 1) / 301, blockCV_dR2=dcv))
b4 = pd.DataFrame(b4)
b4.to_csv(OUT + "B4_incremental_R2.csv", index=False)

# B5: standardized rank-regression coefficients with all (full model) for reference
b5 = []
for fam, outs in outcomes.items():
    for o in outs:
        cols = base + static_rough + resid_diag + conf
        sub = d.dropna(subset=cols + [o]).reset_index(drop=True)
        y = rankdata(sub[o].values); y = (y - y.mean()) / y.std()
        X = rank_mat(sub, cols); X = (X - X.mean(0)) / X.std(0)
        def coefs(ix):
            Xt = np.column_stack([np.ones(len(ix)), X[ix]])
            return np.linalg.lstsq(Xt, y[ix], rcond=None)[0][1:]
        pt = coefs(np.arange(len(sub)))
        bb_ = sub.block.values; ubk = np.unique(bb_); ib = {b: np.where(bb_ == b)[0] for b in ubk}
        vals = np.array([coefs(np.concatenate([ib[b] for b in rng.choice(ubk, len(ubk), replace=True)])) for _ in range(500)])
        lo, hi = np.percentile(vals, [2.5, 97.5], axis=0)
        for k, c in enumerate(cols):
            b5.append(dict(family=fam, outcome=o, n=len(sub), term=c, beta=pt[k], lo=lo[k], hi=hi[k]))
b5 = pd.DataFrame(b5); b5["excl0"] = (b5.lo > 0) | (b5.hi < 0)
b5.to_csv(OUT + "B5_full_rank_regression.csv", index=False)

# figures B
fig, axs = plt.subplots(1, 2, figsize=(15, 9), sharey=True)
order = [k for k in proxies]
for ax, (fam, outs) in zip(axs, outcomes.items()):
    for j, o in enumerate(outs):
        s = b1[b1.outcome == o].set_index("proxy").reindex(order)
        yy = np.arange(len(order)) + (j - (len(outs) - 1) / 2) * 0.25
        ax.errorbar(s.rho, yy, xerr=[s.rho - s.lo, s.hi - s.rho], fmt="o", ms=4, capsize=2, label=o)
    ax.axvline(0, color="k", lw=.5); ax.set_title(fam + ": Spearman (block-boot 95% CI)"); ax.legend(fontsize=8)
axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels(order, fontsize=8); axs[0].invert_yaxis()
fig.suptitle("EXPLORATORY: roughness-related proxies vs benefit (positive = L1 better than OL)")
fig.tight_layout(); fig.savefig(OUT + "B1_proxy_spearman_forest.png", dpi=130); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(15, 8), sharey=True)
order2 = test_prox
for ax, (fam, outs) in zip(axs, outcomes.items()):
    for j, o in enumerate(outs):
        s = b2[b2.outcome == o].set_index("proxy").reindex(order2)
        yy = np.arange(len(order2)) + (j - (len(outs) - 1) / 2) * 0.25
        ax.errorbar(s.pr_A, yy, xerr=[s.pr_A - s.lo_A, s.hi_A - s.pr_A], fmt="o", ms=4, capsize=2, label=o)
    ax.axvline(0, color="k", lw=.5); ax.set_title(fam + ": partial Spearman | rs_smterm, log N"); ax.legend(fontsize=8)
axs[0].set_yticks(range(len(order2))); axs[0].set_yticklabels(order2, fontsize=8); axs[0].invert_yaxis()
fig.suptitle("EXPLORATORY: partial associations beyond SM-term skill and obs count")
fig.tight_layout(); fig.savefig(OUT + "B2_partial_spearman_forest.png", dpi=130); plt.close(fig)

# console summary
pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_rows", 500)
print("\n== gains used\n", gain.round(2))
print("\n== A1 mean anomalies (selected)\n", a1[a1.covariate.isin(["rs_smterm", "sd_smterm", "snr", "sfmc_std", "LAI", "frac_out15", "R_impl", "dyn_ratio", "rs_H", "precip_mm", "H_std_db"])]
      .assign(txt=lambda x: x.apply(lambda r: f"{r.mean_anom:+.3f}[{r.anom_lo:+.3f},{r.anom_hi:+.3f}]", axis=1)).pivot(index="covariate", columns="season", values="txt")[SEAS])
print("\n== A1b n=4 rank assoc\n", a1b.round(2))
print("\n== A2 sig\n", a2[a2.sig].round(2))
print("\n== A3 year\n", pd.DataFrame(varrows).round(2)); print(pd.DataFrame(a3).round(2).to_string())
print("\n== B1\n", b1.pivot(index="proxy", columns="outcome", values="rho").reindex(order).round(2))
print("\n== B1 sig\n", b1[b1.ci_excludes0][["proxy", "outcome", "n", "rho", "lo", "hi"]].round(2))
print("\n== B2 sig (A / B)\n", b2[b2.A_excl0 | b2.B_excl0][["proxy", "outcome", "pr_A", "lo_A", "hi_A", "pr_B", "lo_B", "hi_B"]].round(2))
print("\n== B3\n", pd.DataFrame(b3).round(2))
print("\n== B4\n", b4.round(3))
print("\n== B5 sig\n", b5[b5.excl0].round(2))
print("\nA tile-season map match: n tile-seasons", len(rd))
