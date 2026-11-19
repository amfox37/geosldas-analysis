"""Roughness-sensitive residual diagnostics per tile x season (local data only).

Input : output/cygnss_da_quality_bundle/per_obs/cygnss_ol_innovations_with_quality_all.csv
        (OL run, 2020-2022, 469k obs; d_f = scaled O minus scaled H(x); fcst = scaled H(x);
         fcstvar = ensemble HPH^T in dB^2; QC fields incl. inc angle, SNR, coherency, DEM, slope).
Output: output/benefit_drivers/roughness_diagnostics_by_tile_season.csv
        output/benefit_drivers/roughness_diagnostics_summary.txt  (domain-level numbers quoted in the note)

Caveat: the table has no geom_db / SM-term split (that lives in the Discover-only
output/operator_skill obs table, not present locally). The ensemble spread fcstvar is
used as the SM-term variance proxy; var(fcst) within tile-season is total H variability.
Tile ids are local 909-tile AZ ids; tile lon/lat = median specular-point position of its obs.

Run: /Users/amfox/mamba/envs/regrid/bin/python -I scripts/benefit_drivers/roughness_diagnostics.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "output/cygnss_da_quality_bundle/per_obs/cygnss_ol_innovations_with_quality_all.csv")
OUT = os.path.join(ROOT, "output/benefit_drivers")
MIN_N, MIN_N_SLOPE, CLIP = 30, 60, 15.0
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
          6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}

cols = ["tilenum", "sp_lon", "sp_lat", "obs_date", "sc_num", "obs", "fcst", "d_f", "fcstvar", "obsvar",
        "sp_inc_angle", "ddm_snr", "coherency_state", "coherency_ratio", "srtm_dem_alt", "srtm_slope",
        "modis_land_cover", "coefficient_weighted_opacity", "MWRTM_VEGCLS"]
df = pd.read_csv(SRC, usecols=cols)
df = df[np.isfinite(df.d_f) & np.isfinite(df.fcst)].copy()
n_all = len(df)
df["month"] = pd.to_datetime(df.obs_date).dt.month
df["season"] = df.month.map(SEASON)
df["inc10"] = df.sp_inc_angle / 10.0
df["lcoh"] = np.log10(df.coherency_ratio.clip(lower=1e-3))
# clipped residual keeps the weak-H outliers (|d_f|>15) from dominating moments
df["r"] = df.d_f.clip(-CLIP, CLIP)


def slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < MIN_N_SLOPE or np.std(x[m]) == 0:
        return np.nan
    return stats.linregress(x[m], y[m]).slope


def mad_std(a):
    return 1.4826 * np.median(np.abs(a - np.median(a)))


rows = []
for (t, s), g in df.groupby(["tilenum", "season"]):
    if len(g) < MIN_N:
        continue
    r = g.r.values
    fv = g.fcstvar.mean()
    rows.append(dict(
        tilenum=t, season=s, n=len(g),
        lon=g.sp_lon.median(), lat=g.sp_lat.median(),
        resid_mean=r.mean(), resid_std=r.std(ddof=1), resid_mad_std=mad_std(g.d_f.values),
        frac_weakH_gt15=(g.d_f.abs() > CLIP).mean(),
        slope_inc_per10deg=slope(g.inc10.values, r),
        slope_snr_per_dB=slope(g.ddm_snr.values, r),
        slope_logcoh=slope(g.lcoh.values, r),
        resid_mean_coherent=g.r[g.coherency_ratio >= 0.5].mean() if (g.coherency_ratio >= 0.5).sum() >= 15 else np.nan,
        resid_mean_diffuse=g.r[g.coherency_ratio < 0.5].mean() if (g.coherency_ratio < 0.5).sum() >= 15 else np.nan,
        frac_diffuse_coh_lt05=(g.coherency_ratio < 0.5).mean(),
        coherency_ratio_median=g.coherency_ratio.median(),
        snr_median=g.ddm_snr.median(), inc_mean=g.sp_inc_angle.mean(),
        fcst_std_within=g.fcst.std(ddof=1),       # total H variability (scaled dB)
        mean_fcstvar=fv,                          # SM-term (ensemble) variance proxy
        resid_var_over_fcstvar=g.r.var(ddof=1) / fv if fv > 0 else np.nan,
        implied_R=g.r.var(ddof=1) - fv,
        dem_mean=g.srtm_dem_alt.mean(), srtm_slope_mean=g.srtm_slope.mean(),
        opacity_mean=g.coefficient_weighted_opacity.mean(),
        landcover_mode=g.modis_land_cover.mode().iat[0] if g.modis_land_cover.notna().any() else np.nan,
    ))
tab = pd.DataFrame(rows)
os.makedirs(OUT, exist_ok=True)
tab.to_csv(os.path.join(OUT, "roughness_diagnostics_by_tile_season.csv"), index=False, float_format="%.5g")

# ---------------- domain-level summaries ----------------
L = []
P = lambda *a: L.append(" ".join(str(x) for x in a))
P(f"obs used {n_all}; tile-season rows {len(tab)}; tiles {tab.tilenum.nunique()}")

# 1. seasonal bias that per-tile scaling leaves (tile-demeaned seasonal means)
tab["w"] = tab.n
tw = tab.groupby("tilenum").apply(lambda g: np.average(g.resid_mean, weights=g.n)).rename("tile_mean")
tab = tab.join(tw, on="tilenum")
tab["resid_mean_anom"] = tab.resid_mean - tab.tile_mean
P("\n[1] Tile-demeaned seasonal residual mean (dB), median [P10,P90] across tiles")
for s in ["DJF", "MAM", "JJA", "SON"]:
    a = tab.loc[tab.season == s, "resid_mean_anom"]
    P(f"  {s}: n_tiles={len(a)} median={a.median():+.3f} P10={a.quantile(.1):+.3f} P90={a.quantile(.9):+.3f}")
pool = df.groupby("season").r.agg(["mean", "std", "count"])
P("  pooled all obs by season (clipped resid): ", pool.round(3).to_dict("index"))

# 2. slopes
P("\n[2] Within tile-season slopes of clipped residual, median [P10,P90] (n tile-seasons)")
for c in ["slope_inc_per10deg", "slope_snr_per_dB", "slope_logcoh"]:
    a = tab[c].dropna()
    P(f"  {c}: n={len(a)} median={a.median():+.3f} P10={a.quantile(.1):+.3f} P90={a.quantile(.9):+.3f} frac>0={np.mean(a > 0):.2f}")
    for s in ["DJF", "MAM", "JJA", "SON"]:
        b = tab.loc[tab.season == s, c].dropna()
        P(f"     {s}: median={b.median():+.3f}")

# 3. do residual bias / slopes covary with terrain / vegetation / season? (tile-season Spearman)
P("\n[3] Spearman(tile-season diagnostic, covariate), n tile-seasons")
covs = ["dem_mean", "srtm_slope_mean", "opacity_mean", "frac_diffuse_coh_lt05", "lat", "lon"]
for d in ["resid_mean_anom", "resid_mean", "resid_std", "slope_inc_per10deg", "slope_snr_per_dB", "resid_var_over_fcstvar"]:
    out = []
    for c in covs:
        m = tab[[d, c]].dropna()
        out.append(f"{c}={stats.spearmanr(m[d], m[c])[0]:+.2f}")
    P(f"  {d}: " + "  ".join(out))
# terrain tercile contrast of the incidence slope and residual std
tab["slope_terc"] = pd.qcut(tab.srtm_slope_mean, 3, labels=["flat", "mid", "rough"])
tab["opac_terc"] = pd.qcut(tab.opacity_mean, 3, labels=["low", "mid", "high"])
for k in ["slope_terc", "opac_terc"]:
    P(f"\n[3b] by {k}: medians")
    P(tab.groupby(k, observed=True)[["resid_std", "slope_inc_per10deg", "slope_snr_per_dB", "resid_mean_anom",
                                      "resid_var_over_fcstvar", "frac_diffuse_coh_lt05"]].median().round(3).to_string())
# season x terrain: is the seasonal anomaly bigger over rough/vegetated tiles?
P("\n[3c] median seasonal anomaly of residual mean (dB) by season x slope tercile")
P(tab.pivot_table(index="season", columns="slope_terc", values="resid_mean_anom", aggfunc="median", observed=True).round(3).to_string())
P("\n[3d] median seasonal anomaly of residual mean (dB) by season x opacity tercile")
P(tab.pivot_table(index="season", columns="opac_terc", values="resid_mean_anom", aggfunc="median", observed=True).round(3).to_string())

# 4. coherent vs diffuse
a = tab.dropna(subset=["resid_mean_coherent", "resid_mean_diffuse"])
P(f"\n[4] coherent(coh>=0.5) minus diffuse(coh<0.5) residual mean, within tile-season (n={len(a)}): "
  f"median={np.median(a.resid_mean_coherent - a.resid_mean_diffuse):+.3f} dB; "
  f"std ratio coherent/diffuse pooled = {df.r[df.coherency_ratio>=.5].std()/df.r[df.coherency_ratio<.5].std():.2f}")
P("  by season (pooled obs mean resid): coherent / diffuse")
for s, g in df.groupby("season"):
    P(f"   {s}: {g.r[g.coherency_ratio>=.5].mean():+.3f} / {g.r[g.coherency_ratio<.5].mean():+.3f}  "
      f"frac_diffuse={np.mean(g.coherency_ratio<.5):.3f}")

# 5. residual variance vs SM-term variance
P("\n[5] residual variance relative to SM-term (ensemble) variance, tile-season medians by season")
P(tab.groupby("season")[["resid_std", "mean_fcstvar", "resid_var_over_fcstvar", "fcst_std_within", "implied_R"]].median().round(3).to_string())
P(f"  overall median resid_var/fcstvar = {tab.resid_var_over_fcstvar.median():.2f}")

# 6. spacecraft
P("\n[6] pooled mean clipped residual by spacecraft (dB):")
P(df.groupby("sc_num").r.agg(["mean", "count"]).round(3).T.to_string())
open(os.path.join(OUT, "roughness_diagnostics_summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
