#!/usr/bin/env python3
"""
Stratify the CYGNSS L1 O-F by native `coherency_state`/`coherency_ratio` (Steps 3-5 of
the coherency follow-up to the L1-vs-L3 operator diagnosis). Reads the join produced by
extract_cygl1_coherency_join.py.

Pre-registered prediction (stated independently of the result, per the brief):
mwRTM_get_lr_reflectivity (the L1 operator's reflectivity model) is flat Fresnel
reflectivity with NO roughness/vegetation correction in code (confirmed during the
earlier diagnosis) -- a coherent-scattering formulation. Physics therefore predicts
coherent-flagged returns (coherency_state==1) should show a TIGHTER/better-fit O-F
(smaller |O-F|, smaller sd) than incoherent-flagged returns (coherency_state==0). If the
opposite holds, that's equally informative -- see runs/cygl1_coherency_stratification.md
for the verdict.

Step 3 (obsvar sanity): obsvar is checked for variability before deciding whether to
report a normalized (z = O-F / sqrt(obsvar)) stratification in addition to raw O-F.

Step 4: bucket by coherency_state (all 4 categories, 0/1/2/3) and by coherency_ratio
(quantile bins, capped at 2-3 given sample-size guidance in the brief). Per bucket:
N, mean/median O-F, O-F sd, mean bias, mean fcst (L1 operator forecast in dB) and mean
co-located CYGNSS_SM_6hr fcst (soil-moisture-equivalent context), plus top-5 tiles by
count and % share (flagged if any tile > 20% of a bucket).

Step 5 (within-tile control, the important refinement): restrict to tiles with BOTH
coherency_state==1 (coherent) and coherency_state==0 (not-coherent) obs in the Jan-Feb
window; compute the coherent-minus-incoherent O-F difference PAIRED BY TILE (mean O-F
within each category, per tile), not pooled naively. This differences out static tile
properties (wetness/water-adjacency/roughness) that could make a naive pooled comparison
misleading, since coherent returns preferentially come from smooth/wet/water-adjacent
surfaces.

Output: printed summary tables (also captured to a text side-file for the markdown
write-up) -- no figures required by the brief for this follow-up.
"""
import gzip
import io
import os

import numpy as np
import pandas as pd
from scipy import stats

CSV_PATH = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygnss_l1_coherency_joined_202001_202002.csv.gz"
OUT_TXT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygl1_coherency_stratification_summary.txt"

STATE_NAMES = {0: "not_coherent", 1: "coherent", 2: "mixed", 3: "indeterminate"}


def read_csv_skip_header(path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    n_hdr = sum(1 for line in lines if line.startswith("#"))
    return pd.read_csv(io.StringIO("".join(lines[n_hdr:])))


def top5_tiles(sub):
    n = len(sub)
    vc = sub["tile_id"].value_counts().head(5)
    out = []
    for tid, cnt in vc.items():
        out.append(f"tile {tid}: {cnt} ({100.0*cnt/n:.1f}%)")
    flag = " <== >20% of bucket" if (len(vc) and vc.iloc[0] / n > 0.20) else ""
    return "; ".join(out) + flag


def bucket_summary(df, groupcol, lines):
    lines.append(f"\n=== Bucketed by {groupcol} ===")
    for key, sub in df.groupby(groupcol, dropna=False):
        n = len(sub)
        lines.append(f"\n-- {groupcol}={key} (label={STATE_NAMES.get(key, key) if groupcol=='coherency_state' else key}) --")
        lines.append(f"   N = {n}")
        lines.append(f"   mean O-F   = {sub['omf'].mean():.4f} dB   (= mean bias)")
        lines.append(f"   median O-F = {sub['omf'].median():.4f} dB")
        lines.append(f"   sd O-F     = {sub['omf'].std():.4f} dB")
        lines.append(f"   mean z(O-F)= {sub['z'].mean():.4f}   sd z(O-F) = {sub['z'].std():.4f}")
        lines.append(f"   mean fcst (L1, dB)      = {sub['fcst_db'].mean():.4f}")
        n_l3 = sub['l3_sm6hr_fcst'].notna().sum()
        if n_l3 > 0:
            lines.append(f"   mean co-located L3 SM6hr fcst = {sub['l3_sm6hr_fcst'].mean():.4f} m3/m3  (N with L3 match={n_l3}/{n})")
        else:
            lines.append(f"   mean co-located L3 SM6hr fcst = N/A (no co-located L3 obs in same tile/window, N=0/{n})")
        lines.append(f"   top-5 tiles: {top5_tiles(sub)}")
        lines.append(f"   n_distinct_tiles = {sub['tile_id'].nunique()}")


def main():
    df = read_csv_skip_header(CSV_PATH)
    df["datetime"] = pd.to_datetime(df["datetime"])
    df["omf"] = df["obs_db"] - df["fcst_db"]

    lines = []
    lines.append(f"Loaded {len(df)} joined CYGNSS L1 obs, {df['tile_id'].nunique()} distinct tiles.")
    lines.append(f"Pooled O-F: mean={df['omf'].mean():.4f}  median={df['omf'].median():.4f}  sd={df['omf'].std():.4f}")

    # ---------------- Step 3: obsvar sanity ----------------
    n_distinct_obsvar = df["obsvar"].round(4).nunique()
    lines.append(f"\n=== Step 3: obsvar sanity ===")
    lines.append(f"obsvar: min={df['obsvar'].min():.4f} max={df['obsvar'].max():.4f} "
                 f"mean={df['obsvar'].mean():.4f} sd={df['obsvar'].sd() if hasattr(df['obsvar'],'sd') else df['obsvar'].std():.4f}")
    lines.append(f"distinct obsvar values (rounded 4dp): {n_distinct_obsvar} / {len(df)}")
    obsvar_is_constant = n_distinct_obsvar <= 1
    lines.append(f"obsvar constant across records? {'YES -- raw and normalized stratifications are identical, reporting raw only' if obsvar_is_constant else 'NO -- obsvar varies substantially, reporting BOTH raw O-F and normalized z=O-F/sqrt(obsvar)'}")
    df["z"] = df["omf"] / np.sqrt(df["obsvar"])

    # ---------------- Step 4: bucket by coherency_state ----------------
    bucket_summary(df, "coherency_state", lines)

    states_present = sorted(df["coherency_state"].dropna().unique().tolist())
    states_missing = sorted(set([0, 1, 2, 3]) - set(states_present))
    lines.append(f"\ncoherency_state categories present in Jan-Feb 2020 AZ-domain data: {states_present}")
    lines.append(f"coherency_state categories NOT observed in this window: {states_missing}")

    # ---------------- Step 4: coherency_ratio ----------------
    lines.append("\n=== coherency_ratio, full population (all states pooled) ===")
    rho_omf, p_omf = stats.spearmanr(df["coherency_ratio"], df["omf"])
    rho_absomf, p_absomf = stats.spearmanr(df["coherency_ratio"], df["omf"].abs())
    lines.append(f"Spearman rho(coherency_ratio, O-F)       = {rho_omf:.4f}  (p={p_omf:.3g}, N={len(df)})")
    lines.append(f"Spearman rho(coherency_ratio, |O-F|)     = {rho_absomf:.4f}  (p={p_absomf:.3g}, N={len(df)})")
    lines.append("  (negative rho for |O-F| would mean higher coherency_ratio -> tighter fit, consistent with the prediction)")

    lines.append("\n=== coherency_ratio quantile bins, WITHIN coherent (state==1) subset only ===")
    coh = df[df["coherency_state"] == 1].copy()
    lines.append(f"N coherent (state==1) = {len(coh)}")
    n_bins = 3 if len(coh) >= 3 * 450 else 2
    try:
        coh["ratio_qbin"] = pd.qcut(coh["coherency_ratio"], n_bins, duplicates="drop")
        for qb, sub in coh.groupby("ratio_qbin", observed=True):
            lines.append(f"  bin {qb}: N={len(sub)}  mean O-F={sub['omf'].mean():.4f}  median O-F={sub['omf'].median():.4f}  sd O-F={sub['omf'].std():.4f}  mean ratio={sub['coherency_ratio'].mean():.3f}")
    except Exception as e:
        lines.append(f"  quantile binning failed: {e}")

    lines.append("\n=== coherency_ratio quantile bins, WITHIN not-coherent (state==0) subset only (for comparison) ===")
    notcoh = df[df["coherency_state"] == 0].copy()
    lines.append(f"N not-coherent (state==0) = {len(notcoh)}")
    notcoh["ratio_qbin"] = pd.qcut(notcoh["coherency_ratio"], 3, duplicates="drop")
    for qb, sub in notcoh.groupby("ratio_qbin", observed=True):
        lines.append(f"  bin {qb}: N={len(sub)}  mean O-F={sub['omf'].mean():.4f}  median O-F={sub['omf'].median():.4f}  sd O-F={sub['omf'].std():.4f}  mean ratio={sub['coherency_ratio'].mean():.3f}")

    # ---------------- Step 5: within-tile paired control ----------------
    lines.append("\n=== Step 5: within-tile coherent-minus-incoherent O-F control ===")
    coh_tiles = set(df.loc[df["coherency_state"] == 1, "tile_id"].unique())
    noncoh_tiles = set(df.loc[df["coherency_state"] == 0, "tile_id"].unique())
    both_tiles = sorted(coh_tiles & noncoh_tiles)
    lines.append(f"Tiles with coherency_state==1 obs: {len(coh_tiles)}")
    lines.append(f"Tiles with coherency_state==0 obs: {len(noncoh_tiles)}")
    lines.append(f"Tiles with BOTH (qualify for within-tile control): {len(both_tiles)}")

    diffs = []
    per_tile_rows = []
    for tid in both_tiles:
        sub_c = df[(df["tile_id"] == tid) & (df["coherency_state"] == 1)]
        sub_n = df[(df["tile_id"] == tid) & (df["coherency_state"] == 0)]
        d = sub_c["omf"].mean() - sub_n["omf"].mean()
        d_absomf = sub_c["omf"].abs().mean() - sub_n["omf"].abs().mean()
        d_sd = sub_c["omf"].std() - sub_n["omf"].std()
        diffs.append(d)
        per_tile_rows.append((
            tid, len(sub_c), len(sub_n), sub_c["omf"].mean(), sub_n["omf"].mean(), d,
            sub_c["omf"].abs().mean(), sub_n["omf"].abs().mean(), d_absomf,
            sub_c["omf"].std(), sub_n["omf"].std(), d_sd,
        ))

    if len(diffs) > 0:
        diffs = np.array(diffs)
        pt_df = pd.DataFrame(per_tile_rows, columns=[
            "tile_id", "n_coherent", "n_incoherent", "mean_omf_coherent", "mean_omf_incoherent",
            "diff_coherent_minus_incoherent",
            "mean_absomf_coherent", "mean_absomf_incoherent", "diff_absomf_coherent_minus_incoherent",
            "sd_omf_coherent", "sd_omf_incoherent", "diff_sd_coherent_minus_incoherent",
        ])
        lines.append(f"within-tile (coherent mean O-F) - (incoherent mean O-F) [SIGNED BIAS], N_tiles={len(diffs)}:")
        lines.append(f"  mean diff   = {diffs.mean():.4f} dB")
        lines.append(f"  median diff = {np.median(diffs):.4f} dB")
        lines.append(f"  sd diff     = {diffs.std():.4f} dB")
        lines.append(f"  N tiles with diff < 0 (coherent LESS positive-biased than incoherent): {int((diffs<0).sum())} / {len(diffs)} ({100.0*(diffs<0).mean():.1f}%)")
        if len(diffs) >= 8:
            tstat, pval = stats.wilcoxon(diffs)
            lines.append(f"  Wilcoxon signed-rank test on signed diff (!= 0): stat={tstat:.3f} p={pval:.4g}")

        d_abs = pt_df["diff_absomf_coherent_minus_incoherent"].values
        d_sd_arr = pt_df["diff_sd_coherent_minus_incoherent"].values
        lines.append(f"\nwithin-tile mean(|O-F|)_coherent - mean(|O-F|)_incoherent [ERROR-MAGNITUDE / 'tightness'], N_tiles={len(d_abs)}:")
        lines.append(f"  mean diff = {d_abs.mean():.4f} dB   (negative => coherent obs have SMALLER |O-F| than incoherent in the same tile, i.e. tighter fit, consistent with the prediction)")
        lines.append(f"  median diff = {np.median(d_abs):.4f} dB")
        lines.append(f"  N tiles with coherent tighter (mean|O-F| smaller): {int((d_abs<0).sum())} / {len(d_abs)} ({100.0*(d_abs<0).mean():.1f}%)")
        if len(d_abs) >= 8:
            tstat2, pval2 = stats.wilcoxon(d_abs)
            lines.append(f"  Wilcoxon signed-rank test on |O-F| diff (!= 0): stat={tstat2:.3f} p={pval2:.4g}")

        d_sd_valid = d_sd_arr[~np.isnan(d_sd_arr)]
        n_sd_nan = int(np.isnan(d_sd_arr).sum())
        lines.append(f"\nwithin-tile sd(O-F)_coherent - sd(O-F)_incoherent [SPREAD], N_tiles={len(d_sd_valid)} (dropped {n_sd_nan} tiles with <2 obs in a category, sd undefined):")
        lines.append(f"  mean diff = {d_sd_valid.mean():.4f} dB   (negative => coherent obs have LOWER spread than incoherent in the same tile, consistent with the prediction)")
        lines.append(f"  median diff = {np.median(d_sd_valid):.4f} dB")
        lines.append(f"  N tiles with coherent lower spread: {int((d_sd_valid<0).sum())} / {len(d_sd_valid)} ({100.0*(d_sd_valid<0).mean():.1f}%)")
        if len(d_sd_valid) >= 8:
            tstat3, pval3 = stats.wilcoxon(d_sd_valid)
            lines.append(f"  Wilcoxon signed-rank test on sd diff (!= 0): stat={tstat3:.3f} p={pval3:.4g}")

        lines.append("\nRobustness to per-tile sample size (min N in EACH category):")
        for minn in (1, 3, 5, 10):
            m = (pt_df["n_coherent"] >= minn) & (pt_df["n_incoherent"] >= minn)
            s = pt_df[m]
            if len(s) == 0:
                continue
            lines.append(
                f"  min_n={minn:2d}: N_tiles={len(s):3d}  "
                f"signed_diff mean={s['diff_coherent_minus_incoherent'].mean():+.3f} median={s['diff_coherent_minus_incoherent'].median():+.3f}  "
                f"|O-F|_diff mean={s['diff_absomf_coherent_minus_incoherent'].mean():+.3f}  "
                f"sd_diff mean={s['diff_sd_coherent_minus_incoherent'].mean():+.3f}"
            )

        pt_path = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygl1_coherency_within_tile_pairs.csv"
        pt_df.to_csv(pt_path, index=False)
        lines.append(f"\n  per-tile pairs written to {pt_path}")
    else:
        lines.append("  NO tiles qualify -- cannot compute within-tile control.")

    # ---------------- also |O-F| tightness comparison, pooled (naive) ----------------
    lines.append("\n=== Naive pooled comparison (for contrast with within-tile control) ===")
    for st in [0, 1, 2]:
        sub = df[df["coherency_state"] == st]
        if len(sub) == 0:
            continue
        lines.append(f"  state={st} ({STATE_NAMES[st]}): N={len(sub)}  mean|O-F|={sub['omf'].abs().mean():.4f}  sd(O-F)={sub['omf'].std():.4f}")

    text = "\n".join(lines)
    print(text)
    with open(OUT_TXT, "w") as fo:
        fo.write(text + "\n")
    print(f"\nWrote summary to {OUT_TXT}")


if __name__ == "__main__":
    main()
