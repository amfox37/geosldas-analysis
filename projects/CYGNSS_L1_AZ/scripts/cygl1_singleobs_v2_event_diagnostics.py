#!/usr/bin/env python3
"""
Event-level single-observation diagnostic for the CYGNSS L1 isolated-obs
dataset (Step 1.5 v2, DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs_v2,
job 58019976, errstd=3.0 dB, N=49 genuinely isolated obs, >5deg from every
other kept obs so no localization/HPHt confound).

For each retained obs:
    d_f        = O - F
    d_a        = O - A
    delta_h    = A - F
    gain_proxy = delta_h / d_f

Checks requested:
  1. |O-A| closer to zero than |O-F| (aggregate)
  2. |O-A| < |O-F| for most observations (per-event)
  3. d_f * delta_h > 0 (analysis moves toward the observation)
  4. 0 < gain_proxy < 1 usually; negative = wrong direction, >1 = overshoot
  5. Stronger weighting (higher naive K) -> larger gain_proxy, monotonically
  6. No systematic difference between positive and negative innovations

Event-level output (not annual mean/std) written to CSV; summary stats to txt.
"""
import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_singleobs_v2_gain_consistency_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_singleobs_v2_event_diagnostics.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_singleobs_v2_event_diagnostics_summary.txt")

DFLOOR = 0.1  # dB; below this, gain_proxy is treated as ill-conditioned (near-zero denominator)


def main():
    df = pd.read_csv(IN_CSV, comment="#")
    df["datetime"] = pd.to_datetime(df["datetime"])

    df["d_f"] = df["obs"] - df["fcst"]
    df["d_a"] = df["obs"] - df["ana"]
    df["delta_h"] = df["ana"] - df["fcst"]

    # sanity: matches the already-extracted d / actual_incr columns exactly
    assert np.allclose(df["d_f"], df["d"]), "d_f mismatch vs precomputed d"
    assert np.allclose(df["delta_h"], df["actual_incr"]), "delta_h mismatch vs precomputed actual_incr"

    well_cond = df["d_f"].abs() >= DFLOOR
    df["gain_proxy"] = np.where(well_cond, df["delta_h"] / df["d_f"], np.nan)

    df["abs_d_f"] = df["d_f"].abs()
    df["abs_d_a"] = df["d_a"].abs()
    df["improved"] = df["abs_d_a"] < df["abs_d_f"]
    df["toward_obs"] = df["d_f"] * df["delta_h"] > 0
    df["toward_obs_or_zero"] = df["d_f"] * df["delta_h"] >= 0

    df = df.sort_values(["datetime", "tilenum"]).reset_index(drop=True)

    n = len(df)
    n_wc = int(well_cond.sum())

    lines = []
    lines.append("CYGNSS L1 single-observation event-level diagnostic")
    lines.append(f"experiment: DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs_v2 (job 58019976)")
    lines.append(f"errstd = 3.0 dB, N = {n} genuinely isolated obs (>5deg separation, no HPHt/localization confound)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned): N = {n_wc}/{n}")
    lines.append("")

    # --- 1/2: |O-A| vs |O-F| ---
    mean_abs_df = df["abs_d_f"].mean()
    mean_abs_da = df["abs_d_a"].mean()
    rmse_df = np.sqrt((df["d_f"] ** 2).mean())
    rmse_da = np.sqrt((df["d_a"] ** 2).mean())
    n_improved = int(df["improved"].sum())
    lines.append("1/2. |O-A| vs |O-F|")
    lines.append(f"   mean|d_f| = {mean_abs_df:.4f} dB, mean|d_a| = {mean_abs_da:.4f} dB "
                  f"({'IMPROVED' if mean_abs_da < mean_abs_df else 'NOT improved'} on average, "
                  f"{100*(1-mean_abs_da/mean_abs_df):.1f}% reduction)")
    lines.append(f"   RMSE(O-F) = {rmse_df:.4f} dB, RMSE(O-A) = {rmse_da:.4f} dB "
                  f"({100*(1-rmse_da/rmse_df):.1f}% reduction)")
    lines.append(f"   |O-A| < |O-F| (improved) in {n_improved}/{n} obs ({100*n_improved/n:.1f}%)")
    lines.append("")

    # --- 3: d_f * delta_h > 0 ---
    n_toward = int(df["toward_obs"].sum())
    n_toward_or_zero = int(df["toward_obs_or_zero"].sum())
    lines.append("3. d_f * delta_h > 0 (analysis moves toward the observation)")
    lines.append(f"   toward obs (strict >0): {n_toward}/{n} ({100*n_toward/n:.1f}%)")
    lines.append(f"   toward or zero (>=0): {n_toward_or_zero}/{n} ({100*n_toward_or_zero/n:.1f}%)")
    wrong_dir = df[~df["toward_obs_or_zero"]]
    if len(wrong_dir):
        lines.append(f"   WRONG-DIRECTION events ({len(wrong_dir)}):")
        for _, r in wrong_dir.iterrows():
            lines.append(f"     {r['datetime']} tile {int(r['tilenum'])}: d_f={r['d_f']:.4f}, delta_h={r['delta_h']:.4f}")
    lines.append("")

    # --- 4: gain_proxy distribution ---
    gp = df.loc[well_cond, "gain_proxy"]
    n_neg = int((gp < 0).sum())
    n_in01 = int(((gp >= 0) & (gp <= 1)).sum())
    n_gt1 = int((gp > 1).sum())
    lines.append("4. gain_proxy = delta_h / d_f distribution (well-conditioned events only)")
    lines.append(f"   N = {len(gp)}, mean = {gp.mean():.4f}, median = {gp.median():.4f}, std = {gp.std():.4f}")
    lines.append(f"   min = {gp.min():.4f}, max = {gp.max():.4f}")
    lines.append(f"   quantiles: 10%={gp.quantile(.1):.4f}  25%={gp.quantile(.25):.4f}  "
                  f"75%={gp.quantile(.75):.4f}  90%={gp.quantile(.9):.4f}")
    lines.append(f"   negative (wrong direction): {n_neg}/{len(gp)} ({100*n_neg/len(gp):.1f}%)")
    lines.append(f"   in [0,1] (textbook single-obs range): {n_in01}/{len(gp)} ({100*n_in01/len(gp):.1f}%)")
    lines.append(f"   >1 (overshoot): {n_gt1}/{len(gp)} ({100*n_gt1/len(gp):.1f}%)")
    lines.append("")

    # --- 5: monotonic in naive K (proxy for "stronger weighting") ---
    lines.append("5. Monotonicity of gain_proxy in naive K = fcstvar/(fcstvar+obsvar)")
    lines.append(f"   obsvar range: {df['obsvar'].min():.4f} - {df['obsvar'].max():.4f} "
                  f"({'CONSTANT' if np.isclose(df['obsvar'].min(), df['obsvar'].max()) else 'VARIES'} "
                  f"-- errstd fixed at 3.0dB so obsvar variation, if any, comes from per-tile/species obs-error floors)")
    lines.append(f"   fcstvar range: {df['fcstvar'].min():.4f} - {df['fcstvar'].max():.4f} (varies by tile/cycle ensemble spread)")
    sub = df.loc[well_cond].copy()
    if len(sub) > 3:
        rho, pval = stats.spearmanr(sub["K"], sub["gain_proxy"])
        pear_r, pear_p = stats.pearsonr(sub["K"], sub["gain_proxy"])
        lines.append(f"   Spearman rho(K, gain_proxy) = {rho:.3f} (p={pval:.3f}), "
                      f"Pearson r = {pear_r:.3f} (p={pear_p:.3f})")
        sub["K_tertile"] = pd.qcut(sub["K"], 3, labels=["low K", "mid K", "high K"])
        tert = sub.groupby("K_tertile", observed=True)["gain_proxy"].agg(["count", "mean", "median"])
        lines.append("   gain_proxy by K tertile (expect monotonically increasing mean if weighting is behaving):")
        for lbl, row in tert.iterrows():
            lines.append(f"     {lbl:8s}: N={int(row['count']):3d}  mean={row['mean']:.4f}  median={row['median']:.4f}")
    lines.append("")

    # --- 6: positive vs negative innovation asymmetry ---
    lines.append("6. Positive vs. negative innovation (d_f) asymmetry")
    pos = df[df["d_f"] > 0]
    neg = df[df["d_f"] < 0]
    lines.append(f"   N(d_f>0) = {len(pos)}, N(d_f<0) = {len(neg)}")
    lines.append(f"   improved rate: pos={pos['improved'].mean():.3f}, neg={neg['improved'].mean():.3f}")
    lines.append(f"   toward-obs rate: pos={pos['toward_obs'].mean():.3f}, neg={neg['toward_obs'].mean():.3f}")
    pos_gp = df.loc[well_cond & (df["d_f"] > 0), "gain_proxy"]
    neg_gp = df.loc[well_cond & (df["d_f"] < 0), "gain_proxy"]
    lines.append(f"   gain_proxy: pos mean={pos_gp.mean():.4f} (N={len(pos_gp)}), "
                  f"neg mean={neg_gp.mean():.4f} (N={len(neg_gp)})")
    if len(pos_gp) > 2 and len(neg_gp) > 2:
        u_stat, u_p = stats.mannwhitneyu(pos_gp, neg_gp, alternative="two-sided")
        lines.append(f"   Mann-Whitney U test (pos vs neg gain_proxy): U={u_stat:.1f}, p={u_p:.3f} "
                      f"({'no significant difference' if u_p > 0.05 else 'SIGNIFICANT difference'})")
    lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["datetime", "tilenum", "tile_id", "lon", "lat", "obs", "fcst", "ana",
                  "obsvar", "fcstvar", "K", "d_f", "d_a", "delta_h", "gain_proxy",
                  "improved", "toward_obs"]
    df[event_cols].to_csv(EVENT_CSV, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return df


if __name__ == "__main__":
    main()
