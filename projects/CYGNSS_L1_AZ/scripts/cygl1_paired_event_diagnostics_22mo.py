#!/usr/bin/env python3
"""
22-month extension of cygl1_paired_event_diagnostics.py, restricted to the
DA-intermediate arm only (the only DA arm actually extended past Oct 2020),
May2020-Dec2021. Same six-check event-level diagnostic as the original, full
period + monthly-condensed tables. No local-interaction-count stratification
this time (that covariate's lookup table was only ever built over the
Jan-Oct 2020 thinning run and does not cover the extension months; the
monthly table below serves the same "does it degrade over time" purpose).

Input: cygl1_paired_gain_consistency_22mo_ofa.csv.gz
       (extract_cygl1_paired_gain_consistency_22mo.py)
Outputs: cygl1_paired_event_diagnostics_22mo.csv / _summary.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_paired_gain_consistency_22mo_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_paired_event_diagnostics_22mo.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_paired_event_diagnostics_22mo_summary.txt")

DFLOOR = 0.1  # dB
ARM = "intermediate"


def verbose_block(df, label):
    lines = []
    n = len(df)
    well_cond = df["d_f"].abs() >= DFLOOR
    n_wc = int(well_cond.sum())
    lines.append(f"--- {label} ---  N = {n}, well-conditioned (|d_f|>={DFLOOR}) N = {n_wc}")

    mean_abs_df = df["abs_d_f"].mean()
    mean_abs_da = df["abs_d_a"].mean()
    rmse_df = np.sqrt((df["d_f"] ** 2).mean())
    rmse_da = np.sqrt((df["d_a"] ** 2).mean())
    n_improved = int(df["improved"].sum())
    pct_reduction = 100 * (1 - mean_abs_da / mean_abs_df) if mean_abs_df > 0 else np.nan
    lines.append("1/2. |O-A| vs |O-F|")
    lines.append(f"   mean|d_f| = {mean_abs_df:.4f} dB, mean|d_a| = {mean_abs_da:.4f} dB "
                  f"({'IMPROVED' if mean_abs_da < mean_abs_df else 'NOT improved'}, {pct_reduction:.1f}% reduction)")
    if rmse_df > 0:
        lines.append(f"   RMSE(O-F) = {rmse_df:.4f} dB, RMSE(O-A) = {rmse_da:.4f} dB "
                      f"({100*(1-rmse_da/rmse_df):.1f}% reduction)")
    lines.append(f"   |O-A| < |O-F| (improved) in {n_improved}/{n} obs ({100*n_improved/n:.1f}%)")

    n_toward = int(df["toward_obs"].sum())
    n_toward_or_zero = int(df["toward_obs_or_zero"].sum())
    lines.append("3. d_f * delta_h > 0 (analysis moves toward the observation)")
    lines.append(f"   toward obs (strict >0): {n_toward}/{n} ({100*n_toward/n:.1f}%)")
    lines.append(f"   toward or zero (>=0): {n_toward_or_zero}/{n} ({100*n_toward_or_zero/n:.1f}%)")

    gp = df.loc[well_cond, "gain_proxy"]
    if len(gp) > 0:
        n_neg = int((gp < 0).sum())
        n_in01 = int(((gp >= 0) & (gp <= 1)).sum())
        n_gt1 = int((gp > 1).sum())
        lines.append("4. gain_proxy = delta_h / d_f distribution (well-conditioned events only)")
        lines.append(f"   N = {len(gp)}, mean = {gp.mean():.4f}, median = {gp.median():.4f}, std = {gp.std():.4f}")
        lines.append(f"   negative (wrong direction): {n_neg}/{len(gp)} ({100*n_neg/len(gp):.1f}%)")
        lines.append(f"   in [0,1] (textbook single-obs range): {n_in01}/{len(gp)} ({100*n_in01/len(gp):.1f}%)")
        lines.append(f"   >1 (overshoot): {n_gt1}/{len(gp)} ({100*n_gt1/len(gp):.1f}%)")

    lines.append("5. Monotonicity of gain_proxy in naive K = fcstvar/(fcstvar+obsvar)")
    sub = df.loc[well_cond].copy()
    if len(sub) > 3 and sub["K"].std() > 0:
        rho, pval = stats.spearmanr(sub["K"], sub["gain_proxy"])
        lines.append(f"   Spearman rho(K, gain_proxy) = {rho:.3f} (p={pval:.3f}), N={len(sub)}")
        try:
            sub["K_tertile"] = pd.qcut(sub["K"], 3, labels=["low K", "mid K", "high K"], duplicates="drop")
            tert = sub.groupby("K_tertile", observed=True)["gain_proxy"].agg(["count", "mean", "median"])
            for lbl, row in tert.iterrows():
                lines.append(f"     {lbl:8s}: N={int(row['count']):4d}  mean={row['mean']:.4f}  median={row['median']:.4f}")
        except ValueError:
            lines.append("   (qcut failed -- insufficient distinct K values)")
    else:
        lines.append("   (insufficient N or zero-variance K for monotonicity check)")

    lines.append("6. Positive vs. negative innovation (d_f) asymmetry")
    pos = df[df["d_f"] > 0]
    neg = df[df["d_f"] < 0]
    lines.append(f"   N(d_f>0) = {len(pos)}, N(d_f<0) = {len(neg)}")
    if len(pos) and len(neg):
        lines.append(f"   improved rate: pos={pos['improved'].mean():.3f}, neg={neg['improved'].mean():.3f}")
        lines.append(f"   toward-obs rate: pos={pos['toward_obs'].mean():.3f}, neg={neg['toward_obs'].mean():.3f}")
    pos_gp = df.loc[well_cond & (df["d_f"] > 0), "gain_proxy"]
    neg_gp = df.loc[well_cond & (df["d_f"] < 0), "gain_proxy"]
    if len(pos_gp) and len(neg_gp):
        lines.append(f"   gain_proxy: pos mean={pos_gp.mean():.4f} (N={len(pos_gp)}), "
                      f"neg mean={neg_gp.mean():.4f} (N={len(neg_gp)})")
        if len(pos_gp) > 2 and len(neg_gp) > 2:
            u_stat, u_p = stats.mannwhitneyu(pos_gp, neg_gp, alternative="two-sided")
            lines.append(f"   Mann-Whitney U (pos vs neg gain_proxy): U={u_stat:.1f}, p={u_p:.3f} "
                          f"({'no significant difference' if u_p > 0.05 else 'SIGNIFICANT difference'})")
    lines.append("")
    return lines


def condensed_row(df, label):
    n = len(df)
    if n == 0:
        return f"{label:14s} N=0"
    well_cond = df["d_f"].abs() >= DFLOOR
    gp = df.loc[well_cond, "gain_proxy"]
    mean_abs_df = df["abs_d_f"].mean()
    mean_abs_da = df["abs_d_a"].mean()
    pct_red = 100 * (1 - mean_abs_da / mean_abs_df) if mean_abs_df > 0 else np.nan
    n_improved = int(df["improved"].sum())
    pct_improved = 100 * n_improved / n
    if len(gp) > 0:
        gp_mean = gp.mean()
        gp_median = gp.median()
        pct_neg = 100 * (gp < 0).sum() / len(gp)
        pct_in01 = 100 * ((gp >= 0) & (gp <= 1)).sum() / len(gp)
        pct_over = 100 * (gp > 1).sum() / len(gp)
    else:
        gp_mean = gp_median = pct_neg = pct_in01 = pct_over = np.nan
    return (f"{label:14s} N={n:6d}  mean|d_f|={mean_abs_df:7.4f} mean|d_a|={mean_abs_da:7.4f} "
            f"(red={pct_red:6.1f}%)  %improved={pct_improved:6.1f}  "
            f"gain_proxy: mean={gp_mean:7.4f} med={gp_median:7.4f} "
            f"%neg={pct_neg:5.1f} %in01={pct_in01:5.1f} %over={pct_over:5.1f}  (Nwc={len(gp)})")


def main():
    df = pd.read_csv(IN_CSV, comment="#")
    df["datetime"] = pd.to_datetime(df["datetime"])
    print(f"Loaded {len(df)} events, arm={ARM}")

    df["d_f"] = df["obs"] - df["fcst"]
    df["d_a"] = df["obs"] - df["ana"]
    df["delta_h"] = df["ana"] - df["fcst"]
    assert np.allclose(df["d_f"], df["d"]), "d_f mismatch vs precomputed d"
    assert np.allclose(df["delta_h"], df["actual_incr"]), "delta_h mismatch vs precomputed actual_incr"

    well_cond = df["d_f"].abs() >= DFLOOR
    df["gain_proxy"] = np.where(well_cond, df["delta_h"] / df["d_f"], np.nan)
    df["abs_d_f"] = df["d_f"].abs()
    df["abs_d_a"] = df["d_a"].abs()
    df["improved"] = df["abs_d_a"] < df["abs_d_f"]
    df["toward_obs"] = df["d_f"] * df["delta_h"] > 0
    df["toward_obs_or_zero"] = df["d_f"] * df["delta_h"] >= 0
    df["month"] = df["datetime"].dt.strftime("%Y-%m")

    df = df.sort_values(["datetime", "tilenum"]).reset_index(drop=True)

    lines = []
    lines.append("CYGNSS L1 paired thinning-density experiment: event-level diagnostics, 22-month extension")
    lines.append("arm: intermediate (min_sep=2.40deg) -- the only DA arm extended past Oct 2020")
    lines.append("period: May2020-Dec2021 (evaluation window; Jan-Apr 2020 spin-up excluded)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned)")
    lines.append("")

    lines.append("=" * 100)
    lines.append("FULL-PERIOD (May2020-Dec2021 pooled) -- verbose six-check block")
    lines.append("=" * 100)
    lines.extend(verbose_block(df, f"arm={ARM}"))

    lines.append("=" * 100)
    lines.append("MONTHLY -- condensed (watch for any drift/degradation as cycles accumulate)")
    lines.append("=" * 100)
    for month, g in df.groupby("month"):
        lines.append(condensed_row(g, month))
    lines.append("")

    # yearly split too (2020 H2 vs 2021 full) for an easy at-a-glance comparison
    lines.append("=" * 100)
    lines.append("YEARLY SPLIT -- condensed")
    lines.append("=" * 100)
    df["year"] = df["datetime"].dt.year
    for year, g in df.groupby("year"):
        lines.append(condensed_row(g, str(year)))
    lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["arm", "datetime", "tilenum", "tile_id", "lon", "lat", "obs", "fcst", "ana",
                  "obsvar", "fcstvar", "K", "d_f", "d_a", "delta_h", "gain_proxy",
                  "improved", "toward_obs", "month"]
    df[event_cols].to_csv(EVENT_CSV, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return df


if __name__ == "__main__":
    main()
