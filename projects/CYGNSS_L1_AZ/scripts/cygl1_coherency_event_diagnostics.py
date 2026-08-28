#!/usr/bin/env python3
"""
CYGNSS L1 coherency-screening experiment: event-level single-observation
diagnostic (same six-check framework as cygl1_paired_event_diagnostics.py /
cygl1_singleobs_v2_event_diagnostics.py), computed per arm
(coherency_screened, coherency_randmatch), full-period (Jan-Dec 2020 pooled)
and monthly. No local-interaction-count stratification here (that covariate
belongs to the thinning-density family, not this experiment).

For each retained obs:
    d_f        = O - F
    d_a        = O - A
    delta_h    = A - F
    gain_proxy = delta_h / d_f   (only for well-conditioned |d_f| >= 0.1 dB)

Inputs:
  cygl1_coherency_gain_consistency_ofa.csv.gz  (extract_cygl1_coherency_gain_consistency.py)

Outputs:
  cygl1_coherency_event_diagnostics.csv
  cygl1_coherency_event_diagnostics_summary.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_coherency_gain_consistency_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_coherency_event_diagnostics.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_coherency_event_diagnostics_summary.txt")

DFLOOR = 0.1  # dB
ARMS = ["coherency_screened", "coherency_randmatch"]


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
    else:
        lines.append("   (insufficient N or zero-variance K for monotonicity check)")

    lines.append("6. Positive vs. negative innovation (d_f) asymmetry")
    pos = df[df["d_f"] > 0]
    neg = df[df["d_f"] < 0]
    lines.append(f"   N(d_f>0) = {len(pos)}, N(d_f<0) = {len(neg)}")
    if len(pos) and len(neg):
        lines.append(f"   improved rate: pos={pos['improved'].mean():.3f}, neg={neg['improved'].mean():.3f}")
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
    print(f"Loaded {len(df)} events across arms: {df['arm'].value_counts().to_dict()}")

    df["d_f"] = df["obs"] - df["fcst"]
    df["d_a"] = df["obs"] - df["ana"]
    df["delta_h"] = df["ana"] - df["fcst"]
    well_cond = df["d_f"].abs() >= DFLOOR
    df["gain_proxy"] = np.where(well_cond, df["delta_h"] / df["d_f"], np.nan)
    df["abs_d_f"] = df["d_f"].abs()
    df["abs_d_a"] = df["d_a"].abs()
    df["improved"] = df["abs_d_a"] < df["abs_d_f"]
    df["toward_obs"] = df["d_f"] * df["delta_h"] > 0
    df["toward_obs_or_zero"] = df["d_f"] * df["delta_h"] >= 0
    df["month"] = df["datetime"].dt.strftime("%Y-%m")

    df = df.sort_values(["arm", "datetime", "tilenum"]).reset_index(drop=True)

    lines = []
    lines.append("CYGNSS L1 coherency-screening experiment: event-level diagnostics")
    lines.append("arms: coherency_screened (coherency_ratio>=0.523), coherency_randmatch (random-matched control)")
    lines.append("period: Jan-Dec 2020 (full calendar year, no spin-up excluded)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned)")
    lines.append("")

    lines.append("=" * 100)
    lines.append("FULL-PERIOD (Jan-Dec pooled), per arm -- verbose six-check blocks")
    lines.append("=" * 100)
    for arm in ARMS:
        sub = df[df["arm"] == arm]
        lines.extend(verbose_block(sub, f"arm={arm}"))

    lines.append("=" * 100)
    lines.append("MONTHLY, per arm -- condensed")
    lines.append("=" * 100)
    for arm in ARMS:
        lines.append(f"--- arm={arm} ---")
        sub = df[df["arm"] == arm]
        for month, g in sub.groupby("month"):
            lines.append(condensed_row(g, month))
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
