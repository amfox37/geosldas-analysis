#!/usr/bin/env python3
"""
CYGNSS L1 coherency-screening experiment, multi-seed re-draw: same six-check
event-level diagnostic as cygl1_coherency_event_diagnostics.py, but across
all 4 arms (Arm A = coherency_screened; Arm B seed1/seed2/seed3 =
coherency_randmatch[/_seed2/_seed3]). Full-period (Jan-Dec 2020 pooled) only
-- monthly breakdown is skipped here since the question is single-number
seed-to-seed spread, not within-arm seasonality (already checked for the
original single-seed pair).

Also emits a seed-spread summary block: mean/std/min-max range across the 3
Arm-B seeds for each key metric, and an explicit check of whether Arm A's
value falls outside that range.

Inputs:
  cygl1_coherency_multiseed_gain_consistency_ofa.csv.gz
  (extract_cygl1_coherency_multiseed_gain_consistency.py)

Outputs:
  cygl1_coherency_multiseed_event_diagnostics.csv
  cygl1_coherency_multiseed_event_diagnostics_summary.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_coherency_multiseed_gain_consistency_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_coherency_multiseed_event_diagnostics.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_coherency_multiseed_event_diagnostics_summary.txt")

DFLOOR = 0.1  # dB
ARM_A = "coherency_screened"
ARMS_B = ["coherency_randmatch", "coherency_randmatch_seed2", "coherency_randmatch_seed3"]
ARMS = [ARM_A] + ARMS_B


def compute_metrics(df):
    """Return a dict of scalar metrics for one arm's event dataframe."""
    n = len(df)
    well_cond = df["d_f"].abs() >= DFLOOR
    gp = df.loc[well_cond, "gain_proxy"]

    mean_abs_df = df["abs_d_f"].mean()
    mean_abs_da = df["abs_d_a"].mean()
    pct_red = 100 * (1 - mean_abs_da / mean_abs_df) if mean_abs_df > 0 else np.nan
    n_improved = int(df["improved"].sum())
    pct_improved = 100 * n_improved / n
    n_toward = int(df["toward_obs"].sum())
    pct_toward = 100 * n_toward / n

    if len(gp) > 3 and df.loc[well_cond, "K"].std() > 0:
        rho, pval = stats.spearmanr(df.loc[well_cond, "K"], gp)
    else:
        rho, pval = np.nan, np.nan

    return {
        "N": n,
        "N_wellcond": len(gp),
        "mean_abs_d_f": mean_abs_df,
        "mean_abs_d_a": mean_abs_da,
        "pct_reduction": pct_red,
        "pct_improved": pct_improved,
        "pct_toward_obs": pct_toward,
        "gain_proxy_mean": gp.mean() if len(gp) else np.nan,
        "gain_proxy_median": gp.median() if len(gp) else np.nan,
        "pct_gp_negative": 100 * (gp < 0).sum() / len(gp) if len(gp) else np.nan,
        "pct_gp_in01": 100 * ((gp >= 0) & (gp <= 1)).sum() / len(gp) if len(gp) else np.nan,
        "pct_gp_over1": 100 * (gp > 1).sum() / len(gp) if len(gp) else np.nan,
        "rho_K_gainproxy": rho,
        "rho_pval": pval,
    }


def condensed_row(label, m):
    return (f"{label:28s} N={m['N']:6d}  mean|d_f|={m['mean_abs_d_f']:7.4f} "
            f"mean|d_a|={m['mean_abs_d_a']:7.4f} (red={m['pct_reduction']:6.2f}%)  "
            f"%improved={m['pct_improved']:6.2f}  %toward_obs={m['pct_toward_obs']:6.2f}  "
            f"gain_proxy: mean={m['gain_proxy_mean']:7.4f} med={m['gain_proxy_median']:7.4f} "
            f"%neg={m['pct_gp_negative']:5.1f} %in01={m['pct_gp_in01']:5.1f} %over={m['pct_gp_over1']:5.1f}  "
            f"rho(K,gp)={m['rho_K_gainproxy']:.3f} (p={m['rho_pval']:.3g})  (Nwc={m['N_wellcond']})")


def seed_spread_block(metrics_by_arm, key, label, higher_is_more_skill=True):
    a_val = metrics_by_arm[ARM_A][key]
    b_vals = np.array([metrics_by_arm[a][key] for a in ARMS_B])
    b_mean, b_std = b_vals.mean(), b_vals.std(ddof=1)
    b_min, b_max = b_vals.min(), b_vals.max()
    outside = (a_val < b_min) or (a_val > b_max)
    direction = ""
    if outside:
        direction = "ABOVE B-range" if a_val > b_max else "BELOW B-range"
    lines = [
        f"{label}:",
        f"  Arm A           = {a_val:.4f}",
        f"  Arm B seed1/2/3 = {b_vals[0]:.4f} / {b_vals[1]:.4f} / {b_vals[2]:.4f}",
        f"  Arm B mean±std  = {b_mean:.4f} ± {b_std:.4f}  (range {b_min:.4f} to {b_max:.4f})",
        f"  Arm A vs B-seed range: {'OUTSIDE (' + direction + ')' if outside else 'WITHIN B-seed range (overlaps, not distinguishable from seed noise)'}",
    ]
    return lines


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

    df = df.sort_values(["arm", "datetime", "tilenum"]).reset_index(drop=True)

    metrics_by_arm = {arm: compute_metrics(df[df["arm"] == arm]) for arm in ARMS}

    lines = []
    lines.append("CYGNSS L1 coherency-screening experiment: MULTI-SEED event-level diagnostics")
    lines.append("Arm A = coherency_screened (coherency_ratio>=0.523)")
    lines.append("Arm B seed1/seed2/seed3 = coherency_randmatch[/_seed2/_seed3] (random-matched control, 3 independent draws)")
    lines.append("period: Jan-Dec 2020 (full calendar year, no spin-up excluded), pooled (no monthly breakdown here)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned)")
    lines.append("")
    lines.append("=" * 110)
    lines.append("PER-ARM condensed metrics")
    lines.append("=" * 110)
    for arm in ARMS:
        lines.append(condensed_row(arm, metrics_by_arm[arm]))

    lines.append("")
    lines.append("=" * 110)
    lines.append("ARM A vs ARM-B 3-SEED SPREAD -- is A's value distinguishable from B seed-to-seed noise?")
    lines.append("=" * 110)
    lines.extend(seed_spread_block(metrics_by_arm, "pct_improved", "%improved (|O-A|<|O-F|)"))
    lines.append("")
    lines.extend(seed_spread_block(metrics_by_arm, "pct_reduction", "% reduction in mean|O-F| -> mean|O-A|"))
    lines.append("")
    lines.extend(seed_spread_block(metrics_by_arm, "gain_proxy_mean", "gain_proxy mean"))
    lines.append("")
    lines.extend(seed_spread_block(metrics_by_arm, "gain_proxy_median", "gain_proxy median"))
    lines.append("")
    lines.extend(seed_spread_block(metrics_by_arm, "rho_K_gainproxy", "Spearman rho(K, gain_proxy)"))
    lines.append("")
    lines.extend(seed_spread_block(metrics_by_arm, "pct_toward_obs", "%toward_obs (d_f*delta_h>0)"))

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["arm", "datetime", "tilenum", "tile_id", "lon", "lat", "obs", "fcst", "ana",
                  "obsvar", "fcstvar", "K", "d_f", "d_a", "delta_h", "gain_proxy",
                  "improved", "toward_obs"]
    df[event_cols].to_csv(EVENT_CSV, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return df, metrics_by_arm


if __name__ == "__main__":
    main()
