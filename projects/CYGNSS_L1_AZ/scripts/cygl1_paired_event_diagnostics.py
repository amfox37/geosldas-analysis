#!/usr/bin/env python3
"""
CYGNSS L1 paired thinning-density experiment: event-level single-observation
diagnostic (generalizes cygl1_singleobs_v2_event_diagnostics.py) run against
the combined 3-arm (sparse/intermediate/dense) CYGNSS L1 event table -- NO
isolation filter this time, since the whole point is to see whether the
previously-established clean single-obs analysis behavior (validated on the
genuinely isolated singleobs_v2 set) holds as density increases.

For each retained obs:
    d_f        = O - F
    d_a        = O - A
    delta_h    = A - F
    gain_proxy = delta_h / d_f   (only for well-conditioned |d_f| >= 0.1 dB)

Six checks (same as cygl1_singleobs_v2_event_diagnostics.py):
  1/2. mean|O-A| vs mean|O-F| (+ RMSE), %improved
  3. d_f * delta_h > 0 (toward-obs rate)
  4. gain_proxy distribution: %negative / %in[0,1] / %overshoot
  5. monotonicity of gain_proxy vs naive K = fcstvar/(fcstvar+obsvar)
  6. positive vs negative innovation (d_f) asymmetry

Computed separately per arm, BOTH full-period (May-Oct pooled) AND monthly,
AND stratified by local-interaction-count bin (0, 1, 2, 3+) within each arm
-- the central new analysis: does gain_proxy/improvement degrade as local
obs density (number of other kept obs whose Gaspari-Cohn ellipse overlaps
this obs's, same window) increases?

Inputs:
  cygl1_paired_gain_consistency_ofa.csv.gz   (extract_cygl1_paired_gain_consistency.py)
  cygl1_paired_local_interaction_lookup.csv  (cygl1_paired_local_interaction_join.py)
    -- left-joined on (arm, datetime, tilenum); unmatched rows are logged and
       kept in the full-period/monthly tables (local_interaction_count only
       needed for the stratified tables) rather than silently dropped.

Outputs:
  cygl1_paired_event_diagnostics.csv    (event-level, all 3 arms, +local_interaction_count)
  cygl1_paired_event_diagnostics_summary.txt  (all tables: full-period verbose,
    monthly condensed, count-bin-stratified condensed, per arm)
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_paired_gain_consistency_ofa.csv.gz")
LOOKUP_CSV = os.path.join(OUT_DIR, "cygl1_paired_local_interaction_lookup.csv")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_paired_event_diagnostics.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_paired_event_diagnostics_summary.txt")

DFLOOR = 0.1  # dB; below this, gain_proxy is ill-conditioned (near-zero denominator)
ARMS = ["sparse", "intermediate", "dense"]


def bin_count(c):
    if pd.isna(c):
        return np.nan
    c = int(c)
    if c >= 3:
        return "3+"
    return str(c)


def verbose_block(df, label):
    """Full six-check verbose block, same style as cygl1_singleobs_v2_event_diagnostics.py."""
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
    lines.append(f"   RMSE(O-F) = {rmse_df:.4f} dB, RMSE(O-A) = {rmse_da:.4f} dB "
                  f"({100*(1-rmse_da/rmse_df):.1f}% reduction)" if rmse_df > 0 else "")
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
    """One-line condensed metrics for monthly / stratified tables."""
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

    lookup = pd.read_csv(LOOKUP_CSV)
    lookup["datetime"] = pd.to_datetime(lookup["datetime"])
    lookup_small = lookup[["arm", "datetime", "tilenum", "local_interaction_count"]].drop_duplicates(
        subset=["arm", "datetime", "tilenum"])
    print(f"Loaded {len(lookup)} lookup rows (before dedup: local-interaction join key); "
          f"{len(lookup_small)} unique (arm,datetime,tilenum) keys")

    before = len(df)
    df = df.merge(lookup_small, on=["arm", "datetime", "tilenum"], how="left")
    assert len(df) == before, "merge changed row count -- duplicate join keys in lookup"
    n_unmatched = int(df["local_interaction_count"].isna().sum())
    print(f"Joined local_interaction_count onto {before} events: {n_unmatched} unmatched "
          f"({100*n_unmatched/before:.2f}%)")
    if n_unmatched > 0:
        print("  Unmatched breakdown by arm:")
        print(df[df["local_interaction_count"].isna()].groupby("arm").size())

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
    df["count_bin"] = df["local_interaction_count"].apply(bin_count)

    df = df.sort_values(["arm", "datetime", "tilenum"]).reset_index(drop=True)

    lines = []
    lines.append("CYGNSS L1 paired thinning-density experiment: event-level diagnostics")
    lines.append("arms: sparse (min_sep=5.0deg, isolated), intermediate (min_sep=2.40deg), dense (full unthinned stream)")
    lines.append("period: May-Oct 2020 (evaluation window)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned)")
    lines.append(f"local_interaction_count: joined from cygl1_paired_local_interaction_join.py; "
                  f"{n_unmatched}/{before} events unmatched ({100*n_unmatched/before:.2f}%)")
    lines.append("")

    lines.append("=" * 100)
    lines.append("FULL-PERIOD (May-Oct pooled), per arm -- verbose six-check blocks")
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

    lines.append("=" * 100)
    lines.append("STRATIFIED BY LOCAL-INTERACTION-COUNT BIN (0/1/2/3+), per arm -- condensed")
    lines.append("(sparse is all-0 by construction; shown for completeness)")
    lines.append("=" * 100)
    for arm in ARMS:
        lines.append(f"--- arm={arm} ---")
        sub = df[(df["arm"] == arm) & df["count_bin"].notna()]
        n_dropped = int((df["arm"] == arm).sum()) - len(sub)
        if n_dropped:
            lines.append(f"   ({n_dropped} events dropped from this table: unmatched local_interaction_count)")
        for cbin in ["0", "1", "2", "3+"]:
            g = sub[sub["count_bin"] == cbin]
            lines.append(condensed_row(g, f"count={cbin}"))
        lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["arm", "datetime", "tilenum", "tile_id", "lon", "lat", "obs", "fcst", "ana",
                  "obsvar", "fcstvar", "K", "d_f", "d_a", "delta_h", "gain_proxy",
                  "improved", "toward_obs", "local_interaction_count", "month", "count_bin"]
    df[event_cols].to_csv(EVENT_CSV, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return df


if __name__ == "__main__":
    main()
