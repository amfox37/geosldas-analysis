#!/usr/bin/env python3
"""
CYGNSS L1 paired thinning-density experiment: withheld-monitor matched
comparison. Answers whether DA-intermediate/DA-dense reproduce the
previously-established (unrelated, earlier R-sweep experiment) finding that
assimilating CYGNSS L1 with progressively lower obs-error R monotonically
degraded all 8 monitor-only SMOS/SMAP Tb species' O-F std (worst at
quarter-R; ASCAT/CYGNSS_SM_6hr degraded too but mildly), or whether it stays
clean like the isolated single-obs case.

For the 12 monitor-only species (all except CYGNSS_L1_DDM3X5_CROP_SCALAR),
join each DA arm's (species_name, tile_id, datetime) rows to OL's rows on
the same key (one shared OL run suffices -- monitor-only obs can't feed back
into OL state, so OL's fcst is a valid, arm-independent baseline forecast
trajectory to compare each DA arm's fcst against).

Per event:
    OF_DA = obs - fcst   (DA arm's fcst)
    OF_OL = obs - fcst   (OL's fcst, same species/tile/time)
    diff_abs = |OF_DA| - |OF_OL|     (negative = DA improved vs OL, positive = degraded)
    diff_sq  = OF_DA^2 - OF_OL^2

Report per species per arm: mean/median of both diff metrics, %improved
(|OF_DA| < |OF_OL|), and monthly evolution (May..Oct).

Input: cygl1_paired_withheld_monitor_ofa.csv.gz (extract_cygl1_paired_withheld_monitor.py)
Outputs: cygl1_paired_withheld_monitor_diagnostics.csv (event-level joined table)
         cygl1_paired_withheld_monitor_diagnostics_summary.txt
"""
import os
import gzip
import numpy as np
import pandas as pd

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_paired_withheld_monitor_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_paired_withheld_monitor_diagnostics.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_paired_withheld_monitor_diagnostics_summary.txt")

DA_ARMS = ["sparse", "intermediate", "dense"]

SPECIES_ORDER = [
    "SMOS_fit_Tbh_A", "SMOS_fit_Tbh_D", "SMOS_fit_Tbv_A", "SMOS_fit_Tbv_D",
    "SMAP_L1C_Tbh_A", "SMAP_L1C_Tbh_D", "SMAP_L1C_Tbv_A", "SMAP_L1C_Tbv_D",
    "ASCAT_HSAF_META_SM", "ASCAT_HSAF_METB_SM", "ASCAT_HSAF_METC_SM", "CYGNSS_SM_6hr",
]


def main():
    print(f"Reading {IN_CSV} ...", flush=True)
    with gzip.open(IN_CSV, "rt") as fi:
        df = pd.read_csv(fi, comment="#")
    df["datetime"] = pd.to_datetime(df["datetime"])
    print(f"Loaded {len(df)} rows. Row counts by arm:")
    print(df.groupby("arm").size())

    ol = df[df["arm"] == "OL"][["datetime", "tile_id", "species_name", "obs", "fcst"]].rename(
        columns={"obs": "obs_OL", "fcst": "fcst_OL"})
    print(f"\nOL baseline rows: {len(ol)}")
    ol_dupe = ol.duplicated(subset=["datetime", "tile_id", "species_name"]).sum()
    if ol_dupe:
        print(f"WARNING: {ol_dupe} duplicate OL (datetime,tile_id,species_name) keys -- "
              f"unexpected, should be unique per obs")

    joined_parts = []
    for arm in DA_ARMS:
        da = df[df["arm"] == arm][["datetime", "tile_id", "species_name", "obs", "fcst"]].rename(
            columns={"obs": "obs_DA", "fcst": "fcst_DA"})
        n_da = len(da)
        m = da.merge(ol, on=["datetime", "tile_id", "species_name"], how="left")
        n_unmatched = int(m["fcst_OL"].isna().sum())
        print(f"[{arm}] DA rows={n_da}, matched to OL={n_da - n_unmatched}, "
              f"unmatched={n_unmatched} ({100*n_unmatched/n_da:.2f}%)")
        m["arm"] = arm
        joined_parts.append(m)

    j = pd.concat(joined_parts, ignore_index=True)
    n_before_drop = len(j)
    j = j.dropna(subset=["fcst_OL"]).copy()
    print(f"\nTotal joined rows: {n_before_drop}, after dropping unmatched: {len(j)} "
          f"({100*(n_before_drop-len(j))/n_before_drop:.2f}% dropped, logged above per-arm)")

    # sanity: obs should be essentially identical between DA arm's copy and OL's copy
    # (same physical observation, monitor-only, not perturbed by assimilation)
    obs_diff = (j["obs_DA"] - j["obs_OL"]).abs()
    n_obs_mismatch = int((obs_diff > 1e-6).sum())
    print(f"obs_DA vs obs_OL mismatch (>1e-6): {n_obs_mismatch}/{len(j)} "
          f"({100*n_obs_mismatch/len(j):.3f}%) -- should be ~0 (monitor obs unperturbed by assim)")

    j["OF_DA"] = j["obs_DA"] - j["fcst_DA"]
    j["OF_OL"] = j["obs_OL"] - j["fcst_OL"]
    j["diff_abs"] = j["OF_DA"].abs() - j["OF_OL"].abs()
    j["diff_sq"] = j["OF_DA"] ** 2 - j["OF_OL"] ** 2
    j["improved"] = j["OF_DA"].abs() < j["OF_OL"].abs()
    j["month"] = j["datetime"].dt.strftime("%Y-%m")

    j = j.sort_values(["arm", "species_name", "datetime", "tile_id"]).reset_index(drop=True)

    lines = []
    lines.append("CYGNSS L1 paired thinning-density experiment: withheld-monitor matched comparison")
    lines.append("12 monitor-only species (assim_flag==0 in all 4 arms), OL baseline = OLv8_M36_AZ_paired_monitor")
    lines.append("period: May-Oct 2020")
    lines.append(f"obs_DA vs obs_OL mismatch: {n_obs_mismatch}/{len(j)} ({100*n_obs_mismatch/len(j):.3f}%)")
    lines.append("diff_abs = |OF_DA| - |OF_OL|  (negative = DA improved vs OL, positive = DA degraded)")
    lines.append("diff_sq  = OF_DA^2 - OF_OL^2")
    lines.append("")

    lines.append("=" * 120)
    lines.append("FULL-PERIOD (May-Oct pooled), per species x arm")
    lines.append("=" * 120)
    hdr = (f"{'species':22s} {'arm':14s} {'N':>7s} {'mean_diff_abs':>14s} {'median_diff_abs':>16s} "
           f"{'mean_diff_sq':>13s} {'%improved':>10s} {'flag':>12s}")
    lines.append(hdr)
    full_period_rows = []
    for sp in SPECIES_ORDER:
        for arm in DA_ARMS:
            sub = j[(j["species_name"] == sp) & (j["arm"] == arm)]
            n = len(sub)
            if n == 0:
                lines.append(f"{sp:22s} {arm:14s} {'0':>7s}  (no matched events)")
                continue
            mean_da = sub["diff_abs"].mean()
            med_da = sub["diff_abs"].median()
            mean_dsq = sub["diff_sq"].mean()
            pct_imp = 100 * sub["improved"].mean()
            if mean_da > 0.01:
                flag = "DEGRADED"
            elif mean_da < -0.01:
                flag = "improved"
            else:
                flag = "neutral"
            lines.append(f"{sp:22s} {arm:14s} {n:7d} {mean_da:14.5f} {med_da:16.5f} "
                          f"{mean_dsq:13.5f} {pct_imp:9.2f}% {flag:>12s}")
            full_period_rows.append((sp, arm, n, mean_da, med_da, mean_dsq, pct_imp, flag))
        lines.append("")

    full_period_df = pd.DataFrame(full_period_rows, columns=[
        "species_name", "arm", "N", "mean_diff_abs", "median_diff_abs", "mean_diff_sq", "pct_improved", "flag"])

    lines.append("=" * 120)
    lines.append("MONTHLY EVOLUTION of mean_diff_abs (|OF_DA|-|OF_OL|), per species x arm")
    lines.append("(positive = DA degraded vs OL that month; watch for monotonic growth = R-sweep-style dose-response)")
    lines.append("=" * 120)
    months = ["2020-05", "2020-06", "2020-07", "2020-08", "2020-09", "2020-10"]
    hdr2 = f"{'species':22s} {'arm':14s} " + " ".join(f"{m:>10s}" for m in months)
    lines.append(hdr2)
    for sp in SPECIES_ORDER:
        for arm in DA_ARMS:
            sub = j[(j["species_name"] == sp) & (j["arm"] == arm)]
            vals = []
            for m in months:
                gm = sub[sub["month"] == m]
                vals.append(f"{gm['diff_abs'].mean():10.4f}" if len(gm) else f"{'n/a':>10s}")
            lines.append(f"{sp:22s} {arm:14s} " + " ".join(vals))
        lines.append("")

    lines.append("=" * 120)
    lines.append("SUMMARY: species/arms showing DEGRADATION (mean_diff_abs > 0.01) vs neutral vs improvement")
    lines.append("=" * 120)
    for flag_name in ["DEGRADED", "neutral", "improved"]:
        rows = full_period_df[full_period_df["flag"] == flag_name]
        lines.append(f"{flag_name}: {len(rows)}/{len(full_period_df)}")
        for _, r in rows.iterrows():
            lines.append(f"   {r['species_name']:22s} {r['arm']:14s}  mean_diff_abs={r['mean_diff_abs']:.5f}  "
                          f"%improved={r['pct_improved']:.1f}%")
    lines.append("")

    # Tb-species-only monotonicity-in-density check (dose-response signature from R-sweep finding)
    tb_species = [s for s in SPECIES_ORDER if "Tb" in s]
    lines.append("=" * 120)
    lines.append("DOSE-RESPONSE CHECK: mean_diff_abs vs arm density (sparse -> intermediate -> dense), Tb species only")
    lines.append("(R-sweep finding was monotonic degradation as R decreased / assim strength increased;")
    lines.append(" here density increases analogously -- check whether degradation grows sparse->intermediate->dense)")
    lines.append("=" * 120)
    for sp in tb_species:
        vals = []
        for arm in DA_ARMS:
            row = full_period_df[(full_period_df["species_name"] == sp) & (full_period_df["arm"] == arm)]
            v = row["mean_diff_abs"].values[0] if len(row) else np.nan
            vals.append(v)
        monotonic = all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1) if not np.isnan(vals[i]) and not np.isnan(vals[i+1]))
        lines.append(f"   {sp:22s} sparse={vals[0]:.5f}  intermediate={vals[1]:.5f}  dense={vals[2]:.5f}  "
                      f"monotonic_increase={monotonic}")

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["arm", "datetime", "tile_id", "species_name", "obs_DA", "fcst_DA", "obs_OL", "fcst_OL",
                  "OF_DA", "OF_OL", "diff_abs", "diff_sq", "improved", "month"]
    with gzip.open(EVENT_CSV, "wt") as fo:
        j[event_cols].to_csv(fo, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return j, full_period_df


if __name__ == "__main__":
    main()
