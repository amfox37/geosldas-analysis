#!/usr/bin/env python3
"""
Filter a CYGNSS L1 obs stream tier down to obs with coherency_ratio >=
COHERENCY_THRESHOLD, for the "<tier> + coherency>=0.5, xcorr/ycorr=1.25 (or
0.625 for dense), errstd=2.75, ungated binary" family of experiments
(intermediate-coh05 / dense075-coh05 / dense-coh05 -- a density spectrum, all
directly comparable, differing only in obs density).

coherency_ratio comes from the per-satellite/day QC-pass CSVs used by
build_cygl1_coherency_screening_experiment.py
(cygnss_l1_qc_pass_<date>_cyg<NN>.csv), joined by (sc_num, sample_id, ch_id)
-- NOT via that script's cached file_idx/obs_idx join log (which is keyed to
its own SRC_ROOT/date-range combination and does not line up with a
different tier's row indices). Obs whose join fails (no matching QC-CSV row)
are excluded, same convention as the screening experiment.

Reuses write_thinned_files() from thin_cygl1_nested_density_6mo.py unchanged
for the obs/support netCDF remapping -- only the file_idx/obs_idx/kept_mask
bookkeeping here is new. Only writes files for the dates passed in, so this
can be (and has been) rerun for successive non-overlapping date ranges to
extend an existing coh05 tree without touching earlier dates already written.

Consolidates six prior near-duplicate scripts that differed only in tier
(intermediate/dense075/dense) and BEG_DATE/END_DATE
(filter_cygl1_{dense,dense075,intermediate}_by_coherency.py +
filter_cygl1_{dense075,intermediate}_by_coherency_{jul_dec,2021}.py) --
merged 2026-09-08, now driven by CLI args instead of hardcoded values.

Usage: filter_cygl1_by_coherency.py --tier {intermediate,dense075,dense} \\
           --beg-date 20200101 --end-date 20200630 [--threshold 0.5]
"""
import argparse
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from thin_cygl1_nested_density_6mo import write_thinned_files  # noqa: E402

BASE = "/discover/nobackup/projects/land_da/cygl1_operator_test/"
QC_ROOT = "/gpfsm/dnb06/projects/p284/CYGNSS_operator/artifacts/out_images"

TIERS = {
    "intermediate": dict(
        src=BASE + "CYGNSS_L1_thinned_intermediate_6mo",
        dst=BASE + "CYGNSS_L1_thinned_intermediate_coh05",
        label="nested-superset-of-sparse, min_sep_deg=2.4, xcompact=ycompact=1.25deg, "
              "THEN coherency_ratio>={thresh} filter",
    ),
    "dense075": dict(
        src=BASE + "CYGNSS_L1_thinned_dense075_6mo",
        dst=BASE + "CYGNSS_L1_thinned_dense075_coh05",
        label="nested-superset-of-intermediate, min_sep_deg=0.75, xcompact=ycompact=1.25deg, "
              "THEN coherency_ratio>={thresh} filter",
    ),
    "dense": dict(
        src=BASE + "CYGNSS_L1",
        dst=BASE + "CYGNSS_L1_thinned_dense_coh05",
        label="full/unthinned stream, THEN coherency_ratio>={thresh} filter",
    ),
}


def qc_csv_path(date_str, sc_num):
    nn = f"{sc_num:02d}"
    return os.path.join(
        QC_ROOT,
        f"cygnss_qc_m36_window_counts_{date_str}_cyg{nn}",
        f"cygnss_l1_qc_pass_{date_str}_cyg{nn}.csv",
    )


def load_day_qc(date_str, sc_nums_needed):
    frames = []
    for sc in sc_nums_needed:
        path = qc_csv_path(date_str, int(sc))
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, usecols=["sample_id", "ch_id", "coherency_ratio"])
        df["sc_num"] = int(sc)
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["sc_num", "sample_id", "ch_id", "coherency_ratio"])
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", required=True, choices=list(TIERS))
    ap.add_argument("--beg-date", required=True, help="YYYYMMDD")
    ap.add_argument("--end-date", required=True, help="YYYYMMDD, inclusive")
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    tier = TIERS[args.tier]
    src_root, dst_root = tier["src"], tier["dst"]
    label = tier["label"].format(thresh=args.threshold)

    dates = pd.date_range(args.beg_date, args.end_date, freq="D")
    src_paths = []
    for dt in dates:
        y, m, ymd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%Y%m%d")
        pattern = os.path.join(src_root, f"Y{y}", f"M{m}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{ymd}_all_cyg.nc4")
        matches = glob.glob(pattern)
        if not matches:
            print(f"WARNING: no source file for {ymd} ({pattern})", file=sys.stderr)
            continue
        src_paths.append(matches[0])

    all_rows = []
    n_total = 0
    n_joined = 0
    n_pass = 0
    for fi, src_path in enumerate(src_paths):
        ymd = pd.Timestamp(dates[fi]).strftime("%Y%m%d")
        with nc.Dataset(src_path) as src:
            n_obs = src.dimensions["obs"].size
            sample_id = src.variables["sample_id"][:]
            ch_id = src.variables["ch_id"][:]
            sc_num = src.variables["sc_num"][:]

        obs_df = pd.DataFrame({
            "obs_idx": np.arange(n_obs),
            "sample_id": sample_id,
            "ch_id": ch_id,
            "sc_num": sc_num,
        })
        qc_df = load_day_qc(ymd, sorted(set(sc_num.tolist())))
        merged = obs_df.merge(qc_df, on=["sc_num", "sample_id", "ch_id"], how="left")

        n_total += n_obs
        n_joined += merged["coherency_ratio"].notna().sum()
        keep = merged["coherency_ratio"] >= args.threshold
        n_pass += keep.sum()

        for obs_idx in merged.loc[keep, "obs_idx"]:
            all_rows.append((fi, int(obs_idx)))

        if fi % 30 == 0:
            print(f"  {ymd}: {n_obs} {args.tier} obs -> {int(keep.sum())} coherency-pass", flush=True)

    df = pd.DataFrame(all_rows, columns=["file_idx", "obs_idx"])
    kept_mask = np.ones(len(df), dtype=bool)  # every row in df already passed the filter

    print(f"Total {args.tier} obs ({args.beg_date}-{args.end_date}): {n_total}")
    print(f"Joined to QC coherency_ratio: {n_joined} ({100.0*n_joined/n_total:.1f}%)")
    print(f"Passed coherency_ratio>={args.threshold}: {n_pass} ({100.0*n_pass/n_total:.1f}% of total, "
          f"{100.0*n_pass/n_joined:.1f}% of joined)")

    write_thinned_files(df, kept_mask, src_paths, dates, dst_root, label)
    print(f"Done. Output: {dst_root}")


if __name__ == "__main__":
    main()
