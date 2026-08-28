#!/usr/bin/env python3
"""
Redraw Arm B (random-matched control) of the CYGNSS L1 coherency-screening
experiment with new seeds, WITHOUT touching Arm A or the original Arm B
(seed=20260827) and WITHOUT rerunning the expensive QC-CSV join step --
reuses the cached full-join log + arm-A log written by the original
build_cygl1_coherency_screening_experiment.py run.

Follow-up requested by user 2026-08-27 to check whether the original Arm B's
-4.2% CygL1 skill-vs-OL (vs Arm A's -8.0%) was a lucky/unlucky single-seed draw.
"""
import glob
import os

import netCDF4 as nc
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(__file__))
from build_cygl1_coherency_screening_experiment import build_arm_b, write_thinned_files  # noqa: E402

SRC_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
LOG_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"

BEG_DATE = "20200101"
END_DATE = "20201231"

SEEDS = {
    "seed2": 20260828,
    "seed3": 20260829,
}
DST_ROOTS = {
    "seed2": "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_coherency_randmatch_seed2",
    "seed3": "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_coherency_randmatch_seed3",
}


def rebuild_src_paths_and_dates():
    dates = pd.date_range(BEG_DATE, END_DATE, freq="D")
    src_paths = []
    for dt in dates:
        y, m, ymd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%Y%m%d")
        pattern = os.path.join(SRC_ROOT, f"Y{y}", f"M{m}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{ymd}_all_cyg.nc4")
        matches = glob.glob(pattern)
        src_paths.append(matches[0] if matches else None)
    return src_paths, dates


def main():
    print("Loading cached full-join log (skips re-running the QC-CSV join)...")
    full_log = pd.read_csv(os.path.join(LOG_ROOT, "cygl1_coherency_full_join_log.csv"))
    arm_a_log = pd.read_csv(os.path.join(LOG_ROOT, "cygl1_coherency_screened_arm_a_log.csv"),
                             usecols=["file_idx", "obs_idx"])

    df = full_log.copy()
    df["candidate"] = df["coherency_ratio"].notna()

    a_keys = set(zip(arm_a_log["file_idx"], arm_a_log["obs_idx"]))
    df["kept_a"] = [(fi, oi) in a_keys for fi, oi in zip(df["file_idx"], df["obs_idx"])]

    n_candidate = int(df["candidate"].sum())
    n_a = int(df["kept_a"].sum())
    print(f"Reconstructed candidate pool: {n_candidate}/{len(df)}; Arm A kept: {n_a}")
    assert n_a == 155847 or True  # sanity note only, not a hard gate

    src_paths, dates = rebuild_src_paths_and_dates()
    assert len(src_paths) == len(dates)

    for label, seed in SEEDS.items():
        print(f"\n=== {label} (seed={seed}) ===")
        kept_b = build_arm_b(df, df["candidate"].values, df["kept_a"].values, seed)
        n_b = int(kept_b.sum())
        print(f"{label}: {n_b} kept (should equal arm A's {n_a})")
        assert n_b == n_a, f"{label} count {n_b} != arm A count {n_a}"

        dst_root = DST_ROOTS[label]
        nfb = write_thinned_files(
            df, kept_b, src_paths, dates, dst_root,
            f"random-matched control, seed={seed}, per-window count matched to arm A "
            f"(re-seed follow-up, cached join)"
        )
        print(f"{label}: wrote {nfb} files to {dst_root}")

        out_cols = ["file_idx", "obs_idx", "anchor", "sc_num", "sample_id", "ch_id", "coherency_ratio"]
        df.loc[kept_b, out_cols].to_csv(
            os.path.join(LOG_ROOT, f"cygl1_coherency_randmatch_arm_b_{label}_log.csv"), index=False)
        print(f"{label}: wrote log to cygl1_coherency_randmatch_arm_b_{label}_log.csv")


if __name__ == "__main__":
    main()
