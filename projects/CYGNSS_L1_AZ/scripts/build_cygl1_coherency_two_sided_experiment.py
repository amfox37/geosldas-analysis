#!/usr/bin/env python3
"""
Build Arm C (two-sided coherency screen) of the CYGNSS L1 coherency-screening
experiment, per cygl1_coherency_screening_experiment_spec.md sec 3/7's
explicitly-flagged follow-up: "is the top end independently bad, or was
cutting it just discarding good data".

Arm C: keep obs with COHERENCY_LO <= coherency_ratio <= COHERENCY_HI, i.e. the
middle band -- same lower bound as Arm A (0.523, UNCHANGED, so Arm A and Arm C
are directly comparable / differ only in the added top cut), plus a new upper
cut at P80 of the same candidate population Arm A/B were built from.

Reuses the cached full-join log written by the original
build_cygl1_coherency_screening_experiment.py run (full year 2020, AZ domain,
all satellites, join-succeeded candidate pool) -- does NOT rerun the expensive
QC-CSV join. Same reuse pattern as build_cygl1_coherency_randmatch_reseed.py.

No new random-matched control is built for Arm C: the existing Arm B (3
seeds) already isolates "which obs" from "how many obs" at Arm A's count; Arm
C's job is to show whether trimming the top quintile helps or hurts relative
to Arm A's own already-banked numbers (spec sec 7: "record the one-sided
statistics from the same join alongside it so the trade is visible").
"""
import glob
import os

import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(__file__))
from build_cygl1_coherency_screening_experiment import write_thinned_files  # noqa: E402

SRC_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
LOG_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
DST_C = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_coherency_twosided"

BEG_DATE = "20200101"
END_DATE = "20201231"

COHERENCY_LO = 0.523  # unchanged from Arm A -- keeps the trade-off apples-to-apples
# COHERENCY_HI computed below as P80 of the candidate pool (same population Arm A/B used)


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
    df = pd.read_csv(os.path.join(LOG_ROOT, "cygl1_coherency_full_join_log.csv"))
    df["candidate"] = df["coherency_ratio"].notna()

    candidate = df["candidate"].values
    cr = df["coherency_ratio"]
    coherency_hi = float(cr[df["candidate"]].quantile(0.80))
    print(f"COHERENCY_LO = {COHERENCY_LO} (unchanged from Arm A)")
    print(f"COHERENCY_HI = {coherency_hi:.6f} (P80 of candidate pool, computed fresh here)")

    kept_c = candidate & (cr >= COHERENCY_LO).values & (cr <= coherency_hi).values
    n_candidate = int(candidate.sum())
    n_c = int(kept_c.sum())
    n_a_expected = int((candidate & (cr >= COHERENCY_LO).values).sum())
    print(f"Candidate pool: {n_candidate}")
    print(f"Arm A count (for reference, >= {COHERENCY_LO} only): {n_a_expected} "
          f"({100.0*n_a_expected/n_candidate:.1f}%)")
    print(f"Arm C count ({COHERENCY_LO} <= x <= {coherency_hi:.4f}): {n_c} "
          f"({100.0*n_c/n_candidate:.1f}% of candidate pool)")

    src_paths, dates = rebuild_src_paths_and_dates()
    assert len(src_paths) == len(dates)

    nfc = write_thinned_files(
        df, kept_c, src_paths, dates, DST_C,
        f"two-sided coherency screen, {COHERENCY_LO}<=coherency_ratio<={coherency_hi:.4f} "
        f"(P20-P80 of full-year-2020 candidate pool, lower bound unchanged from arm A)"
    )
    print(f"Arm C: wrote {nfc} files to {DST_C}")

    out_cols = ["file_idx", "obs_idx", "anchor", "sc_num", "sample_id", "ch_id", "coherency_ratio"]
    df.loc[kept_c, out_cols].to_csv(
        os.path.join(LOG_ROOT, "cygl1_coherency_twosided_arm_c_log.csv"), index=False)
    print(f"Wrote log to cygl1_coherency_twosided_arm_c_log.csv")
    print(f"\nCOHERENCY_LO={COHERENCY_LO}, COHERENCY_HI={coherency_hi:.6f}, "
          f"date range={BEG_DATE}-{END_DATE}")


if __name__ == "__main__":
    main()

# ====================== EOF =========================================================
