#!/usr/bin/env python3
"""
Build the two coherency-screening-experiment CYGNSS L1 obs arms, per
cygl1_operator_test/cygl1_coherency_screening_experiment_spec.md:

  Arm A (coherency-screened): keep obs with coherency_ratio >= COHERENCY_THRESHOLD
  Arm B (random-matched control): same kept-count per assimilation window as arm A,
    drawn at random (fixed seed) from the same join-succeeded candidate pool.

coherency_ratio comes from the CYGNSS_operator repo's per-satellite/day QC-pass
CSVs (cygnss_l1_qc_pass_<date>_cyg<NN>.csv), joined onto the merged obs GEOSldas
actually ingests (CYGNSS_L1/Y*/M*/cygnss_l1_ddm3x5_crop_scalar_m36_<date>_all_cyg.nc4)
by (sc_num->cygNN, sample_id, ch_id). Obs whose join fails (no matching QC-CSV
row -- missing day/satellite QC file, or no row for that sample_id/ch_id) are
EXCLUDED from the candidate pool for BOTH arms, not just missing from one --
see spec sec 2.

Window/anchor computation (true centered 3-hr assimilation window, handles
day-boundary-straddling correctly) reuses the same method already validated in
thin_cygl1_nested_density_6mo.py -- do not reinvent.
"""
import argparse
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

SRC_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
DST_A = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_coherency_screened"
DST_B = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_coherency_randmatch"
QC_ROOT = "/gpfsm/dnb06/projects/p284/CYGNSS_operator/artifacts/out_images"

BEG_DATE = "20200101"
END_DATE = "20201231"  # inclusive -- full calendar year 2020, no spin-up excluded (user direction)

DTSTEP_ASSIM_SEC = 10800  # LANDASSIM_DT
LANDASSIM_T0_SEC = 0

COHERENCY_THRESHOLD = 0.523  # ~P20 of coherency_ratio; drop bottom quintile, keep top -- see spec sec 3
RANDOM_SEED = 20260827        # arbitrary, fixed; single seed per user decision (spec sec 4)

LOG_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"


def qc_csv_path(date_str, sc_num):
    nn = f"{sc_num:02d}"
    return os.path.join(
        QC_ROOT,
        f"cygnss_qc_m36_window_counts_{date_str}_cyg{nn}",
        f"cygnss_l1_qc_pass_{date_str}_cyg{nn}.csv",
    )


def load_day_qc(date_str, sc_nums_needed):
    """Load and concat this date's QC-pass CSVs for just the satellites present
    in that day's obs file. Returns a DataFrame with columns [sc_num, sample_id,
    ch_id, coherency_ratio], or an empty frame if none found."""
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


def load_all_obs():
    dates = pd.date_range(BEG_DATE, END_DATE, freq="D")
    src_paths = []
    for dt in dates:
        y, m, ymd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%Y%m%d")
        pattern = os.path.join(SRC_ROOT, f"Y{y}", f"M{m}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{ymd}_all_cyg.nc4")
        matches = glob.glob(pattern)
        if not matches:
            print(f"WARNING: no source file for {ymd} ({pattern})", file=sys.stderr)
            src_paths.append(None)
            continue
        src_paths.append(matches[0])

    all_frames = []
    n_total_good = 0
    n_join_ok = 0
    for fi, (dt, src_path) in enumerate(zip(dates, src_paths)):
        if src_path is None:
            continue
        ymd = dt.strftime("%Y%m%d")
        with nc.Dataset(src_path) as src:
            status = src.variables["status"][:]
            lon = src.variables["sp_lon"][:]
            lat = src.variables["sp_lat"][:]
            year = src.variables["year"][:]
            day = src.variables["day"][:]
            secs = src.variables["ddm_timestamp_utc_sec"][:]
            sample_id = src.variables["sample_id"][:]
            ch_id = src.variables["ch_id"][:]
            sc_num = src.variables["sc_num"][:]

        good = np.where(status == 1)[0]
        n_total_good += len(good)
        if len(good) == 0:
            continue

        base = pd.to_datetime(year[good].astype(str), format="%Y").values.astype("datetime64[s]")
        abs_ts = (
            base.astype("int64")
            + (day[good].astype(np.int64) - 1) * 86400
            + np.round(secs[good]).astype(np.int64)
        )
        anchor = np.round((abs_ts - LANDASSIM_T0_SEC) / DTSTEP_ASSIM_SEC).astype(np.int64) \
            * DTSTEP_ASSIM_SEC + LANDASSIM_T0_SEC

        day_df = pd.DataFrame({
            "file_idx": fi,
            "obs_idx": good,
            "lon": lon[good].astype(float),
            "lat": lat[good].astype(float),
            "anchor": anchor,
            "sample_id": sample_id[good].astype(int),
            "ch_id": ch_id[good].astype(int),
            "sc_num": sc_num[good].astype(int),
        })

        qc_df = load_day_qc(ymd, np.unique(day_df["sc_num"].values))
        merged = day_df.merge(qc_df, on=["sc_num", "sample_id", "ch_id"], how="left")
        # guard against a QC CSV containing duplicate (sc_num,sample_id,ch_id) rows,
        # which would silently fan out the merge -- assert 1:1 cardinality instead.
        assert len(merged) == len(day_df), (
            f"{ymd}: join fanned out ({len(day_df)} obs -> {len(merged)} rows) -- "
            "duplicate (sc_num,sample_id,ch_id) in a QC CSV?"
        )
        n_join_ok += merged["coherency_ratio"].notna().sum()
        all_frames.append(merged)

        if fi % 30 == 0:
            print(f"  {ymd}: {len(day_df)} good obs, "
                  f"{merged['coherency_ratio'].notna().sum()} joined", flush=True)

    df = pd.concat(all_frames, ignore_index=True)
    print(f"\nLoaded {len(df)} good obs across {len(src_paths)} days ({BEG_DATE}-{END_DATE})")
    print(f"Join match rate: {n_join_ok}/{n_total_good} = {100.0*n_join_ok/n_total_good:.2f}%")
    return df, src_paths, dates


def build_arm_a(df):
    """Coherency screen: keep obs with coherency_ratio >= COHERENCY_THRESHOLD,
    restricted to the join-succeeded candidate pool."""
    candidate = df["coherency_ratio"].notna()
    kept = candidate & (df["coherency_ratio"] >= COHERENCY_THRESHOLD)
    return candidate, kept


def build_arm_b(df, candidate_mask, kept_mask_a, seed):
    """Random-matched control: for each window (anchor), draw exactly as many
    obs at random (no replacement) from that window's candidate pool as arm A
    kept in that same window."""
    rng = np.random.default_rng(seed)
    kept_b = np.zeros(len(df), dtype=bool)

    df_a = pd.DataFrame({"anchor": df["anchor"], "candidate": candidate_mask, "kept_a": kept_mask_a})
    for anchor, grp in df_a[df_a["candidate"]].groupby("anchor"):
        k = int(grp["kept_a"].sum())
        if k == 0:
            continue
        pool_idx = grp.index.values
        assert k <= len(pool_idx), f"anchor {anchor}: arm A kept {k} > candidate pool {len(pool_idx)}"
        chosen = rng.choice(pool_idx, size=k, replace=False)
        kept_b[chosen] = True
    return kept_b


def write_thinned_files(df, kept_mask, src_paths, dates, dst_root, label):
    n_files_written = 0
    for fi, src_path in enumerate(src_paths):
        if src_path is None:
            continue
        kept_obs_idx = sorted(df.loc[(df.file_idx == fi) & kept_mask, "obs_idx"].tolist())
        with nc.Dataset(src_path) as src:
            tile_start = src.variables["tile_start"][:]
            tile_count = src.variables["tile_count"][:]

            kept = np.array(kept_obs_idx, dtype=int)
            n_kept = len(kept)

            new_tile_start = np.zeros(n_kept, dtype=tile_start.dtype)
            new_tile_count = tile_count[kept].copy() if n_kept else tile_count[kept]
            support_take = []
            cursor = 0
            for k, i in enumerate(kept):
                new_tile_start[k] = cursor
                s0, s1 = int(tile_start[i]), int(tile_start[i]) + int(tile_count[i])
                support_take.append(np.arange(s0, s1))
                cursor += (s1 - s0)
            support_take = np.concatenate(support_take) if support_take else np.array([], dtype=int)

            dst_path = os.path.join(
                dst_root,
                pd.Timestamp(dates[fi]).strftime("Y%Y"),
                pd.Timestamp(dates[fi]).strftime("M%m"),
                os.path.basename(src_path),
            )
            os.makedirs(os.path.dirname(dst_path), exist_ok=True)
            with nc.Dataset(dst_path, "w", format=src.file_format) as dst:
                dst.createDimension("obs", n_kept)
                dst.createDimension("support", len(support_take))

                for name, var in src.variables.items():
                    if var.dimensions == ("obs",):
                        if name == "tile_start":
                            data = new_tile_start
                        elif name == "tile_count":
                            data = new_tile_count
                        else:
                            data = var[:][kept]
                        newvar = dst.createVariable(name, var.dtype, ("obs",))
                        newvar[:] = data
                    elif var.dimensions == ("support",):
                        data = var[:][support_take]
                        newvar = dst.createVariable(name, var.dtype, ("support",))
                        newvar[:] = data
                    else:
                        raise ValueError(f"unexpected dims for {name}: {var.dimensions}")
                    for attname in var.ncattrs():
                        newvar.setncattr(attname, var.getncattr(attname))

                for attname in src.ncattrs():
                    dst.setncattr(attname, src.getncattr(attname))
                dst.setncattr("thinning_applied", label)
                dst.setncattr("thinning_script",
                               "geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/"
                               "build_cygl1_coherency_screening_experiment.py")
                dst.setncattr("thinning_source_file", src_path)
                dst.setncattr("thinning_n_obs_original", int((df.file_idx == fi).sum()))
                dst.setncattr("thinning_n_obs_kept", n_kept)
        n_files_written += 1
        if fi % 30 == 0:
            print(f"  {pd.Timestamp(dates[fi]).strftime('%Y%m%d')}: "
                  f"{(df.file_idx==fi).sum()} -> {n_kept} kept  ({dst_root})", flush=True)
    return n_files_written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beg-date", type=str, default=None, help="override BEG_DATE (YYYYMMDD)")
    ap.add_argument("--end-date", type=str, default=None, help="override END_DATE (YYYYMMDD, inclusive)")
    ap.add_argument("--dry-run", action="store_true",
                    help="compute + log arm A/B membership only, do not write thinned nc4 files")
    args = ap.parse_args()

    global BEG_DATE, END_DATE
    if args.beg_date:
        BEG_DATE = args.beg_date
    if args.end_date:
        END_DATE = args.end_date

    df, src_paths, dates = load_all_obs()

    candidate_mask, kept_a = build_arm_a(df)
    n_candidate = int(candidate_mask.sum())
    n_a = int(kept_a.sum())
    print(f"\nCandidate pool (join succeeded): {n_candidate}/{len(df)}")
    print(f"Arm A (coherency_ratio >= {COHERENCY_THRESHOLD}): {n_a} kept "
          f"({100.0*n_a/n_candidate:.1f}% of candidate pool)")

    kept_b = build_arm_b(df, candidate_mask, kept_a, RANDOM_SEED)
    n_b = int(kept_b.sum())
    print(f"Arm B (random-matched, seed={RANDOM_SEED}): {n_b} kept "
          f"(should equal arm A's {n_a})")
    assert n_b == n_a, f"arm B count {n_b} != arm A count {n_a}"

    if args.dry_run:
        print("\n--dry-run: skipping file writes.")
        return

    print("\n--- Writing arm A files ---")
    nfa = write_thinned_files(df, kept_a.values, src_paths, dates, DST_A,
                               f"coherency-screened, coherency_ratio>={COHERENCY_THRESHOLD}, "
                               f"threshold~P20 of full-year-2020 distribution")
    print(f"Arm A: wrote {nfa} files to {DST_A}")

    print("\n--- Writing arm B files ---")
    nfb = write_thinned_files(df, kept_b, src_paths, dates, DST_B,
                               f"random-matched control, seed={RANDOM_SEED}, "
                               f"per-window count matched to arm A")
    print(f"Arm B: wrote {nfb} files to {DST_B}")

    # logs
    os.makedirs(LOG_ROOT, exist_ok=True)
    log_cols = ["file_idx", "obs_idx", "lon", "lat", "anchor", "sc_num", "sample_id",
                "ch_id", "coherency_ratio"]
    df[candidate_mask & kept_a][log_cols].to_csv(
        os.path.join(LOG_ROOT, "cygl1_coherency_screened_arm_a_log.csv"), index=False)
    df.loc[kept_b, log_cols].to_csv(
        os.path.join(LOG_ROOT, "cygl1_coherency_randmatch_arm_b_log.csv"), index=False)
    df[["file_idx", "obs_idx", "anchor", "sc_num", "sample_id", "ch_id",
        "coherency_ratio"]].to_csv(
        os.path.join(LOG_ROOT, "cygl1_coherency_full_join_log.csv"), index=False)
    print(f"\nWrote logs to {LOG_ROOT}")
    print(f"Threshold={COHERENCY_THRESHOLD}, seed={RANDOM_SEED}, "
          f"date range={BEG_DATE}-{END_DATE}")


if __name__ == "__main__":
    main()
