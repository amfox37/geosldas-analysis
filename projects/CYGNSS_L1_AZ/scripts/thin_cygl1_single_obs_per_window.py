#!/usr/bin/env python3
"""
Thin CYGNSS L1 preprocessed obs files to at most one observation per land-
analysis window, spaced far enough apart that their Gaspari-Cohn
compact-support ellipses (xcompact=ycompact=2.5deg in the actual R-sweep nml,
nml_HSAF_ascat_fov12p5_AZ_scaled_cygl1assim) cannot overlap. This isolates
each kept observation so the naive single-obs K*d identity
(pred_incr = fcstvar/(fcstvar+obsvar) * (obs-fcst)) is actually a valid test
of the analysis update, unlike the full obs stream (see memory
`cygl1-single-obs-diagnosis-roadmap`, Step 1, where ~56 CYGNSS L1 obs/cycle
fall inside any given obs's 2.5deg neighborhood).

CORRECTED 2026-08-21 (v2): the original version of this script grouped obs
into midnight-aligned bins ((ts%86400)//10800, i.e. [00:00,03:00), etc.) for
the isolation/thinning decision. GEOSldas's actual CYGNSS L1 assimilation
window, confirmed directly from source
(clsm_ensupd_read_obs.F90:read_obs_cygnss_l1_scalar, per its own comment) is
instead CENTERED on each analysis time:

    (date_time - dtstep_assim/2, date_time + dtstep_assim/2]

with dtstep_assim=LANDASSIM_DT=10800s (3h) and LANDASSIM_T0=000000 for every
experiment in this project (confirmed in the exeinp template's commented
defaults). That's a 1.5h offset from the old midnight-aligned bins, which
caused the isolation/thinning step to group obs by the WRONG partition of
the day -- contributing (along with legitimate owner-tile dedup, scaling
rejection, and model-based QC, all separate/expected) to a large gap between
obs kept by the (old) thinning script and obs actually reaching
assim_flag==1 in the run (104 kept -> only 47 actually assimilated in the
2026-08-21 v1 singleobs run, job 58019133 -- see
cygl1-single-obs-diagnosis-roadmap memory).

This version computes each obs's TRUE window anchor as the nearest absolute
multiple of dtstep_assim to its (year, day-of-year, ddm_timestamp_utc_sec)
timestamp -- an absolute, continuous quantity, so day-boundary-straddling
windows (e.g. the 00Z window, (22:30 prev day, 01:30]) are handled correctly
by loading and grouping obs across ALL source days together before deciding
what to keep, rather than thinning each day's file independently as v1 did.

Input schema (schema_version 0.5, "cygnss_tile_coefficient_preprocessor_netcdf"):
  dims: obs, support
  obs-dim vars used here: status, sp_lon, sp_lat, year, day,
                          ddm_timestamp_utc_sec, tile_start, tile_count
                          (+ all other obs-dim vars, copied through unchanged
                          for the kept rows)
  support-dim vars: tile_index0, tile_ig, tile_jg, coefficient,
                    coefficient_weight, active_pixel_count_by_tile
  Each obs i's support tiles are support-dim rows
  [tile_start(i) : tile_start(i)+tile_count(i)]. Thinning drops rows, so
  tile_start/tile_count must be remapped to the new, compacted support array.

Usage: python3 thin_cygl1_single_obs_per_window.py
"""
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

SRC_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
DST_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_singleobs"

BEG_DATE = "20200101"
END_DATE = "20200110"  # inclusive, 10 days

DTSTEP_ASSIM_SEC = 10800  # LANDASSIM_DT
LANDASSIM_T0_SEC = 0      # LANDASSIM_T0 (seconds of day), always 000000 in this project
MIN_SEP_DEG = 5.0  # 2 x xcompact (2.5 deg) -- Gaspari-Cohn ellipses cannot overlap

LOG_ROWS = []


def load_one_file(src_path):
    with nc.Dataset(src_path) as src:
        status = src.variables["status"][:]
        lon = src.variables["sp_lon"][:]
        lat = src.variables["sp_lat"][:]
        year = src.variables["year"][:]
        day = src.variables["day"][:]  # day-of-year, 1-based
        secs = src.variables["ddm_timestamp_utc_sec"][:]

        good = np.where(status == 1)[0]

        # Absolute epoch seconds, via pandas: Jan 1 of `year` + (day-of-year - 1) days + secs.
        base = pd.to_datetime(year[good].astype(str), format="%Y").values.astype("datetime64[s]")
        abs_ts = (
            base.astype("int64")
            + (day[good].astype(np.int64) - 1) * 86400
            + np.round(secs[good]).astype(np.int64)
        )

        # Nearest synoptic analysis time (absolute seconds), i.e. the true
        # GEOSldas window this obs falls into: (anchor-dtstep/2, anchor+dtstep/2].
        anchor = np.round((abs_ts - LANDASSIM_T0_SEC) / DTSTEP_ASSIM_SEC).astype(np.int64) \
            * DTSTEP_ASSIM_SEC + LANDASSIM_T0_SEC

    return good, lon[good], lat[good], anchor


def main():
    dates = pd.date_range(BEG_DATE, END_DATE, freq="D")
    src_paths = []
    for dt in dates:
        y, m, ymd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%Y%m%d")
        pattern = os.path.join(SRC_ROOT, f"Y{y}", f"M{m}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{ymd}_all_cyg.nc4")
        matches = glob.glob(pattern)
        if not matches:
            print(f"WARNING: no source file for {ymd} ({pattern})", file=sys.stderr)
            continue
        src_paths.append(matches[0])

    # Load every obs (across all days) with its TRUE absolute window anchor.
    all_rows = []  # (file_idx, obs_idx_in_good, lon, lat, anchor)
    for fi, src_path in enumerate(src_paths):
        good, lon, lat, anchor = load_one_file(src_path)
        for k in range(len(good)):
            all_rows.append((fi, int(good[k]), float(lon[k]), float(lat[k]), int(anchor[k])))

    df = pd.DataFrame(all_rows, columns=["file_idx", "obs_idx", "lon", "lat", "anchor"])
    total_orig = len(df)
    print(f"Loaded {total_orig} good obs across {len(src_paths)} days", flush=True)

    # Greedily thin within each TRUE window (anchor), across file boundaries.
    kept_mask = np.zeros(len(df), dtype=bool)
    for anchor, grp in df.groupby("anchor"):
        kept_idx = []
        for idx, row in grp.iterrows():
            ok = True
            for j in kept_idx:
                d = np.hypot(row.lon - df.loc[j, "lon"], row.lat - df.loc[j, "lat"])
                if d < MIN_SEP_DEG:
                    ok = False
                    break
            if ok:
                kept_idx.append(idx)
        kept_mask[kept_idx] = True
        for idx in kept_idx:
            r = df.loc[idx]
            LOG_ROWS.append((
                os.path.basename(src_paths[int(r.file_idx)]), int(anchor), int(r.obs_idx),
                float(r.lon), float(r.lat),
            ))

    df["kept"] = kept_mask
    total_kept = int(kept_mask.sum())
    print(f"Total: {total_orig} -> {total_kept} kept across {len(src_paths)} days "
          f"(true centered-window grouping)", flush=True)

    # Write out per-day files, keeping only the kept obs, remapping support-tile indices.
    for fi, src_path in enumerate(src_paths):
        kept_obs_idx = sorted(df.loc[(df.file_idx == fi) & df.kept, "obs_idx"].tolist())
        with nc.Dataset(src_path) as src:
            tile_start = src.variables["tile_start"][:]
            tile_count = src.variables["tile_count"][:]
            n_support = src.dimensions["support"].size

            kept = np.array(kept_obs_idx, dtype=int)
            n_kept = len(kept)

            new_tile_start = np.zeros(n_kept, dtype=tile_start.dtype)
            new_tile_count = tile_count[kept].copy()
            support_take = []
            cursor = 0
            for k, i in enumerate(kept):
                new_tile_start[k] = cursor
                s0, s1 = int(tile_start[i]), int(tile_start[i]) + int(tile_count[i])
                support_take.append(np.arange(s0, s1))
                cursor += (s1 - s0)
            support_take = np.concatenate(support_take) if support_take else np.array([], dtype=int)

            dst_path = os.path.join(
                DST_ROOT,
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
                dst.setncattr("thinning_applied",
                               "single-obs-per-true-centered-window (v2), min_sep_deg=%.2f" % MIN_SEP_DEG)
                dst.setncattr("thinning_script",
                               "geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/thin_cygl1_single_obs_per_window.py")
                dst.setncattr("thinning_source_file", src_path)
                dst.setncattr("thinning_n_obs_original", int((df.file_idx == fi).sum()))
                dst.setncattr("thinning_n_obs_kept", n_kept)

        print(f"{pd.Timestamp(dates[fi]).strftime('%Y%m%d')}: "
              f"{(df.file_idx==fi).sum()} -> {n_kept} kept  ({dst_path})")

    log_df = pd.DataFrame(LOG_ROWS, columns=["source_file", "window_anchor_epoch_sec", "obs_idx", "sp_lon", "sp_lat"])
    log_path = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygl1_singleobs_thinning_log.csv"
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    log_df.to_csv(log_path, index=False)
    print(f"Wrote thinning log: {log_path} ({len(log_df)} kept obs)")


if __name__ == "__main__":
    main()
