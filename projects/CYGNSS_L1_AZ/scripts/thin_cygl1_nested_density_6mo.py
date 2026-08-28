#!/usr/bin/env python3
"""
Build the DA-sparse and DA-intermediate CYGNSS L1 obs streams for the 6-month
(2020-04-01 to 2020-10-31) paired thinning-density experiment, nested so that
DA-sparse obs subset-of DA-intermediate obs subset-of the original full stream
(changes attributable to added observations, not a different sample).

DA-sparse: identical algorithm to thin_cygl1_single_obs_per_window.py (single
obs per true centered 3-hr assimilation window, min pairwise separation
MIN_SEP_SPARSE_DEG=5.0deg -- unchanged, already validated in Steps 1.5/1.5v2).

DA-intermediate: starts from DA-sparse's kept set (force-kept, guarantees the
nesting property) and greedily admits additional candidate obs (fixed time
order) subject to a cap on each obs's own "local-interaction count" -- the
number of OTHER kept obs (same window) whose Gaspari-Cohn ellipse
(dlon/xcompact)^2+(dlat/ycompact)^2 < 1 contains it, using the ACTUAL
xcompact=ycompact=1.25deg from nml_HSAF_ascat_fov12p5_AZ_scaled_cygl1assim
(corrected 2026-08-21 from an earlier, wrong 2.5deg figure in project memory).

Run in CALIBRATE mode first (default) to sweep candidate caps and report the
resulting median/p90 local-interaction-count distribution, before locking in
an actual cap and writing thinned obs files (--build-intermediate CAP).
"""
import argparse
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

SRC_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
DST_SPARSE = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_sparse_6mo"
DST_INTERMEDIATE = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_intermediate_6mo"

BEG_DATE = "20200101"  # extended 2026-08-21 from 20200401: actual restart source available
                       # (LS_OLv8_M36_v2, global, proven working) is Jan-1-anchored, not Apr-1
                       # (an Apr-1 restart from OLv8_M36_all_sensors_AZ's own AZ-cropped output
                       # hit a real ldas_setup bug -- domain-cropped .til.domain file mishandled
                       # as if it were netCDF -- so all four paired experiments now restart from
                       # the standard global source and simply get 3 extra months of spin-up)
END_DATE = "20201031"  # inclusive; overridable via --beg-date/--end-date (windows/anchors are
                       # purely per-3hr-timestamp with no cross-month state, so any date range
                       # can be (re)generated independently -- used 2026-08-25 to extend
                       # DA-intermediate's obs 14 more months without touching Jan-Oct 2020)

DTSTEP_ASSIM_SEC = 10800  # LANDASSIM_DT
LANDASSIM_T0_SEC = 0

XCOMPACT_DEG = 1.25  # corrected 2026-08-21; was wrongly documented as 2.5 in earlier memory
YCOMPACT_DEG = 1.25
MIN_SEP_SPARSE_DEG = 5.0  # 2x the OLDER (wrong) 2.5deg figure; still > 2x1.25=2.5, so DA-sparse
                          # remains fully isolated under the corrected xcompact -- no change needed

LOG_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"


def load_all_obs():
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

    all_rows = []
    for fi, src_path in enumerate(src_paths):
        with nc.Dataset(src_path) as src:
            status = src.variables["status"][:]
            lon = src.variables["sp_lon"][:]
            lat = src.variables["sp_lat"][:]
            year = src.variables["year"][:]
            day = src.variables["day"][:]
            secs = src.variables["ddm_timestamp_utc_sec"][:]

            good = np.where(status == 1)[0]
            base = pd.to_datetime(year[good].astype(str), format="%Y").values.astype("datetime64[s]")
            abs_ts = (
                base.astype("int64")
                + (day[good].astype(np.int64) - 1) * 86400
                + np.round(secs[good]).astype(np.int64)
            )
            anchor = np.round((abs_ts - LANDASSIM_T0_SEC) / DTSTEP_ASSIM_SEC).astype(np.int64) \
                * DTSTEP_ASSIM_SEC + LANDASSIM_T0_SEC

            for k in range(len(good)):
                all_rows.append((fi, int(good[k]), float(lon[k]), float(lat[k]), int(anchor[k])))

    df = pd.DataFrame(all_rows, columns=["file_idx", "obs_idx", "lon", "lat", "anchor"])
    print(f"Loaded {len(df)} good obs across {len(src_paths)} days ({BEG_DATE}-{END_DATE})", flush=True)
    return df, src_paths, dates


def build_sparse(df):
    """Single obs per true centered window, min_sep=MIN_SEP_SPARSE_DEG -- identical algorithm
    to thin_cygl1_single_obs_per_window.py."""
    kept_mask = np.zeros(len(df), dtype=bool)
    for anchor, grp in df.groupby("anchor"):
        kept_idx = []
        for idx, row in grp.iterrows():
            ok = True
            for j in kept_idx:
                d = np.hypot(row.lon - df.loc[j, "lon"], row.lat - df.loc[j, "lat"])
                if d < MIN_SEP_SPARSE_DEG:
                    ok = False
                    break
            if ok:
                kept_idx.append(idx)
        kept_mask[kept_idx] = True
    return kept_mask


def local_interaction_counts(df, kept_mask):
    """For each kept obs, count OTHER kept obs in the same window that could jointly influence
    some shared target tile -- i.e. whose own xcompact/ycompact ellipse (centered on that obs)
    overlaps this obs's ellipse. For two identical axis-aligned ellipses of semi-axes
    (xcompact,ycompact), that overlap condition is (dlon/(2*xcompact))^2+(dlat/(2*ycompact))^2<1
    (Minkowski-sum / "sum of radii" boundary), NOT distance<xcompact (which is the much
    stricter, wrong condition -- one obs's CENTER falling inside the other's ellipse -- fixed
    2026-08-21 after the first version of this function used the wrong radius and made every
    min_sep>=1.5deg trivially show zero interaction)."""
    counts = np.zeros(len(df), dtype=int)
    kept_df = df[kept_mask]
    for anchor, grp in kept_df.groupby("anchor"):
        lons = grp["lon"].values
        lats = grp["lat"].values
        idxs = grp.index.values
        n = len(idxs)
        for a in range(n):
            c = 0
            for b in range(n):
                if a == b:
                    continue
                d2 = ((lons[a] - lons[b]) / (2 * XCOMPACT_DEG)) ** 2 + ((lats[a] - lats[b]) / (2 * YCOMPACT_DEG)) ** 2
                if d2 < 1.0:
                    c += 1
            counts[idxs[a]] = c
    return counts


def build_intermediate_candidate(df, sparse_mask, min_sep_deg):
    """Greedy nested relaxation: start from sparse_mask (force-kept), walk remaining obs in
    fixed (file, obs_idx) order, admit if it is >= min_sep_deg from every already-kept obs in
    its window. Unlike a hard per-obs neighbor-count cap (which saturates to a degenerate
    uniform distribution -- every kept obs pinned exactly at the cap, since a 7-month domain has
    enough raw candidates to fill any cap almost everywhere), a distance threshold produces a
    genuine, non-degenerate local-interaction-count distribution: whether a kept pair ends up
    inside each other's TRUE Gaspari-Cohn ellipse (radius xcompact=1.25deg) depends on local
    obs geometry, not just the enforced exclusion radius, as long as min_sep_deg is allowed to
    be smaller than 2xcompact=2.5deg (the sparse threshold, 5.0deg, enforces zero overlap;
    anything below 2.5deg starts letting some -- not all -- pairs fall inside each other's true
    ellipse)."""
    kept_mask = sparse_mask.copy()
    remaining = df[~sparse_mask].sort_values(["file_idx", "obs_idx"])

    window_members = {}
    for anchor, grp in df[kept_mask].groupby("anchor"):
        window_members[anchor] = list(zip(grp["lon"].values, grp["lat"].values))

    for idx, row in remaining.iterrows():
        anchor = row["anchor"]
        members = window_members.get(anchor, [])
        ok = True
        for (mlon, mlat) in members:
            if np.hypot(row.lon - mlon, row.lat - mlat) < min_sep_deg:
                ok = False
                break
        if not ok:
            continue
        kept_mask[idx] = True
        window_members.setdefault(anchor, []).append((row.lon, row.lat))

    return kept_mask


def write_thinned_files(df, kept_mask, src_paths, dates, dst_root, label):
    for fi, src_path in enumerate(src_paths):
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
                               "geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/thin_cygl1_nested_density_6mo.py")
                dst.setncattr("thinning_source_file", src_path)
                dst.setncattr("thinning_n_obs_original", int((df.file_idx == fi).sum()))
                dst.setncattr("thinning_n_obs_kept", n_kept)
        if fi % 30 == 0:
            print(f"  {pd.Timestamp(dates[fi]).strftime('%Y%m%d')}: "
                  f"{(df.file_idx==fi).sum()} -> {n_kept} kept  ({dst_root})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true", help="sweep min_sep thresholds, report only, write nothing")
    ap.add_argument("--build-sparse", action="store_true")
    ap.add_argument("--build-intermediate", type=float, default=None,
                     help="min_sep_deg to lock in and write files")
    ap.add_argument("--beg-date", type=str, default=None, help="override BEG_DATE (YYYYMMDD)")
    ap.add_argument("--end-date", type=str, default=None, help="override END_DATE (YYYYMMDD, inclusive)")
    args = ap.parse_args()

    global BEG_DATE, END_DATE
    if args.beg_date:
        BEG_DATE = args.beg_date
    if args.end_date:
        END_DATE = args.end_date

    df, src_paths, dates = load_all_obs()
    sparse_mask = build_sparse(df)
    print(f"DA-sparse: {sparse_mask.sum()} kept (of {len(df)})", flush=True)
    sparse_counts = local_interaction_counts(df, sparse_mask)
    sc = sparse_counts[sparse_mask]
    print(f"  DA-sparse local-interaction counts: median={np.median(sc):.2f}, "
          f"p90={np.percentile(sc,90):.2f}, max={sc.max() if len(sc) else 'n/a'}", flush=True)

    if args.calibrate or (not args.build_sparse and args.build_intermediate is None):
        print("\n=== Calibration sweep for DA-intermediate (min_sep_deg thresholds) ===", flush=True)
        rows = []
        for min_sep in [2.49, 2.47, 2.45, 2.42, 2.40, 2.37, 2.35, 2.32, 2.30, 2.25]:
            inter_mask = build_intermediate_candidate(df, sparse_mask, min_sep)
            counts = local_interaction_counts(df, inter_mask)
            c = counts[inter_mask]
            n_kept = int(inter_mask.sum())
            n_added = n_kept - int(sparse_mask.sum())
            median_c, p90_c, p99_c = np.median(c), np.percentile(c, 90), np.percentile(c, 99)
            frac_zero = float((c == 0).mean())
            rows.append((min_sep, n_kept, n_added, median_c, p90_c, p99_c, c.max(), frac_zero))
            print(f"min_sep={min_sep:.2f}deg: kept={n_kept} (+{n_added} vs sparse), "
                  f"median_local_count={median_c:.2f}, p90={p90_c:.2f}, p99={p99_c:.2f}, "
                  f"max={c.max()}, frac_zero_neighbors={frac_zero:.3f}",
                  flush=True)
        cal_df = pd.DataFrame(rows, columns=["min_sep_deg", "n_kept", "n_added_vs_sparse",
                                              "median_local_count", "p90_local_count",
                                              "p99_local_count", "max_local_count", "frac_zero_neighbors"])
        out_path = os.path.join(LOG_ROOT, "cygl1_intermediate_thinning_calibration.csv")
        cal_df.to_csv(out_path, index=False)
        print(f"\nWrote calibration table: {out_path}")
        return

    if args.build_sparse:
        write_thinned_files(df, sparse_mask, src_paths, dates, DST_SPARSE,
                             "single-obs-per-true-centered-window, min_sep_deg=%.2f" % MIN_SEP_SPARSE_DEG)
        log_df = df[sparse_mask][["file_idx", "obs_idx", "lon", "lat", "anchor"]].copy()
        log_df.to_csv(os.path.join(LOG_ROOT, "cygl1_sparse_6mo_thinning_log.csv"), index=False)
        print(f"DA-sparse: wrote {sparse_mask.sum()} obs to {DST_SPARSE}")

    if args.build_intermediate is not None:
        min_sep = args.build_intermediate
        inter_mask = build_intermediate_candidate(df, sparse_mask, min_sep)
        assert np.all(inter_mask[sparse_mask]), "nesting violated: sparse not subset of intermediate"
        write_thinned_files(df, inter_mask, src_paths, dates, DST_INTERMEDIATE,
                             f"nested-superset-of-sparse, min_sep_deg={min_sep}, "
                             f"xcompact=ycompact={XCOMPACT_DEG}deg")
        log_df = df[inter_mask][["file_idx", "obs_idx", "lon", "lat", "anchor"]].copy()
        log_df["in_sparse"] = sparse_mask[inter_mask]
        log_df.to_csv(os.path.join(LOG_ROOT, "cygl1_intermediate_6mo_thinning_log.csv"), index=False)
        print(f"DA-intermediate (min_sep={min_sep}): wrote {inter_mask.sum()} obs to {DST_INTERMEDIATE}")


if __name__ == "__main__":
    main()
