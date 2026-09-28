#!/usr/bin/env python3
"""
Replay the two per-tile CYGNSS L1 obs-selection rules in GEOSldas and flag where they disagree.

  reader   (clsm_ensupd_read_obs.F90, read_obs_CYGNSS_L1_scalar): per owner tile, the non-nodata
           obs with the smallest sp_nearest_tile_distance_km among obs INSIDE the assimilation
           window (t - dt/2, t + dt/2]. Supplies the observed value.
  operator (cygnss_preprocessed_obs.F90, cygnss_preproc_load + cygnss_preproc_find_obs): per owner
           tile, the obs with the smallest sp_nearest_tile_distance_km among ALL obs in every daily
           file the window touches -- no time-window or nodata filter. Supplies sp_inc_angle and the
           support-tile coefficients used to compute the model prediction.

When they pick different obs, the prediction is computed with another obs's incidence angle and
footprint than the obs it is compared against. Output: one row per (tile, cycle) where the reader
keeps an obs, with a mismatch flag and the time / incidence-angle / footprint differences.

Usage:
  cygl1_operator_selection_mismatch.py --obs-dir .../CYGNSS_L1 --tag full \\
      --tile-ref <any tavg24_1d_lnd_Nt file of the experiment domain> --start 20200101 --end 20220101
"""
import argparse
import glob
import os

import netCDF4 as nc
import numpy as np
import pandas as pd

OUT_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/operator_selection_mismatch"
NODATA, TOL = -9999.0, 1e-3
VARS = ["sp_nearest_tile_ig", "sp_nearest_tile_jg", "observed_y_db", "year", "day",
        "ddm_timestamp_utc_sec", "sp_nearest_tile_distance_km", "sp_inc_angle", "tile_start", "tile_count"]


def read_day(obs_dir, day, ij_to_row):
    f = os.path.join(obs_dir, f"Y{day:%Y}", f"M{day:%m}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{day:%Y%m%d}_all_cyg.nc4")
    if not os.path.exists(f):
        return None
    with nc.Dataset(f) as d:
        if d.dimensions["obs"].size == 0:
            return None
        df = pd.DataFrame({v: np.asarray(d[v][:]) for v in VARS})
        tig, tjg = np.asarray(d["tile_ig"][:]), np.asarray(d["tile_jg"][:])
    # footprint signature: the sorted support tiles (ig,jg) of each obs
    df["support"] = [tuple(sorted(zip(tig[s:s + c], tjg[s:s + c]))) for s, c in zip(df.tile_start, df.tile_count)]
    df["row"] = [ij_to_row.get((i, j), -1) for i, j in zip(df.sp_nearest_tile_ig, df.sp_nearest_tile_jg)]
    df["t"] = (pd.to_datetime(df.year.astype(str), format="%Y") + pd.to_timedelta(df.day - 1, "D")
               + pd.to_timedelta(np.rint(df.ddm_timestamp_utc_sec), "s"))
    df["file_day"] = day
    df["order"] = np.arange(len(df))  # file order; ties go to the first index in both Fortran loops
    return df[df.row >= 0]


def first_min(g):
    return g.sort_values(["sp_nearest_tile_distance_km", "file_day", "order"]).iloc[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--obs-dir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--tile-ref", required=True, help="tile-space nc4 with IG/JG (row order = OFA tilenum-1)")
    ap.add_argument("--start", default="20200101")
    ap.add_argument("--end", default="20220101")
    ap.add_argument("--dt-assim", type=int, default=10800, help="assimilation window length [s]")
    args = ap.parse_args()

    with nc.Dataset(args.tile_ref) as d:
        ig, jg = np.asarray(d["IG"][:]), np.asarray(d["JG"][:])
    ij_to_row = {(int(a), int(b)): r for r, (a, b) in enumerate(zip(ig, jg))}

    start, end = pd.Timestamp(args.start), pd.Timestamp(args.end)
    half = pd.Timedelta(seconds=args.dt_assim // 2)
    days = {}
    for day in pd.date_range(start - pd.Timedelta(days=1), end + pd.Timedelta(days=1), freq="D"):
        days[day] = read_day(args.obs_dir, day, ij_to_row)

    rows = []
    for t in pd.date_range(start + pd.Timedelta(hours=3), end, freq=f"{args.dt_assim}s"):
        lo, up = t - half, t + half
        files = [days.get(d) for d in pd.date_range(lo.normalize(), up.normalize(), freq="D")]
        files = [f for f in files if f is not None]
        if not files:
            continue
        cand = pd.concat(files, ignore_index=True)
        rd = cand[(cand.t > lo) & (cand.t <= up) & (np.abs(cand.observed_y_db - NODATA) > TOL)]
        if rd.empty:
            continue
        reader = rd.groupby("row", group_keys=False).apply(first_min, include_groups=False)
        oper = cand[cand.row.isin(reader.index)].groupby("row", group_keys=False).apply(first_min, include_groups=False)
        n_day = cand[cand.row.isin(reader.index)].groupby("row").size()
        n_win = rd.groupby("row").size()
        for r in reader.index:
            a, b = reader.loc[r], oper.loc[r]
            same = (a.file_day == b.file_day) and (a.order == b.order)
            rows.append((r, t, same, int(n_win[r]), int(n_day[r]),
                         abs((b.t - a.t).total_seconds()) / 3600.0,
                         float(b.sp_inc_angle - a.sp_inc_angle), a.support == b.support,
                         float(a.observed_y_db)))
    out = pd.DataFrame(rows, columns=["tile", "time", "same_obs", "n_in_window", "n_in_day_files",
                                      "dt_hours", "d_inc_angle", "same_support", "obs_db"])
    os.makedirs(OUT_ROOT, exist_ok=True)
    p = os.path.join(OUT_ROOT, f"selection_{args.tag}_{args.start}_{args.end}.parquet")
    out.to_parquet(p)
    mm = ~out.same_obs
    print(f"[{args.tag}] reader-selected tile-cycles: {len(out)}; operator uses a DIFFERENT obs: "
          f"{mm.sum()} ({100 * mm.mean():.1f}%); of those, different footprint: {(~out.same_support[mm]).mean() * 100:.1f}%, "
          f"median |dt| {out.dt_hours[mm].median():.1f} h, median |d_inc_angle| {out.d_inc_angle[mm].abs().median():.1f} deg")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
