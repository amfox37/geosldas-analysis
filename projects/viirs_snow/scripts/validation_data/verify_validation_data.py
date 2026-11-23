#!/usr/bin/env python3
"""Verify the SNOTEL / GHCN / IMS products in validation_data/ and write validation_report.json.

Checks: date coverage (first, last, missing days), row/station counts, value ranges, file counts and sizes, and spot
checks against the products that already existed in the repo (SNOTEL parquet, GHCN OL/DA validation time series, IMS
raw files). With --delete-ims-gz the raw IMS .gz cache is removed once the yearly NetCDFs are verified.

Run: python -I verify_validation_data.py --root /discover/nobackup/projects/land_da/VIIRS_testing/validation_data
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

REPO = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects")
BEGIN, END = pd.Timestamp("2020-08-01"), pd.Timestamp("2026-07-31")


def dir_stats(p: Path) -> dict:
    files = [f for f in p.rglob("*") if f.is_file()]
    return {"files": len(files), "bytes": int(sum(f.stat().st_size for f in files))}


def snotel(root: Path) -> dict:
    pq = next((root / "snotel").glob("snotel_daily_wteq_snwd_*.parquet"))
    d = pd.read_parquet(pq)
    dates = pd.date_range("2025-10-01", END)
    per_day = d.groupby("date").size()
    rep = {"file": str(pq), "rows": len(d), "stations": int(d.stationTriplet.nunique()),
           "first": str(d.date.min().date()), "last": str(d.date.max().date()),
           "missing_dates": [str(x.date()) for x in dates.difference(per_day.index)],
           "stations_per_day_min_max": [int(per_day.min()), int(per_day.max())],
           "WTEQ_in_min_max": [float(d.WTEQ.min()), float(d.WTEQ.max())],
           "SNWD_in_min_max": [float(d.SNWD.min()), float(d.SNWD.max())],
           "WTEQ_nan_frac": float(d.WTEQ.isna().mean()), "SNWD_nan_frac": float(d.SNWD.isna().mean())}
    rep["extension_report"] = json.loads((root / "snotel" / "snotel_extension_report.json").read_text())
    return rep


def ghcn(root: Path) -> dict:
    g = root / "ghcn"
    parts = [pd.read_parquet(p) for p in sorted(g.glob("ghcn_snwd_20*.parquet"))]
    d = pd.concat(parts, ignore_index=True)
    all_dates = pd.date_range(BEGIN, END)
    per_day = d.groupby("date").station_id.nunique()
    rep = {"rows": len(d), "stations": int(d.station_id.nunique()), "first": str(d.date.min().date()),
           "last": str(d.date.max().date()), "missing_dates": [str(x.date()) for x in all_dates.difference(per_day.index)],
           "stations_per_day_min_max": [int(per_day.min()), int(per_day.max())],
           "snwd_mm_min_max": [float(d.snwd_mm.min()), float(d.snwd_mm.max())],
           "qflag_nonblank_rows": int((d.qflag.fillna("") != "").sum()),
           "duplicates_station_date": int(d.duplicated(["station_id", "date"]).sum())}
    # spot check against the existing OL/DA validation time series (obs_snwd_mm, 10,119 selected stations)
    nc = netCDF4.Dataset(REPO / "GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_snwd_raw_timeseries_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.nc")
    st = np.asarray(nc["station"][:]).astype(str)
    t0 = pd.Timestamp("2000-01-01")
    tt = t0 + pd.to_timedelta(np.asarray(nc["time"][:]), unit="D")
    sel = np.where((tt >= BEGIN) & (tt <= pd.Timestamp("2024-12-31")))[0]
    old = np.ma.filled(nc["obs_snwd_mm"][sel[0]:sel[-1] + 1, :], np.nan)
    o = pd.DataFrame(old, index=tt[sel], columns=st).stack(dropna=True).rename("old").reset_index()
    o.columns = ["date", "station_id", "old"]
    m = o.merge(d[["station_id", "date", "snwd_mm"]], on=["station_id", "date"], how="outer", indicator=True)
    both = m[m["_merge"] == "both"]
    # new rows for stations in the old selection only (rows of other stations are expected to be "right_only")
    in_sel = m["station_id"].isin(set(st))
    rep["spot_check_vs_existing_validation_timeseries"] = {
        "period": "2020-08-01..2024-12-31", "old_station_days": int(len(o)),
        "present_in_both": int(len(both)),
        "equal_in_both": int((np.abs(both["old"] - both["snwd_mm"]) < 0.5).sum()),
        "only_in_old": int((m["_merge"] == "left_only").sum()),
        "only_in_new_for_old_stations": int(((m["_merge"] == "right_only") & in_sel).sum())}
    rep["files"] = dir_stats(g)
    rep["downstream_selection_rules_used_by_existing_validation"] = {
        k: str(getattr(nc, k)) for k in ("max_distance_km", "max_abs_elev_diff_m", "min_valid_days",
                                         "min_avg_snow_days_per_year", "model_snwd_definition")}
    return rep


def ims(root: Path, delete_gz: bool) -> dict:
    out = {}
    ok_all = True
    for y in (2025, 2026):
        f = root / "ims" / f"ims_snowcover_24km_{y}.nc4"
        if not f.exists():
            out[str(y)] = "missing"; ok_all = False; continue
        d = netCDF4.Dataset(f)
        n_exp = 365 if y == 2025 else 212
        v = d["ims_snowcover"]
        nan_days, vals = [], set()
        for ti in range(v.shape[0]):
            a = np.ma.filled(v[ti], np.nan)
            if np.isnan(a).all():
                nan_days.append(ti)
            else:
                vals |= set(np.unique(a[~np.isnan(a)]).astype(int).tolist())
        tm = np.asarray(d["time"][:])
        out[str(y)] = {"file": str(f), "n_time": int(v.shape[0]), "expected": n_exp, "all_nan_days_index": nan_days,
                       "values_present": sorted(vals), "time_first_last_days_since_Jan1": [float(tm[0]), float(tm[-1])],
                       "bytes": f.stat().st_size}
        ok_all &= (v.shape[0] == n_exp and not nan_days and vals <= {0, 1, 2, 3, 4})
    out["verified"] = bool(ok_all)
    out["raw_gz_files_before_delete"] = len(list((root / "ims" / "gz").glob("*.gz"))) if (root / "ims" / "gz").exists() else 0
    out["fetch_report"] = json.loads((root / "ims" / "ims_fetch_report.json").read_text()) if (root / "ims" / "ims_fetch_report.json").exists() else None
    if delete_gz and ok_all and (root / "ims" / "gz").exists():
        for p in (root / "ims" / "gz").glob("ims20*_00UTC_24km_v1.3.asc.gz"):
            p.unlink()
        (root / "ims" / "gz").rmdir()
        out["raw_gz_deleted"] = True
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--delete-ims-gz", action="store_true")
    a = ap.parse_args()
    root = Path(a.root)
    rep = {"snotel": snotel(root), "ghcn": ghcn(root), "ims": ims(root, a.delete_ims_gz)}
    rep["totals"] = dir_stats(root)
    (root / "validation_report.json").write_text(json.dumps(rep, indent=2, default=str))
    print(json.dumps({k: (v if k == "totals" else "see validation_report.json") for k, v in rep.items()}, indent=2))


if __name__ == "__main__":
    main()
