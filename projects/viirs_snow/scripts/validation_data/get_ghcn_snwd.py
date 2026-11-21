#!/usr/bin/env python3
"""Build GHCN-Daily snow depth (SNWD) for 2020-08-01..2026-07-31 as one consolidated parquet per calendar year.

Script version of projects/GHCN_snwd/notebooks/ghcn_snwd_global_1998_present_build.ipynb, restricted to the period
needed for the VIIRS transition analysis. Same rules as the notebook:
  * source: NOAA/NCEI GHCN-Daily `by_year/<YYYY>.csv.gz` (streamed in chunks);
  * keep ELEMENT == "SNWD" (snow depth, mm);
  * drop missing (DATA_VALUE == -9999);
  * drop any value with a non-blank QFLAG (failed NOAA QA) (DROP_QFLAGGED = True in the notebook);
  * date from YYYYMMDD, snwd_mm as float32; flags and OBS_TIME kept as strings.
Differences: output is ONE parquet per year (file quota), rows outside [--begin, --end] are dropped, and raw .csv.gz
downloads are kept in <out>/tmp only until the year's parquet is verified, then deleted.

Also writes station metadata (stations with SNWD in the NOAA inventory overlapping the period), station coverage
(first/last date, n_days, max depth within the period), the ghcnd-*.txt source files, a build log and a manifest.

Run (login node, GEOSpyD py3.13):
  python get_ghcn_snwd.py --out /discover/nobackup/projects/land_da/VIIRS_testing/validation_data/ghcn
"""
from __future__ import annotations

import argparse
import gzip
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

BASE = "https://www.ncei.noaa.gov/pub/data/ghcn/daily"
UA = {"User-Agent": "snow-validation-data/1.0"}
COLS = ["station_id", "date_str", "element", "data_value", "mflag", "qflag", "sflag", "obs_time"]
DT = {"station_id": "string", "date_str": "string", "element": "string", "data_value": "Int32",
      "mflag": "string", "qflag": "string", "sflag": "string", "obs_time": "string"}


def download(url: str, dest: Path, tries: int = 4) -> Path:
    tmp = dest.with_suffix(dest.suffix + ".part")
    for i in range(tries):
        try:
            with requests.get(url, stream=True, timeout=180, headers=UA) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for ch in r.iter_content(1 << 20):
                        if ch:
                            f.write(ch)
            tmp.replace(dest)
            return dest
        except Exception:
            tmp.unlink(missing_ok=True)
            if i == tries - 1:
                raise
            time.sleep(2 ** (i + 1))


def read_stations(p: Path) -> pd.DataFrame:
    spec = [(0, 11), (12, 20), (21, 30), (31, 37), (38, 40), (41, 71), (72, 75), (76, 79), (80, 85)]
    names = ["station_id", "lat", "lon", "elev", "state", "name", "gsn_flag", "hcn_crn_flag", "wmo_id"]
    df = pd.read_fwf(p, colspecs=spec, names=names, dtype=str)
    for c in ("lat", "lon", "elev"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("station_id", "state", "name", "gsn_flag", "hcn_crn_flag", "wmo_id"):
        df[c] = df[c].fillna("").str.strip()
    return df


def read_inventory(p: Path) -> pd.DataFrame:
    return pd.read_fwf(p, colspecs=[(0, 11), (12, 20), (21, 30), (31, 35), (36, 40), (41, 45)],
                       names=["station_id", "lat", "lon", "element", "first_year", "last_year"],
                       dtype={"station_id": str, "lat": float, "lon": float, "element": str,
                              "first_year": int, "last_year": int})


def process_year(year: int, gz: Path, begin: pd.Timestamp, end: pd.Timestamp, out_pq: Path) -> dict:
    n_snwd = n_missing = n_qflag = 0
    parts = []
    for chunk in pd.read_csv(gz, names=COLS, header=None, dtype=DT, chunksize=2_000_000, low_memory=False):
        c = chunk.loc[chunk["element"] == "SNWD"].copy()
        if c.empty:
            continue
        n_snwd += len(c)
        miss = c["data_value"].isna() | (c["data_value"] == -9999)
        n_missing += int(miss.sum())
        c = c.loc[~miss]
        for col in ("mflag", "qflag", "sflag", "obs_time"):
            c[col] = c[col].fillna("").astype("string").str.strip()
        flagged = c["qflag"] != ""
        n_qflag += int(flagged.sum())
        c = c.loc[~flagged]
        c["date"] = pd.to_datetime(c["date_str"], format="%Y%m%d", errors="coerce")
        c = c.loc[c["date"].notna() & (c["date"] >= begin) & (c["date"] <= end)]
        if c.empty:
            continue
        c["snwd_mm"] = c["data_value"].astype("float32")
        parts.append(c[["station_id", "date", "snwd_mm", "mflag", "qflag", "sflag", "obs_time"]])
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["station_id", "date", "snwd_mm", "mflag", "qflag", "sflag", "obs_time"])
    df = df.sort_values(["station_id", "date"]).reset_index(drop=True)
    df.to_parquet(out_pq, index=False)
    return {"year": year, "rows": int(len(df)), "stations": int(df["station_id"].nunique()),
            "snwd_rows_raw": n_snwd, "dropped_missing": n_missing, "dropped_qflagged": n_qflag,
            "date_min": str(df["date"].min().date()) if len(df) else None,
            "date_max": str(df["date"].max().date()) if len(df) else None,
            "n_dates": int(df["date"].nunique()) if len(df) else 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--begin", default="2020-08-01")
    ap.add_argument("--end", default="2026-07-31")
    ap.add_argument("--keep-raw", action="store_true", help="keep the downloaded by_year gz files")
    a = ap.parse_args()
    out = Path(a.out); tmp = out / "tmp"
    out.mkdir(parents=True, exist_ok=True); tmp.mkdir(exist_ok=True)
    b, e = pd.Timestamp(a.begin), pd.Timestamp(a.end)

    for f in ("ghcnd-stations.txt", "ghcnd-inventory.txt"):
        download(f"{BASE}/{f}", out / f)
    download(f"{BASE}/readme.txt", out / "ghcnd-readme.txt")
    stations, inv = read_stations(out / "ghcnd-stations.txt"), read_inventory(out / "ghcnd-inventory.txt")
    sel = inv.query("element == 'SNWD'").loc[lambda d: (d["last_year"] >= b.year)]
    meta = sel.merge(stations[["station_id", "elev", "state", "name", "gsn_flag", "hcn_crn_flag", "wmo_id"]],
                     on="station_id", how="left")
    meta.to_parquet(out / "ghcn_snwd_station_metadata.parquet", index=False)

    log, cov = [], []
    for y in range(b.year, e.year + 1):
        gz = download(f"{BASE}/by_year/{y}.csv.gz", tmp / f"{y}.csv.gz")
        with gzip.open(gz, "rb") as f:
            f.read(1)  # gzip sanity
        pq = out / f"ghcn_snwd_{y}.parquet"
        r = process_year(y, gz, b, e, pq)
        log.append(r); print(r, flush=True)
        d = pd.read_parquet(pq, columns=["station_id", "date", "snwd_mm"])
        if len(d):
            cov.append(d.groupby("station_id").agg(first_date=("date", "min"), last_date=("date", "max"),
                                                    n_days=("date", "count"), max_snwd_mm=("snwd_mm", "max")).reset_index())
        if not a.keep_raw and r["rows"] > 0:
            gz.unlink()
    c = pd.concat(cov).groupby("station_id").agg(first_date=("first_date", "min"), last_date=("last_date", "max"),
                                                  n_days=("n_days", "sum"), max_snwd_mm=("max_snwd_mm", "max")).reset_index()
    c = c.merge(stations[["station_id", "lat", "lon", "elev", "state", "name", "gsn_flag", "hcn_crn_flag", "wmo_id"]],
                on="station_id", how="left").sort_values("n_days", ascending=False)
    c.to_parquet(out / "ghcn_snwd_station_coverage.parquet", index=False)
    pd.DataFrame(log).to_csv(out / "ghcn_snwd_build_log.csv", index=False)
    manifest = {"dataset": "NOAA/NCEI GHCN-Daily SNWD subset (viirs_snow transition analysis)", "element": "SNWD",
                "units": "mm", "begin": a.begin, "end": a.end, "drop_missing": True, "drop_qflagged": True,
                "source": BASE + "/by_year/<YYYY>.csv.gz", "built_utc": datetime.now(timezone.utc).isoformat(),
                "columns": ["station_id", "date", "snwd_mm", "mflag", "qflag", "sflag", "obs_time"],
                "note": "qflag is always blank (flagged rows dropped); station coordinates/elevation are in ghcn_snwd_station_coverage.parquet",
                "years": log}
    (out / "ghcn_snwd_manifest.json").write_text(json.dumps(manifest, indent=2))
    print("done; stations with data:", len(c))


if __name__ == "__main__":
    main()
