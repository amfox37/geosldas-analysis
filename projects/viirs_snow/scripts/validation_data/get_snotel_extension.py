#!/usr/bin/env python3
"""Extend the SNOTEL daily WTEQ (SWE) and SNWD (snow depth) record to 2026-07-31.

Logic follows projects/SNOTEL/scripts/download_snotel_swe_snwd.py (AWDB SOAP for the station list and metadata,
AWDB REST v1 for daily data), with two differences to respect the file quota and the servers:
  * the REST data call is made for batches of stations (comma-separated triplets), not one call per station;
  * output is ONE consolidated parquet (no per-station CSVs).

It fetches OVERLAP_BEGIN..END, checks the overlap days against the existing consolidated parquet, and writes only
dates >= OUT_BEGIN. Units are the AWDB native units (inches), as in the existing parquet.

Run (login node, GEOSpyD py3.13):
  python get_snotel_extension.py --out /discover/nobackup/projects/land_da/VIIRS_testing/validation_data/snotel
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests
import zeep

WSDL_URL = "https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL"
REST = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data"
ELEMENTS = ["WTEQ", "SNWD"]
EXISTING = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet")
UA = {"User-Agent": "snow-validation-data/1.0"}


def retry(fn, tries=4, backoff=2.0):
    for i in range(tries):
        try:
            return fn()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(backoff * 2 ** i)


def station_list_and_meta():
    client = zeep.Client(WSDL_URL)
    trip = list(retry(lambda: client.service.getStations(networkCds="SNTL", logicalAnd=False)))
    rows = []
    for i in range(0, len(trip), 200):
        chunk = trip[i:i + 200]
        for m in retry(lambda: client.service.getStationMetadataMultiple(stationTriplets=chunk)):
            rows.append({k: getattr(m, k, None) for k in
                         ("stationTriplet", "name", "latitude", "longitude", "elevation", "state", "countyName", "huc",
                          "beginDate", "endDate")})
    meta = pd.DataFrame(rows).dropna(subset=["stationTriplet"])
    return trip, meta


def fetch_batch(session, triplets, begin, end):
    params = {"stationTriplets": ",".join(triplets), "beginDate": begin, "endDate": end,
              "elements": ",".join(ELEMENTS), "duration": "DAILY"}
    r = retry(lambda: session.get(REST, params=params, timeout=120, headers=UA))
    r.raise_for_status()
    out = []
    for st in r.json():
        trip = st.get("stationTriplet")
        for blk in st.get("data", []) or []:
            se = blk.get("stationElement", {}) or {}
            el, du = se.get("elementCode"), se.get("durationName")
            if el in ELEMENTS and du == "DAILY":
                v = pd.DataFrame(blk.get("values", []) or [])
                if v.empty or "value" not in v.columns:
                    continue  # block with no valid values in the window (all missing)
                v = v[["date", "value"]].rename(columns={"value": el})
                v["stationTriplet"] = trip
                out.append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--overlap-begin", default="2025-09-24")
    ap.add_argument("--out-begin", default="2025-10-01")
    ap.add_argument("--end", default="2026-07-31")
    ap.add_argument("--batch", type=int, default=25)
    ap.add_argument("--sleep", type=float, default=1.0)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    trip, meta = station_list_and_meta()
    print(f"{len(trip)} SNTL stations from AWDB")
    session = requests.Session()
    parts = []
    for i in range(0, len(trip), a.batch):
        parts += fetch_batch(session, trip[i:i + a.batch], a.overlap_begin, a.end)
        time.sleep(a.sleep)
        if (i // a.batch) % 10 == 0:
            print(f"  fetched {min(i + a.batch, len(trip))}/{len(trip)} stations", flush=True)
    raw = pd.concat(parts, ignore_index=True)
    raw["date"] = pd.to_datetime(raw["date"])
    # merge elements per (station, date)
    df = None
    for el in ELEMENTS:
        sub = raw[["stationTriplet", "date", el]].dropna(subset=[el]) if el in raw.columns else None
        if sub is None or sub.empty:
            continue
        df = sub if df is None else df.merge(sub, on=["stationTriplet", "date"], how="outer")
    for el in ELEMENTS:
        if el not in df.columns:
            df[el] = float("nan")
    df = df.sort_values(["stationTriplet", "date"]).reset_index(drop=True)

    # ---- overlap check against the existing parquet
    ex = pd.read_parquet(EXISTING, columns=["date", "WTEQ", "SNWD", "stationTriplet"])
    ex = ex[ex["date"] >= a.overlap_begin]
    ov = df[df["date"] < a.out_begin]
    cmp = ex.merge(ov, on=["stationTriplet", "date"], how="outer", suffixes=("_old", "_new"), indicator=True)
    both = cmp[cmp["_merge"] == "both"]
    rep = {"overlap_days": f"{a.overlap_begin}..{a.out_begin}(excl)",
           "rows_old": int(len(ex)), "rows_new": int(len(ov)), "rows_both": int(len(both)),
           "only_old": int((cmp["_merge"] == "left_only").sum()), "only_new": int((cmp["_merge"] == "right_only").sum())}
    bad = pd.Series(False, index=both.index)
    for el in ELEMENTS:
        x, y = both[f"{el}_old"], both[f"{el}_new"]
        same = (x.isna() & y.isna()) | ((x - y).abs() < 1e-6)
        rep[f"{el}_mismatch_in_overlap"] = int((~same).sum())
        bad |= ~same
    # keep the differing rows (value revisions, or rows present in only one product) for inspection
    diff = pd.concat([both[bad].assign(kind="value_differs"),
                      cmp[cmp["_merge"] != "both"].assign(kind=lambda d: d["_merge"].astype(str))])
    diff.drop(columns="_merge").to_csv(out / "snotel_overlap_differences.csv", index=False)

    new = df[df["date"] >= a.out_begin].copy()
    new = new.merge(meta[["stationTriplet", "name", "latitude", "longitude", "elevation"]], on="stationTriplet", how="left")
    fn = out / f"snotel_daily_wteq_snwd_{a.out_begin.replace('-', '')}_{a.end.replace('-', '')}.parquet"
    new.to_parquet(fn, index=False)
    meta.to_csv(out / "snotel_stations_metadata.csv", index=False)

    old_st = set(pd.read_parquet(EXISTING, columns=["stationTriplet"])["stationTriplet"].unique())
    new_st = set(new["stationTriplet"].unique())
    nd = new.groupby("date").size()
    rep.update({"output": str(fn), "rows": int(len(new)), "stations": len(new_st),
                "date_min": str(new["date"].min().date()), "date_max": str(new["date"].max().date()),
                "n_dates": int(nd.shape[0]), "expected_dates": int((pd.Timestamp(a.end) - pd.Timestamp(a.out_begin)).days + 1),
                "stations_new_since_existing_parquet": sorted(new_st - old_st),
                "stations_in_existing_but_no_new_data": sorted(old_st - new_st),
                "stations_per_day_min_max": [int(nd.min()), int(nd.max())],
                "WTEQ_range_in": [float(new["WTEQ"].min()), float(new["WTEQ"].max())],
                "SNWD_range_in": [float(new["SNWD"].min()), float(new["SNWD"].max())],
                "units": "inches (AWDB native), as in the existing parquet"})
    (out / "snotel_extension_report.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
