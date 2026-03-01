#!/usr/bin/env python3
"""
Bulk download SNOTEL daily SWE (WTEQ) and Snow Depth (SNWD) + station metadata.

- Station list + metadata: AWDB SOAP (WSDL)
- Time series: AWDB REST v1 data endpoint

Outputs:
  out_dir/
    snotel_station_triplets.txt
    stations_metadata.csv
    stations/
      <stationTriplet>.csv
    all_stations_daily_wteq_snwd.parquet

Refs:
- AWDB Web Service User Guide (SOAP methods like getStations, getStationMetadataMultiple) and WSDL
  https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL  (see NRCS guide)
- REST endpoint:
  https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data

Example useage:
    python download_snotel_swe_snwd.py --out ./snotel_test_10

    python download_snotel_swe_snwd.py --out ./snotel_test_10 --sleep 0.5
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, TypeVar

import pandas as pd
import requests
from tqdm import tqdm
import zeep

WSDL_URL = "https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL"
REST_DATA_ENDPOINT = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data"

NETWORK = "SNTL"                 # SNOTEL network code
ELEMENTS = ["WTEQ", "SNWD"]       # SWE and Snow Depth
DURATION = "DAILY"               # daily series
T = TypeVar("T")


def with_retries(fn: Callable[[], T], retries: int = 3, backoff: float = 1.0) -> T:
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            return fn()
        except Exception as e:
            last_err = e
            if attempt == retries - 1:
                break
            time.sleep(backoff * (2 ** attempt))
    assert last_err is not None
    raise last_err


def get_snotel_station_triplets(client: zeep.Client) -> List[str]:
    # SOAP method described in AWDB user guide: getStations(networkCds=..., logicalAnd=...)
    triplets = with_retries(
        lambda: client.service.getStations(networkCds=NETWORK, logicalAnd=False)
    )
    return list(triplets)


def chunked(seq: List[str], n: int) -> List[List[str]]:
    return [seq[i : i + n] for i in range(0, len(seq), n)]


def get_station_metadata_multiple(client: zeep.Client, triplets: List[str]) -> pd.DataFrame:
    """
    SOAP getStationMetadataMultiple returns StationMetadata objects.
    We keep: stationTriplet, name, latitude, longitude, elevation (and units if present).
    """
    # Some SOAP servers dislike huge lists; keep chunks moderate.
    rows: List[Dict[str, Any]] = []
    for chunk in tqdm(chunked(triplets, 200), desc="Fetching station metadata (SOAP)"):
        metas = with_retries(
            lambda: client.service.getStationMetadataMultiple(stationTriplets=chunk)
        )
        for m in metas:
            # zeep objects behave like attributes
            rows.append(
                {
                    "stationTriplet": getattr(m, "stationTriplet", None),
                    "name": getattr(m, "name", None),
                    "latitude": getattr(m, "latitude", None),
                    "longitude": getattr(m, "longitude", None),
                    "elevation": getattr(m, "elevation", None),
                    "state": getattr(m, "state", None),
                    "countyName": getattr(m, "countyName", None),
                    "huc": getattr(m, "huc", None),
                }
            )
    df = pd.DataFrame(rows).dropna(subset=["stationTriplet"])
    return df


def fetch_station_daily_elements(
    session: requests.Session,
    station_triplet: str,
    begin_date: str,
    end_date: str,
    elements: List[str],
    duration: str = DURATION,
    timeout: int = 60,
) -> pd.DataFrame:
    """
    REST call that requests multiple elements for one station.
    Returns df with columns: date, WTEQ, SNWD (as available).

    NOTE: AWDB returns a list with one station object (usually).
    Each element/duration is a block under ["data"] with ["values"].
    """
    params = {
        "stationTriplets": station_triplet,
        "beginDate": begin_date,
        "endDate": end_date,
        "elements": ",".join(elements),
        "duration": duration,
    }
    r = with_retries(
        lambda: session.get(REST_DATA_ENDPOINT, params=params, timeout=timeout)
    )
    r.raise_for_status()
    js = r.json()

    if not js:
        return pd.DataFrame(columns=["date"] + elements)

    blocks = js[0].get("data", [])
    series_by_elem: Dict[str, pd.DataFrame] = {}

    for blk in blocks:
        se = blk.get("stationElement", {}) or {}
        elem = se.get("elementCode")
        dur = se.get("durationName")
        if elem in elements and dur == duration:
            vals = blk.get("values", []) or []
            df = pd.DataFrame(vals)
            if df.empty:
                continue
            df = df[["date", "value"]].copy()
            df["date"] = pd.to_datetime(df["date"])
            df = df.rename(columns={"value": elem})
            series_by_elem[elem] = df

    # Merge on date (outer) so missing days/elements remain NaN
    out = None
    for elem in elements:
        df = series_by_elem.get(elem)
        if df is None:
            continue
        out = df if out is None else out.merge(df, on="date", how="outer")

    if out is None:
        return pd.DataFrame(columns=["date"] + elements)

    out = out.sort_values("date")
    # Ensure all element columns exist
    for elem in elements:
        if elem not in out.columns:
            out[elem] = pd.NA
    return out[["date"] + elements]


def safe_filename(station_triplet: str) -> str:
    return station_triplet.replace(":", "_")


def load_existing_station_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--begin", default="1999-10-01", help="YYYY-MM-DD (default: 1999-10-01)")
    ap.add_argument("--end", default="2025-09-30", help="YYYY-MM-DD (default: 2025-09-30)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--limit", type=int, default=None, help="Number of stations to download (default: all)")
    ap.add_argument("--sleep", type=float, default=0.2, help="Seconds to sleep between REST requests (default: 0.2)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing per-station CSVs")
    args = ap.parse_args()

    out_dir = Path(args.out)
    stations_dir = out_dir / "stations"
    out_dir.mkdir(parents=True, exist_ok=True)
    stations_dir.mkdir(parents=True, exist_ok=True)

    # SOAP client
    client = zeep.Client(WSDL_URL)

    print("Fetching SNOTEL station triplets (SOAP getStations)…")
    triplets = get_snotel_station_triplets(client)
    print(f"Found {len(triplets)} SNOTEL stations total")

    # Keep test subset
    if args.limit is not None:
        triplets = triplets[: args.limit]
    (out_dir / "snotel_station_triplets.txt").write_text("\n".join(triplets) + "\n")

    print("Fetching station metadata (SOAP getStationMetadataMultiple)…")
    meta = get_station_metadata_multiple(client, triplets)
    meta.to_csv(out_dir / "stations_metadata.csv", index=False)

    # REST time series
    session = requests.Session()
    ok, empty, failed, reused = 0, 0, 0, 0
    all_rows: List[pd.DataFrame] = []

    for st in tqdm(triplets, desc="Downloading DAILY WTEQ+SNWD"):
        fn = stations_dir / f"{safe_filename(st)}.csv"
        if fn.exists() and not args.overwrite:
            reused_df = load_existing_station_csv(fn)
            if not reused_df.empty:
                if "stationTriplet" not in reused_df.columns:
                    reused_df["stationTriplet"] = st
                all_rows.append(reused_df)
                reused += 1
            continue

        try:
            df = fetch_station_daily_elements(
                session=session,
                station_triplet=st,
                begin_date=args.begin,
                end_date=args.end,
                elements=ELEMENTS,
                duration=DURATION,
            )

            # Add stationTriplet for combined output
            df["stationTriplet"] = st

            if df.empty:
                empty += 1
                df.to_csv(fn, index=False)
            else:
                df.to_csv(fn, index=False)
                ok += 1
                all_rows.append(df)

        except Exception as e:
            failed += 1
            (stations_dir / f"{safe_filename(st)}.ERROR.txt").write_text(f"{type(e).__name__}: {e}\n")

        time.sleep(args.sleep)

    if all_rows:
        big = pd.concat(all_rows, ignore_index=True)
        # Attach lat/lon/elev/name (left join)
        big = big.merge(meta[["stationTriplet", "name", "latitude", "longitude", "elevation"]], on="stationTriplet", how="left")
        big.to_parquet(out_dir / "all_stations_daily_wteq_snwd.parquet", index=False)

    print(f"Done. ok={ok}, empty={empty}, failed={failed}")
    print(f"Reused existing station CSVs: {reused}")
    print(f"Wrote: {out_dir / 'stations_metadata.csv'}")
    print(f"Wrote per-station CSVs under: {stations_dir}")
    if all_rows:
        print(f"Wrote combined parquet: {out_dir / 'all_stations_daily_wteq_snwd.parquet'}")


if __name__ == "__main__":
    main()
