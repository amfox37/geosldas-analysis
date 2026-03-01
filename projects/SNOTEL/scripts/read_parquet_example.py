#!/usr/bin/env python3
"""Example: read and query SNOTEL combined parquet output."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--parquet",
        default="SNOTEL/output/all_stations_daily_wteq_snwd.parquet",
        help="Path to parquet file",
    )
    ap.add_argument("--station", default=None, help="Optional stationTriplet filter")
    ap.add_argument("--start", default=None, help="Optional start date YYYY-MM-DD")
    ap.add_argument("--end", default=None, help="Optional end date YYYY-MM-DD")
    args = ap.parse_args()

    parquet_path = Path(args.parquet)
    if not parquet_path.exists():
        raise FileNotFoundError(f"Parquet not found: {parquet_path}")

    df = pd.read_parquet(parquet_path)
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"], errors="coerce")

    print("Rows:", len(df))
    print("Columns:", list(df.columns))
    print("\nHead:")
    print(df.head(10).to_string(index=False))

    # Optional station filter
    if args.station:
        df = df[df["stationTriplet"] == args.station]
        print(f"\nAfter station filter ({args.station}): {len(df)} rows")

    # Optional date filters
    if args.start:
        df = df[df["date"] >= pd.to_datetime(args.start)]
    if args.end:
        df = df[df["date"] <= pd.to_datetime(args.end)]
    if args.start or args.end:
        print(
            f"After date filter ({args.start or '-inf'} to {args.end or '+inf'}): {len(df)} rows"
        )

    if {"WTEQ", "SNWD"}.issubset(df.columns):
        print("\nSummary stats for WTEQ and SNWD:")
        print(df[["WTEQ", "SNWD"]].describe().to_string())


if __name__ == "__main__":
    main()
