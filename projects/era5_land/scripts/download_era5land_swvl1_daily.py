#!/usr/bin/env python3
"""
Download ERA5-Land hourly surface soil moisture (layer 1) and create daily files.

What this script does
---------------------
For each month in the requested date range, it:
  1. Downloads hourly ERA5-Land data for volumetric_soil_water_layer_1
  2. Saves the hourly monthly file
  3. Computes daily means from the hourly file
  4. Saves:
       - one monthly daily-mean NetCDF, and/or
       - one NetCDF per day

Variable:
  volumetric_soil_water_layer_1  (ERA5-Land layer 1 soil moisture, ~0-7 cm)

Requirements
------------
  pip install "cdsapi>=0.7.7" xarray netCDF4 numpy pandas

CDS API setup
-------------
1. Create/login to your CDS account
2. Accept the Terms of Use for ERA5-Land on the dataset page
3. Put your API token in ~/.cdsapirc, for example:

   url: https://cds.climate.copernicus.eu/api
   key: <PERSONAL-ACCESS-TOKEN>

Notes
-----
- ERA5-Land is served from dataset: "reanalysis-era5-land"
- This script requests NetCDF output
- Daily means are computed from the hourly data
- By default, days are UTC days based on the downloaded "time" coordinate

Example
-------
python download_era5land_swvl1_daily.py \
    --start 2020-01-01 \
    --end 2020-03-31 \
    --outdir ./era5land_swvl1

Optional subsetting
-------------------
Use --area "north,west,south,east", e.g.
    --area "60,-130,20,-60"

This follows CDS ordering: North, West, South, East.
"""

from __future__ import annotations

import argparse
import calendar
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

import pandas as pd
import xarray as xr

try:
    import cdsapi
except ImportError as exc:
    raise SystemExit(
        "cdsapi is not installed. Install with: pip install 'cdsapi>=0.7.7'"
    ) from exc


DATASET = "reanalysis-era5-land"
VARIABLE = "volumetric_soil_water_layer_1"
TIME_STRINGS = [f"{hour:02d}:00" for hour in range(24)]


@dataclass(frozen=True)
class MonthWindow:
    year: int
    month: int
    first_day: pd.Timestamp
    last_day: pd.Timestamp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download ERA5-Land hourly swvl1 and create daily files."
    )
    parser.add_argument(
        "--start",
        required=True,
        help="Start date, inclusive, in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--end",
        required=True,
        help="End date, inclusive, in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for hourly and daily files.",
    )
    parser.add_argument(
        "--area",
        default=None,
        help=(
            "Optional spatial subset in CDS order: 'north,west,south,east'. "
            "Example: '60,-130,20,-60'"
        ),
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Do not download if monthly hourly file already exists.",
    )
    parser.add_argument(
        "--keep-hourly",
        action="store_true",
        help="Keep monthly hourly files after creating daily products.",
    )
    parser.add_argument(
        "--write-monthly-daily",
        action="store_true",
        help="Write one monthly daily-mean NetCDF per month.",
    )
    parser.add_argument(
        "--write-daily-files",
        action="store_true",
        help="Write one NetCDF file per day.",
    )
    parser.add_argument(
        "--engine",
        default=None,
        help=(
            "Optional xarray engine for reading/writing NetCDF, e.g. 'netcdf4'. "
            "Leave unset to let xarray decide."
        ),
    )

    args = parser.parse_args()

    if not args.write_monthly_daily and not args.write_daily_files:
        parser.error("Choose at least one output mode: --write-monthly-daily and/or --write-daily-files")

    try:
        start = pd.Timestamp(args.start)
        end = pd.Timestamp(args.end)
    except Exception as exc:
        parser.error(f"Could not parse start/end dates: {exc}")

    if end < start:
        parser.error("--end must be on or after --start")

    if args.area is not None:
        try:
            parts = [float(x.strip()) for x in args.area.split(",")]
            if len(parts) != 4:
                raise ValueError("Need exactly four comma-separated values")
        except Exception as exc:
            parser.error(f"Invalid --area value: {exc}")

    return args


def month_windows(start: pd.Timestamp, end: pd.Timestamp) -> Iterator[MonthWindow]:
    current = pd.Timestamp(year=start.year, month=start.month, day=1)

    while current <= end:
        year = current.year
        month = current.month
        last_dom = calendar.monthrange(year, month)[1]
        month_start = pd.Timestamp(year=year, month=month, day=1)
        month_end = pd.Timestamp(year=year, month=month, day=last_dom)

        first_day = max(month_start, start.normalize())
        last_day = min(month_end, end.normalize())

        yield MonthWindow(
            year=year,
            month=month,
            first_day=first_day,
            last_day=last_day,
        )

        current = month_end + pd.Timedelta(days=1)


def build_request(window: MonthWindow, area: Optional[str] = None) -> dict:
    request = {
        "variable": [VARIABLE],
        "year": [f"{window.year:04d}"],
        "month": [f"{window.month:02d}"],
        "day": [f"{day:02d}" for day in range(window.first_day.day, window.last_day.day + 1)],
        "time": TIME_STRINGS,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }

    if area is not None:
        request["area"] = [float(x.strip()) for x in area.split(",")]

    return request


def download_month(
    client: cdsapi.Client,
    window: MonthWindow,
    hourly_path: Path,
    area: Optional[str],
    skip_download: bool,
) -> None:
    if hourly_path.exists() and skip_download:
        print(f"[skip] Hourly file exists: {hourly_path}")
        return

    request = build_request(window, area=area)
    print(
        f"[download] {window.year:04d}-{window.month:02d} "
        f"{window.first_day.date()} to {window.last_day.date()} -> {hourly_path}"
    )
    client.retrieve(DATASET, request, str(hourly_path))


def ensure_time_sorted(ds: xr.Dataset) -> xr.Dataset:
    """
    Ensure the dataset has a sortable 'time' coordinate.
    ERA5-Land NetCDFs from CDS often use 'valid_time' instead of 'time'.
    """
    if "time" in ds.coords or "time" in ds.dims:
        pass
    elif "valid_time" in ds.coords or "valid_time" in ds.dims:
        ds = ds.rename({"valid_time": "time"})
    else:
        raise ValueError(
            f"Dataset does not contain a usable time coordinate. "
            f"Coords: {list(ds.coords)} | Dims: {list(ds.dims)}"
        )

    return ds.sortby("time")


def compute_daily_means(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert hourly ERA5-Land data to daily means.
    Handles ERA5-Land naming quirks:
      - time may be called 'valid_time'
      - soil moisture layer 1 may be called 'swvl1'
    """
    ds = ensure_time_sorted(ds)

    # ERA5-Land NetCDF usually stores the variable as 'swvl1'
    candidate_names = [
        VARIABLE,   # 'volumetric_soil_water_layer_1'
        "swvl1",
    ]
    found = [v for v in candidate_names if v in ds.data_vars]
    if not found:
        raise ValueError(
            f"Expected one of {candidate_names} not found. "
            f"Available variables: {list(ds.data_vars)}"
        )

    varname = found[0]
    ds = ds[[varname]]

    daily = ds.resample(time="1D").mean(keep_attrs=True)

    # Rename to the long descriptive name if you want consistent downstream naming
    if varname != VARIABLE:
        daily = daily.rename({varname: VARIABLE})

    daily.attrs = dict(ds.attrs)
    prev_history = daily.attrs.get("history", "")
    new_history = "Daily means computed from hourly ERA5-Land data by resample(time='1D').mean()"
    daily.attrs["history"] = f"{prev_history} | {new_history}".strip(" |")

    daily[VARIABLE].attrs = dict(ds[varname].attrs)
    daily[VARIABLE].attrs["cell_methods"] = "time: mean within days"

    return daily

def write_monthly_daily_file(
    daily: xr.Dataset,
    outpath: Path,
    engine: Optional[str] = None,
) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        var: {"zlib": True, "complevel": 4}
        for var in daily.data_vars
    }
    print(f"[write] Monthly daily file: {outpath}")
    daily.to_netcdf(outpath, mode="w", format="NETCDF4", engine=engine, encoding=encoding)


def write_one_file_per_day(
    daily: xr.Dataset,
    outdir: Path,
    engine: Optional[str] = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for t in daily.time.values:
        one_day = daily.sel(time=[t])
        date_str = pd.Timestamp(t).strftime("%Y%m%d")
        outpath = outdir / f"era5land_swvl1_daily_{date_str}.nc"
        encoding = {
            var: {"zlib": True, "complevel": 4}
            for var in one_day.data_vars
        }
        print(f"[write] Daily file: {outpath}")
        one_day.to_netcdf(outpath, mode="w", format="NETCDF4", engine=engine, encoding=encoding)


def main() -> int:
    args = parse_args()

    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)
    outdir = Path(args.outdir)
    hourly_dir = outdir / "hourly"
    monthly_daily_dir = outdir / "daily_monthly"
    daily_files_dir = outdir / "daily_files"

    hourly_dir.mkdir(parents=True, exist_ok=True)

    try:
        client = cdsapi.Client()
    except Exception as exc:
        print(
            "Failed to create CDS API client. Check your ~/.cdsapirc and network access.\n"
            f"Error: {exc}",
            file=sys.stderr,
        )
        return 2

    for window in month_windows(start, end):
        yyyymm = f"{window.year:04d}{window.month:02d}"
        hourly_path = hourly_dir / f"era5land_swvl1_hourly_{yyyymm}.nc"

        download_month(
            client=client,
            window=window,
            hourly_path=hourly_path,
            area=args.area,
            skip_download=args.skip_download,
        )

        if not hourly_path.exists():
            print(f"Expected hourly file not found after download: {hourly_path}", file=sys.stderr)
            return 3

        print(f"[open] {hourly_path}")
        ds = xr.open_dataset(hourly_path, engine=args.engine)

        try:
            daily = compute_daily_means(ds)

            if args.write_monthly_daily:
                monthly_daily_path = monthly_daily_dir / f"era5land_swvl1_daily_{yyyymm}.nc"
                write_monthly_daily_file(daily, monthly_daily_path, engine=args.engine)

            if args.write_daily_files:
                write_one_file_per_day(daily, daily_files_dir, engine=args.engine)

        finally:
            ds.close()

        if not args.keep_hourly:
            try:
                print(f"[cleanup] Removing hourly file: {hourly_path}")
                hourly_path.unlink()
            except Exception as exc:
                print(f"Warning: could not remove {hourly_path}: {exc}", file=sys.stderr)

    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())