#!/usr/bin/env python3
"""
Download monthly ERA5-Land fields from the Copernicus Climate Data Store (CDS),
unzip CDS NetCDF packages if needed, and merge the yearly files into one
NetCDF.

The package ``cdsapi`` is the official Python download client for CDS. It is
not a science-processing library; it is only used here to submit ERA5-Land
requests to Copernicus and save the returned files. To use it, you need a CDS
account, you need to accept the ERA5-Land terms of use on the CDS website, and
you need a personal CDS API token in ``~/.cdsapirc``.

This script writes yearly ERA5-Land monthly files and one merged analysis file.

By default it downloads the monthly variables used in the GEOS-LDAS comparison:

    snow_depth                         -> sde
    snow_depth_water_equivalent        -> sd
    snow_cover                         -> snowc
    soil_temperature_level_1           -> stl1
    volumetric_soil_water_layer_1      -> swvl1
    volumetric_soil_water_layer_2      -> swvl2
    volumetric_soil_water_layer_3      -> swvl3

Requirements
------------
Install the Python packages:

    pip install "cdsapi>=0.7.7" xarray netCDF4 numpy

CDS API setup
-------------
1. Create/login to a Copernicus Climate Data Store account.
2. Accept the Terms of Use for ERA5-Land monthly means.
3. Install ``cdsapi`` with the command above.
4. Put your API token in ~/.cdsapirc, for example:

       url: https://cds.climate.copernicus.eu/api
       key: <PERSONAL-ACCESS-TOKEN>

Example
-------
Download global monthly ERA5-Land fields for 2000-2025 and merge them:

    python download_era5land_monthly.py \
        --start-year 2000 \
        --end-year 2025 \
        --outdir ./era5land_monthly \
        --merged-name ERA5_land_monthly_subset.nc

Optional spatial subsetting
---------------------------
Use CDS area order: north,west,south,east

    --area "60,-130,20,-60"

Notes
-----
- The CDS request uses dataset "reanalysis-era5-land-monthly-means".
- ERA5-Land monthly means use product_type "monthly_averaged_reanalysis".
- CDS may return a ZIP archive even when the target filename ends in ".nc".
  This script detects ZIP files, extracts the inner NetCDF, then merges the
  unzipped yearly files.
- Variables are requested one at a time and then merged. This is slower than
  one large CDS request, but it avoids ERA5-Land/CDS quirks where mixed
  variable requests can silently return only part of the requested list.
- Longitudes are normalized to 0-360 before merging, because different
  ERA5-Land variables can be returned with different longitude conventions for
  the same requested area.
- Existing yearly downloads and unzipped files are reused unless --overwrite is
  supplied.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import xarray as xr

try:
    import cdsapi
except ImportError as exc:
    raise SystemExit(
        "cdsapi is not installed. Install with: pip install 'cdsapi>=0.7.7'"
    ) from exc


DATASET = "reanalysis-era5-land-monthly-means"

DEFAULT_VARIABLES = [
    "snow_depth",
    "snow_depth_water_equivalent",
    "snow_cover",
    "soil_temperature_level_1",
    "volumetric_soil_water_layer_1",
    "volumetric_soil_water_layer_2",
    "volumetric_soil_water_layer_3",
]

MONTHS = [f"{month:02d}" for month in range(1, 13)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and merge monthly ERA5-Land fields from CDS."
    )
    parser.add_argument(
        "--start-year",
        type=int,
        required=True,
        help="First year to download, e.g. 2000.",
    )
    parser.add_argument(
        "--end-year",
        type=int,
        required=True,
        help="Last year to download, inclusive, e.g. 2025.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory for yearly downloads, unzipped files, and merged output.",
    )
    parser.add_argument(
        "--merged-name",
        default="ERA5_land_monthly_subset.nc",
        help="Filename for the merged all-years NetCDF.",
    )
    parser.add_argument(
        "--area",
        default=None,
        help=(
            "Optional spatial subset in CDS order: 'north,west,south,east'. "
            "Example: '60,-130,20,-60'. Omit for global."
        ),
    )
    parser.add_argument(
        "--variables",
        default=",".join(DEFAULT_VARIABLES),
        help=(
            "Comma-separated CDS variable names. Defaults to the snow, soil "
            "temperature, and soil moisture variables used by the GEOS-LDAS "
            "ERA5-Land comparison."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload/recreate yearly and merged files even if they already exist.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=5,
        help="Maximum CDS retry attempts per year.",
    )
    parser.add_argument(
        "--compression",
        type=int,
        default=4,
        help="zlib compression level for the merged NetCDF.",
    )
    return parser.parse_args()


def parse_area(area_text: str | None) -> list[float] | None:
    if area_text is None:
        return None
    parts = [p.strip() for p in area_text.split(",")]
    if len(parts) != 4:
        raise ValueError("--area must contain four comma-separated values")
    return [float(p) for p in parts]


def magic_bytes(path: Path, nbytes: int = 8) -> bytes:
    with path.open("rb") as stream:
        return stream.read(nbytes)


def looks_like_zip(path: Path) -> bool:
    return magic_bytes(path, 4) == b"PK\x03\x04"


def looks_like_netcdf(path: Path) -> bool:
    header = magic_bytes(path, 8)
    return header.startswith(b"CDF") or header.startswith(b"\x89HDF")


def variable_slug(variable: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", variable).strip("_")


def download_year_variable(
    client: cdsapi.Client,
    year: int,
    variable: str,
    outdir: Path,
    area: list[float] | None,
    overwrite: bool,
    max_attempts: int,
) -> Path:
    target = outdir / f"era5l_monthly_{year}_{variable_slug(variable)}.nc"
    if target.exists() and not overwrite:
        print(f"[skip] {target} exists")
        return target

    request = {
        "product_type": "monthly_averaged_reanalysis",
        "variable": [variable],
        "year": f"{year}",
        "month": MONTHS,
        "time": "00:00",
        "format": "netcdf",
    }
    if area is not None:
        request["area"] = area

    for attempt in range(1, max_attempts + 1):
        try:
            print(f"[CDS] requesting {year} {variable} -> {target}")
            client.retrieve(DATASET, request, str(target))
            print(f"[ok] wrote {target}")
            return target
        except Exception as exc:
            if attempt == max_attempts:
                raise
            wait_seconds = 30 * attempt
            print(
                f"[warn] {year} {variable} request failed: {exc}\n"
                f"       retrying in {wait_seconds}s "
                f"({attempt + 1}/{max_attempts})"
            )
            time.sleep(wait_seconds)

    raise RuntimeError(f"failed to download {year} {variable}")


def extract_or_copy_netcdf(path: Path, overwrite: bool) -> Path:
    out = path.with_name(path.stem + "_unzipped.nc")
    if out.exists() and not overwrite:
        print(f"[skip] {out} exists")
        return out

    if looks_like_zip(path):
        print(f"[unzip] {path.name} -> {out.name}")
        with zipfile.ZipFile(path, "r") as archive:
            members = [name for name in archive.namelist() if name.lower().endswith(".nc")]
            if not members:
                raise RuntimeError(f"No NetCDF member found inside {path}")
            with archive.open(members[0]) as source, out.open("wb") as target:
                shutil.copyfileobj(source, target)
    elif looks_like_netcdf(path):
        print(f"[copy] {path.name} is already NetCDF -> {out.name}")
        shutil.copy2(path, out)
    else:
        raise RuntimeError(
            f"{path} is neither ZIP nor NetCDF; first bytes are {magic_bytes(path, 16)!r}"
        )

    if not looks_like_netcdf(out):
        raise RuntimeError(f"{out} was created but is not a readable NetCDF/HDF5 file")
    return out


def time_name(ds: xr.Dataset) -> str:
    for name in ("valid_time", "time"):
        if name in ds.coords or name in ds.dims:
            return name
    raise RuntimeError("Dataset has no recognizable time coordinate: expected valid_time or time")


def normalize_longitude(ds: xr.Dataset) -> xr.Dataset:
    if "longitude" not in ds.coords:
        return ds
    lon = ds["longitude"]
    if float(lon.min()) < 0:
        ds = ds.assign_coords(longitude=(lon % 360))
    return ds.sortby("longitude")


def merge_variable_files_for_year(
    paths: list[Path],
    year: int,
    outdir: Path,
    compression: int,
    overwrite: bool,
) -> Path:
    out_file = outdir / f"era5l_monthly_{year}_combined.nc"
    if out_file.exists() and not overwrite:
        print(f"[skip] yearly merged file exists: {out_file}")
        return out_file

    print(f"[merge-year] {year}: opening {len(paths)} variable files")
    datasets = []
    try:
        for path in sorted(paths):
            ds_var = xr.open_dataset(path, chunks=None)
            ds_var = normalize_longitude(ds_var)
            datasets.append(ds_var.sortby(time_name(ds_var)))

        merged = xr.merge(datasets, compat="override", join="exact")
        tname = time_name(merged)
        merged = merged.sortby(tname)

        encoding = {
            variable: {"zlib": True, "complevel": compression}
            for variable in merged.data_vars
        }
        for coord in (tname, "latitude", "longitude"):
            if coord in merged.coords:
                encoding[coord] = {"zlib": True, "complevel": compression}

        print(f"[write] {out_file}")
        merged.to_netcdf(out_file, encoding=encoding)
        return out_file
    finally:
        for ds in datasets:
            ds.close()


def merge_yearly_files(paths: list[Path], out_file: Path, compression: int, overwrite: bool) -> None:
    if out_file.exists() and not overwrite:
        print(f"[skip] merged file exists: {out_file}")
        return

    print(f"[merge] opening {len(paths)} yearly NetCDF files")
    datasets = []
    try:
        datasets = []
        for path in sorted(paths):
            ds_year = xr.open_dataset(path, chunks=None)
            tname = time_name(ds_year)
            datasets.append(ds_year.sortby(tname))

        tname = time_name(datasets[0])
        ds = xr.concat(
            datasets,
            dim=tname,
            data_vars="minimal",
            coords="minimal",
            compat="override",
        )
        ds = ds.sortby(tname)

        _, unique_indices = np.unique(ds[tname].values, return_index=True)
        if len(unique_indices) < ds.sizes[tname]:
            dropped = ds.sizes[tname] - len(unique_indices)
            print(f"[info] dropping {dropped} duplicate time records")
            ds = ds.isel({tname: np.sort(unique_indices)})

        encoding = {
            variable: {"zlib": True, "complevel": compression}
            for variable in ds.data_vars
        }
        for coord in (tname, "latitude", "longitude"):
            if coord in ds.coords:
                encoding[coord] = {"zlib": True, "complevel": compression}

        out_file.parent.mkdir(parents=True, exist_ok=True)
        print(f"[write] {out_file}")
        ds.to_netcdf(out_file, encoding=encoding)
        print("[done]")
    finally:
        for ds in datasets:
            ds.close()


def main() -> None:
    args = parse_args()
    if args.start_year > args.end_year:
        raise SystemExit("--start-year must be <= --end-year")

    variables = [v.strip() for v in args.variables.split(",") if v.strip()]
    if not variables:
        raise SystemExit("At least one variable is required")

    area = parse_area(args.area)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[info] dataset: {DATASET}")
    print(f"[info] years: {args.start_year}-{args.end_year}")
    print(f"[info] variables: {', '.join(variables)}")
    print(f"[info] outdir: {args.outdir}")
    if area is not None:
        print(f"[info] area: {area}")

    client = cdsapi.Client(timeout=600)
    yearly_combined = []
    for year in range(args.start_year, args.end_year + 1):
        downloaded = [
            download_year_variable(
                client=client,
                year=year,
                variable=variable,
                outdir=args.outdir,
                area=area,
                overwrite=args.overwrite,
                max_attempts=args.max_attempts,
            )
            for variable in variables
        ]
        unzipped = [extract_or_copy_netcdf(path, args.overwrite) for path in downloaded]
        yearly_combined.append(
            merge_variable_files_for_year(
                paths=unzipped,
                year=year,
                outdir=args.outdir,
                compression=args.compression,
                overwrite=args.overwrite,
            )
        )

    merge_yearly_files(
        paths=yearly_combined,
        out_file=args.outdir / args.merged_name,
        compression=args.compression,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[stopped]", file=sys.stderr)
        raise SystemExit(130)
