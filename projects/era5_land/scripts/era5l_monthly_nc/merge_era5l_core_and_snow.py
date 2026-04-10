#!/usr/bin/env python3
"""
Merge ERA5-Land monthly yearly files from two streams:
  1) core variables (e.g., stl1, snowc, swvl1-3)
  2) snow variables (sd, sde)

Produces:
  - optional yearly combined files
  - one final merged NetCDF across all years

This avoids in-place append workflows that can silently corrupt large files.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import xarray as xr


def log(msg: str) -> None:
    print(msg, flush=True)


def year_from_name(path: Path) -> str | None:
    m = re.search(r"(\d{4})", path.name)
    return m.group(1) if m else None


def discover_year_files(directory: Path, suffix: str) -> dict[str, Path]:
    patt = str(directory / f"era5l_monthly_*{suffix}")
    out: dict[str, Path] = {}
    for p in sorted(glob.glob(patt)):
        pp = Path(p)
        y = year_from_name(pp)
        if y:
            out[y] = pp
    return out


def assert_coord_equal(a: xr.Dataset, b: xr.Dataset, coord: str) -> None:
    if coord in a.coords and coord in b.coords:
        if not np.array_equal(a[coord].values, b[coord].values):
            raise ValueError(f"Coordinate mismatch for '{coord}'")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge ERA5-Land core + snow yearly files into one complete dataset")
    p.add_argument("--core-dir", type=Path, required=True, help="Directory containing yearly core files")
    p.add_argument("--snow-dir", type=Path, required=True, help="Directory containing yearly snow files")
    p.add_argument(
        "--suffix",
        type=str,
        default="_unzipped.nc",
        help="Yearly filename suffix (default: _unzipped.nc; example names: era5l_monthly_2000_unzipped.nc)",
    )
    p.add_argument(
        "--core-vars",
        type=str,
        default="stl1,snowc,swvl1,swvl2,swvl3",
        help="Comma-separated core variables to keep",
    )
    p.add_argument(
        "--snow-vars",
        type=str,
        default="sd,sde",
        help="Comma-separated snow variables to keep",
    )
    p.add_argument(
        "--out-yearly-dir",
        type=Path,
        default=None,
        help="Optional output directory for yearly combined files",
    )
    p.add_argument(
        "--out-file",
        type=Path,
        required=True,
        help="Output NetCDF path for all-years merged dataset",
    )
    p.add_argument("--compression", type=int, default=4, help="zlib compression level")
    p.add_argument(
        "--engine",
        type=str,
        default=None,
        help="xarray engine for reading/writing (optional)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    core_vars = [v.strip() for v in args.core_vars.split(",") if v.strip()]
    snow_vars = [v.strip() for v in args.snow_vars.split(",") if v.strip()]

    core_files = discover_year_files(args.core_dir, args.suffix)
    snow_files = discover_year_files(args.snow_dir, args.suffix)

    if not core_files:
        raise FileNotFoundError(f"No core files found in {args.core_dir} with suffix '{args.suffix}'")
    if not snow_files:
        raise FileNotFoundError(f"No snow files found in {args.snow_dir} with suffix '{args.suffix}'")

    years = sorted(set(core_files).intersection(snow_files))
    if not years:
        raise RuntimeError("No overlapping years between core and snow files")

    log(f"[info] overlapping years: {years[0]}..{years[-1]} ({len(years)} years)")

    if args.out_yearly_dir is not None:
        args.out_yearly_dir.mkdir(parents=True, exist_ok=True)

    yearly_outputs: list[Path] = []

    for y in years:
        cpath = core_files[y]
        spath = snow_files[y]
        log(f"[merge-year] {y}:\n  core={cpath}\n  snow={spath}")

        cds = xr.open_dataset(cpath, engine=args.engine)
        sds = xr.open_dataset(spath, engine=args.engine)

        missing_core = [v for v in core_vars if v not in cds.data_vars]
        missing_snow = [v for v in snow_vars if v not in sds.data_vars]
        if missing_core:
            cds.close(); sds.close()
            raise KeyError(f"{cpath.name}: missing core vars {missing_core}")
        if missing_snow:
            cds.close(); sds.close()
            raise KeyError(f"{spath.name}: missing snow vars {missing_snow}")

        # Strict coordinate compatibility checks.
        for coord in ("valid_time", "time", "latitude", "lat", "longitude", "lon"):
            if coord in cds.coords or coord in sds.coords:
                try:
                    assert_coord_equal(cds, sds, coord)
                except ValueError:
                    if coord in ("valid_time", "time"):
                        cv = cds[coord].values if coord in cds.coords else np.array([])
                        sv = sds[coord].values if coord in sds.coords else np.array([])
                        msg = (
                            f"{y}: Coordinate mismatch for '{coord}' "
                            f"(core={len(cv)} timestamps, snow={len(sv)} timestamps)"
                        )
                        cds.close(); sds.close()
                        raise ValueError(msg)
                    cds.close(); sds.close()
                    raise

        merged = xr.merge([cds[core_vars], sds[snow_vars]], compat="override", join="exact")

        # Normalize time axis ordering and remove accidental duplicates.
        tname = "valid_time" if "valid_time" in merged.coords else ("time" if "time" in merged.coords else None)
        if tname is not None:
            merged = merged.sortby(tname)
            _, uniq_idx = np.unique(merged[tname].values, return_index=True)
            if len(uniq_idx) < merged.sizes[tname]:
                dropped = merged.sizes[tname] - len(uniq_idx)
                log(f"  [info] dropping {dropped} duplicate timesteps in {y}")
                merged = merged.isel({tname: np.sort(uniq_idx)})

        if args.out_yearly_dir is not None:
            yout = args.out_yearly_dir / f"era5l_monthly_{y}_combined.nc"
            enc = {v: {"zlib": True, "complevel": args.compression} for v in merged.data_vars}
            for c in ("latitude", "longitude", "valid_time", "time"):
                if c in merged.coords:
                    enc[c] = {"zlib": True, "complevel": args.compression}
            merged.to_netcdf(yout, encoding=enc)
            yearly_outputs.append(yout)
            log(f"  [write] {yout}")

        merged.close()
        cds.close()
        sds.close()

    if args.out_yearly_dir is None:
        # Build merged-all directly from source yearly files to avoid extra disk use.
        combine_paths = []
        for y in years:
            combine_paths.append(str(core_files[y]))
            combine_paths.append(str(snow_files[y]))
        # We cannot directly open alternating core/snow files and expect complete vars;
        # instead generate an in-memory list of yearly merged datasets lazily.
        dsets = []
        for y in years:
            cds = xr.open_dataset(core_files[y], engine=args.engine)[core_vars]
            sds = xr.open_dataset(snow_files[y], engine=args.engine)[snow_vars]
            dsets.append(xr.merge([cds, sds], compat="override", join="exact"))
        ds_all = xr.concat(dsets, dim="valid_time")
        for d in dsets:
            d.close()
    else:
        yearly_paths = sorted(map(str, yearly_outputs))
        log(f"[merge-all] concatenating {len(yearly_paths)} yearly combined files")
        ds_all = xr.open_mfdataset(
            yearly_paths,
            combine="by_coords",
            coords="minimal",
            compat="override",
            join="exact",
            engine=args.engine,
        )

    tname = "valid_time" if "valid_time" in ds_all.coords else ("time" if "time" in ds_all.coords else None)
    if tname is None:
        raise RuntimeError("Merged dataset has no time coordinate ('valid_time' or 'time')")

    ds_all = ds_all.sortby(tname)
    _, uniq_idx = np.unique(ds_all[tname].values, return_index=True)
    if len(uniq_idx) < ds_all.sizes[tname]:
        dropped = ds_all.sizes[tname] - len(uniq_idx)
        log(f"[info] dropping {dropped} duplicate timesteps in merged-all")
        ds_all = ds_all.isel({tname: np.sort(uniq_idx)})

    args.out_file.parent.mkdir(parents=True, exist_ok=True)
    enc = {v: {"zlib": True, "complevel": args.compression} for v in ds_all.data_vars}
    for c in ("latitude", "longitude", "valid_time", "time"):
        if c in ds_all.coords:
            enc[c] = {"zlib": True, "complevel": args.compression}

    log(f"[write] final merged file: {args.out_file}")
    ds_all.to_netcdf(args.out_file, encoding=enc)
    ds_all.close()

    log("[done]")


if __name__ == "__main__":
    main()
