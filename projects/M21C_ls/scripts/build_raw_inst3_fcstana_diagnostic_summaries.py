#!/usr/bin/env python3
"""Build monthly diagnostic ANA-FCST summaries from raw inst3 files.

This script reads raw/submonthly GEOS-LDAS inst3_1d_lndfcstana_Nt files and
summarizes analysis-minus-forecast diagnostic state differences over each
calendar month. These fields are burden/activity metrics in their native state
units (m3 m-3, K), not cumulative water increments.
"""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

xr.set_options(keep_attrs=True)

FILL_THRESHOLD = 1.0e10
DEFAULT_ROOT = (
    "/discover/nobackup/projects/land_da/M21C_land_sweeper/"
    "LS_DAv8_M36_v3/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg"
)
DEFAULT_OUTPUT = (
    "/discover/nobackup/projects/land_da/geosldas-analysis/projects/M21C_ls/"
    "output/raw_increment_summaries/inst3_fcstana_raw_monthly_diagnostic_200006_202405.nc"
)


def month_starts(start: str, end: str) -> pd.DatetimeIndex:
    start_ts = pd.Timestamp(start).replace(day=1)
    end_ts = pd.Timestamp(end).replace(day=1)
    return pd.date_range(start_ts, end_ts, freq="MS")


def raw_files_for_month(root: Path, prefix: str, month: pd.Timestamp) -> list[Path]:
    month_dir = root / f"Y{month:%Y}" / f"M{month:%m}"
    if not month_dir.exists():
        return []
    _, last_day = calendar.monthrange(month.year, month.month)
    files: list[Path] = []
    for day in range(1, last_day + 1):
        stem = f"{prefix}.inst3_1d_lndfcstana_Nt.{month:%Y%m}{day:02d}"
        daily = month_dir / f"{stem}.nc4"
        if daily.exists():
            files.append(daily)
        else:
            files.extend(sorted(f for f in month_dir.glob(f"{stem}*.nc4") if ".monthly." not in f.name))
    return [f for f in files if ".monthly." not in f.name]


def open_month(files: list[Path], state_vars: list[str], engine: str) -> xr.Dataset:
    needed = []
    for state_var in state_vars:
        needed.extend([f"{state_var}_ANA", f"{state_var}_FCST"])

    def preprocess(ds: xr.Dataset) -> xr.Dataset:
        return ds[[v for v in needed if v in ds.variables]]

    ds = xr.open_mfdataset(
        [str(f) for f in files],
        combine="nested",
        concat_dim="time",
        data_vars="minimal",
        coords="minimal",
        preprocess=preprocess,
        engine=engine,
        parallel=False,
    )
    for name in needed:
        if name in ds:
            ds[name] = ds[name].where(ds[name] < FILL_THRESHOLD, other=np.float32(np.nan))
    return ds


def _read_latlon(path: Path, engine: str) -> tuple[np.ndarray, np.ndarray]:
    with xr.open_dataset(str(path), engine=engine) as ds0:
        lat = ds0["lat"].values
        lon = ds0["lon"].values
    if lat.ndim > 1:
        lat = lat[0]
    if lon.ndim > 1:
        lon = lon[0]
    return lat, lon


def summarize_month(files: list[Path], month: pd.Timestamp, state_vars: list[str], engine: str) -> xr.Dataset:
    lat, lon = _read_latlon(files[0], engine)
    with open_month(files, state_vars=state_vars, engine=engine) as ds:
        out = xr.Dataset()
        first_count = None

        for state_var in state_vars:
            ana_name = f"{state_var}_ANA"
            fcst_name = f"{state_var}_FCST"
            if ana_name not in ds or fcst_name not in ds:
                print(f"[warn] {month:%Y-%m}: missing {ana_name} or {fcst_name}; skipping {state_var}")
                continue

            diff = ds[ana_name] - ds[fcst_name]
            count = diff.notnull().sum("time")
            if first_count is None:
                first_count = count

            prefix = state_var.upper()
            out[f"{prefix}_INC_MEAN"] = diff.mean("time", skipna=True)
            out[f"{prefix}_INC_ABS_MEAN"] = abs(diff).mean("time", skipna=True)
            out[f"{prefix}_INC_RMS"] = np.sqrt((diff ** 2).mean("time", skipna=True))

            units = ds[ana_name].attrs.get("units", "")
            out[f"{prefix}_INC_MEAN"].attrs.update(
                long_name=f"monthly mean {state_var} analysis-minus-forecast diagnostic difference",
                units=units,
            )
            out[f"{prefix}_INC_ABS_MEAN"].attrs.update(
                long_name=f"monthly mean absolute {state_var} analysis-minus-forecast diagnostic difference",
                units=units,
            )
            out[f"{prefix}_INC_RMS"].attrs.update(
                long_name=f"monthly RMS {state_var} analysis-minus-forecast diagnostic difference",
                units=units,
            )

        if not out.data_vars:
            raise KeyError(f"No requested ANA/FCST variable pairs found in {files[0]}")

        if first_count is not None:
            out["time_step_count"] = first_count
            out["time_step_count"].attrs.update(
                long_name="number of finite raw inst3 ANA-FCST time steps",
                units="1",
            )
        out = out.assign_coords(lat=("tile", lat), lon=("tile", lon))
        out = out.load()

    out = out.expand_dims(time=[np.datetime64(f"{month:%Y-%m}-01")])
    for v in list(out.data_vars):
        if v != "time_step_count":
            out[v] = out[v].astype("float32")
    return out


def write_output(ds: xr.Dataset, output: Path, engine: str, complevel: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    time_chunk = min(int(ds.sizes.get("time", 1)), 12)
    tile_chunk = min(int(ds.sizes.get("tile", 1)), 100_000)
    encoding = {}
    for name, da in ds.data_vars.items():
        enc = {"zlib": True, "complevel": complevel}
        if da.dims == ("time", "tile"):
            enc["chunksizes"] = (time_chunk, tile_chunk)
        encoding[name] = enc
    ds.to_netcdf(output, engine=engine, encoding=encoding)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize raw inst3 ANA-FCST diagnostic differences into monthly metrics."
    )
    parser.add_argument("--root", type=Path, default=Path(DEFAULT_ROOT))
    parser.add_argument("--file-prefix", default="LS_DAv8_M36")
    parser.add_argument("--start", default="2000-06-01")
    parser.add_argument("--end", default="2024-05-31")
    parser.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT))
    parser.add_argument("--engine", choices=["h5netcdf", "netcdf4"], default="h5netcdf")
    parser.add_argument(
        "--state-vars",
        nargs="+",
        default=["SFMC", "RZMC", "PRMC", "TSURF", "TSOIL1"],
        help="State variable prefixes with matching *_ANA and *_FCST variables.",
    )
    parser.add_argument("--complevel", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    monthly = []
    for month in month_starts(args.start, args.end):
        files = raw_files_for_month(args.root, args.file_prefix, month)
        if not files:
            print(f"[warn] no raw inst3_1d_lndfcstana_Nt files for {month:%Y-%m}")
            continue
        print(f"{month:%Y-%m}: {len(files)} files")
        monthly.append(summarize_month(files, month, args.state_vars, args.engine))

    if not monthly:
        raise FileNotFoundError(f"No raw inst3_1d_lndfcstana_Nt files found under {args.root}")

    out = xr.concat(monthly, dim="time")
    out.attrs.update(
        source_root=str(args.root),
        file_prefix=args.file_prefix,
        state_vars=", ".join(args.state_vars),
        note=(
            "Monthly diagnostics from raw/submonthly inst3_1d_lndfcstana_Nt files. "
            "ANA-FCST differences are diagnostic state corrections, not cumulative water increments."
        ),
    )
    write_output(out, args.output, args.engine, args.complevel)
    print(f"Wrote {args.output}")
    print("Variables:", list(out.data_vars))
    print("Dims:", dict(out.sizes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
