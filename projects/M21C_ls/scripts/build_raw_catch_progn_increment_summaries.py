#!/usr/bin/env python3
"""Build cumulative monthly summaries from raw catch_progn_incr files.

This script reads raw/submonthly GEOS-LDAS catch_progn_incr files and sums the
native prognostic increments over each calendar month. It intentionally does
not read the catch_progn_incr.monthly.YYYYMM.nc4 files because those are
time means in the M21C output and are not suitable for cumulative water-budget
interpretation.
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
    "output/raw_increment_summaries/catch_progn_raw_monthly_cumulative_200006_202405.nc"
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
        stem = f"{prefix}.catch_progn_incr.{month:%Y%m}{day:02d}"
        daily = month_dir / f"{stem}.nc4"
        if daily.exists():
            files.append(daily)
        else:
            files.extend(sorted(f for f in month_dir.glob(f"{stem}*.nc4") if ".monthly." not in f.name))
    return [f for f in files if ".monthly." not in f.name]


def open_month(files: list[Path], engine: str) -> xr.Dataset:
    needed = [
        "WESNN1_INCR",
        "WESNN2_INCR",
        "WESNN3_INCR",
        "SRFEXC_INCR",
        "RZEXC_INCR",
        "CATDEF_INCR",
    ]

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


def sum_or_nan(da: xr.DataArray) -> xr.DataArray:
    count = da.notnull().sum("time")
    out = da.fillna(0.0).sum("time")
    return out.where(count > 0)


def _read_latlon(path: Path, engine: str) -> tuple[np.ndarray, np.ndarray]:
    with xr.open_dataset(str(path), engine=engine) as ds0:
        lat = ds0["lat"].values
        lon = ds0["lon"].values
    if lat.ndim > 1:
        lat = lat[0]
    if lon.ndim > 1:
        lon = lon[0]
    return lat, lon


def summarize_month(files: list[Path], month: pd.Timestamp, engine: str, threshold: float) -> xr.Dataset:
    lat, lon = _read_latlon(files[0], engine)
    with open_month(files, engine=engine) as ds:
        missing = [v for v in ["WESNN1_INCR", "WESNN2_INCR", "WESNN3_INCR"] if v not in ds]
        if missing:
            raise KeyError(f"Missing required snow increment variables in {files[0]}: {missing}")

        snow_step = ds["WESNN1_INCR"] + ds["WESNN2_INCR"] + ds["WESNN3_INCR"]
        out = xr.Dataset()
        out["snow_net"] = sum_or_nan(snow_step)
        out["snow_abs_netpack"] = sum_or_nan(abs(snow_step))
        out["snow_abs_layers"] = sum_or_nan(
            abs(ds["WESNN1_INCR"]) + abs(ds["WESNN2_INCR"]) + abs(ds["WESNN3_INCR"])
        )
        out["snow_event_count"] = ((abs(snow_step) > threshold) & snow_step.notnull()).sum("time")
        out["time_step_count"] = snow_step.notnull().sum("time")

        # CATDEF is a deficit: positive CATDEF_INCR means drier, opposite the
        # water-content sign of SRFEXC_INCR/RZEXC_INCR. Do not sum catdef_net
        # with the other two net fields without flipping its sign.
        soil_vars = ["SRFEXC_INCR", "RZEXC_INCR", "CATDEF_INCR"]
        if all(v in ds for v in soil_vars):
            srf = ds["SRFEXC_INCR"]
            rz = ds["RZEXC_INCR"]
            catdef = ds["CATDEF_INCR"]
            out["srfexc_net"] = sum_or_nan(srf)
            out["rzexc_net"] = sum_or_nan(rz)
            out["catdef_net"] = sum_or_nan(catdef)
            out["soil_water_net_approx"] = sum_or_nan(srf + rz - catdef)
            out["soil_water_abs_activity"] = sum_or_nan(abs(srf) + abs(rz) + abs(catdef))

            out["srfexc_net"].attrs.update(
                long_name="cumulative surface-excess increment, raw catch_progn_incr",
                units="kg m-2",
            )
            out["rzexc_net"].attrs.update(
                long_name="cumulative root-zone-excess increment, raw catch_progn_incr",
                units="kg m-2",
            )
            out["catdef_net"].attrs.update(
                long_name="cumulative catchment-deficit increment, raw catch_progn_incr",
                units="kg m-2",
                note=(
                    "CATDEF is a deficit variable: positive values mean drier, "
                    "opposite the water-content sign of srfexc_net/rzexc_net. "
                    "Flip sign before adding to those fields; see soil_water_net_approx."
                ),
            )

        out = out.assign_coords(lat=("tile", lat), lon=("tile", lon))
        out = out.load()

    out = out.expand_dims(time=[np.datetime64(f"{month:%Y-%m}-01")])
    out["snow_net"].attrs.update(
        long_name="cumulative snowpack increment, raw catch_progn_incr",
        units="kg m-2",
    )
    out["snow_abs_netpack"].attrs.update(
        long_name="cumulative absolute net snowpack increment activity",
        units="kg m-2",
    )
    out["snow_abs_layers"].attrs.update(
        long_name="cumulative absolute layer-wise snow increment activity",
        units="kg m-2",
    )
    out["snow_event_count"].attrs.update(
        long_name=f"number of raw time steps with abs(net snow increment) > {threshold:g} kg m-2",
        units="1",
    )
    out["time_step_count"].attrs.update(
        long_name="number of finite raw catch_progn_incr time steps",
        units="1",
    )
    if "soil_water_net_approx" in out:
        out["soil_water_net_approx"].attrs.update(
            long_name="approximate cumulative soil-water increment, SRFEXC + RZEXC - CATDEF",
            units="kg m-2",
            note=(
                "CATDEF is a deficit variable, so its sign is flipped for water-content "
                "interpretation. Treat as approximate because Catchment prognostics are "
                "not independent additive reservoirs."
            ),
        )
        out["soil_water_abs_activity"].attrs.update(
            long_name="cumulative absolute native soil prognostic increment activity",
            units="kg m-2",
        )
    skip = {"snow_event_count", "time_step_count"}
    for v in list(out.data_vars):
        if v not in skip:
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
        description="Sum raw catch_progn_incr files into monthly cumulative increment fields."
    )
    parser.add_argument("--root", type=Path, default=Path(DEFAULT_ROOT))
    parser.add_argument("--file-prefix", default="LS_DAv8_M36")
    parser.add_argument("--start", default="2000-06-01")
    parser.add_argument("--end", default="2024-05-31")
    parser.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT))
    parser.add_argument("--engine", choices=["h5netcdf", "netcdf4"], default="h5netcdf")
    parser.add_argument("--snow-event-threshold", type=float, default=0.001)
    parser.add_argument("--complevel", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    monthly = []
    for month in month_starts(args.start, args.end):
        files = raw_files_for_month(args.root, args.file_prefix, month)
        if not files:
            print(f"[warn] no raw catch_progn_incr files for {month:%Y-%m}")
            continue
        print(f"{month:%Y-%m}: {len(files)} files")
        monthly.append(summarize_month(files, month, args.engine, args.snow_event_threshold))

    if not monthly:
        raise FileNotFoundError(f"No raw catch_progn_incr files found under {args.root}")

    out = xr.concat(monthly, dim="time")
    out.attrs.update(
        source_root=str(args.root),
        file_prefix=args.file_prefix,
        note=(
            "Monthly cumulative sums from raw/submonthly catch_progn_incr files. "
            "Do not confuse with catch_progn_incr.monthly files, which are time means."
        ),
    )
    write_output(out, args.output, args.engine, args.complevel)
    print(f"Wrote {args.output}")
    print("Variables:", list(out.data_vars))
    print("Dims:", dict(out.sizes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
