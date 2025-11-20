"""
Shared utilities to reproduce the MATLAB obs-scaling workflow in Python.

This module mirrors the helper routines that the MATLAB scripts rely on:
parsing `ldas_obsparam` tables, reading binary ObsFcstAna files, date/time
book-keeping, and NetCDF writing on an Earth-fixed lat/lon grid.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Iterable, Iterator, List, Sequence

import math
import shlex
import struct

import netCDF4 as nc
import numpy as np


# --------------------------------------------------------------------------------------
# Data containers


@dataclass
class DateTime:
    """Container that mimics the MATLAB struct used in the driver script."""

    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: int
    dofyr: int | None = None
    pentad: int | None = None


@dataclass
class ObsParam:
    """Subset of observation-parameter metadata."""

    descr: str
    species: int
    orbit: int
    pol: int
    n_ang: int
    ang: np.ndarray
    freq: float
    fov: float
    fov_units: str
    assim: str
    scale: str
    getinnov: str
    rtm_id: int
    bias_npar: float
    bias_trel: float
    bias_tcut: float
    nodata: float
    varname: str
    units: str
    path: str
    name: str
    maskpath: str
    maskname: str
    scalepath: str
    scalename: str
    flistpath: str
    flistname: str
    errstd: float
    std_normal_max: float
    zeromean: str
    coarsen_pert: str
    xcorr: float
    ycorr: float
    adapt: float


@dataclass
class ObsFcstAnaRecord:
    """Arrays read from a single ObsFcstAna binary file."""

    date_time: DateTime
    obs_assim: np.ndarray
    obs_species: np.ndarray
    obs_tilenum: np.ndarray
    obs_lon: np.ndarray
    obs_lat: np.ndarray
    obs_obs: np.ndarray
    obs_obsvar: np.ndarray
    obs_fcst: np.ndarray
    obs_fcstvar: np.ndarray
    obs_ana: np.ndarray
    obs_anavar: np.ndarray


# --------------------------------------------------------------------------------------
# Date helpers


def is_leap_year(year: int) -> bool:
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def days_in_month(year: int, month: int) -> int:
    if month < 1 or month > 12:
        raise ValueError(f"Invalid month: {month}")
    leap = is_leap_year(year)
    days = [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return days[month - 1]


def pentad_of_year(day_of_year: int, year: int) -> int:
    if is_leap_year(year) and day_of_year >= 59:
        return math.floor((day_of_year - 2) / 5) + 1
    return math.floor((day_of_year - 1) / 5) + 1


def get_dofyr_pentad(date_time: DateTime) -> DateTime:
    dofyr = date_time.day
    for month in range(1, date_time.month):
        dofyr += days_in_month(date_time.year, month)
    pentad = pentad_of_year(dofyr, date_time.year)
    return replace(date_time, dofyr=dofyr, pentad=pentad)


def augment_date_time(dtstep: int, date_time_old: DateTime) -> DateTime:
    """
    Add/subtract seconds to a DateTime and recompute DOY/pentad.

    Follows the MATLAB implementation, updating one day at a time to avoid
    accumulated floating-point drift when the window is large.
    """

    if math.isnan(dtstep):
        raise ValueError("dtstep cannot be NaN")

    date_time = replace(date_time_old)

    if dtstep == 0:
        return get_dofyr_pentad(date_time)

    seconds_remaining = dtstep
    while seconds_remaining != 0:
        step = 86400
        if seconds_remaining > 0:
            step = min(step, seconds_remaining)
        else:
            step = max(-step, seconds_remaining)
        seconds_remaining -= step

        secs_in_day = date_time.hour * 3600 + date_time.minute * 60 + date_time.second
        secs_in_day += step

        if step > 0 and secs_in_day >= 86400:
            secs_in_day -= 86400
            date_time = _increment_calendar(date_time, forward=True)
        elif step < 0 and secs_in_day < 0:
            secs_in_day += 86400
            date_time = _increment_calendar(date_time, forward=False)

        date_time = replace(
            date_time,
            hour=int(secs_in_day // 3600),
            minute=int((secs_in_day % 3600) // 60),
            second=int(secs_in_day % 60),
        )

    return get_dofyr_pentad(date_time)


def _increment_calendar(date_time: DateTime, forward: bool) -> DateTime:
    if forward:
        if date_time.day == days_in_month(date_time.year, date_time.month):
            if date_time.month == 12:
                return replace(date_time, year=date_time.year + 1, month=1, day=1)
            return replace(date_time, month=date_time.month + 1, day=1)
        return replace(date_time, day=date_time.day + 1)

    # backward
    if date_time.day == 1:
        if date_time.month == 1:
            return replace(date_time, year=date_time.year - 1, month=12, day=31)
        prev_month = date_time.month - 1
        return replace(date_time, month=prev_month, day=days_in_month(date_time.year, prev_month))
    return replace(date_time, day=date_time.day - 1)


# --------------------------------------------------------------------------------------
# File readers


def read_obs_param(fname: str | Path) -> list[ObsParam]:
    path = Path(fname)
    if not path.exists():
        raise FileNotFoundError(f"ldas_obsparam file not found: {path}")

    with path.open() as f:
        lexer = shlex.shlex(f, posix=True)
        lexer.whitespace_split = True
        lexer.commenters = ""
        tokens = list(lexer)

    token_iter = iter(tokens)
    n_entries = int(next(token_iter))
    params: list[ObsParam] = []
    for _ in range(n_entries):
        descr = _strip_quotes(next(token_iter))
        species = int(float(next(token_iter)))
        orbit = int(float(next(token_iter)))
        pol = int(float(next(token_iter)))
        n_ang = int(float(next(token_iter)))
        ang_vals = np.array([float(next(token_iter)) for _ in range(n_ang)], dtype=float)
        freq = float(next(token_iter))
        fov = float(next(token_iter))
        fov_units = _strip_quotes(next(token_iter))
        assim = _strip_quotes(next(token_iter))
        scale = _strip_quotes(next(token_iter))
        getinnov = _strip_quotes(next(token_iter))
        rtm_id = int(float(next(token_iter)))
        bias_npar = float(next(token_iter))
        bias_trel = float(next(token_iter))
        bias_tcut = float(next(token_iter))
        nodata = float(next(token_iter))
        varname = _strip_quotes(next(token_iter))
        units = _strip_quotes(next(token_iter))
        path_str = _strip_quotes(next(token_iter))
        name = _strip_quotes(next(token_iter))
        maskpath = _strip_quotes(next(token_iter))
        maskname = _strip_quotes(next(token_iter))
        scalepath = _strip_quotes(next(token_iter))
        scalename = _strip_quotes(next(token_iter))
        flistpath = _strip_quotes(next(token_iter))
        flistname = _strip_quotes(next(token_iter))
        errstd = float(next(token_iter))
        std_normal_max = float(next(token_iter))
        zeromean = _strip_quotes(next(token_iter))
        coarsen_pert = _strip_quotes(next(token_iter))
        xcorr = float(next(token_iter))
        ycorr = float(next(token_iter))
        adapt = float(next(token_iter))

        params.append(
            ObsParam(
                descr=descr,
                species=species,
                orbit=orbit,
                pol=pol,
                n_ang=n_ang,
                ang=ang_vals,
                freq=freq,
                fov=fov,
                fov_units=fov_units,
                assim=assim,
                scale=scale,
                getinnov=getinnov,
                rtm_id=rtm_id,
                bias_npar=bias_npar,
                bias_trel=bias_trel,
                bias_tcut=bias_tcut,
                nodata=nodata,
                varname=varname,
                units=units,
                path=path_str,
                name=name,
                maskpath=maskpath,
                maskname=maskname,
                scalepath=scalepath,
                scalename=scalename,
                flistpath=flistpath,
                flistname=flistname,
                errstd=errstd,
                std_normal_max=std_normal_max,
                zeromean=zeromean,
                coarsen_pert=coarsen_pert,
                xcorr=xcorr,
                ycorr=ycorr,
                adapt=adapt,
            )
        )

    return params


def _strip_quotes(token: str) -> str:
    if len(token) >= 2 and token[0] == token[-1] and token[0] in {"'", '"'}:
        return token[1:-1]
    return token


def read_obs_fcst_ana(fname: str | Path, is_ldassa: bool | None = None) -> ObsFcstAnaRecord | None:
    path = Path(fname)
    if not path.exists():
        return None

    attempts = []
    if is_ldassa is None:
        attempts = ["<", ">"]
    else:
        attempts = [">" if is_ldassa else "<"]

    last_error: Exception | None = None
    for endian in attempts:
        try:
            return _read_obs_fcst_ana_struct(path, endian)
        except (EOFError, OSError, ValueError, struct.error) as exc:
            last_error = exc
            continue
    if last_error is not None:
        raise last_error
    return None


def _read_obs_fcst_ana_struct(path: Path, endian: str) -> ObsFcstAnaRecord:
    dtype_int = np.dtype(f"{endian}i4")
    dtype_float = np.dtype(f"{endian}f4")

    with path.open("rb") as fp:
        def read_int() -> int:
            data = fp.read(4)
            if len(data) != 4:
                raise EOFError("Unexpected EOF while reading int record")
            return struct.unpack(f"{endian}i", data)[0]

        def read_int_array(count: int) -> np.ndarray:
            data = fp.read(count * 4)
            if len(data) != count * 4:
                raise EOFError("Unexpected EOF while reading int array")
            return np.frombuffer(data, dtype=dtype_int).copy()

        def read_float_array(count: int) -> np.ndarray:
            data = fp.read(count * 4)
            if len(data) != count * 4:
                raise EOFError("Unexpected EOF while reading float array")
            return np.frombuffer(data, dtype=dtype_float).copy()

        _ = read_int()
        n_obs = read_int()
        timestamp = np.frombuffer(fp.read(32), dtype=dtype_int)
        if timestamp.size != 8:
            raise EOFError("Unexpected EOF while reading timestamp")
        year, month, day, hour, minute, second, dofyr, pentad = [int(x) for x in timestamp]
        _ = read_int()

        date_time = DateTime(year, month, day, hour, minute, second, dofyr, pentad)

        _ = read_int()
        tmp_assim = read_int_array(n_obs)
        _ = read_int()
        obs_assim = (tmp_assim != 0).astype(bool)

        _ = read_int()
        obs_species = read_int_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_tilenum = read_int_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_lon = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_lat = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_obs = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_obsvar = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_fcst = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_fcstvar = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_ana = read_float_array(n_obs)
        _ = read_int()

        _ = read_int()
        obs_anavar = read_float_array(n_obs)
        _ = read_int()

    nodata = -9999.0
    obs_obsvar = np.where(obs_obsvar == nodata, np.nan, obs_obsvar)
    obs_fcst = np.where(obs_fcst == nodata, np.nan, obs_fcst)
    obs_fcstvar = np.where(obs_fcstvar == nodata, np.nan, obs_fcstvar)
    obs_ana = np.where(obs_ana == nodata, np.nan, obs_ana)
    obs_anavar = np.where(obs_anavar == nodata, np.nan, obs_anavar)

    return ObsFcstAnaRecord(
        date_time=date_time,
        obs_assim=obs_assim,
        obs_species=obs_species,
        obs_tilenum=obs_tilenum,
        obs_lon=obs_lon,
        obs_lat=obs_lat,
        obs_obs=obs_obs,
        obs_obsvar=obs_obsvar,
        obs_fcst=obs_fcst,
        obs_fcstvar=obs_fcstvar,
        obs_ana=obs_ana,
        obs_anavar=obs_anavar,
    )


# --------------------------------------------------------------------------------------
# NetCDF writer


def write_netcdf_latlon_grid(
    fname: str | Path,
    colind: np.ndarray,
    rowind: np.ndarray,
    ll_lons: np.ndarray,
    ll_lats: np.ndarray,
    data: np.ndarray,
    pentad: Sequence[int] | int,
    start_time: Sequence[DateTime],
    end_time: Sequence[DateTime],
    overwrite: bool,
    n_out_fields: int,
    ll_lon: float,
    ll_lat: float,
    d_lon: float,
    d_lat: float,
) -> None:
    path = Path(fname)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists and overwrite is False")
    path.parent.mkdir(parents=True, exist_ok=True)

    data = np.asarray(data, dtype=float)
    if data.shape[0] != n_out_fields:
        raise ValueError("First dimension of data must equal n_out_fields")
    if data.ndim not in (2, 3):
        raise ValueError("data must have shape (Nf, Ncells) or (Nf, Ncells, Npentad)")

    pentad_vals = np.atleast_1d(pentad).astype(int)
    n_pentad = pentad_vals.size
    n_lon = len(ll_lons)
    n_lat = len(ll_lats)
    n_cells = data.shape[1]
    if len(colind) != n_cells or len(rowind) != n_cells:
        raise ValueError("colind/rowind length must match number of grid cells in data")
    has_time = data.ndim == 3

    fill_value = -999.0
    grid_data = np.full((n_out_fields, n_lat, n_lon, n_pentad), fill_value, dtype=float)
    colind = np.asarray(colind, dtype=int)
    rowind = np.asarray(rowind, dtype=int)
    if has_time:
        for cell in range(n_cells):
            grid_data[:, rowind[cell], colind[cell], :] = data[:, cell, :]
    else:
        for cell in range(n_cells):
            grid_data[:, rowind[cell], colind[cell], 0] = data[:, cell]

    def _to_days_since_1950(items: Sequence[DateTime]) -> np.ndarray:
        epoch = datetime(1950, 1, 1)
        values = []
        for item in items:
            dt_obj = datetime(item.year, item.month, item.day, item.hour, item.minute, item.second)
            values.append((dt_obj - epoch).total_seconds() / 86400.0)
        return np.array(values, dtype=float)

    start_days = _to_days_since_1950(start_time)
    end_days = _to_days_since_1950(end_time)

    grid_legacy = np.transpose(grid_data, (0, 3, 2, 1))  # fields, pentad, lon, lat

    dataset = nc.Dataset(path, "w", format="NETCDF4")
    try:
        dataset.createDimension("pentad", n_pentad)
        dataset.createDimension("lon", n_lon)
        dataset.createDimension("lat", n_lat)

        dataset.createVariable("version", "i4")[:] = 0

        def _scalar(name: str, value: float, attrs: dict[str, str]) -> None:
            var = dataset.createVariable(name, "f8")
            var[:] = value
            for key, val in attrs.items():
                var.setncattr(key, val)

        _scalar(
            "ll_lon",
            ll_lon,
            {
                "standard_name": "longitude of lower left corner",
                "long_name": "longitude of lower left corner",
                "units": "degrees_east",
                "axis": "X",
            },
        )
        _scalar(
            "ll_lat",
            ll_lat,
            {
                "standard_name": "latitude of lower left corner",
                "long_name": "latitude of lower left corner",
                "units": "degrees_north",
                "axis": "Y",
            },
        )
        _scalar(
            "d_lon",
            d_lon,
            {
                "standard_name": "longitude grid spacing",
                "long_name": "longitude grid spacing",
                "units": "degrees",
                "axis": "X",
            },
        )
        _scalar(
            "d_lat",
            d_lat,
            {
                "standard_name": "latitude grid spacing",
                "long_name": "latitude grid spacing",
                "units": "degrees",
                "axis": "Y",
            },
        )

        pentad_var = dataset.createVariable("pentad", "i4", ("pentad",))
        pentad_var[:] = pentad_vals
        pentad_var.setncattr("standard_name", "pentad")
        pentad_var.setncattr("long_name", "pentad")
        pentad_var.setncattr("units", "1")
        pentad_var.setncattr("axis", "T")

        dataset.createVariable("start_time", "f8", ("pentad",))[:] = start_days
        dataset.createVariable("end_time", "f8", ("pentad",))[:] = end_days

        def _var(name: str, dimensions: tuple[str, ...], attrs: dict[str, str]) -> nc.Variable:
            compression = dict(zlib=True, complevel=5)
            var = dataset.createVariable(name, "f8", dimensions, **compression)
            for key, val in attrs.items():
                var.setncattr(key, val)
            return var

        stat_attrs = {
            "o_mean": {
                "standard_name": "observation mean",
                "long_name": "Observation mean for pentad calculated over all years for window length",
                "units": "Degree of saturation (0-1)",
            },
            "o_std": {
                "standard_name": "observation standard deviation",
                "long_name": "Observation standard deviation for pentad calculated over all years for window length",
                "units": "Degree of saturation (0-1)",
            },
            "m_mean": {
                "standard_name": "model mean",
                "long_name": "Model mean for pentad calculated over all years for window length",
                "units": "Surface soil moisture (m^3 m^-3)",
            },
            "m_std": {
                "standard_name": "model standard deviation",
                "long_name": "Model standard deviation for pentad calculated over all years for window length",
                "units": "Surface soil moisture (m^3 m^-3)",
            },
            "n_data": {
                "standard_name": "number of data points",
                "long_name": "Number of data points for pentad calculated over all years for window length",
                "units": "1",
            },
            "m_min": {
                "standard_name": "model minimum",
                "long_name": "Model minimum calculated over all years",
                "units": "Surface soil moisture (m^3 m^-3)",
            },
            "m_max": {
                "standard_name": "model maximum",
                "long_name": "Model maximum calculated over all years",
                "units": "Surface soil moisture (m^3 m^-3)",
            },
        }

        om = _var("o_mean", ("pentad", "lon", "lat"), stat_attrs["o_mean"])
        ov = _var("o_std", ("pentad", "lon", "lat"), stat_attrs["o_std"])
        mm = _var("m_mean", ("pentad", "lon", "lat"), stat_attrs["m_mean"])
        mv = _var("m_std", ("pentad", "lon", "lat"), stat_attrs["m_std"])
        ndata_var = _var("n_data", ("pentad", "lon", "lat"), stat_attrs["n_data"])
        mmin = _var("m_min", ("lon", "lat"), stat_attrs["m_min"])
        mmax = _var("m_max", ("lon", "lat"), stat_attrs["m_max"])

        om[:, :, :] = grid_legacy[0]
        ov[:, :, :] = grid_legacy[1]
        mm[:, :, :] = grid_legacy[2]
        mv[:, :, :] = grid_legacy[3]
        ndata_var[:, :, :] = grid_legacy[4]

        with np.errstate(all="ignore"):
            mmax[:, :] = np.nanmax(grid_legacy[6], axis=0)

            min_data = np.array(grid_legacy[5], dtype=float)
            min_data[min_data < -9998] = np.nan
            min_vals = np.nanmin(min_data, axis=0)
        min_vals = np.where(np.isnan(min_vals), -9999.0, min_vals)
        mmin[:, :] = min_vals
    finally:
        dataset.close()


# --------------------------------------------------------------------------------------
# Helpers for sliding window accumulation


def nan_add(accumulator: float, values: Iterable[float]) -> float:
    vals = [v for v in values if not np.isnan(v)]
    if not vals:
        return accumulator
    if np.isnan(accumulator):
        accumulator = 0.0
    return accumulator + float(np.sum(vals))


def nan_count(accumulator: float, values: Iterable[float]) -> float:
    if np.isnan(accumulator):
        accumulator = 0.0
    return accumulator + sum(0 if np.isnan(v) else 1 for v in values)


def nan_min(accumulator: float, values: Iterable[float]) -> float:
    valid = [v for v in values if not np.isnan(v)]
    if not valid:
        return accumulator
    if np.isnan(accumulator):
        return min(valid)
    return min(accumulator, min(valid))


def nan_max(accumulator: float, values: Iterable[float]) -> float:
    valid = [v for v in values if not np.isnan(v)]
    if not valid:
        return accumulator
    if np.isnan(accumulator):
        return max(valid)
    return max(accumulator, max(valid))


__all__ = [
    "DateTime",
    "ObsParam",
    "ObsFcstAnaRecord",
    "augment_date_time",
    "days_in_month",
    "get_dofyr_pentad",
    "is_leap_year",
    "pentad_of_year",
    "read_obs_fcst_ana",
    "read_obs_param",
    "write_netcdf_latlon_grid",
    "nan_add",
    "nan_count",
    "nan_min",
    "nan_max",
]
