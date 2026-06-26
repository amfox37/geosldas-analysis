#!/usr/bin/env python3
"""Aggregate monthly ObsFcstAna sums into temporal_stats NetCDF files.

This is intended for the M21C Land Sweeper paper figures. It reads monthly
``*.ens_avg.ldas_ObsFcstAna.YYYYMM*_stats*.nc4`` files and writes the
``temporal_stats_{OL,DA}_YYYYMMDD_YYYYMMDD.nc4`` products expected by the
paper plotting notebooks.
"""

from __future__ import annotations

import argparse
import calendar
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


VAR_LIST = [
    "obs_obs",
    "obs_obsvar",
    "obs_fcst",
    "obs_fcstvar",
    "obs_ana",
    "obs_anavar",
]

OUTPUT_VARS = [
    "O_mean",
    "O_stdv",
    "F_mean",
    "F_stdv",
    "A_mean",
    "A_stdv",
    "OmF_mean",
    "OmF_stdv",
    "OmA_mean",
    "OmA_stdv",
    "OmF_norm_mean",
    "OmF_norm_stdv",
    "N_data",
]


@dataclass(frozen=True)
class Period:
    start: date
    end: date

    @property
    def start_ymd(self) -> str:
        return self.start.strftime("%Y%m%d")

    @property
    def end_ymd(self) -> str:
        return self.end.strftime("%Y%m%d")


def parse_date(value: str) -> date:
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            pass
    raise argparse.ArgumentTypeError(
        f"Expected date as YYYY-MM-DD or YYYYMMDD, got {value!r}"
    )


def iter_month_starts(period: Period):
    current = date(period.start.year, period.start.month, 1)
    last = date(period.end.year, period.end.month, 1)
    while current <= last:
        yield current
        if current.month == 12:
            current = date(current.year + 1, 1, 1)
        else:
            current = date(current.year, current.month + 1, 1)


def validate_full_month_period(period: Period) -> None:
    if period.start.day != 1:
        raise ValueError(f"Period start must be first of month: {period.start}")

    last_day = calendar.monthrange(period.end.year, period.end.month)[1]
    if period.end.day != last_day:
        raise ValueError(f"Period end must be last of month: {period.end}")

    if period.end < period.start:
        raise ValueError(f"Period end precedes start: {period.start} to {period.end}")


def monthly_suffix_candidates(exp_tag: str, suffix: str | None) -> list[str]:
    if suffix:
        return [suffix]

    if exp_tag.upper() == "OL":
        return ["_stats_CROSSMASKED.nc4", "_stats.nc4"]
    return ["_stats.nc4", "_stats_CROSSMASKED.nc4"]


def find_monthly_file(
    monthly_root: Path, expid: str, month: date, suffixes: list[str]
) -> Path:
    stem = f"{expid}.ens_avg.ldas_ObsFcstAna.{month:%Y%m}"
    month_dir = monthly_root / f"Y{month:%Y}" / f"M{month:%m}"

    for suffix in suffixes:
        candidate = month_dir / f"{stem}{suffix}"
        if candidate.exists():
            return candidate

    tried = ", ".join(str(month_dir / f"{stem}{suffix}") for suffix in suffixes)
    raise FileNotFoundError(f"Missing monthly stats file for {month:%Y-%m}: {tried}")


def as_float_array(value, fill_value=np.nan) -> np.ndarray:
    if hasattr(value, "filled"):
        value = value.filled(fill_value)
    return np.asarray(value, dtype=np.float64)


def read_monthly_sums(path: Path) -> tuple[
    np.ndarray,
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    data_sum = {}
    data2_sum = {}
    with Dataset(path, "r") as nc:
        n_data = as_float_array(nc.variables["N_data"][:], fill_value=0)
        oxf_sum = as_float_array(nc.variables["obsxfcst_sum"][:])
        oxa_sum = as_float_array(nc.variables["obsxana_sum"][:])
        fxa_sum = as_float_array(nc.variables["fcstxana_sum"][:])
        for var in VAR_LIST:
            data_sum[var] = as_float_array(nc.variables[f"{var}_sum"][:])
            data2_sum[var] = as_float_array(nc.variables[f"{var}2_sum"][:])
    return n_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum


def aggregate_monthlies(
    monthly_root: Path, expid: str, exp_tag: str, period: Period, suffix: str | None
) -> tuple[dict[str, np.ndarray], list[Path]]:
    suffixes = monthly_suffix_candidates(exp_tag, suffix)
    monthly_files = []

    n_data = None
    oxf_sum = None
    oxa_sum = None
    fxa_sum = None
    data_sum = {}
    data2_sum = {}

    for month in iter_month_starts(period):
        path = find_monthly_file(monthly_root, expid, month, suffixes)
        monthly_files.append(path)
        mn_data, mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum = read_monthly_sums(path)

        if n_data is None:
            n_data = np.zeros_like(mn_data, dtype=np.float64)
            oxf_sum = np.zeros_like(moxf_sum, dtype=np.float64)
            oxa_sum = np.zeros_like(moxa_sum, dtype=np.float64)
            fxa_sum = np.zeros_like(mfxa_sum, dtype=np.float64)
            data_sum = {var: np.zeros_like(mdata_sum[var], dtype=np.float64) for var in VAR_LIST}
            data2_sum = {var: np.zeros_like(mdata2_sum[var], dtype=np.float64) for var in VAR_LIST}

        n_data += mn_data
        oxf_sum += moxf_sum
        oxa_sum += moxa_sum
        fxa_sum += mfxa_sum
        for var in VAR_LIST:
            data_sum[var] += mdata_sum[var]
            data2_sum[var] += mdata2_sum[var]

    if n_data is None or oxf_sum is None or oxa_sum is None or fxa_sum is None:
        raise RuntimeError("No monthly files were read")

    return compute_temporal_stats(n_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum), monthly_files


def safe_divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    out = np.full(numerator.shape, np.nan, dtype=np.float64)
    np.divide(numerator, denominator, out=out, where=denominator != 0)
    return out


def compute_temporal_stats(
    n_data: np.ndarray,
    data_sum: dict[str, np.ndarray],
    data2_sum: dict[str, np.ndarray],
    oxf_sum: np.ndarray,
    oxa_sum: np.ndarray,
    fxa_sum: np.ndarray,
) -> dict[str, np.ndarray]:
    data_mean = {}
    data2_mean = {}
    data_var = {}

    for var in VAR_LIST:
        data_mean[var] = safe_divide(data_sum[var], n_data)
        data2_mean[var] = safe_divide(data2_sum[var], n_data)
        data_var[var] = data2_mean[var] - data_mean[var] ** 2

    oxf_mean = safe_divide(oxf_sum, n_data)
    oxa_mean = safe_divide(oxa_sum, n_data)

    omf_mean = data_mean["obs_obs"] - data_mean["obs_fcst"]
    oma_mean = data_mean["obs_obs"] - data_mean["obs_ana"]

    omf_var = (
        data_var["obs_obs"]
        + data_var["obs_fcst"]
        - 2.0 * (oxf_mean - data_mean["obs_obs"] * data_mean["obs_fcst"])
    )
    oma_var = (
        data_var["obs_obs"]
        + data_var["obs_ana"]
        - 2.0 * (oxa_mean - data_mean["obs_obs"] * data_mean["obs_ana"])
    )

    norm_denom = data_mean["obs_obsvar"] + data_mean["obs_fcstvar"]

    with np.errstate(invalid="ignore", divide="ignore"):
        stats = {
            "O_mean": data_mean["obs_obs"],
            "O_stdv": np.sqrt(data_var["obs_obs"]),
            "F_mean": data_mean["obs_fcst"],
            "F_stdv": np.sqrt(data_var["obs_fcst"]),
            "A_mean": data_mean["obs_ana"],
            "A_stdv": np.sqrt(data_var["obs_ana"]),
            "OmF_mean": omf_mean,
            "OmF_stdv": np.sqrt(omf_var),
            "OmA_mean": oma_mean,
            "OmA_stdv": np.sqrt(oma_var),
            "OmF_norm_mean": omf_mean / np.sqrt(norm_denom),
            "OmF_norm_stdv": np.sqrt(omf_var / norm_denom),
            "N_data": n_data,
        }

    return stats


def write_temporal_stats(path: Path, stats: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sample = stats["N_data"]

    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("tile", sample.shape[0])
        nc.createDimension("species", sample.shape[1])

        for name in OUTPUT_VARS:
            var = nc.createVariable(name, "f4", ("tile", "species"), zlib=True, complevel=4)
            var[:, :] = stats[name].astype(np.float32)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build M21C temporal_stats files from monthly ObsFcstAna sum files."
    )
    parser.add_argument("--monthly-root", required=True, type=Path)
    parser.add_argument("--expid", required=True, help="Experiment ID, e.g. LS_DAv8_M36")
    parser.add_argument("--exp-tag", required=True, choices=["OL", "DA"])
    parser.add_argument("--start", required=True, type=parse_date)
    parser.add_argument("--end", required=True, type=parse_date)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--monthly-suffix",
        help="Force a monthly suffix, e.g. _stats.nc4 or _stats_CROSSMASKED.nc4",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing temporal_stats output file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    period = Period(args.start, args.end)
    validate_full_month_period(period)

    output = args.output_dir / f"temporal_stats_{args.exp_tag}_{period.start_ymd}_{period.end_ymd}.nc4"
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"{output} exists; pass --overwrite to replace it")

    stats, monthly_files = aggregate_monthlies(
        args.monthly_root, args.expid, args.exp_tag, period, args.monthly_suffix
    )
    write_temporal_stats(output, stats)

    print(f"Wrote {output}")
    print(f"Read {len(monthly_files)} monthly files from {monthly_files[0]} to {monthly_files[-1]}")


if __name__ == "__main__":
    main()
