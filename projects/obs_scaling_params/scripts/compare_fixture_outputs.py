"""Compare Python and MATLAB obs-scaling fixture NetCDF outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np


EXP_RUN = "LS_DAv8_M36_as_test2"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
DEFAULT_VARIABLES = ("o_mean", "o_std", "m_mean", "m_std", "n_data", "m_min", "m_max")
DEFAULT_COORDS = ("ll_lon", "ll_lat", "d_lon", "d_lat", "pentad", "start_time", "end_time")
MATCH_TOKEN = "zscore_stats_"


@dataclass
class VariableDiff:
    name: str
    shape_python: tuple[int, ...]
    shape_matlab: tuple[int, ...]
    finite_python: int
    finite_matlab: int
    finite_mask_mismatch: int
    compared_points: int
    max_abs_diff: float
    mean_abs_diff: float
    status: str


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    python_dir, matlab_dir = resolve_dirs(args)
    variables = parse_csv(args.variables) or list(DEFAULT_VARIABLES)
    coords = parse_csv(args.coords) if args.compare_coords else []

    print("Obs-scaling fixture comparison")
    print(f"  python_dir: {python_dir}")
    print(f"  matlab_dir: {matlab_dir}")
    print(f"  variables:  {', '.join(variables)}")
    if coords:
        print(f"  coords:     {', '.join(coords)}")
    print(f"  tolerance:  atol={args.atol:g}, rtol={args.rtol:g}")

    python_files = index_outputs(python_dir)
    matlab_files = index_outputs(matlab_dir)
    if not python_files:
        raise FileNotFoundError(f"No NetCDF outputs found in Python directory: {python_dir}")
    if not matlab_files:
        raise FileNotFoundError(f"No NetCDF outputs found in MATLAB directory: {matlab_dir}")

    python_keys = set(python_files)
    matlab_keys = set(matlab_files)
    common_keys = sorted(python_keys & matlab_keys)
    if args.key_contains:
        common_keys = [key for key in common_keys if args.key_contains in key]
    only_python = sorted(python_keys - matlab_keys)
    only_matlab = sorted(matlab_keys - python_keys)

    print(f"  matched files: {len(common_keys)}")
    if only_python:
        print(f"  only Python:   {len(only_python)}")
        for key in only_python[: args.max_missing_list]:
            print(f"    {python_files[key].name}")
    if only_matlab:
        print(f"  only MATLAB:   {len(only_matlab)}")
        for key in only_matlab[: args.max_missing_list]:
            print(f"    {matlab_files[key].name}")

    if args.require_same_files and (only_python or only_matlab):
        print("ERROR: output file sets differ")
        return 2
    if not common_keys:
        print("ERROR: no matching output files to compare")
        return 2

    overall_failed = False
    for key in common_keys:
        python_path = python_files[key]
        matlab_path = matlab_files[key]
        print()
        print(f"File key: {key}")
        print(f"  Python: {python_path.name}")
        print(f"  MATLAB: {matlab_path.name}")
        diffs = compare_file(
            python_path,
            matlab_path,
            variables=variables,
            coords=coords,
            atol=args.atol,
            rtol=args.rtol,
        )
        print_table(diffs)
        if any(diff.status != "OK" for diff in diffs):
            overall_failed = True

    if overall_failed:
        print()
        print("Comparison failed tolerance or mask checks.")
        return 1

    print()
    print("Comparison passed.")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture-root",
        type=Path,
        default=Path("/tmp/obs_scaling_fixture"),
        help="Root containing EXP_RUN/output/DOMAIN/stats.",
    )
    parser.add_argument("--python-dir", type=Path, help="Explicit Python output directory.")
    parser.add_argument("--matlab-dir", type=Path, help="Explicit MATLAB output directory.")
    parser.add_argument("--python-out-dir", default="python_fixture_zscore_1deg")
    parser.add_argument("--matlab-out-dir", default="matlab_fixture_zscore_1deg")
    parser.add_argument(
        "--variables",
        default=",".join(DEFAULT_VARIABLES),
        help="Comma-separated data variables to compare.",
    )
    parser.add_argument(
        "--coords",
        default=",".join(DEFAULT_COORDS),
        help="Comma-separated coordinate/metadata variables to compare when --compare-coords is set.",
    )
    parser.add_argument("--compare-coords", action="store_true")
    parser.add_argument("--atol", type=float, default=1.0e-6)
    parser.add_argument("--rtol", type=float, default=0.0)
    parser.add_argument("--require-same-files", action="store_true")
    parser.add_argument(
        "--key-contains",
        default="",
        help="Optional substring filter for matched file keys, for example DOY121.",
    )
    parser.add_argument("--max-missing-list", type=int, default=8)
    return parser.parse_args(argv)


def resolve_dirs(args: argparse.Namespace) -> tuple[Path, Path]:
    stats_dir = args.fixture_root / EXP_RUN / "output" / DOMAIN / "stats"
    python_dir = args.python_dir or stats_dir / args.python_out_dir
    matlab_dir = args.matlab_dir or stats_dir / args.matlab_out_dir
    return python_dir, matlab_dir


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def index_outputs(directory: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in sorted(directory.glob("*.nc4")):
        out[comparison_key(path)] = path
    return out


def comparison_key(path: Path) -> str:
    name = path.name
    if MATCH_TOKEN in name:
        return name.split(MATCH_TOKEN, 1)[1]
    return name


def compare_file(
    python_path: Path,
    matlab_path: Path,
    *,
    variables: list[str],
    coords: list[str],
    atol: float,
    rtol: float,
) -> list[VariableDiff]:
    diffs: list[VariableDiff] = []
    with nc.Dataset(python_path) as py_ds, nc.Dataset(matlab_path) as ml_ds:
        py_dims = {name: len(dim) for name, dim in py_ds.dimensions.items()}
        ml_dims = {name: len(dim) for name, dim in ml_ds.dimensions.items()}
        if py_dims != ml_dims:
            print(f"  dimension mismatch: python={py_dims} matlab={ml_dims}")

        for name in coords + variables:
            diffs.append(compare_variable(py_ds, ml_ds, name, atol=atol, rtol=rtol))
    return diffs


def compare_variable(py_ds: nc.Dataset, ml_ds: nc.Dataset, name: str, *, atol: float, rtol: float) -> VariableDiff:
    if name not in py_ds.variables or name not in ml_ds.variables:
        return VariableDiff(
            name=name,
            shape_python=variable_shape(py_ds, name),
            shape_matlab=variable_shape(ml_ds, name),
            finite_python=0,
            finite_matlab=0,
            finite_mask_mismatch=0,
            compared_points=0,
            max_abs_diff=np.nan,
            mean_abs_diff=np.nan,
            status="MISSING",
        )

    py_data = read_numeric(py_ds.variables[name])
    ml_data = read_numeric(ml_ds.variables[name])
    if py_data.shape != ml_data.shape:
        return VariableDiff(
            name=name,
            shape_python=py_data.shape,
            shape_matlab=ml_data.shape,
            finite_python=int(np.isfinite(py_data).sum()),
            finite_matlab=int(np.isfinite(ml_data).sum()),
            finite_mask_mismatch=-1,
            compared_points=0,
            max_abs_diff=np.nan,
            mean_abs_diff=np.nan,
            status="SHAPE",
        )

    py_finite = np.isfinite(py_data)
    ml_finite = np.isfinite(ml_data)
    common = py_finite & ml_finite
    mask_mismatch = int(np.count_nonzero(py_finite != ml_finite))

    if np.any(common):
        abs_diff = np.abs(py_data[common] - ml_data[common])
        max_abs = float(np.max(abs_diff))
        mean_abs = float(np.mean(abs_diff))
        tolerance = atol + rtol * np.abs(ml_data[common])
        numeric_ok = bool(np.all(abs_diff <= tolerance))
    else:
        max_abs = 0.0 if py_data.size == 0 or not (np.any(py_finite) or np.any(ml_finite)) else np.nan
        mean_abs = max_abs
        numeric_ok = True

    status = "OK" if mask_mismatch == 0 and numeric_ok else "DIFF"
    return VariableDiff(
        name=name,
        shape_python=py_data.shape,
        shape_matlab=ml_data.shape,
        finite_python=int(np.count_nonzero(py_finite)),
        finite_matlab=int(np.count_nonzero(ml_finite)),
        finite_mask_mismatch=mask_mismatch,
        compared_points=int(np.count_nonzero(common)),
        max_abs_diff=max_abs,
        mean_abs_diff=mean_abs,
        status=status,
    )


def variable_shape(dataset: nc.Dataset, name: str) -> tuple[int, ...]:
    if name not in dataset.variables:
        return ()
    return tuple(dataset.variables[name].shape)


def read_numeric(variable: nc.Variable) -> np.ndarray:
    data = np.ma.filled(variable[...], np.nan)
    arr = np.asarray(data)
    if arr.dtype.kind in {"S", "U", "O"}:
        raise TypeError(f"Cannot compare non-numeric variable {variable.name!r}")
    return arr.astype(float, copy=False)


def print_table(diffs: list[VariableDiff]) -> None:
    header = (
        "  variable        status  finite_py  finite_ml  mask_diff  "
        "points      max_abs      mean_abs"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for diff in diffs:
        print(
            f"  {diff.name:<14} {diff.status:<7}"
            f" {diff.finite_python:>9}"
            f" {diff.finite_matlab:>10}"
            f" {diff.finite_mask_mismatch:>10}"
            f" {diff.compared_points:>7}"
            f" {format_float(diff.max_abs_diff):>12}"
            f" {format_float(diff.mean_abs_diff):>13}"
        )


def format_float(value: float) -> str:
    if np.isnan(value):
        return "nan"
    return f"{value:.4e}"


if __name__ == "__main__":
    raise SystemExit(main())
