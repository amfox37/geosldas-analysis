#!/usr/bin/env python
"""Compare sparse Python TC output with a full-grid MATLAB TC MAT file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


DEFAULT_R2 = ("R2_TC_SMOS", "R2_TC_mod", "R2_TC_ASC")
DEFAULT_SIGMA2 = ("sigma2_SMOS", "sigma2_mod", "sigma2_ASC")
DEFAULT_CORRELATION = ("R_mod_SMOS", "R_ASC_SMOS", "R_mod_ASC")
DEFAULT_COVARIANCE = ("C_SMOS_mod", "C_SMOS_ASC", "C_mod_ASC")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    mappings = (
        ("R2_TC", parse_names(args.matlab_r2), "r2"),
        ("sigma2", parse_names(args.matlab_sigma2), "source"),
        (
            "correlation",
            parse_names(args.matlab_correlation),
            "pair",
        ),
        ("covariance", parse_names(args.matlab_covariance), "pair"),
    )
    matlab_fields = ("N_sm",) + tuple(
        name for _, names, _ in mappings for name in names
    )
    matlab = load_mat_fields(args.matlab, matlab_fields)

    failed = False
    with np.load(args.python) as python:
        grid_size = int(python["grid_size"])
        idx = np.asarray(python["idx0"], dtype=np.int64)
        python_count = expand_sparse(
            idx, python["N_sm"], grid_size, fill=0, dtype=np.int64
        )
        matlab_count = flatten(matlab["N_sm"]).astype(np.int64)
        require_size(matlab_count, grid_size, "N_sm")
        count_mismatch = int(np.count_nonzero(python_count != matlab_count))
        print(
            f"N_sm: python_nonzero={np.count_nonzero(python_count):,} "
            f"matlab_nonzero={np.count_nonzero(matlab_count):,} "
            f"mismatched={count_mismatch:,}"
        )
        failed |= count_mismatch > 0

        for python_field, matlab_names, _ in mappings:
            values = np.asarray(python[python_field], dtype=np.float64)
            if values.shape != (idx.size, 3):
                raise ValueError(
                    f"Python {python_field} has shape {values.shape}; "
                    f"expected {(idx.size, 3)}"
                )
            for column, matlab_name in enumerate(matlab_names):
                python_full = expand_sparse(
                    idx,
                    values[:, column],
                    grid_size,
                    fill=np.nan,
                    dtype=np.float64,
                )
                matlab_full = flatten(matlab[matlab_name]).astype(np.float64)
                require_size(matlab_full, grid_size, matlab_name)
                python_valid = np.isfinite(python_full)
                matlab_valid = np.isfinite(matlab_full)
                mask_mismatch = int(
                    np.count_nonzero(python_valid != matlab_valid)
                )
                common = python_valid & matlab_valid
                difference = np.abs(python_full[common] - matlab_full[common])
                max_abs = float(np.max(difference)) if difference.size else np.nan
                mean_abs = float(np.mean(difference)) if difference.size else np.nan
                numeric_ok = bool(
                    np.allclose(
                        python_full[common],
                        matlab_full[common],
                        rtol=args.rtol,
                        atol=args.atol,
                    )
                )
                print(
                    f"{python_field}[{column}] vs {matlab_name}: "
                    f"python_valid={python_valid.sum():,} "
                    f"matlab_valid={matlab_valid.sum():,} "
                    f"mask_mismatched={mask_mismatch:,} "
                    f"mean_abs={mean_abs:.6g} max_abs={max_abs:.6g} "
                    f"within_tolerance={numeric_ok}"
                )
                failed |= mask_mismatch > 0 or not numeric_ok

    if failed:
        raise SystemExit(1)
    print("Parity check passed.")


def expand_sparse(idx, values, size, *, fill, dtype):
    idx = np.asarray(idx, dtype=np.int64)
    values = np.asarray(values, dtype=dtype)
    if idx.ndim != 1 or values.ndim != 1 or idx.size != values.size:
        raise ValueError("Sparse index and value arrays must be equal-length vectors")
    if np.any(idx < 0) or np.any(idx >= size):
        raise ValueError("Sparse index is outside grid_size")
    output = np.full(size, fill, dtype=dtype)
    output[idx] = values
    return output


def flatten(values) -> np.ndarray:
    return np.asarray(values).reshape(-1, order="F")


def require_size(values: np.ndarray, expected: int, field: str) -> None:
    if values.size != expected:
        raise ValueError(
            f"MATLAB {field} has {values.size} values; expected grid_size={expected}"
        )


def load_mat_fields(path: Path, fields: tuple[str, ...]) -> dict[str, np.ndarray]:
    try:
        from scipy.io import loadmat

        data = loadmat(path, variable_names=list(fields))
        return {field: np.asarray(data[field]) for field in fields}
    except NotImplementedError:
        import h5py

        with h5py.File(path) as data:
            return {field: np.asarray(data[field]) for field in fields}


def parse_names(value: str) -> tuple[str, str, str]:
    names = tuple(name.strip() for name in value.split(","))
    if len(names) != 3 or any(not name for name in names):
        raise ValueError("MATLAB field mappings require exactly three comma-separated names")
    return names


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path, required=True, help="Step-4 TC NPZ.")
    parser.add_argument("--matlab", type=Path, required=True, help="MATLAB TC MAT.")
    parser.add_argument("--matlab-r2", default=",".join(DEFAULT_R2))
    parser.add_argument("--matlab-sigma2", default=",".join(DEFAULT_SIGMA2))
    parser.add_argument(
        "--matlab-correlation", default=",".join(DEFAULT_CORRELATION)
    )
    parser.add_argument("--matlab-covariance", default=",".join(DEFAULT_COVARIANCE))
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
