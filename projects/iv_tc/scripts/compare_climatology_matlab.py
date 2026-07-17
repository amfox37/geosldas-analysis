#!/usr/bin/env python
"""Compare a sparse Python step-3 climatology with a MATLAB step-3 MAT file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


FIELDS = ("mod_sm_clim", "obs_sm_clim")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    matlab = load_mat_fields(args.matlab, ("N_sm_clim", *FIELDS))
    with np.load(args.python) as python:
        grid_size = int(python["grid_size"])
        idx = np.asarray(python["idx0"], dtype=np.int64)
        n_pentad = python["mod_sm_clim"].shape[1]

        failed = False
        matlab_count = orient(matlab["N_sm_clim"], grid_size, n_pentad).astype(np.int64)
        require_shape(matlab_count, grid_size, n_pentad, "N_sm_clim")
        count = expand_sparse(
            idx, python["N_sm_clim"], grid_size, n_pentad, fill=0, dtype=np.int64
        )
        count_mismatch = int(np.count_nonzero(count != matlab_count))
        print(
            f"N_sm_clim: python_nonzero={np.count_nonzero(count):,} "
            f"matlab_nonzero={np.count_nonzero(matlab_count):,} "
            f"mismatched={count_mismatch:,}"
        )
        failed |= count_mismatch > 0

        for field in FIELDS:
            python_full = expand_sparse(
                idx, python[field], grid_size, n_pentad, fill=np.nan, dtype=np.float64
            )
            matlab_full = orient(matlab[field], grid_size, n_pentad).astype(np.float64)
            require_shape(matlab_full, grid_size, n_pentad, field)
            python_valid = np.isfinite(python_full)
            matlab_valid = np.isfinite(matlab_full)
            mask_mismatch = int(np.count_nonzero(python_valid != matlab_valid))
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
                f"{field}: python_valid={python_valid.sum():,} "
                f"matlab_valid={matlab_valid.sum():,} "
                f"mask_mismatched={mask_mismatch:,} "
                f"mean_abs={mean_abs:.6g} max_abs={max_abs:.6g} "
                f"within_tolerance={numeric_ok}"
            )
            failed |= mask_mismatch > 0 or not numeric_ok

    if failed:
        raise SystemExit(1)
    print("Parity check passed.")


def expand_sparse(idx, values, size, n_pentad, *, fill, dtype):
    """Expand one sparse Python field onto MATLAB's full linear grid."""

    idx = np.asarray(idx, dtype=np.int64)
    values = np.asarray(values, dtype=dtype)
    if idx.ndim != 1 or values.shape[0] != idx.size or values.shape[1] != n_pentad:
        raise ValueError("Sparse index and value arrays must have matching row counts")
    if np.any(idx < 0) or np.any(idx >= size):
        raise ValueError("Sparse index is outside grid_size")
    output = np.full((size, n_pentad), fill, dtype=dtype)
    output[idx] = values
    return output


def orient(values: np.ndarray, grid_size: int, n_pentad: int) -> np.ndarray:
    """Return a MATLAB climatology array as (grid_size, n_pentad).

    MATLAB stores the array logically as (grid_size, n_pentad), but its
    column-major layout means the v7.3/HDF5 reader (h5py) returns the raw
    storage order transposed, i.e. (n_pentad, grid_size).
    """

    values = np.asarray(values)
    if values.shape == (grid_size, n_pentad):
        return values
    if values.shape == (n_pentad, grid_size):
        return values.T
    raise ValueError(
        f"Unexpected MATLAB array shape {values.shape}; "
        f"expected ({grid_size}, {n_pentad}) or its transpose"
    )


def require_shape(values: np.ndarray, grid_size: int, n_pentad: int, field: str) -> None:
    if values.shape != (grid_size, n_pentad):
        raise ValueError(
            f"MATLAB {field} has shape {values.shape}; "
            f"expected ({grid_size}, {n_pentad})"
        )


def load_mat_fields(path: Path, fields: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Read ordinary or v7.3 MATLAB arrays without changing their values."""

    try:
        from scipy.io import loadmat

        data = loadmat(path, variable_names=list(fields))
        return {field: np.asarray(data[field]) for field in fields}
    except NotImplementedError:
        import h5py

        with h5py.File(path) as data:
            return {field: np.asarray(data[field]) for field in fields}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path, required=True, help="Step-3 NPZ file.")
    parser.add_argument("--matlab", type=Path, required=True, help="MATLAB step-3 MAT file.")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
