#!/usr/bin/env python
"""Compare a sparse Python step-4 result with a MATLAB step-4 MAT file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


FIELDS = (
    "R2_ivd_mod",
    "R2_ivd_obs",
    "R2_ivs_mod",
    "R2_ivs_obs",
    "R_mod_obs",
)
UNDEFINED_CORRELATION = -9999.0


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    matlab = load_mat_fields(args.matlab, ("N_sm", *FIELDS))
    with np.load(args.python) as python:
        grid_size = int(python["grid_size"])
        idx = np.asarray(python["idx0"], dtype=np.int64)
        count = expand_sparse(idx, python["N_sm"], grid_size, fill=0, dtype=np.int64)

        failed = False
        matlab_count = flatten(matlab["N_sm"]).astype(np.int64)
        require_size(matlab_count, grid_size, "N_sm")
        count_mismatch = int(np.count_nonzero(count != matlab_count))
        print(
            f"N_sm: python_nonzero={np.count_nonzero(count):,} "
            f"matlab_nonzero={np.count_nonzero(matlab_count):,} "
            f"mismatched={count_mismatch:,}"
        )
        failed |= count_mismatch > 0

        for field in FIELDS:
            fill = UNDEFINED_CORRELATION if field == "R_mod_obs" else np.nan
            python_full = expand_sparse(
                idx, python[field], grid_size, fill=fill, dtype=np.float64
            )
            matlab_full = flatten(matlab[field]).astype(np.float64)
            require_size(matlab_full, grid_size, field)
            python_valid = valid_mask(field, python_full)
            matlab_valid = valid_mask(field, matlab_full)
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


def expand_sparse(idx, values, size, *, fill, dtype):
    """Expand one sparse Python field onto MATLAB's full linear grid."""

    idx = np.asarray(idx, dtype=np.int64)
    values = np.asarray(values, dtype=dtype)
    if idx.ndim != 1 or values.ndim != 1 or idx.size != values.size:
        raise ValueError("Sparse index and value arrays must be equal-length vectors")
    if np.any(idx < 0) or np.any(idx >= size):
        raise ValueError("Sparse index is outside grid_size")
    output = np.full(size, fill, dtype=dtype)
    output[idx] = values
    return output


def valid_mask(field: str, values: np.ndarray) -> np.ndarray:
    mask = np.isfinite(values)
    if field == "R_mod_obs":
        mask &= values != UNDEFINED_CORRELATION
    return mask


def flatten(values) -> np.ndarray:
    return np.asarray(values).reshape(-1, order="F")


def require_size(values: np.ndarray, expected: int, field: str) -> None:
    if values.size != expected:
        raise ValueError(
            f"MATLAB {field} has {values.size} values; expected grid_size={expected}"
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
    parser.add_argument("--python", type=Path, required=True, help="Step-4 NPZ file.")
    parser.add_argument("--matlab", type=Path, required=True, help="MATLAB step-4 MAT file.")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
