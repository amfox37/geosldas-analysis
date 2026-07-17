#!/usr/bin/env python
"""Compute sparse MATLAB-compatible step-4 IVd/IVs statistics."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.climatology import (  # noqa: E402
    climatology_output_path,
    daily_pair_paths,
    load_pentad_climatology_npz,
)
from iv_tc.independent_validation import (  # noqa: E402
    compute_independent_validation_from_paths,
    independent_validation_npz_is_valid,
    iv_output_path,
    save_independent_validation_npz,
)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    start_date = parse_date(args.start_date)
    end_date = parse_date(args.end_date)
    climatology_root = args.climatology_root or args.pair_root
    output_root = args.output_root or args.pair_root

    written = 0
    existing = 0
    failed = 0
    for sensor in args.sensor:
        for run_name in args.run:
            output_path = iv_output_path(
                output_root,
                sensor,
                run_name,
                start_date,
                end_date,
                args.lag_days,
                args.min_count,
            )
            if not args.overwrite and independent_validation_npz_is_valid(output_path):
                existing += 1
                print(f"[exists] {sensor} {run_name} -> {output_path}")
                continue

            clim_path = climatology_output_path(
                climatology_root,
                sensor,
                run_name,
                start_date,
                end_date,
                args.window_days,
            )
            if not clim_path.exists():
                failed += 1
                print(f"[missing] {sensor} {run_name}: climatology {clim_path}")
                continue

            paths, missing_dates = daily_pair_paths(
                args.pair_root, sensor, run_name, start_date, end_date
            )
            if not paths:
                failed += 1
                print(
                    f"[missing] {sensor} {run_name}: no daily pairs in "
                    f"{start_date} through {end_date}"
                )
                continue

            try:
                climatology = load_pentad_climatology_npz(clim_path)
                result = compute_independent_validation_from_paths(
                    paths,
                    climatology,
                    lag_days=args.lag_days,
                    min_count=args.min_count,
                    start_date=start_date,
                    end_date=end_date,
                )
            except (OSError, ValueError) as exc:
                failed += 1
                print(f"[failed] {sensor} {run_name}: {exc}")
                continue

            save_independent_validation_npz(result, output_path)
            written += 1
            finite_ivd = int(np.isfinite(result.r2_ivd_model).sum())
            print(
                f"[written] {sensor} {run_name}: source_days={result.source_days} "
                f"missing={len(missing_dates)} paired_days={result.paired_days}/"
                f"{result.candidate_days} cells={result.idx.size} "
                f"finite_R2_ivd_mod={finite_ivd} Nmin={result.min_count} -> "
                f"{output_path}"
            )

    print(f"Summary: written={written}, exists={existing}, failed={failed}")
    if failed:
        raise SystemExit(1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", required=True, help="First date, inclusive.")
    parser.add_argument("--end-date", required=True, help="Last date, inclusive.")
    parser.add_argument("--sensor", action="append", required=True)
    parser.add_argument("--run", action="append", required=True)
    parser.add_argument(
        "--pair-root",
        type=Path,
        required=True,
        help="Root containing step2_pairs/{sensor}/{run}/YYYYMMDD.npz.",
    )
    parser.add_argument(
        "--climatology-root",
        type=Path,
        help="Root containing step3_climatology; defaults to --pair-root.",
    )
    parser.add_argument(
        "--output-root", type=Path, help="Output root; defaults to --pair-root."
    )
    parser.add_argument("--lag-days", type=int, default=2)
    parser.add_argument("--min-count", type=int, default=100)
    parser.add_argument(
        "--window-days",
        type=int,
        default=31,
        help="Window token used to resolve the step-3 climatology filename.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def parse_date(value: str):
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            pass
    raise ValueError(f"Date must be YYYY-MM-DD or YYYYMMDD, got {value!r}")


if __name__ == "__main__":
    main()
