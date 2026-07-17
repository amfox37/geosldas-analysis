"""Compute sparse pentad climatologies from IV/TC step-2 daily pairs."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.climatology import (  # noqa: E402
    climatology_npz_is_valid,
    climatology_output_path,
    compute_pentad_climatology_from_paths,
    daily_pair_paths,
    matlab_default_min_count,
    save_pentad_climatology_npz,
)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    start_date = parse_date(args.start_date)
    end_date = parse_date(args.end_date)
    min_count = (
        matlab_default_min_count(start_date, end_date)
        if args.min_count is None
        else args.min_count
    )
    output_root = args.output_root or args.pair_root

    written = 0
    existing = 0
    failed = 0
    for sensor in args.sensor:
        for run_name in args.run:
            output_path = climatology_output_path(
                output_root,
                sensor,
                run_name,
                start_date,
                end_date,
                args.window_days,
            )
            if not args.overwrite and climatology_npz_is_valid(output_path):
                existing += 1
                print(f"[exists] {sensor} {run_name} -> {output_path}")
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
                climatology = compute_pentad_climatology_from_paths(
                    paths,
                    min_count=min_count,
                    window_days=args.window_days,
                )
            except (OSError, ValueError) as exc:
                failed += 1
                print(f"[failed] {sensor} {run_name}: {exc}")
                continue

            save_pentad_climatology_npz(climatology, output_path)
            written += 1
            finite = int(np.isfinite(climatology.obs).sum())
            print(
                f"[written] {sensor} {run_name}: days={climatology.source_days} "
                f"missing={len(missing_dates)} cells={climatology.idx.size} "
                f"finite_cell_pentads={finite} Nday_min={min_count} -> {output_path}"
            )

    print(f"Summary: written={written}, exists={existing}, failed={failed}")
    if failed:
        raise SystemExit(1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", required=True, help="First date, inclusive.")
    parser.add_argument("--end-date", required=True, help="Last date, inclusive.")
    parser.add_argument("--sensor", action="append", required=True, help="Step-2 sensor token.")
    parser.add_argument("--run", action="append", required=True, help="Step-2 run name.")
    parser.add_argument(
        "--pair-root",
        type=Path,
        required=True,
        help="Root containing step2_pairs/{sensor}/{run}/YYYYMMDD.npz.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Output root; defaults to --pair-root.",
    )
    parser.add_argument("--window-days", type=int, default=31, help="Odd circular window length.")
    parser.add_argument(
        "--min-count",
        type=int,
        help="Minimum samples per cell/pentad; default matches MATLAB's 4 * year span.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace valid climatology files.")
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
