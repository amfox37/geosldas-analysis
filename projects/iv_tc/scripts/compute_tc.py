#!/usr/bin/env python
"""Compute sparse MATLAB-compatible triple-collocation statistics."""

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
from iv_tc.generation import canonical_sensor  # noqa: E402
from iv_tc.triple_collocation import (  # noqa: E402
    compute_triple_collocation_from_paths,
    save_triple_collocation_npz,
    tc_output_path,
    triple_collocation_npz_is_valid,
)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    start_date = parse_date(args.start_date)
    end_date = parse_date(args.end_date)
    climatology_start_date = (
        parse_date(args.climatology_start_date)
        if args.climatology_start_date
        else start_date
    )
    climatology_end_date = (
        parse_date(args.climatology_end_date)
        if args.climatology_end_date
        else end_date
    )
    climatology_root = args.climatology_root or args.pair_root
    output_root = args.output_root or args.pair_root
    primary_name = args.primary_name or canonical_sensor(args.primary_sensor)
    secondary_name = args.secondary_name or canonical_sensor(args.secondary_sensor)
    secondary_scale = _secondary_scale(args.secondary_sensor, args.secondary_scale)

    written = 0
    existing = 0
    failed = 0
    for run_name in args.run:
        source_names = (primary_name, args.model_name, secondary_name)
        output_path = tc_output_path(
            output_root,
            source_names,
            run_name,
            start_date,
            end_date,
            season=args.season,
            min_count=args.min_count,
        )
        if not args.overwrite and triple_collocation_npz_is_valid(output_path):
            existing += 1
            print(f"[exists] {' + '.join(source_names)} {run_name} -> {output_path}")
            continue

        primary_clim_path = climatology_output_path(
            climatology_root,
            args.primary_sensor,
            run_name,
            climatology_start_date,
            climatology_end_date,
            args.window_days,
        )
        secondary_clim_path = climatology_output_path(
            climatology_root,
            args.secondary_sensor,
            run_name,
            climatology_start_date,
            climatology_end_date,
            args.window_days,
        )
        missing_climatologies = [
            path
            for path in (primary_clim_path, secondary_clim_path)
            if not path.exists()
        ]
        if missing_climatologies:
            failed += 1
            print(
                f"[missing] {run_name}: climatology "
                + ", ".join(str(path) for path in missing_climatologies)
            )
            continue

        primary_paths, primary_missing = daily_pair_paths(
            args.pair_root,
            args.primary_sensor,
            run_name,
            start_date,
            end_date,
        )
        secondary_paths, secondary_missing = daily_pair_paths(
            args.pair_root,
            args.secondary_sensor,
            run_name,
            start_date,
            end_date,
        )
        try:
            primary_climatology = load_pentad_climatology_npz(primary_clim_path)
            secondary_climatology = load_pentad_climatology_npz(
                secondary_clim_path
            )
            result = compute_triple_collocation_from_paths(
                primary_paths,
                primary_climatology,
                secondary_paths,
                secondary_climatology,
                primary_name=primary_name,
                secondary_name=secondary_name,
                model_name=args.model_name,
                primary_scale=args.primary_scale,
                secondary_scale=secondary_scale,
                model_scale=args.model_scale,
                model_values_source=args.model_values_source,
                model_climatology_source=args.model_climatology_source,
                start_date=start_date,
                end_date=end_date,
                season=args.season,
                min_count=args.min_count,
                r2_floor=args.r2_floor,
                strict_missing=not args.allow_missing,
                model_agreement_tolerance=args.model_agreement_tolerance,
            )
        except (OSError, ValueError) as exc:
            failed += 1
            print(f"[failed] {' + '.join(source_names)} {run_name}: {exc}")
            continue

        save_triple_collocation_npz(result, output_path)
        written += 1
        finite_r2 = tuple(int(np.isfinite(result.r2[:, i]).sum()) for i in range(3))
        print(
            f"[written] {' + '.join(source_names)} {run_name}: "
            f"paired_days={result.paired_days}/{result.candidate_days} "
            f"missing={result.missing_days} cells={result.idx.size} "
            f"finite_R2={finite_r2} Nmin={result.min_count} "
            f"scales={result.source_scales} -> {output_path}"
        )
        if primary_missing or secondary_missing:
            print(
                f"  input gaps: primary={len(primary_missing)} "
                f"secondary={len(secondary_missing)}"
            )

    print(f"Summary: written={written}, exists={existing}, failed={failed}")
    if failed:
        raise SystemExit(1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", required=True, help="First date, inclusive.")
    parser.add_argument("--end-date", required=True, help="Last date, inclusive.")
    parser.add_argument("--primary-sensor", required=True)
    parser.add_argument("--secondary-sensor", required=True)
    parser.add_argument(
        "--climatology-start-date",
        help="Step-3 filename start date; defaults to --start-date.",
    )
    parser.add_argument(
        "--climatology-end-date",
        help="Step-3 filename end date; defaults to --end-date.",
    )
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
    parser.add_argument("--primary-name")
    parser.add_argument("--secondary-name")
    parser.add_argument("--model-name", default="model")
    parser.add_argument("--primary-scale", type=float, default=1.0)
    parser.add_argument(
        "--secondary-scale",
        type=float,
        help="Defaults to 0.005 for ASCAT and 1 otherwise.",
    )
    parser.add_argument("--model-scale", type=float, default=1.0)
    parser.add_argument(
        "--model-values-source", choices=("primary", "secondary"), default="primary"
    )
    parser.add_argument(
        "--model-climatology-source",
        choices=("primary", "secondary"),
        default="secondary",
        help="MATLAB's two-file TC scripts use the secondary climatology.",
    )
    parser.add_argument(
        "--season", choices=("annual", "summer", "winter"), default="annual"
    )
    parser.add_argument("--min-count", type=int, default=20)
    parser.add_argument("--r2-floor", type=float, default=0.001)
    parser.add_argument("--window-days", type=int, default=31)
    parser.add_argument(
        "--model-agreement-tolerance",
        type=float,
        help="Fail if overlapping model values differ by more than this amount.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing days; MATLAB parity mode fails on any missing pair.",
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


def _secondary_scale(sensor: str, requested: float | None) -> float:
    if requested is not None:
        return requested
    return 0.005 if canonical_sensor(sensor).startswith("ascat_") else 1.0


if __name__ == "__main__":
    main()
