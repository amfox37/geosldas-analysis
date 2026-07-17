"""Generate daily IV/TC step-2 obs/model pair files.

Outputs are written as:
    output_root/step2_pairs/{sensor}/{run_name}/YYYYMMDD.npz
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.config import ProductRoots, parse_run_spec  # noqa: E402
from iv_tc.generation import date_range, generate_daily_pairs  # noqa: E402


DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "test_data" / "outputs"


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    dates = date_range(parse_date(args.start_date), parse_date(args.end_date))
    runs = [parse_run_spec(value) for value in args.run]
    roots = ProductRoots(
        ascat_h121_root=args.ascat_h121_root,
        ascat_h119_h120_root=args.ascat_h119_h120_root,
        smosic_root=args.smosic_root,
        smap_l3_root=args.smap_l3_root,
        domain=args.domain,
        collection=args.collection,
    )

    results = generate_daily_pairs(
        dates=dates,
        sensors=args.sensor,
        runs=runs,
        roots=roots,
        output_root=args.output_root,
        skip_existing=not args.overwrite,
        model_variable=args.model_variable,
        h119_h120_method=args.h119_h120_method,
        h119_h120_fill_nearest=args.h119_h120_fill_nearest,
    )

    counts = Counter(result.status for result in results)
    for result in results:
        detail = f" n={result.n_points}" if result.n_points else ""
        message = f" ({result.message})" if result.message else ""
        print(
            f"[{result.status}] {result.sensor} {result.run} {result.date:%Y-%m-%d}"
            f"{detail} -> {result.output_path}{message}"
        )

    summary = ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
    print(f"Summary: {summary}; total={len(results)}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", required=True, help="First date, YYYY-MM-DD or YYYYMMDD.")
    parser.add_argument("--end-date", required=True, help="Last date, YYYY-MM-DD or YYYYMMDD.")
    parser.add_argument(
        "--sensor",
        action="append",
        required=True,
        help=(
            "Sensor to process. Repeat for multiple sensors. "
            "Use smosic, smap_l3, ascat_h121, or ascat_h119_h120."
        ),
    )
    parser.add_argument(
        "--run",
        action="append",
        required=True,
        help="GEOSldas run as NAME=/path/to/run/root. Repeat for multiple runs.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Root for step2_pairs output.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Regenerate existing valid pair files.")
    parser.add_argument("--model-variable", default="SFMC", help="Model tile-space variable to pair.")
    parser.add_argument(
        "--h119-h120-method",
        default="linear",
        choices=("linear", "nearest", "cubic"),
        help="SciPy griddata method for ASCAT H119/H120.",
    )
    parser.add_argument(
        "--h119-h120-fill-nearest",
        action="store_true",
        help=(
            "Fill gaps outside the interpolation convex hull with nearest-neighbor "
            "values. Off by default to match the legacy MATLAB griddata(...,'natural') "
            "behavior, which leaves those cells as NaN rather than extrapolating."
        ),
    )
    parser.add_argument(
        "--ascat-h121-root",
        type=Path,
        default=ProductRoots.ascat_h121_root,
        help="Root containing ASCAT_SSM_CDR/H121.",
    )
    parser.add_argument(
        "--ascat-h119-h120-root",
        type=Path,
        default=ProductRoots.ascat_h119_h120_root,
        help="Root containing H119_H120_processed and Auxiliary.",
    )
    parser.add_argument(
        "--smosic-root",
        type=Path,
        default=ProductRoots.smosic_root,
        help="Root containing smos_ic_sm_m36_YYYYMMDD.nc files.",
    )
    parser.add_argument(
        "--smap-l3-root",
        type=Path,
        default=ProductRoots.smap_l3_root,
        help="Root containing the flat SPL3SMP_v009/YYYYY tree.",
    )
    parser.add_argument("--domain", default=ProductRoots.domain, help="GEOSldas domain name.")
    parser.add_argument("--collection", default=ProductRoots.collection, help="GEOSldas output collection tag.")
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
