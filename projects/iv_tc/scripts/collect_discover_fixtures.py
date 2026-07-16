"""Collect tiny Discover fixtures for IV/TC adapter development.

Default behavior is a dry run. Add --copy to copy existing files into the
portable fixture tree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.config import ProductRoots, parse_run_spec  # noqa: E402
from iv_tc.fixtures import (  # noqa: E402
    collect_fixture_files,
    copy_fixture_files,
    parse_date,
    print_fixture_plan,
)


DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "test_data" / "inputs" / "discover_sample"


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    dates = [parse_date(value) for value in args.date]
    runs = [parse_run_spec(value) for value in args.run]
    roots = ProductRoots(
        ascat_h121_root=args.ascat_h121_root,
        ascat_h119_h120_root=args.ascat_h119_h120_root,
        smosic_root=args.smosic_root,
        smap_l3_root=args.smap_l3_root,
        domain=args.domain,
        collection=args.collection,
    )

    files = collect_fixture_files(dates, args.output_root, roots, runs)
    print_fixture_plan(files)

    found = sum(1 for item in files if item.exists)
    missing = len(files) - found
    print("")
    print(f"Summary: {found} found, {missing} missing, {len(files)} planned")

    if not args.copy:
        print("Dry run only. Add --copy to copy existing files.")
        return

    copied, skipped = copy_fixture_files(files)
    print(f"Copied {copied} files; skipped {skipped} missing files.")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--date",
        action="append",
        required=True,
        help="Date to collect, YYYY-MM-DD or YYYYMMDD. Repeat for multiple dates.",
    )
    parser.add_argument(
        "--run",
        action="append",
        default=[],
        help="GEOSldas run as NAME=/path/to/run/root. Repeat for multiple runs.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Portable fixture output root.",
    )
    parser.add_argument(
        "--copy",
        action="store_true",
        help="Copy existing files. Without this flag, only print the plan.",
    )
    parser.add_argument(
        "--ascat-h121-root",
        "--ascat-cdr-root",
        dest="ascat_h121_root",
        type=Path,
        default=ProductRoots.ascat_h121_root,
        help="Root containing ASCAT_SSM_CDR/H121 and ASCAT_SSM_CDR/flists.",
    )
    parser.add_argument(
        "--ascat-h119-h120-root",
        dest="ascat_h119_h120_root",
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
    parser.add_argument(
        "--domain",
        default=ProductRoots.domain,
        help="GEOSldas domain name.",
    )
    parser.add_argument(
        "--collection",
        default=ProductRoots.collection,
        help="GEOSldas output collection tag.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
