#!/usr/bin/env python3
"""Generate GEOS-LDAS tile-space z-score scaling parameters for SMAP Tb."""

from __future__ import annotations

import argparse
from math import ceil
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.io import read_obs_param  # noqa: E402
from obs_scaling.tb_tile_stats import generate_tb_scaling_params  # noqa: E402
from scripts.run_scaling_params import (  # noqa: E402
    compute_year_bounds,
    parse_int_csv,
    parse_year_month,
)


def parse_float_csv(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("At least one incidence angle is required")
    return values


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate SMAP Tb scaling climatology in GEOS-LDAS tile-space binary format."
    )
    parser.add_argument("--exp-path", type=Path, required=True)
    parser.add_argument("--exp-run", required=True)
    parser.add_argument("--domain", default="SMAP_EASEv2_M09_GLOBAL")
    parser.add_argument("--start", required=True, help="Start month, YYYY-MM")
    parser.add_argument("--end", required=True, help="End month, YYYY-MM")
    parser.add_argument("--orbit", choices=("A", "D"), default="D")
    parser.add_argument("--angles", default="40", help="Comma-separated incidence angles")
    parser.add_argument("--description-contains", default="SMAP_L1C")
    parser.add_argument("--window-days", type=int, default=75)
    parser.add_argument("--ndata-min", type=int, default=20)
    parser.add_argument("--dt-assim", type=int, default=3 * 60 * 60)
    parser.add_argument("--t0-assim", type=int, default=0)
    parser.add_argument("--obsfcstana-format", choices=("auto", "bin", "nc4"), default="bin")
    parser.add_argument("--run-months", default="")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--obsparam", type=Path, help="Explicit ldas_obsparam path")
    parser.add_argument("--no-convert-grid", action="store_true", help="Keep every model tile")
    parser.add_argument("--each-doy", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    start_year, start_month = parse_year_month(args.start)
    end_year, end_month = parse_year_month(args.end)
    run_months = (
        parse_int_csv(args.run_months)
        if args.run_months
        else list(range(1, 13)) + list(range(1, ceil(args.window_days / 30) + 1))
    )
    earliest_year, latest_year = compute_year_bounds(
        run_months, start_year, start_month, end_year, end_month
    )
    obsparam = args.obsparam or (
        args.exp_path / args.exp_run / "output" / args.domain
        / f"rc_out/Y{start_year:04d}/M{start_month:02d}"
        / f"{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt"
    )
    obs_params = read_obs_param(obsparam)
    orbit = 1 if args.orbit == "A" else 2
    print("SMAP Tb tile-scaling configuration")
    print("  experiment:", args.exp_path / args.exp_run)
    print("  domain:", args.domain)
    print("  period:", args.start, "to", args.end)
    print("  orbit:", args.orbit)
    print("  angles:", parse_float_csv(args.angles))
    print("  obsparam:", obsparam)
    written = generate_tb_scaling_params(
        run_months=run_months,
        exp_path=args.exp_path,
        exp_run=args.exp_run,
        domain=args.domain,
        start_year=earliest_year,
        end_year=latest_year,
        dt_assim=args.dt_assim,
        t0_assim=args.t0_assim,
        obs_params=obs_params,
        description_contains=args.description_contains,
        orbit=orbit,
        angles=parse_float_csv(args.angles),
        window_days=args.window_days,
        ndata_min=args.ndata_min,
        prefix=args.prefix,
        convert_grid=None if args.no_convert_grid else "EASEv2_M36",
        obsfcstana_format=args.obsfcstana_format,
        print_each_doy=args.each_doy,
    )
    print(f"Wrote {len(written)} files")


if __name__ == "__main__":
    main()
