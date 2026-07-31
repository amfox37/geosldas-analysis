#!/usr/bin/env python3
"""Generate owner-tile z-score scaling parameters for CYGNSS L1."""

from __future__ import annotations

import argparse
from math import ceil
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.io import read_obs_param  # noqa: E402
from obs_scaling.owner_tile_stats import (  # noqa: E402
    CYGNSS_L1_DESCRIPTION,
    generate_cygnss_l1_scaling_params,
)
from scripts.run_scaling_params import (  # noqa: E402
    compute_year_bounds,
    parse_int_csv,
    parse_year_month,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a CYGNSS L1 scaling climatology on the global M36 owner grid."
    )
    parser.add_argument("--exp-path", type=Path, required=True)
    parser.add_argument("--exp-run", required=True)
    parser.add_argument("--domain", default="SMAP_EASEv2_M36_GLOBAL")
    parser.add_argument("--start", required=True, help="Start month, YYYY-MM")
    parser.add_argument("--end", required=True, help="End month, YYYY-MM")
    parser.add_argument("--description", default=CYGNSS_L1_DESCRIPTION)
    parser.add_argument("--window-days", type=int, default=75)
    parser.add_argument("--ndata-min", type=int, default=20)
    parser.add_argument("--std-epsilon", type=float, default=1e-6)
    parser.add_argument("--dt-assim", type=int, default=3 * 60 * 60)
    parser.add_argument("--t0-assim", type=int, default=0)
    parser.add_argument("--run-months", default="")
    parser.add_argument("--prefix", default="CYGNSS_L1_owner_zscore_")
    parser.add_argument("--out-dir", default="cygnss_l1_z_score_clim")
    parser.add_argument("--obsparam", type=Path, help="Explicit ldas_obsparam path")
    parser.add_argument("--tilecoord", type=Path, help="Explicit ldas_tilecoord.bin path")
    parser.add_argument("--tilegrids", type=Path, help="Explicit ldas_tilegrids.bin path")
    parser.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    parser.set_defaults(overwrite=True)
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
        args.exp_path
        / args.exp_run
        / "output"
        / args.domain
        / f"rc_out/Y{start_year:04d}/M{start_month:02d}"
        / f"{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt"
    )

    print("CYGNSS L1 owner-tile scaling configuration")
    print("  experiment:", args.exp_path / args.exp_run)
    print("  domain:", args.domain)
    print("  period:", args.start, "to", args.end)
    print("  species:", args.description)
    print("  obsparam:", obsparam)
    print("  window_days:", args.window_days)
    print("  ndata_min:", args.ndata_min)

    output = generate_cygnss_l1_scaling_params(
        run_months=run_months,
        exp_path=args.exp_path,
        exp_run=args.exp_run,
        domain=args.domain,
        start_year=earliest_year,
        end_year=latest_year,
        dt_assim=args.dt_assim,
        t0_assim=args.t0_assim,
        obs_params=read_obs_param(obsparam),
        window_days=args.window_days,
        ndata_min=args.ndata_min,
        prefix=args.prefix,
        description=args.description,
        out_dir=args.out_dir,
        std_epsilon=args.std_epsilon,
        tilecoord_path=args.tilecoord,
        tilegrids_path=args.tilegrids,
        overwrite=args.overwrite,
    )
    print("Wrote", output)


if __name__ == "__main__":
    main()
