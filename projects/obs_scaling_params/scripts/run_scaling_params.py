"""Run z-score observation scaling climatology generation.

This is the project driver for the Python port of the legacy MATLAB
`get_model_and_obs_clim_stats_latlon_grid.m` workflow. Defaults match the
land-sweeper ASCAT configuration used during MATLAB/Python parity checks.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from math import ceil
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.clim_stats import get_model_and_obs_clim_stats_latlon_grid
from obs_scaling.io import read_obs_param


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    start_year, start_month = parse_year_month(args.start)
    end_year, end_month = parse_year_month(args.end)
    species_names = parse_csv(args.species)
    if args.run_months:
        run_months = parse_int_csv(args.run_months)
    else:
        run_months = list(range(1, 13)) + list(range(1, ceil(args.window_days / 30) + 1))
    earliest_year, latest_year = compute_year_bounds(
        run_months, start_year, start_month, end_year, end_month
    )

    obs_param_fname = (
        args.exp_path
        / args.exp_run
        / "output"
        / args.domain
        / f"rc_out/Y{start_year:04d}/M{start_month:02d}"
        / f"{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt"
    )

    obs_params = read_obs_param(obs_param_fname)
    species = collect_species_ids(obs_params, species_names)

    print("Observation scaling configuration")
    print("  exp_path:", args.exp_path)
    print("  exp_run:", args.exp_run)
    print("  domain:", args.domain)
    print("  period:", args.start, "to", args.end)
    print("  species:", ", ".join(species_names))
    print("  species ids:", species)
    print("  combine_species:", args.combine_species)
    print("  grid_resolution:", args.grid_resolution)
    print("  window_days:", args.window_days)
    print("  ndata_min:", args.ndata_min)
    print("  obsfcstana_format:", args.obsfcstana_format)
    print("  out_dir:", args.out_dir)

    get_model_and_obs_clim_stats_latlon_grid(
        species_names=species_names,
        run_months=run_months,
        exp_path=str(args.exp_path),
        exp_run=args.exp_run,
        domain=args.domain,
        start_year=earliest_year,
        end_year=latest_year,
        dt_assim=args.dt_assim,
        t0_assim=args.t0_assim,
        species=species,
        combine_species_stats=args.combine_species,
        resol=args.grid_resolution,
        w_days=args.window_days,
        ndata_min=args.ndata_min,
        prefix=args.prefix,
        print_each_DOY=args.print_each_doy,
        print_each_pentad=args.print_each_pentad,
        print_all_pentads=args.print_all_pentads,
        out_dir=args.out_dir,
        enable_dedup=args.enable_dedup,
        obsfcstana_format=args.obsfcstana_format,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate ObsFcstAna-based z-score scaling climatology files."
    )
    parser.add_argument(
        "--exp-path",
        type=Path,
        default=Path("/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36"),
        help="Parent directory containing the experiment run directory.",
    )
    parser.add_argument("--exp-run", default="LS_OLv8_M36", help="Experiment/run id.")
    parser.add_argument("--domain", default="SMAP_EASEv2_M36_GLOBAL", help="GEOS-LDAS output domain.")
    parser.add_argument("--start", default="2007-06", help="Start month, YYYY-MM.")
    parser.add_argument("--end", default="2024-05", help="End month, YYYY-MM.")
    parser.add_argument(
        "--species",
        default="ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM",
        help="Comma-separated observation descriptions from ldas_obsparam.",
    )
    parser.add_argument("--grid-resolution", type=float, default=0.25, help="Lat/lon grid spacing in degrees.")
    parser.add_argument("--window-days", type=int, default=75, help="Moving-window length in days.")
    parser.add_argument("--ndata-min", type=int, default=20, help="Minimum samples required per grid cell/window.")
    parser.add_argument("--dt-assim", type=int, default=3 * 60 * 60, help="Assimilation timestep in seconds.")
    parser.add_argument("--t0-assim", type=int, default=0, help="Assimilation reference offset in seconds.")
    parser.add_argument(
        "--obsfcstana-format",
        choices=("auto", "bin", "nc4"),
        default="auto",
        help="ObsFcstAna input format. auto tries .bin first, then .nc4.",
    )
    parser.add_argument(
        "--run-months",
        default="",
        help="Optional comma-separated month override. Defaults to annual climatology months plus window padding.",
    )
    parser.add_argument("--prefix", default="M36_python_dedup_zscore_stats_", help="Output filename prefix.")
    parser.add_argument(
        "--out-dir",
        default="python_z_score_dedup_clim_quarter_degree",
        help="Output subdirectory under EXP/DOMAIN/stats.",
    )
    parser.add_argument("--no-combine-species", dest="combine_species", action="store_false")
    parser.add_argument("--no-dedup", dest="enable_dedup", action="store_false")
    parser.add_argument("--no-each-doy", dest="print_each_doy", action="store_false")
    parser.add_argument("--each-pentad", dest="print_each_pentad", action="store_true")
    parser.add_argument("--no-all-pentads", dest="print_all_pentads", action="store_false")
    parser.set_defaults(
        combine_species=True,
        enable_dedup=True,
        print_each_doy=True,
        print_each_pentad=False,
        print_all_pentads=True,
    )
    return parser.parse_args(argv)


def parse_year_month(value: str) -> tuple[int, int]:
    try:
        dt = datetime.strptime(value, "%Y-%m")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Expected YYYY-MM, got {value!r}") from exc
    return dt.year, dt.month


def parse_csv(value: str) -> list[str]:
    out = [item.strip() for item in value.split(",") if item.strip()]
    if not out:
        raise argparse.ArgumentTypeError("At least one species description is required")
    return out


def parse_int_csv(value: str) -> list[int]:
    months = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not months:
        raise argparse.ArgumentTypeError("At least one month is required")
    invalid = [month for month in months if month < 1 or month > 12]
    if invalid:
        raise argparse.ArgumentTypeError(f"Months must be 1-12, got {invalid}")
    return months


def compute_year_bounds(
    run_months: Iterable[int],
    start_year: int,
    start_month: int,
    end_year: int,
    end_month: int,
) -> tuple[list[int], list[int]]:
    earliest = []
    latest = []
    start_ref = datetime(start_year, start_month, 1)
    end_ref = datetime(end_year, end_month, 1)
    for month in run_months:
        current_start = datetime(start_year, month, 1)
        earliest.append(start_year + 1 if current_start < start_ref else start_year)
        current_end = datetime(end_year, month, 1)
        latest.append(end_year - 1 if current_end > end_ref else end_year)
    return earliest, latest


def collect_species_ids(obs_params, species_names: list[str]) -> list[int]:
    ids: list[int] = []
    for name in species_names:
        matches = [param.species for param in obs_params if param.descr == name]
        if not matches:
            raise ValueError(f"Species {name!r} not found in ldas_obsparam file")
        ids.extend(matches)
    return sorted(set(ids))


if __name__ == "__main__":
    main()
