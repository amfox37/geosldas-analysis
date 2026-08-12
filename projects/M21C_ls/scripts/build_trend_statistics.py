#!/usr/bin/env python3
"""Build M21C tile-level trend statistics from the audited monthly loader."""

from __future__ import annotations

import argparse
from pathlib import Path

from trend_breakpoint_series import DEFAULT_DATA_DIR, MonthlySeriesLoader
from trend_statistics import DEFAULT_CONFIG, compute_trend_statistics, load_trend_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output" / "trends_breakpoints"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variable", required=True)
    parser.add_argument("--mask", default="valid_land")
    parser.add_argument(
        "--series",
        choices=["auto", "delta", "ol", "da", "value"],
        default="auto",
        help="auto selects delta for paired variables and value for DA-only diagnostics",
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tile-start", type=int, default=0)
    parser.add_argument("--tile-stop", type=int)
    parser.add_argument("--n-jobs", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    settings = load_trend_config(args.config)
    with MonthlySeriesLoader(args.data_dir) as loader:
        tile_stop = args.tile_stop if args.tile_stop is not None else loader.n_tile
        if args.tile_start < 0 or tile_stop > loader.n_tile or tile_stop <= args.tile_start:
            raise ValueError(
                f"Invalid tile range [{args.tile_start}, {tile_stop}) for {loader.n_tile} tiles"
            )
        is_subset = args.tile_start != 0 or tile_stop != loader.n_tile
        if is_subset and args.output is None:
            raise ValueError("Subset trend runs require an explicit --output path")
        fields = loader.load_tile_series(args.variable, mask=args.mask)
        series_name = args.series
        if series_name == "auto":
            series_name = "delta" if "delta" in fields else "value"
        if series_name not in fields:
            raise KeyError(
                f"Series {series_name!r} is unavailable; choices are {list(fields.data_vars)}"
            )
        tile_slice = slice(args.tile_start, tile_stop)
        values = fields[series_name].isel(tile=tile_slice).load()
        eligible = fields["eligible"].isel(tile=tile_slice).load()

    result = compute_trend_statistics(
        values,
        eligible,
        config=settings,
        n_jobs=args.n_jobs,
    )
    result.attrs.update(
        mask=args.mask,
        selected_series=series_name,
        tile_start=int(args.tile_start),
        tile_stop=int(tile_stop),
        total_input_tiles=int(loader.n_tile),
        fdr_scope=(
            "diagnostic selected tile subset; not global FDR"
            if is_subset
            else "all finite tiles eligible under the selected mask"
        ),
        trend_config=str(args.config),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output or (
        args.output_dir
        / f"{args.variable}_{series_name}_{args.mask}_trend_statistics.nc"
    )
    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in result.data_vars
        if result[name].dtype.kind in "fiu"
    }
    result.to_netcdf(output_path, encoding=encoding)
    print(f"Wrote {output_path}")
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
