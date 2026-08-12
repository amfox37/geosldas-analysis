#!/usr/bin/env python3
"""Build M21C tile-level trend statistics from the audited monthly loader."""

from __future__ import annotations

import argparse
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from trend_breakpoint_series import (
    DEFAULT_DATA_DIR,
    DEFAULT_INPUT_CONTRACT,
    DEFAULT_VARIABLE_SELECTION,
    MonthlySeriesLoader,
)
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
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument(
        "--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tile-start", type=int, default=0)
    parser.add_argument("--tile-stop", type=int)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--run-id", default="standalone")
    parser.add_argument("--run-role", choices=["primary", "sensitivity", "standalone"], default="standalone")
    parser.add_argument("--run-matrix", type=Path)
    return parser.parse_args()


def git_commit() -> str:
    """Return the source commit when available without making git a runtime dependency."""

    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT.parents[1], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> int:
    args = parse_args()
    settings = load_trend_config(args.config)
    with MonthlySeriesLoader(
        args.data_dir,
        input_contract=args.input_contract,
        variable_selection=args.variable_selection,
    ) as loader:
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
        input_contract=str(args.input_contract),
        variable_selection=str(args.variable_selection),
        data_directory=str(args.data_dir.resolve()),
        run_id=args.run_id,
        run_role=args.run_role,
        run_matrix=str(args.run_matrix) if args.run_matrix else "",
        generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        source_git_commit=git_commit(),
    )
    output_path = args.output or (
        args.output_dir
        / f"{args.variable}_{series_name}_{args.mask}_trend_statistics.nc"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in result.data_vars
        if result[name].dtype.kind in "fiu"
    }
    temporary_path = output_path.with_name(
        f".{output_path.stem}.incomplete{output_path.suffix}"
    )
    temporary_path.unlink(missing_ok=True)
    try:
        result.to_netcdf(temporary_path, encoding=encoding)
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)
    print(f"Wrote {output_path}")
    month_counts = result["n_calendar_months_used"].values
    unique_counts, frequencies = np.unique(month_counts, return_counts=True)
    print(
        "Calendar months retained (count: tiles): "
        + ", ".join(
            f"{int(count)}: {int(frequency)}"
            for count, frequency in zip(unique_counts, frequencies)
        )
    )
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
