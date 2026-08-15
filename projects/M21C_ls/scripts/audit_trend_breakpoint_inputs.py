#!/usr/bin/env python3
"""Audit the existing M21C monthly products before trend analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from m21c_periods import load_period_frames
from trend_breakpoint_series import (
    find_legacy_increment_variables,
    load_derived_variables,
    load_variable_selection,
    read_tile_area,
)
from trend_statistics import DEFAULT_CONFIG as DEFAULT_TREND_CONFIG
from trend_statistics import load_trend_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = Path(
    "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2"
)
DEFAULT_DIAGNOSTICS_DIR = PROJECT_ROOT / "output" / "monthly_synthesis_diagnostics"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output" / "trends_breakpoints"
DEFAULT_SELECTION = PROJECT_ROOT / "config" / "trend_breakpoint_variable_selection.json"
DEFAULT_INPUT_CONTRACT = PROJECT_ROOT / "config" / "trend_breakpoint_inputs.json"


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_check(rows: list[dict], check: str, passed: bool, detail: str) -> None:
    rows.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})


def audit(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    checks: list[dict] = []
    manifest: list[dict] = []

    registry, fine, validation, sensors = load_period_frames(args.period_registry)
    add_check(
        checks,
        "period_registry",
        len(fine) == 9 and len(validation) == 3,
        f"P periods={len(fine)}, V periods={len(validation)}, sensors={len(sensors)}",
    )
    p7 = fine.set_index("period_id").loc["P7"]
    add_check(
        checks,
        "p7_segment_constraint",
        p7["n_months"] == 15
        and not bool(p7["reliable_for_slope"])
        and bool(p7["changepoint_detection_exempt"]),
        "P7 must be 15 months, slope-unreliable, and changepoint-score exempt",
    )
    trend_config = load_trend_config(args.trend_config)
    add_check(
        checks,
        "trend_statistics_config",
        trend_config["seasonal_adjustment"]
        == "trend_preserving_calendar_month_anomaly"
        and float(trend_config["modified_mk"]["minimum_variance_factor"]) >= 1.0,
        str(args.trend_config),
    )

    coverage_path = args.diagnostics_dir / "monthly_synthesis_time_coverage.csv"
    inventory_path = args.diagnostics_dir / "monthly_synthesis_input_inventory.csv"
    variable_registry_path = args.diagnostics_dir / "monthly_synthesis_variable_registry.csv"
    for path in (coverage_path, inventory_path, variable_registry_path):
        add_check(checks, f"existing_product:{path.name}", path.exists(), str(path))
    if any(not path.exists() for path in (coverage_path, inventory_path, variable_registry_path)):
        return pd.DataFrame(checks), pd.DataFrame(manifest)

    coverage = pd.read_csv(coverage_path).set_index("dataset")
    inventory = pd.read_csv(inventory_path)
    variable_registry = pd.read_csv(variable_registry_path)
    selection = load_variable_selection(args.variable_selection)
    derived_variables = load_derived_variables(args.variable_selection)
    input_contract = json.loads(args.input_contract.read_text())
    dataset_contracts = input_contract["datasets"]
    expected_start = pd.Timestamp(input_contract["record_start"])
    expected_end = pd.Timestamp(input_contract["record_end"])
    expected_months = int(input_contract["n_months"])
    expected_tiles = int(input_contract["n_tiles"])
    tilecoord_path = args.data_dir / input_contract["tilecoord_file"]

    expected_time = pd.date_range(expected_start, expected_end, freq="MS")
    add_check(
        checks,
        "input_contract_time",
        len(expected_time) == expected_months,
        f"{expected_start.date()} to {expected_end.date()}, months={expected_months}",
    )
    add_check(checks, "tilecoord_exists", tilecoord_path.exists(), str(tilecoord_path))
    add_check(
        checks,
        "tile_area_units",
        input_contract.get("tile_area_units") == "km2",
        str(input_contract.get("tile_area_units", "missing")),
    )
    if tilecoord_path.exists():
        tile_area = read_tile_area(tilecoord_path)
        positive_area = tile_area[np.isfinite(tile_area) & (tile_area > 0)]
        tile_area_sum = float(positive_area.sum())
        add_check(
            checks,
            "tile_area_count",
            tile_area.size == expected_tiles,
            f"areas={tile_area.size}, expected={expected_tiles}",
        )
        add_check(
            checks,
            "tile_area_positive",
            positive_area.size == expected_tiles,
            f"positive finite areas={positive_area.size}",
        )
        area_min = float(input_contract["tile_area_sum_min"])
        area_max = float(input_contract["tile_area_sum_max"])
        add_check(
            checks,
            "tile_area_sum",
            area_min <= tile_area_sum <= area_max,
            f"sum={tile_area_sum:.6g} {input_contract['tile_area_units']}; "
            f"expected {area_min:.6g}-{area_max:.6g}",
        )
    actual_variables: dict[str, set[str]] = {}
    actual_units: dict[tuple[str, str], str] = {}
    actual_dims: dict[tuple[str, str], str] = {}

    for dataset, spec in dataset_contracts.items():
        filename = spec["filename"]
        source_token = spec["source_token"]
        path = args.data_dir / filename
        if not path.exists():
            add_check(checks, f"dataset_exists:{dataset}", False, str(path))
            continue

        stat = path.stat()
        with xr.open_dataset(path, decode_times=True) as ds:
            times = pd.DatetimeIndex(ds["time"].values)
            n_time = int(ds.sizes.get("time", 0))
            n_tile = int(ds.sizes.get("tile", 0))
            source_root = str(ds.attrs.get("source_root", ""))
            actual_variables[dataset] = set(ds.data_vars)
            for variable in ds.data_vars:
                actual_units[(dataset, variable)] = str(ds[variable].attrs.get("units", ""))
                actual_dims[(dataset, variable)] = ",".join(ds[variable].dims)

        manifest.append(
            {
                "dataset": dataset,
                "path": str(path),
                "size_bytes": stat.st_size,
                "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                "source_root": source_root,
                "n_time": n_time,
                "n_tile": n_tile,
                "start": times[0].date().isoformat() if len(times) else "",
                "end": times[-1].date().isoformat() if len(times) else "",
                "sha256": sha256sum(path) if args.sha256 else "not_computed",
            }
        )
        add_check(checks, f"dataset_exists:{dataset}", True, str(path))
        add_check(
            checks,
            f"dataset_shape:{dataset}",
            n_time == expected_months and n_tile == expected_tiles,
            f"time={n_time}, tile={n_tile}",
        )
        add_check(
            checks,
            f"dataset_time:{dataset}",
            times.equals(expected_time),
            f"{times[0].date() if len(times) else 'missing'} to {times[-1].date() if len(times) else 'missing'}",
        )
        add_check(
            checks,
            f"dataset_source:{dataset}",
            source_token in source_root,
            source_root or "missing source_root",
        )

        if dataset in coverage.index:
            row = coverage.loc[dataset]
            coverage_ok = (
                int(row["n_time"]) == n_time
                and int(row["n_tile"]) == n_tile
                and pd.Timestamp(row["start"]) == times[0]
                and pd.Timestamp(row["end"]) == times[-1]
            )
            add_check(checks, f"coverage_matches:{dataset}", coverage_ok, str(row.to_dict()))
        else:
            add_check(checks, f"coverage_matches:{dataset}", False, "dataset absent from coverage table")

    inventory_index = inventory.set_index(["file", "variable"])
    registry_index = variable_registry.set_index(["source", "variable"])
    for selected in selection.itertuples(index=False):
        selected_dataset_keys = [
            name for name in (selected.ol_dataset, selected.da_dataset) if name
        ]
        for dataset in selected_dataset_keys:
            in_contract = dataset in dataset_contracts
            add_check(
                checks,
                f"input_contract:{selected.analysis_variable}:{dataset}",
                in_contract,
                dataset,
            )
            if not in_contract:
                continue
            definition = derived_variables.get(selected.analysis_variable)
            variables = (
                definition["source_variables"]
                if definition is not None
                else [selected.source_variable]
            )
            for variable in variables:
                suffix = f":{variable}" if len(variables) > 1 else ""
                exists = dataset in actual_variables and variable in actual_variables[dataset]
                add_check(
                    checks,
                    f"selected_variable:{selected.analysis_variable}:{dataset}{suffix}",
                    exists,
                    f"{dataset}/{variable}",
                )
                if not exists:
                    continue

                filename = dataset_contracts[dataset]["filename"]
                inventory_key = (filename, variable)
                in_inventory = inventory_key in inventory_index.index
                detail = "missing from input inventory"
                if in_inventory:
                    inv = inventory_index.loc[inventory_key]
                    dims_match = str(inv["dims"]) == actual_dims[(dataset, variable)]
                    units_match = str(inv["units"]) == actual_units[(dataset, variable)]
                    in_inventory = bool(dims_match and units_match)
                    detail = (
                        f"dims={actual_dims[(dataset, variable)]}, "
                        f"units={actual_units[(dataset, variable)]}"
                    )
                add_check(
                    checks,
                    f"inventory_contract:{selected.analysis_variable}:{dataset}{suffix}",
                    in_inventory,
                    detail,
                )

                registry_key = (dataset_contracts[dataset]["family"], variable)
                add_check(
                    checks,
                    f"variable_registry:{selected.analysis_variable}:{dataset}{suffix}",
                    registry_key in registry_index.index,
                    f"source={registry_key[0]}, variable={variable}",
                )

    selected_datasets = set(selection["ol_dataset"]) | set(selection["da_dataset"])
    legacy_variables = find_legacy_increment_variables(selection["source_variable"])
    add_check(
        checks,
        "legacy_increment_excluded",
        "monthly_increment_legacy" not in selected_datasets
        and not legacy_variables,
        f"legacy suffix matches={sorted(legacy_variables)}",
    )

    return pd.DataFrame(checks), pd.DataFrame(manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--diagnostics-dir", type=Path, default=DEFAULT_DIAGNOSTICS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--period-registry",
        type=Path,
        default=PROJECT_ROOT / "config" / "observing_system_registry.json",
    )
    parser.add_argument("--variable-selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument("--trend-config", type=Path, default=DEFAULT_TREND_CONFIG)
    parser.add_argument("--sha256", action="store_true", help="Hash all monthly input files")
    parser.add_argument("--no-write", action="store_true", help="Validate without writing reports")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    checks, manifest = audit(args)
    if not args.no_write:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        checks.to_csv(args.output_dir / "trend_breakpoint_input_audit.csv", index=False)
        manifest.to_csv(args.output_dir / "trend_breakpoint_input_manifest.csv", index=False)

    failed = checks[checks["status"] == "FAIL"]
    print(f"Audit checks: {len(checks) - len(failed)} passed, {len(failed)} failed")
    if len(failed):
        print(failed.to_string(index=False))
        return 1
    print(f"Validated {len(manifest)} monthly datasets and the shared period registry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
