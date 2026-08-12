#!/usr/bin/env python3
"""Audit the existing M21C monthly products before trend analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import xarray as xr

from m21c_periods import load_period_frames


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = Path(
    "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2"
)
DEFAULT_DIAGNOSTICS_DIR = PROJECT_ROOT / "output" / "monthly_synthesis_diagnostics"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output" / "trends_breakpoints"
DEFAULT_SELECTION = PROJECT_ROOT / "config" / "trend_breakpoint_variable_selection.json"

EXPECTED_START = pd.Timestamp("2000-06-01")
EXPECTED_END = pd.Timestamp("2024-05-01")
EXPECTED_MONTHS = 288
EXPECTED_TILES = 112_573

DATASETS = {
    "ol_land": ("OLv8_land_variables_2000_2024_compressed.nc", "LS_OLv8_M36_v2"),
    "da_land": ("DAv8_land_variables_2000_2024_compressed.nc", "LS_DAv8_M36_v3"),
    "catch_raw_cumulative": (
        "catch_progn_raw_monthly_cumulative_200006_202405.nc",
        "LS_DAv8_M36_v3",
    ),
    "inst3_raw_diagnostic": (
        "inst3_fcstana_raw_monthly_diagnostic_200006_202405.nc",
        "LS_DAv8_M36_v3",
    ),
    "ol_flux_core": ("OLv8_flux_core_2000_2024_compressed.nc", "LS_OLv8_M36_v2"),
    "da_flux_core": ("DAv8_flux_core_2000_2024_compressed.nc", "LS_DAv8_M36_v3"),
    "ol_latent_components": (
        "OLv8_latent_components_2000_2024_compressed.nc",
        "LS_OLv8_M36_v2",
    ),
    "da_latent_components": (
        "DAv8_latent_components_2000_2024_compressed.nc",
        "LS_DAv8_M36_v3",
    ),
    "ol_water_budget": ("OLv8_water_budget_2000_2024_compressed.nc", "LS_OLv8_M36_v2"),
    "da_water_budget": ("DAv8_water_budget_2000_2024_compressed.nc", "LS_DAv8_M36_v3"),
    "ol_energy_context": ("OLv8_energy_context_2000_2024_compressed.nc", "LS_OLv8_M36_v2"),
    "da_energy_context": ("DAv8_energy_context_2000_2024_compressed.nc", "LS_DAv8_M36_v3"),
}


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_selection(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text())
    frame = pd.DataFrame(payload["variables"], columns=payload["columns"])
    if frame["analysis_variable"].duplicated().any():
        raise ValueError("Duplicate analysis_variable entries in variable selection")
    return frame


def source_family(dataset: str) -> str:
    if dataset.startswith("ol_") or dataset.startswith("da_"):
        return dataset.split("_", 1)[1]
    return dataset


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
    selection = load_selection(args.variable_selection)

    expected_time = pd.date_range(EXPECTED_START, EXPECTED_END, freq="MS")
    actual_variables: dict[str, set[str]] = {}
    actual_units: dict[tuple[str, str], str] = {}
    actual_dims: dict[tuple[str, str], str] = {}

    for dataset, (filename, source_token) in DATASETS.items():
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
            n_time == EXPECTED_MONTHS and n_tile == EXPECTED_TILES,
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
        datasets = [name for name in (selected.ol_dataset, selected.da_dataset) if name]
        for dataset in datasets:
            variable = selected.source_variable
            exists = dataset in actual_variables and variable in actual_variables[dataset]
            add_check(
                checks,
                f"selected_variable:{selected.analysis_variable}:{dataset}",
                exists,
                f"{dataset}/{variable}",
            )
            if not exists:
                continue

            filename = DATASETS[dataset][0]
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
                f"inventory_contract:{selected.analysis_variable}:{dataset}",
                in_inventory,
                detail,
            )

            registry_key = (source_family(dataset), variable)
            add_check(
                checks,
                f"variable_registry:{selected.analysis_variable}:{dataset}",
                registry_key in registry_index.index,
                f"source={registry_key[0]}, variable={variable}",
            )

    selected_datasets = set(selection["ol_dataset"]) | set(selection["da_dataset"])
    add_check(
        checks,
        "legacy_increment_excluded",
        "monthly_increment_legacy" not in selected_datasets
        and not selection["source_variable"].isin(["SFMC_INC", "RZMC_INC", "SNOWMASS_INCR"]).any(),
        "legacy LS_monthly_increments variables must not be selected",
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
