#!/usr/bin/env python3
"""Shared contracts and audits for M21C Phase 1 trend production."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from trend_breakpoint_series import (
    DEFAULT_INPUT_CONTRACT,
    DEFAULT_VARIABLE_SELECTION,
    load_variable_selection,
)
from trend_statistics import DEFAULT_CONFIG, load_trend_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_MATRIX = PROJECT_ROOT / "config" / "phase1_trend_runs.json"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output" / "trends_breakpoints"
ALLOWED_MASKS = {
    "valid_land",
    "snow_possible",
    "seasonal_snow",
    "warm_static",
    "warm_snowfree_monthly",
    "locally_snowy_monthly",
}
REQUIRED_OUTPUT_VARIABLES = {
    "slope",
    "intercept",
    "slope_ci_low",
    "slope_ci_high",
    "slope_ci_low_nominal",
    "slope_ci_high_nominal",
    "p_value",
    "p_value_original_mk",
    "q_value",
    "significant_fdr",
    "ci_excludes_zero",
    "fdr_ci_disagreement",
    "mk_variance_factor",
    "lag1_residual_pearson_autocorrelation",
    "lag1_rank_autocorrelation",
    "n_positive_autocorrelation_lags",
    "n_eligible",
    "n_valid",
    "n_trend",
    "valid_fraction",
    "span_years",
    "status",
    "calendar_month_used",
    "n_calendar_months_used",
}


def load_phase1_runs(
    path: str | Path = DEFAULT_RUN_MATRIX,
    *,
    variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Load and validate the production run matrix against the variable registry."""

    payload = json.loads(Path(path).read_text())
    if int(payload.get("phase", -1)) != 1:
        raise ValueError("Phase 1 run matrix must declare phase=1")
    runs = [dict(run) for run in payload.get("runs", [])]
    if len(runs) != int(payload.get("expected_run_count", -1)):
        raise ValueError("Phase 1 run count does not match expected_run_count")

    selection = load_variable_selection(variable_selection)
    phase1 = selection.loc[selection["phase"] == 1]
    run_ids: set[str] = set()
    signatures: set[tuple[str, str, str]] = set()
    primary_variables: list[str] = []
    for run in runs:
        missing = {
            key
            for key in ("run_id", "variable", "series", "mask", "role", "rationale")
            if not str(run.get(key, "")).strip()
        }
        if missing:
            raise ValueError(f"Phase 1 run is missing fields: {sorted(missing)}")
        run_id = str(run["run_id"])
        variable = str(run["variable"])
        series = str(run["series"])
        mask = str(run["mask"])
        role = str(run["role"])
        if run_id in run_ids:
            raise ValueError(f"Duplicate Phase 1 run_id: {run_id}")
        run_ids.add(run_id)
        signature = (variable, series, mask)
        if signature in signatures:
            raise ValueError(f"Duplicate Phase 1 run signature: {signature}")
        signatures.add(signature)
        if variable not in phase1.index:
            raise ValueError(f"Run {run_id} references a non-Phase-1 variable: {variable}")
        row = phase1.loc[variable]
        expected_series = "delta" if str(row["ol_dataset"]) else "value"
        if series != expected_series:
            raise ValueError(
                f"Run {run_id} uses series={series}; expected {expected_series} for {variable}"
            )
        if mask not in ALLOWED_MASKS:
            raise ValueError(f"Run {run_id} uses unknown mask: {mask}")
        if role not in {"primary", "sensitivity"}:
            raise ValueError(f"Run {run_id} has invalid role: {role}")
        if role == "primary":
            primary_variables.append(variable)

    expected_primary = int(payload.get("expected_primary_variable_count", -1))
    if len(primary_variables) != expected_primary:
        raise ValueError("Primary run count does not match expected_primary_variable_count")
    duplicate_primary = pd.Series(primary_variables).duplicated().any()
    if duplicate_primary or set(primary_variables) != set(phase1.index):
        raise ValueError("Every Phase 1 variable must have exactly one primary run")
    return payload, runs


def output_filename(run: dict[str, str]) -> str:
    """Return the stable production filename for one matrix row."""

    return f"{run['variable']}_{run['series']}_{run['mask']}_trend_statistics.nc"


def audit_trend_output(
    path: str | Path,
    run: dict[str, str],
    *,
    expected_tiles: int,
    expected_total_input_tiles: int,
    require_global_fdr: bool,
    expected_configuration: str | None = None,
) -> dict[str, Any]:
    """Validate one trend NetCDF and return a compact machine-readable summary."""

    path = Path(path)
    summary: dict[str, Any] = {
        "run_id": run["run_id"],
        "variable": run["variable"],
        "series": run["series"],
        "mask": run["mask"],
        "role": run["role"],
        "path": str(path),
        "status": "FAIL",
        "detail": "",
        "n_tiles": 0,
        "n_success": 0,
        "n_significant_fdr": 0,
        "n_fdr_ci_disagreement": 0,
        "fraction_significant_with_ci_disagreement": np.nan,
        "fraction_significant_fdr": np.nan,
        "median_slope_success": np.nan,
        "median_valid_fraction_success": np.nan,
        "median_mk_variance_factor_success": np.nan,
        "area_success_km2": np.nan,
        "area_significant_fdr_km2": np.nan,
        "fraction_area_significant_fdr": np.nan,
        "n_significant_positive_slope": 0,
        "n_significant_negative_slope": 0,
        "n_significant_zero_slope": 0,
        "slope_success_p05": np.nan,
        "slope_success_p95": np.nan,
        "minimum_calendar_months_used_success": -1,
        "maximum_calendar_months_used_success": -1,
        "source_git_commit": "",
        "parallel_execution": "",
        **{f"n_status_{code}": 0 for code in range(5)},
    }
    if not path.exists():
        summary["detail"] = "missing output"
        return summary

    errors: list[str] = []
    try:
        with xr.open_dataset(path) as dataset:
            n_tiles = int(dataset.sizes.get("tile", -1))
            summary["n_tiles"] = n_tiles
            if n_tiles != expected_tiles:
                errors.append(f"tile count {n_tiles}, expected {expected_tiles}")
            missing = REQUIRED_OUTPUT_VARIABLES.difference(dataset.data_vars)
            if missing:
                errors.append(f"missing variables {sorted(missing)}")
            for coordinate in ("lat", "lon", "tile_area"):
                if coordinate not in dataset.coords:
                    errors.append(f"missing coordinate {coordinate}")
            expected_attrs = {
                "analysis_variable": run["variable"],
                "source_series": run["series"],
                "mask": run["mask"],
                "run_id": run["run_id"],
                "run_role": run["role"],
            }
            for name, expected in expected_attrs.items():
                actual = str(dataset.attrs.get(name, ""))
                if actual != expected:
                    errors.append(f"attribute {name}={actual!r}, expected {expected!r}")
            source_git_commit = str(dataset.attrs.get("source_git_commit", ""))
            summary["source_git_commit"] = source_git_commit
            if not source_git_commit or source_git_commit == "unknown":
                errors.append("missing source_git_commit provenance")
            summary["parallel_execution"] = str(
                dataset.attrs.get("parallel_execution", "")
            )
            if expected_configuration is not None:
                actual_configuration = str(dataset.attrs.get("configuration", ""))
                if actual_configuration != expected_configuration:
                    errors.append("embedded trend configuration differs from audited config")
            total_input_tiles = int(dataset.attrs.get("total_input_tiles", -1))
            if total_input_tiles != expected_total_input_tiles:
                errors.append(
                    f"total_input_tiles {total_input_tiles}, expected {expected_total_input_tiles}"
                )
            fdr_scope = str(dataset.attrs.get("fdr_scope", ""))
            if require_global_fdr and not fdr_scope.startswith("all finite tiles"):
                errors.append(f"non-global FDR scope: {fdr_scope!r}")
            if not require_global_fdr and "not global FDR" not in fdr_scope:
                errors.append(f"diagnostic subset lacks non-global FDR label: {fdr_scope!r}")

            if REQUIRED_OUTPUT_VARIABLES.issubset(dataset.data_vars):
                status = np.asarray(dataset["status"].values)
                valid_status = np.isin(status, np.arange(5))
                if not np.all(valid_status):
                    errors.append("status contains values outside 0..4")
                success = status == 0
                for code in range(5):
                    summary[f"n_status_{code}"] = int((status == code).sum())
                n_success = int(success.sum())
                significant = np.asarray(dataset["significant_fdr"].values, dtype=bool)
                n_significant = int(significant.sum())
                summary.update(
                    n_success=n_success,
                    n_significant_fdr=n_significant,
                    fraction_significant_fdr=(
                        n_significant / n_success if n_success else np.nan
                    ),
                )
                p_value = np.asarray(dataset["p_value"].values, dtype="float64")
                q_value = np.asarray(dataset["q_value"].values, dtype="float64")
                if not np.array_equal(np.isfinite(p_value), success):
                    errors.append("finite p-value support differs from status==success")
                if not np.array_equal(np.isfinite(q_value), success):
                    errors.append("finite q-value support differs from status==success")
                if np.any((p_value[success] < 0) | (p_value[success] > 1)):
                    errors.append("p_value outside [0, 1]")
                if np.any((q_value[success] < 0) | (q_value[success] > 1)):
                    errors.append("q_value outside [0, 1]")
                fdr_alpha = float(dataset["significant_fdr"].attrs.get("alpha", 0.05))
                if not np.array_equal(significant, np.isfinite(q_value) & (q_value <= fdr_alpha)):
                    errors.append("significant_fdr is inconsistent with q_value and alpha")
                ci_excludes = np.asarray(dataset["ci_excludes_zero"].values, dtype=bool)
                disagreement = np.asarray(
                    dataset["fdr_ci_disagreement"].values, dtype=bool
                )
                expected_disagreement = significant & ~ci_excludes
                if not np.array_equal(disagreement, expected_disagreement):
                    errors.append("fdr_ci_disagreement does not match FDR and CI fields")
                summary["n_fdr_ci_disagreement"] = int(disagreement.sum())
                summary["fraction_significant_with_ci_disagreement"] = (
                    int(disagreement.sum()) / n_significant
                    if n_significant
                    else 0.0
                )
                n_month = np.asarray(dataset["n_calendar_months_used"].values)
                if np.any((n_month < 0) | (n_month > 12)):
                    errors.append("n_calendar_months_used outside 0..12")
                if n_success:
                    slope = np.asarray(dataset["slope"].values, dtype="float64")
                    valid_fraction = np.asarray(
                        dataset["valid_fraction"].values, dtype="float64"
                    )
                    variance_factor = np.asarray(
                        dataset["mk_variance_factor"].values, dtype="float64"
                    )
                    tile_area = np.asarray(dataset["tile_area"].values, dtype="float64")
                    area_success = float(np.nansum(tile_area[success]))
                    area_significant = float(np.nansum(tile_area[significant]))
                    significant_positive = significant & (slope > 0)
                    significant_negative = significant & (slope < 0)
                    significant_zero = significant & (slope == 0)
                    summary.update(
                        median_slope_success=float(np.nanmedian(slope[success])),
                        median_valid_fraction_success=float(
                            np.nanmedian(valid_fraction[success])
                        ),
                        median_mk_variance_factor_success=float(
                            np.nanmedian(variance_factor[success])
                        ),
                        area_success_km2=area_success,
                        area_significant_fdr_km2=area_significant,
                        fraction_area_significant_fdr=(
                            area_significant / area_success if area_success > 0 else np.nan
                        ),
                        n_significant_positive_slope=int(significant_positive.sum()),
                        n_significant_negative_slope=int(significant_negative.sum()),
                        n_significant_zero_slope=int(significant_zero.sum()),
                        slope_success_p05=float(np.nanquantile(slope[success], 0.05)),
                        slope_success_p95=float(np.nanquantile(slope[success], 0.95)),
                        minimum_calendar_months_used_success=int(n_month[success].min()),
                        maximum_calendar_months_used_success=int(n_month[success].max()),
                    )
    except Exception as exc:
        errors.append(f"cannot read output: {type(exc).__name__}: {exc}")

    if errors:
        summary["detail"] = "; ".join(errors)
    else:
        summary["status"] = "PASS"
        summary["detail"] = "complete global field" if require_global_fdr else "diagnostic subset"
    return summary


def audit_phase1_outputs(
    *,
    run_matrix: str | Path = DEFAULT_RUN_MATRIX,
    variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
    input_contract: str | Path = DEFAULT_INPUT_CONTRACT,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    trend_config: str | Path = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Audit every production output declared by the Phase 1 run matrix."""

    _, runs = load_phase1_runs(run_matrix, variable_selection=variable_selection)
    contract = json.loads(Path(input_contract).read_text())
    expected_tiles = int(contract["n_tiles"])
    expected_configuration = json.dumps(
        load_trend_config(trend_config), sort_keys=True
    )
    output_dir = Path(output_dir)
    rows = [
        audit_trend_output(
            output_dir / output_filename(run),
            run,
            expected_tiles=expected_tiles,
            expected_total_input_tiles=expected_tiles,
            require_global_fdr=True,
            expected_configuration=expected_configuration,
        )
        for run in runs
    ]
    report = pd.DataFrame(rows)
    commits = report.loc[report["status"] == "PASS", "source_git_commit"].unique()
    if commits.size > 1:
        report.loc[:, "status"] = "FAIL"
        report.loc[:, "detail"] = report["detail"] + "; mixed source_git_commit values"
    return report
