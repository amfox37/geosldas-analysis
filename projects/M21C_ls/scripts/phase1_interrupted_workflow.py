#!/usr/bin/env python3
"""Shared assembly and audit helpers for M21C Phase 1 interrupted series."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from interrupted_time_series import load_interrupted_config
from phase1_trend_workflow import DEFAULT_OUTPUT_DIR, DEFAULT_RUN_MATRIX, load_trend_runs
from trend_breakpoint_series import DEFAULT_VARIABLE_SELECTION
from trend_statistics import benjamini_hochberg


DEFAULT_COEFFICIENTS = DEFAULT_OUTPUT_DIR / "phase1_interrupted_series_coefficients.csv"
DEFAULT_MODELS = DEFAULT_OUTPUT_DIR / "phase1_interrupted_series_models.csv"
DEFAULT_MONTHLY = DEFAULT_OUTPUT_DIR / "phase1_interrupted_series_monthly.nc"
DEFAULT_AUDIT = DEFAULT_OUTPUT_DIR / "phase1_interrupted_series_audit.csv"

SERIES_METADATA_COLUMNS = [
    "series_id",
    "run_id",
    "variable",
    "source_series",
    "mask",
    "role",
    "series_role",
    "units",
    "rationale",
]


def phase1_series_catalogue(
    run_matrix: str | Path = DEFAULT_RUN_MATRIX,
    *,
    variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
) -> pd.DataFrame:
    """Expand matrix rows into the OL, DA, delta, and DA-only domain series."""

    _, runs = load_trend_runs(run_matrix, variable_selection=variable_selection)
    rows: list[dict[str, str]] = []
    for run in runs:
        source_series = ["ol", "da", "delta"] if run["series"] == "delta" else ["value"]
        for source in source_series:
            rows.append(
                {
                    "series_id": f"{run['run_id']}__{source}",
                    "run_id": run["run_id"],
                    "variable": run["variable"],
                    "source_series": source,
                    "mask": run["mask"],
                    "role": run["role"],
                    "series_role": (
                        "primary_estimand" if source == run["series"] else "paired_control"
                    ),
                    "units": "",
                    "rationale": run["rationale"],
                }
            )
    catalogue = pd.DataFrame(rows, columns=SERIES_METADATA_COLUMNS)
    if catalogue["series_id"].duplicated().any():
        raise ValueError("Expanded interrupted-series catalogue contains duplicate IDs")
    return catalogue


def assign_transition_fdr(
    coefficients: pd.DataFrame,
    *,
    alpha: float,
) -> pd.DataFrame:
    """Apply BH FDR separately to each scientifically comparable boundary."""

    output = coefficients.copy()
    output["fdr_family"] = ""
    independent = pd.Series(True, index=output.index, dtype=bool)
    period_slope = output["coefficient_type"] == "period_slope"
    if "independently_estimated" in output:
        independent.loc[period_slope] = (
            output.loc[period_slope, "independently_estimated"] == True  # noqa: E712
        )
    eligible = output["coefficient_type"].isin(
        ["level_change", "slope_change", "period_slope"]
    ) & independent
    output.loc[eligible, "fdr_family"] = output.loc[eligible, "coefficient"]
    output["q_value"] = np.nan
    output["significant_fdr"] = False
    for family, indices in output.loc[eligible].groupby("fdr_family").groups.items():
        if not family:
            continue
        q_values, significant = benjamini_hochberg(
            output.loc[indices, "p_value"].to_numpy(dtype="float64"), alpha=alpha
        )
        output.loc[indices, "q_value"] = q_values
        output.loc[indices, "significant_fdr"] = significant
    return output


def expected_fdr_family_size(catalogue: pd.DataFrame) -> int:
    """Return the number of domain series in every transition-level FDR family."""

    return int(len(catalogue))


def audit_interrupted_outputs(
    *,
    coefficients_path: str | Path = DEFAULT_COEFFICIENTS,
    models_path: str | Path = DEFAULT_MODELS,
    monthly_path: str | Path = DEFAULT_MONTHLY,
    run_matrix: str | Path = DEFAULT_RUN_MATRIX,
    variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
    config_path: str | Path,
) -> pd.DataFrame:
    """Audit the complete three-file interrupted-series product."""

    catalogue = phase1_series_catalogue(
        run_matrix, variable_selection=variable_selection
    )
    settings = load_interrupted_config(config_path)
    errors: list[str] = []
    coefficient_path = Path(coefficients_path)
    model_path = Path(models_path)
    monthly_file = Path(monthly_path)
    for path in (coefficient_path, model_path, monthly_file):
        if not path.exists():
            errors.append(f"missing {path}")
    if errors:
        return pd.DataFrame([{"status": "FAIL", "detail": "; ".join(errors)}])

    coefficients = pd.read_csv(coefficient_path)
    models = pd.read_csv(model_path)
    expected_ids = set(catalogue["series_id"])
    if set(models["series_id"]) != expected_ids:
        errors.append("model series IDs differ from the expanded run matrix")
    if models["series_id"].duplicated().any():
        errors.append("duplicate model series IDs")
    if not np.all(models["design_rank"] == 28):
        errors.append("not every model has design rank 28")
    if not np.all(models["transformed_design_rank"] == 28):
        errors.append("not every transformed model has design rank 28")
    if not np.all(models["bootstrap_replicates"] == int(settings["bootstrap"]["replicates"])):
        errors.append("bootstrap replicate count differs from configuration")

    expected_levels = {f"level_change_P{index}" for index in range(2, 10)}
    expected_slopes = {
        f"slope_change_P{index}" for index in (2, 3, 4, 5, 6, 8, 9)
    }
    expected_periods = {f"period_slope_P{index}" for index in range(1, 10)}
    if "estimate_units" not in coefficients:
        errors.append("coefficient table lacks estimate_units")
    else:
        slope_rows = coefficients["coefficient_type"].isin(
            ["baseline_slope", "slope_change", "period_slope"]
        )
        expected_units = np.where(
            slope_rows,
            coefficients["units"].astype(str) + " yr-1",
            coefficients["units"].astype(str),
        )
        if not np.array_equal(coefficients["estimate_units"].astype(str), expected_units):
            errors.append("estimate_units does not distinguish levels from annual slopes")
    for series_id, frame in coefficients.groupby("series_id"):
        names = set(frame["coefficient"])
        if not expected_levels.issubset(names):
            errors.append(f"{series_id} lacks one or more level changes")
        if not expected_slopes.issubset(names):
            errors.append(f"{series_id} lacks one or more slope changes")
        if "slope_change_P7" in names:
            errors.append(f"{series_id} incorrectly estimates a P7 slope change")
        if not expected_periods.issubset(names):
            errors.append(f"{series_id} lacks one or more period slopes")

    fdr_rows = coefficients[coefficients["fdr_family"].fillna("") != ""]
    if fdr_rows.empty:
        errors.append("no FDR families were written")
    else:
        family_sizes = fdr_rows.groupby("fdr_family").size()
        if not np.all(family_sizes == expected_fdr_family_size(catalogue)):
            errors.append("one or more FDR families do not contain every domain series")
        q = fdr_rows["q_value"].to_numpy(dtype="float64")
        if not np.all(np.isfinite(q) & (q >= 0) & (q <= 1)):
            errors.append("FDR q-values are incomplete or outside [0, 1]")
        expected_significant = q <= float(settings["fdr"]["alpha"])
        actual_significant = fdr_rows["significant_fdr"]
        if actual_significant.dtype != bool:
            actual_significant = actual_significant.astype(str).str.lower().map(
                {"true": True, "false": False}
            )
        if not np.array_equal(
            actual_significant.to_numpy(dtype=bool), expected_significant
        ):
            errors.append("significant_fdr is inconsistent with q_value")

    with xr.open_dataset(monthly_file) as monthly:
        if monthly.sizes.get("series") != len(catalogue):
            errors.append("monthly file series count differs from the run matrix")
        if monthly.sizes.get("time") != 288:
            errors.append("monthly file does not contain 288 months")
        if set(monthly["series_id"].values.astype(str)) != expected_ids:
            errors.append("monthly file series IDs differ from the run matrix")
        if str(monthly.attrs.get("configuration", "")) != json.dumps(
            settings, sort_keys=True
        ):
            errors.append("monthly embedded configuration differs from audited config")
        source_commit = str(monthly.attrs.get("source_git_commit", ""))
        if not source_commit or source_commit == "unknown":
            errors.append("monthly file lacks source git commit")
        observed = monthly["observed"]
        for run_id, group in catalogue.groupby("run_id"):
            sources = set(group["source_series"])
            if sources != {"ol", "da", "delta"}:
                continue
            ol = observed.sel(series=f"{run_id}__ol").values
            da = observed.sel(series=f"{run_id}__da").values
            delta = observed.sel(series=f"{run_id}__delta").values
            finite = np.isfinite(ol) & np.isfinite(da) & np.isfinite(delta)
            if not np.array_equal(np.isfinite(ol), finite) or not np.array_equal(
                np.isfinite(da), finite
            ):
                errors.append(f"{run_id} paired domain means have unequal support")
            if not np.allclose(delta[finite], da[finite] - ol[finite], atol=1e-10):
                errors.append(f"{run_id} violates delta = DA - OL")

    return pd.DataFrame(
        [
            {
                "status": "PASS" if not errors else "FAIL",
                "detail": "; ".join(dict.fromkeys(errors)),
                "n_series": len(catalogue),
                "n_models": len(models),
                "n_coefficient_rows": len(coefficients),
                "n_fdr_families": int(fdr_rows["fdr_family"].nunique()),
                "n_significant_fdr": int(expected_significant.sum()),
            }
        ]
    )
