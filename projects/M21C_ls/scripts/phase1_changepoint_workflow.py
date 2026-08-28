#!/usr/bin/env python3
"""Shared output contracts and audit for M21C Phase 1 changepoints."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from changepoint_detection import load_changepoint_config
from phase1_interrupted_workflow import DEFAULT_MONTHLY, phase1_series_catalogue
from phase1_trend_workflow import DEFAULT_OUTPUT_DIR, DEFAULT_RUN_MATRIX
from trend_breakpoint_series import DEFAULT_VARIABLE_SELECTION


DEFAULT_DETECTIONS = DEFAULT_OUTPUT_DIR / "phase1_changepoint_detections.csv"
DEFAULT_PENALTY_GRID = DEFAULT_OUTPUT_DIR / "phase1_changepoint_penalty_grid.csv"
DEFAULT_COMPARISON = DEFAULT_OUTPUT_DIR / "phase1_changepoint_boundary_comparison.csv"
DEFAULT_MODELS = DEFAULT_OUTPUT_DIR / "phase1_changepoint_models.csv"
DEFAULT_MONTHLY_OUTPUT = DEFAULT_OUTPUT_DIR / "phase1_changepoint_monthly.nc"
DEFAULT_AUDIT = DEFAULT_OUTPUT_DIR / "phase1_changepoint_audit.csv"


def audit_changepoint_outputs(
    *,
    input_monthly: str | Path = DEFAULT_MONTHLY,
    detections_path: str | Path = DEFAULT_DETECTIONS,
    penalty_grid_path: str | Path = DEFAULT_PENALTY_GRID,
    comparison_path: str | Path = DEFAULT_COMPARISON,
    models_path: str | Path = DEFAULT_MODELS,
    monthly_path: str | Path = DEFAULT_MONTHLY_OUTPUT,
    run_matrix: str | Path = DEFAULT_RUN_MATRIX,
    variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
    config_path: str | Path,
) -> pd.DataFrame:
    """Audit the complete independently detected changepoint product."""

    settings = load_changepoint_config(config_path)
    catalogue = phase1_series_catalogue(
        run_matrix, variable_selection=variable_selection
    )
    expected_ids = set(catalogue["series_id"])
    paths = [
        Path(input_monthly),
        Path(detections_path),
        Path(penalty_grid_path),
        Path(comparison_path),
        Path(models_path),
        Path(monthly_path),
    ]
    errors = [f"missing {path}" for path in paths if not path.exists()]
    if errors:
        return pd.DataFrame([{"status": "FAIL", "detail": "; ".join(errors)}])

    detections = pd.read_csv(detections_path)
    penalty_grid = pd.read_csv(penalty_grid_path)
    comparison = pd.read_csv(comparison_path)
    models = pd.read_csv(models_path)
    if set(models["series_id"]) != expected_ids or models["series_id"].duplicated().any():
        errors.append("model series IDs differ from the run-matrix catalogue")
    if not np.all(models["n_time"] == 288) or not np.all(models["n_valid"] == 288):
        errors.append("not every model contains 288 finite months")
    expected_configuration = json.dumps(settings, sort_keys=True)
    if not np.all(models["configuration"] == expected_configuration):
        errors.append("model configuration differs from the audited config")

    accepted = detections[detections["accepted_detection"].astype(bool)]
    if not accepted.empty:
        if not np.all(accepted["method"] == settings["methods"]["primary"]):
            errors.append("accepted detections include a non-primary method")
        if not np.all(accepted["method_consensus"].astype(bool)):
            errors.append("accepted detections lack cross-method consensus")
        if not np.all(
            accepted["penalty_stability_fraction"]
            >= float(settings["acceptance"]["minimum_penalty_stability_fraction"])
        ):
            errors.append("accepted detections fail the penalty-stability floor")
        floor = int(settings["minimum_segment_months"])
        if not np.all(
            (accepted["break_index"] >= floor)
            & (accepted["break_index"] <= 288 - floor)
        ):
            errors.append("accepted detection violates the segment-length floor")

    allowed_methods = set(settings["methods"].values()).intersection(
        {"pelt_ar1_linear", "pelt_prewhitened_linear"}
    )
    allowed_multipliers = set(map(float, settings["penalty"]["multipliers"]))
    if not set(penalty_grid["method"]).issubset(allowed_methods):
        errors.append("penalty grid contains an unknown method")
    if not set(penalty_grid["penalty_multiplier"]).issubset(allowed_multipliers):
        errors.append("penalty grid contains an unknown multiplier")

    if len(comparison) != len(catalogue) * 8:
        errors.append("boundary comparison does not contain 8 rows per series")
    if set(comparison["series_id"]) != expected_ids:
        errors.append("boundary comparison series IDs differ from catalogue")
    if set(comparison["period_id"]) != {f"P{index}" for index in range(2, 10)}:
        errors.append("boundary comparison does not cover P2-P9")
    scored_by_period = comparison.groupby("period_id")["scored"].unique().to_dict()
    for period_id in (f"P{index}" for index in range(2, 10)):
        expected = period_id != "P7"
        values = scored_by_period.get(period_id, [])
        if len(values) != 1 or bool(values[0]) != expected:
            errors.append(f"{period_id} scoring exemption is incorrect")
    primary_match = comparison["matched_within_primary_tolerance"].astype(bool)
    sensitivity_match = comparison["matched_within_sensitivity_tolerance"].astype(bool)
    if not np.all(~primary_match | sensitivity_match):
        errors.append("a primary boundary match is absent from sensitivity matches")
    primary_tolerance = int(
        settings["known_boundary_comparison"]["primary_tolerance_months"]
    )
    sensitivity_tolerance = int(
        settings["known_boundary_comparison"]["sensitivity_tolerance_months"]
    )
    distance = comparison["absolute_offset_months"].to_numpy(dtype="float64")
    if not np.array_equal(primary_match.to_numpy(), np.isfinite(distance) & (distance <= primary_tolerance)):
        errors.append("primary boundary-match flags disagree with offsets")
    if not np.array_equal(
        sensitivity_match.to_numpy(), np.isfinite(distance) & (distance <= sensitivity_tolerance)
    ):
        errors.append("sensitivity boundary-match flags disagree with offsets")

    with xr.open_dataset(input_monthly) as source, xr.open_dataset(monthly_path) as monthly:
        if monthly.sizes.get("series") != len(catalogue) or monthly.sizes.get("time") != 288:
            errors.append("changepoint monthly output has incorrect dimensions")
        if set(monthly["series_id"].values.astype(str)) != expected_ids:
            errors.append("changepoint monthly series IDs differ from catalogue")
        if not np.array_equal(monthly.time.values, source.time.values):
            errors.append("changepoint monthly time differs from input")
        if not np.allclose(monthly["observed"].values, source["observed"].values, equal_nan=True):
            errors.append("changepoint observed series differ from interrupted-series input")
        if str(monthly.attrs.get("configuration", "")) != expected_configuration:
            errors.append("monthly configuration differs from audited config")
        source_commit = str(monthly.attrs.get("source_git_commit", ""))
        if not source_commit or source_commit == "unknown":
            errors.append("monthly output lacks source git provenance")
        if str(monthly.attrs.get("input_source_git_commit", "")) != str(
            source.attrs.get("source_git_commit", "")
        ):
            errors.append("monthly input provenance differs from source file")

    return pd.DataFrame(
        [
            {
                "status": "PASS" if not errors else "FAIL",
                "detail": "; ".join(dict.fromkeys(errors)),
                "n_series": len(models),
                "n_detection_rows": len(detections),
                "n_accepted_detections": len(accepted),
                "n_penalty_grid_rows": len(penalty_grid),
                "n_boundary_rows": len(comparison),
                "n_scored_primary_matches": int(
                    (comparison["scored"].astype(bool) & primary_match).sum()
                ),
                "n_scored_sensitivity_matches": int(
                    (comparison["scored"].astype(bool) & sensitivity_match).sum()
                ),
            }
        ]
    )
