#!/usr/bin/env python3
"""Fixed-seed false-positive and recovery validation for M21C changepoints."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from changepoint_detection import (
    DEFAULT_CONFIG,
    detect_changepoints,
    known_boundary_table,
    load_changepoint_config,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "output" / "trends_breakpoints" / "changepoint_validation.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-series", type=int, default=24)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--parallel-prefer", choices=["threads", "processes"], default="processes")
    return parser.parse_args()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _simulate(
    *,
    seed: int,
    phi: float,
    level_changes: dict[int, float] | None = None,
    slope_changes: dict[int, float] | None = None,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    innovation = rng.normal(size=408)
    error = np.zeros(408, dtype="float64")
    for index in range(1, error.size):
        error[index] = phi * error[index - 1] + innovation[index]
    month = np.arange(288, dtype="float64")
    values = (
        10.0
        + 0.02 * month / 12.0
        + 2.0 * np.sin(2.0 * np.pi * month / 12.0)
        + error[-288:]
    )
    for boundary, effect in (level_changes or {}).items():
        values[boundary:] += effect
    for boundary, endpoint_effect in (slope_changes or {}).items():
        values += np.maximum(0.0, month - boundary) * endpoint_effect / (288 - boundary)
    return values


def _evaluate_one(
    scenario: str,
    series_index: int,
    values: np.ndarray,
    targets: tuple[int, ...],
    settings: dict[str, Any],
) -> dict[str, Any]:
    time = pd.date_range("2000-06-01", periods=288, freq="MS")
    result = detect_changepoints(values, time, config=settings)
    accepted = result.detections.loc[
        result.detections["accepted_detection"], "break_index"
    ].astype(int).tolist()
    primary_tolerance = int(
        settings["known_boundary_comparison"]["primary_tolerance_months"]
    )
    sensitivity_tolerance = int(
        settings["known_boundary_comparison"]["sensitivity_tolerance_months"]
    )
    target_distances = [
        min((abs(target - detected) for detected in accepted), default=np.inf)
        for target in targets
    ]
    matched_detections = {
        detected
        for detected in accepted
        if any(abs(target - detected) <= sensitivity_tolerance for target in targets)
    }
    return {
        "scenario": scenario,
        "series_index": series_index,
        "n_accepted": len(accepted),
        "n_targets": len(targets),
        "n_targets_primary": sum(distance <= primary_tolerance for distance in target_distances),
        "n_targets_sensitivity": sum(
            distance <= sensitivity_tolerance for distance in target_distances
        ),
        "n_spurious": len(set(accepted).difference(matched_detections)),
        "median_target_error_months": (
            float(np.median(target_distances)) if targets and np.all(np.isfinite(target_distances)) else np.nan
        ),
        "estimated_global_ar1": result.summary["estimated_global_ar1"],
        "estimated_local_detrended_ar1": result.summary[
            "estimated_local_detrended_ar1"
        ],
    }


def main() -> int:
    args = parse_args()
    if args.n_series < 8:
        raise ValueError("--n-series must be at least 8")
    if args.n_jobs == 0:
        raise ValueError("--n-jobs cannot be zero")
    settings = load_changepoint_config(args.config)
    validation = settings["validation"]
    seed = int(validation["random_seed"])
    phi = float(validation["simulation_ar1"])
    level_effect = float(validation["planted_level_change_innovation_sd"])
    slope_effect = float(validation["planted_slope_endpoint_change_innovation_sd"])
    time = pd.date_range("2000-06-01", periods=288, freq="MS")
    boundaries = known_boundary_table(time).set_index("period_id")
    p3 = int(boundaries.loc["P3", "boundary_index"])
    p6 = int(boundaries.loc["P6", "boundary_index"])
    edge = int(settings["minimum_segment_months"] + 6)
    scenarios = {
        "white_noise_null": {"phi": 0.0, "targets": (), "levels": {}, "slopes": {}},
        "ar1_null": {"phi": phi, "targets": (), "levels": {}, "slopes": {}},
        "ar1_level_P6": {
            "phi": phi,
            "targets": (p6,),
            "levels": {p6: level_effect},
            "slopes": {},
        },
        "ar1_slope_P6": {
            "phi": phi,
            "targets": (p6,),
            "levels": {},
            "slopes": {p6: slope_effect},
        },
        "ar1_two_levels_P3_P6": {
            "phi": phi,
            "targets": (p3, p6),
            "levels": {p3: level_effect, p6: -level_effect},
            "slopes": {},
        },
        "ar1_level_near_start": {
            "phi": phi,
            "targets": (edge,),
            "levels": {edge: level_effect},
            "slopes": {},
        },
    }
    tasks: list[tuple[str, int, np.ndarray, tuple[int, ...]]] = []
    for scenario_index, (scenario, spec) in enumerate(scenarios.items()):
        for series_index in range(args.n_series):
            values = _simulate(
                seed=seed + scenario_index * 10000 + series_index,
                phi=float(spec["phi"]),
                level_changes=spec["levels"],
                slope_changes=spec["slopes"],
            )
            tasks.append((scenario, series_index, values, tuple(spec["targets"])))
    print(f"Running {len(tasks)} synthetic changepoint cases", flush=True)
    rows = Parallel(n_jobs=args.n_jobs, prefer=args.parallel_prefer)(
        delayed(_evaluate_one)(scenario, index, values, targets, settings)
        for scenario, index, values, targets in tasks
    )
    raw = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for scenario, frame in raw.groupby("scenario", sort=False):
        n_targets = int(frame["n_targets"].iloc[0])
        denominator = len(frame) * n_targets
        summary_rows.append(
            {
                "scenario": scenario,
                "n_series": len(frame),
                "n_targets_per_series": n_targets,
                "mean_accepted_breaks": float(frame["n_accepted"].mean()),
                "mean_spurious_breaks": float(frame["n_spurious"].mean()),
                "primary_recovery_fraction": (
                    float(frame["n_targets_primary"].sum() / denominator)
                    if denominator
                    else np.nan
                ),
                "sensitivity_recovery_fraction": (
                    float(frame["n_targets_sensitivity"].sum() / denominator)
                    if denominator
                    else np.nan
                ),
                "median_target_error_months": (
                    float(frame["median_target_error_months"].median())
                    if frame["median_target_error_months"].notna().any()
                    else np.nan
                ),
                "mean_estimated_global_ar1": float(frame["estimated_global_ar1"].mean()),
                "mean_estimated_local_detrended_ar1": float(
                    frame["estimated_local_detrended_ar1"].mean()
                ),
            }
        )
    report = pd.DataFrame(summary_rows)
    report["validation_seed"] = seed
    report["configuration"] = json.dumps(settings, sort_keys=True)
    report["source_git_commit"] = git_commit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.incomplete")
    report.to_csv(temporary, index=False)
    temporary.replace(args.output)

    null = report[report["n_targets_per_series"] == 0]
    effects = report[report["n_targets_per_series"] > 0]
    null_pass = bool(
        np.all(
            null["mean_accepted_breaks"]
            <= float(validation["maximum_mean_accepted_false_breaks"])
        )
    )
    guarded_scenarios = set(validation["recovery_guard_scenarios"])
    guarded_effects = effects[effects["scenario"].isin(guarded_scenarios)]
    if set(guarded_effects["scenario"]) != guarded_scenarios:
        raise ValueError("One or more configured recovery guard scenarios are absent")
    recovery_pass = bool(
        np.all(
            guarded_effects["sensitivity_recovery_fraction"]
            >= float(validation["minimum_known_break_recovery_fraction"])
        )
    )
    print(report.drop(columns=["configuration"]).to_string(index=False), flush=True)
    print(
        f"Validation guards: null={'PASS' if null_pass else 'FAIL'}, "
        f"recovery={'PASS' if recovery_pass else 'FAIL'}",
        flush=True,
    )
    print(f"Wrote {args.output}", flush=True)
    return 0 if null_pass and recovery_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
