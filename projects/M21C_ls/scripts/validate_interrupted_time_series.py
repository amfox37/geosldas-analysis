#!/usr/bin/env python3
"""Fixed-seed calibration of the M21C interrupted-series inference."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from interrupted_time_series import (
    DEFAULT_CONFIG,
    build_interrupted_design,
    fit_interrupted_series,
    load_interrupted_config,
)
from trend_statistics import benjamini_hochberg


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "output" / "trends_breakpoints" / "interrupted_series_validation.csv"
)
TRANSITION_TYPES = {"level_change", "slope_change"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-series", type=int, default=24)
    parser.add_argument("--bootstrap-replicates", type=int, default=499)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--parallel-prefer", choices=["threads", "processes"], default="processes")
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--maximum-null-fdr-discoveries", type=int, default=2)
    parser.add_argument("--minimum-target-coverage", type=float, default=0.80)
    parser.add_argument("--maximum-target-relative-bias", type=float, default=0.50)
    return parser.parse_args()


def _simulate(
    design: np.ndarray,
    beta: np.ndarray,
    *,
    phi: float,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    innovation = rng.normal(size=design.shape[0] + 120)
    error = np.empty(innovation.size, dtype="float64")
    state = 0.0
    for index, value in enumerate(innovation):
        state = phi * state + value
        error[index] = state
    return design @ beta + error[120:]


def _fit_scenario(
    scenario: str,
    series_index: int,
    values: np.ndarray,
    time: pd.DatetimeIndex,
    settings: dict[str, Any],
    seed: int,
) -> pd.DataFrame:
    result = fit_interrupted_series(
        values,
        time,
        config=settings,
        random_seed=seed,
    )
    frame = result.coefficients[
        result.coefficients["coefficient_type"].isin(TRANSITION_TYPES)
    ].copy()
    frame["scenario"] = scenario
    frame["series_index"] = series_index
    return frame


def main() -> int:
    args = parse_args()
    if args.n_series < 8:
        raise ValueError("--n-series must be at least 8")
    if args.bootstrap_replicates < 99:
        raise ValueError("--bootstrap-replicates must be at least 99")
    settings = load_interrupted_config(args.config)
    settings["bootstrap"] = {
        **settings["bootstrap"],
        "replicates": int(args.bootstrap_replicates),
    }
    time = pd.date_range("2000-06-01", "2024-05-01", freq="MS")
    design, metadata, _ = build_interrupted_design(time, config=settings)
    coefficient_index = {
        name: index for index, name in enumerate(metadata["coefficient"])
    }
    scenarios = {
        "ar1_no_transition": {},
        "ar1_level_change_P6": {"level_change_P6": 2.0},
        "ar1_slope_change_P6": {"slope_change_P6": 2.0},
    }
    tasks: list[tuple[str, int, np.ndarray, int]] = []
    truth: dict[tuple[str, str], float] = {}
    task_number = 0
    for scenario_index, (scenario, effects) in enumerate(scenarios.items()):
        beta = np.zeros(design.shape[1], dtype="float64")
        beta[coefficient_index["intercept"]] = 10.0
        for coefficient, value in effects.items():
            beta[coefficient_index[coefficient]] = value
        for coefficient in metadata.loc[
            metadata["coefficient_type"].isin(TRANSITION_TYPES), "coefficient"
        ]:
            truth[(scenario, coefficient)] = effects.get(coefficient, 0.0)
        for series_index in range(args.n_series):
            simulation_seed = args.seed + scenario_index * 10000 + series_index
            values = _simulate(design, beta, phi=0.8, seed=simulation_seed)
            tasks.append((scenario, series_index, values, args.seed + 100000 + task_number))
            task_number += 1

    print(
        f"Fitting {len(tasks)} synthetic series with "
        f"{args.bootstrap_replicates} bootstrap replicates each",
        flush=True,
    )
    frames = Parallel(n_jobs=args.n_jobs, prefer=args.parallel_prefer)(
        delayed(_fit_scenario)(scenario, index, values, time, settings, seed)
        for scenario, index, values, seed in tasks
    )
    estimates = pd.concat(frames, ignore_index=True)
    estimates["truth"] = [
        truth[(scenario, coefficient)]
        for scenario, coefficient in zip(estimates["scenario"], estimates["coefficient"])
    ]
    estimates["covered"] = (
        (estimates["ci_low_bootstrap"] <= estimates["truth"])
        & (estimates["ci_high_bootstrap"] >= estimates["truth"])
    )
    estimates["significant_pointwise"] = estimates["p_value"] <= 0.05
    estimates["q_value"] = np.nan
    estimates["significant_fdr"] = False
    for (_, _), indices in estimates.groupby(["scenario", "coefficient"]).groups.items():
        q_value, significant = benjamini_hochberg(
            estimates.loc[indices, "p_value"].to_numpy(), alpha=0.05
        )
        estimates.loc[indices, "q_value"] = q_value
        estimates.loc[indices, "significant_fdr"] = significant

    rows: list[dict[str, Any]] = []
    for (scenario, coefficient), frame in estimates.groupby(["scenario", "coefficient"]):
        true_value = float(frame["truth"].iloc[0])
        error = frame["estimate"] - true_value
        rows.append(
            {
                "scenario": scenario,
                "coefficient": coefficient,
                "coefficient_type": frame["coefficient_type"].iloc[0],
                "truth": true_value,
                "n_series": len(frame),
                "mean_estimate": float(frame["estimate"].mean()),
                "bias": float(error.mean()),
                "rmse": float(np.sqrt(np.mean(error**2))),
                "bootstrap_ci_coverage": float(frame["covered"].mean()),
                "fraction_pointwise_significant": float(
                    frame["significant_pointwise"].mean()
                ),
                "n_fdr_significant": int(frame["significant_fdr"].sum()),
                "fraction_fdr_significant": float(frame["significant_fdr"].mean()),
            }
        )
    report = pd.DataFrame(rows).sort_values(["scenario", "coefficient"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.incomplete")
    report.to_csv(temporary, index=False)
    temporary.replace(args.output)

    null = report[report["scenario"] == "ar1_no_transition"]
    targets = report[
        ((report["scenario"] == "ar1_level_change_P6") & (report["coefficient"] == "level_change_P6"))
        | ((report["scenario"] == "ar1_slope_change_P6") & (report["coefficient"] == "slope_change_P6"))
    ]
    print(
        f"Null: {int(null['n_fdr_significant'].sum())} BH discoveries across "
        f"{len(null)} transition families",
        flush=True,
    )
    print(targets.to_string(index=False), flush=True)
    null_pass = int(null["n_fdr_significant"].sum()) <= int(
        args.maximum_null_fdr_discoveries
    )
    target_pass = bool(
        np.all(targets["bootstrap_ci_coverage"] >= args.minimum_target_coverage)
        and np.all(np.sign(targets["mean_estimate"]) == np.sign(targets["truth"]))
        and np.all(
            np.abs(targets["bias"] / targets["truth"])
            <= args.maximum_target_relative_bias
        )
    )
    print(
        "Calibration guards: "
        f"null={'PASS' if null_pass else 'FAIL'}, "
        f"known effects={'PASS' if target_pass else 'FAIL'}",
        flush=True,
    )
    print(f"Wrote {args.output}", flush=True)
    return 0 if null_pass and target_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
