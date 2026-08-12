#!/usr/bin/env python3
"""Validate the M21C trend engine on synthetic autocorrelated monthly series."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from trend_statistics import DEFAULT_CONFIG, compute_trend_statistics, load_trend_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "output" / "trends_breakpoints" / "trend_statistics_validation.csv"
)

SCENARIOS = [
    {"scenario": "white_noise_no_trend", "phi": 0.0, "slope": 0.0},
    {"scenario": "ar1_no_trend", "phi": 0.8, "slope": 0.0},
    {"scenario": "ar1_positive_trend", "phi": 0.8, "slope": 0.02},
    {"scenario": "ar1_negative_trend", "phi": 0.8, "slope": -0.02},
]


def ar1_noise(
    rng: np.random.Generator,
    n_time: int,
    n_series: int,
    phi: float,
    innovation_sigma: float,
) -> np.ndarray:
    innovations = rng.normal(scale=innovation_sigma, size=(n_time, n_series))
    values = np.zeros_like(innovations)
    values[0] = innovations[0]
    for index in range(1, n_time):
        values[index] = phi * values[index - 1] + innovations[index]
    return values


def run_validation(
    *,
    n_series: int = 100,
    seed: int = 20260812,
    innovation_sigma: float = 0.08,
    n_jobs: int = 1,
    config_path: str | Path = DEFAULT_CONFIG,
) -> pd.DataFrame:
    """Return scenario-level false-positive, power, bias, and coverage metrics."""

    config = load_trend_config(config_path)
    time = pd.date_range("2000-06-01", periods=288, freq="MS")
    time_years = np.arange(time.size, dtype="float64") / 12.0
    seasonal = 0.2 * np.sin(2.0 * np.pi * np.arange(time.size) / 12.0)
    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    for scenario in SCENARIOS:
        true_slope = float(scenario["slope"])
        values = (
            seasonal[:, None]
            + true_slope * time_years[:, None]
            + ar1_noise(
                rng,
                time.size,
                n_series,
                float(scenario["phi"]),
                innovation_sigma,
            )
        )
        data = xr.DataArray(
            values,
            dims=("time", "tile"),
            coords={"time": time, "tile": np.arange(n_series)},
            name="synthetic",
            attrs={"units": "synthetic unit"},
        )
        result = compute_trend_statistics(data, config=config, n_jobs=n_jobs)
        success = result["status"].values == 0
        estimates = result["slope"].values[success]
        ci_cover = (
            (result["slope_ci_low"].values[success] <= true_slope)
            & (result["slope_ci_high"].values[success] >= true_slope)
        )
        rows.append(
            {
                **scenario,
                "n_series": n_series,
                "n_success": int(success.sum()),
                "ordinary_mk_p_lt_0.05_fraction": float(
                    np.mean(result["p_value_original_mk"].values[success] < 0.05)
                ),
                "modified_mk_p_lt_0.05_fraction": float(
                    np.mean(result["p_value"].values[success] < 0.05)
                ),
                "bh_fdr_detection_fraction": float(
                    np.mean(result["significant_fdr"].values[success])
                ),
                "median_slope": float(np.median(estimates)),
                "slope_bias": float(np.mean(estimates - true_slope)),
                "slope_rmse": float(np.sqrt(np.mean((estimates - true_slope) ** 2))),
                "theil_sen_ci_coverage": float(np.mean(ci_cover)),
                "median_lag1_autocorrelation": float(
                    np.nanmedian(result["lag1_autocorrelation"].values[success])
                ),
                "median_mk_variance_factor": float(
                    np.nanmedian(result["mk_variance_factor"].values[success])
                ),
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-series", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--innovation-sigma", type=float, default=0.08)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = run_validation(
        n_series=args.n_series,
        seed=args.seed,
        innovation_sigma=args.innovation_sigma,
        n_jobs=args.n_jobs,
        config_path=args.config,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, index=False)
    print(summary.to_string(index=False))
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
