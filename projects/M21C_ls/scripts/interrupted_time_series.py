#!/usr/bin/env python3
"""Known-transition interrupted time-series statistics for M21C monthly series."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from m21c_periods import DEFAULT_REGISTRY, load_period_frames
from trend_statistics import monthly_time_axis


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "interrupted_time_series.json"


@dataclass
class InterruptedSeriesResult:
    """Complete result for one monthly response series."""

    coefficients: pd.DataFrame
    period_slopes: pd.DataFrame
    summary: dict[str, Any]
    fitted: np.ndarray
    residual: np.ndarray


def load_interrupted_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load and validate the interrupted-series configuration."""

    config = json.loads(Path(path).read_text())
    if config["model"] != "segmented_prais_winsten_ar1_with_innovation_bootstrap":
        raise ValueError("Unsupported interrupted-series model")
    if not bool(config["monthly_fixed_effects"]):
        raise ValueError("This implementation requires monthly fixed effects")
    reference_month = int(config["reference_calendar_month"])
    if reference_month not in range(1, 13):
        raise ValueError("reference_calendar_month must be in 1..12")
    if int(config["minimum_valid_months"]) < 24:
        raise ValueError("minimum_valid_months must be at least 24")
    valid_fraction = float(config["minimum_valid_fraction"])
    if not 0 < valid_fraction <= 1:
        raise ValueError("minimum_valid_fraction must be in (0, 1]")
    confidence = float(config["confidence_level"])
    if not 0 < confidence < 1:
        raise ValueError("confidence_level must be in (0, 1)")
    bootstrap = config["bootstrap"]
    if bootstrap["method"] != "fitted_ar1_resampled_innovations":
        raise ValueError("Unsupported interrupted-series bootstrap method")
    if int(bootstrap["replicates"]) < 99:
        raise ValueError("bootstrap replicates must be at least 99")
    if int(bootstrap["burn_in_months"]) < 0:
        raise ValueError("bootstrap burn_in_months cannot be negative")
    if bootstrap["ar1_simulation"] != "upper_large_sample_confidence_bound":
        raise ValueError("Unsupported bootstrap AR(1) simulation rule")
    if not 0 < float(bootstrap["ar1_confidence_level"]) < 1:
        raise ValueError("bootstrap ar1_confidence_level must be in (0, 1)")
    autocorrelation = config["autocorrelation"]
    if autocorrelation["estimator"] != "iterative_prais_winsten_ar1":
        raise ValueError("Only iterative Prais-Winsten AR(1) estimation is implemented")
    if int(autocorrelation["maximum_iterations"]) < 1:
        raise ValueError("maximum_iterations must be positive")
    if float(autocorrelation["convergence_tolerance"]) <= 0:
        raise ValueError("convergence_tolerance must be positive")
    maximum_rho = float(autocorrelation["maximum_absolute_ar1"])
    if not 0 < maximum_rho < 1:
        raise ValueError("maximum_absolute_ar1 must be in (0, 1)")
    if autocorrelation["covariance"] != "newey_west_bartlett":
        raise ValueError("Only Newey-West/Bartlett covariance is implemented")
    if int(autocorrelation["maximum_lag_months"]) < 1:
        raise ValueError("maximum_lag_months must be positive")
    return config


def build_interrupted_design(
    time: Any,
    *,
    config: dict[str, Any] | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Build the P1-P9 segmented design and coefficient metadata."""

    settings = config if config is not None else load_interrupted_config()
    time_index = pd.DatetimeIndex(np.asarray(time))
    time_years, months = monthly_time_axis(time_index)
    _, fine, _, _ = load_period_frames(registry_path)
    if time_index[0] != fine["start"].iloc[0] or time_index[-1].to_period("M") != fine[
        "end"
    ].iloc[-1].to_period("M"):
        raise ValueError("Series time coordinate does not cover the P1-P9 registry")

    columns: list[np.ndarray] = []
    rows: list[dict[str, Any]] = []

    def append(
        name: str,
        values: np.ndarray,
        coefficient_type: str,
        *,
        period_id: str = "",
        boundary: pd.Timestamp | None = None,
        rationale: str = "",
    ) -> None:
        columns.append(np.asarray(values, dtype="float64"))
        rows.append(
            {
                "coefficient": name,
                "coefficient_type": coefficient_type,
                "period_id": period_id,
                "boundary": boundary,
                "boundary_rationale": rationale,
            }
        )

    append("intercept", np.ones(time_index.size), "nuisance")
    append("baseline_slope_P1", time_years, "baseline_slope", period_id="P1")

    reference_month = int(settings["reference_calendar_month"])
    for month in range(1, 13):
        if month == reference_month:
            continue
        append(
            f"calendar_month_{month:02d}",
            months == month,
            "monthly_fixed_effect",
        )

    start_month_number = time_index[0].year * 12 + time_index[0].month - 1
    for row in fine.iloc[1:].itertuples():
        boundary = pd.Timestamp(row.start)
        boundary_month_number = boundary.year * 12 + boundary.month - 1
        boundary_year = (boundary_month_number - start_month_number) / 12.0
        after = time_index >= boundary
        append(
            f"level_change_{row.period_id}",
            after,
            "level_change",
            period_id=row.period_id,
            boundary=boundary,
            rationale=str(row.boundary_rationale),
        )
        if bool(row.reliable_for_slope):
            append(
                f"slope_change_{row.period_id}",
                np.where(after, time_years - boundary_year, 0.0),
                "slope_change",
                period_id=row.period_id,
                boundary=boundary,
                rationale=str(row.boundary_rationale),
            )

    design = np.column_stack(columns)
    metadata = pd.DataFrame(rows)
    if "slope_change_P7" in set(metadata["coefficient"]):
        raise RuntimeError("P7 must not receive an independent slope-change term")
    if set(metadata.loc[metadata["coefficient_type"] == "level_change", "period_id"]) != {
        f"P{index}" for index in range(2, 10)
    }:
        raise RuntimeError("The interrupted design must include all P2-P9 level changes")
    return design, metadata, fine


def newey_west_covariance(
    design: np.ndarray,
    residual: np.ndarray,
    valid: np.ndarray,
    *,
    maximum_lag: int,
    finite_sample_correction: bool,
) -> np.ndarray:
    """Return a gap-aware Newey-West covariance with Bartlett lag weights."""

    design = np.asarray(design, dtype="float64")
    residual = np.asarray(residual, dtype="float64")
    valid = np.asarray(valid, dtype=bool)
    x = design[valid]
    u = residual[valid]
    n, n_parameter = x.shape
    bread = np.linalg.inv(x.T @ x)
    score = x * u[:, None]
    meat = score.T @ score

    for lag in range(1, min(maximum_lag, design.shape[0] - 1) + 1):
        pair = valid[lag:] & valid[:-lag]
        if not np.any(pair):
            continue
        score_later = design[lag:][pair] * residual[lag:][pair, None]
        score_earlier = design[:-lag][pair] * residual[:-lag][pair, None]
        cross = score_later.T @ score_earlier
        weight = 1.0 - lag / (maximum_lag + 1.0)
        meat += weight * (cross + cross.T)

    covariance = bread @ meat @ bread
    if finite_sample_correction and n > n_parameter:
        covariance *= n / (n - n_parameter)
    return (covariance + covariance.T) / 2.0


def _estimate_ar1(residual: np.ndarray, valid: np.ndarray, maximum: float) -> float:
    """Estimate lag-1 persistence from actual adjacent valid calendar months."""

    pair = valid[1:] & valid[:-1]
    earlier = residual[:-1][pair]
    later = residual[1:][pair]
    denominator = float(earlier @ earlier)
    if earlier.size < 3 or denominator <= 0:
        return 0.0
    rho = float(earlier @ later / denominator)
    return float(np.clip(rho, -maximum, maximum))


def _prais_winsten_transform(
    values: np.ndarray,
    design: np.ndarray,
    valid: np.ndarray,
    rho: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform each contiguous valid segment while retaining its first row."""

    transformed_values = np.full(values.shape, np.nan, dtype="float64")
    transformed_design = np.full(design.shape, np.nan, dtype="float64")
    first_scale = math.sqrt(max(0.0, 1.0 - rho * rho))
    for index in np.flatnonzero(valid):
        if index == 0 or not valid[index - 1]:
            transformed_values[index] = first_scale * values[index]
            transformed_design[index] = first_scale * design[index]
        else:
            transformed_values[index] = values[index] - rho * values[index - 1]
            transformed_design[index] = design[index] - rho * design[index - 1]
    return transformed_values, transformed_design


def _fit_prais_winsten(
    values: np.ndarray,
    design: np.ndarray,
    valid: np.ndarray,
    autocorrelation: dict[str, Any],
) -> tuple[np.ndarray, float, int, np.ndarray, np.ndarray, np.ndarray]:
    """Fit iterative Prais-Winsten and return transformed fit components."""

    beta, _, _, _ = np.linalg.lstsq(design[valid], values[valid], rcond=None)
    rho = 0.0
    transformed_values = values.copy()
    transformed_design = design.copy()
    n_iteration = 0
    for n_iteration in range(1, int(autocorrelation["maximum_iterations"]) + 1):
        original_residual = values - design @ beta
        updated_rho = _estimate_ar1(
            original_residual,
            valid,
            float(autocorrelation["maximum_absolute_ar1"]),
        )
        transformed_values, transformed_design = _prais_winsten_transform(
            values, design, valid, updated_rho
        )
        updated_beta, _, _, _ = np.linalg.lstsq(
            transformed_design[valid], transformed_values[valid], rcond=None
        )
        converged = (
            abs(updated_rho - rho) < float(autocorrelation["convergence_tolerance"])
            and np.max(np.abs(updated_beta - beta))
            < float(autocorrelation["convergence_tolerance"])
        )
        rho = updated_rho
        beta = updated_beta
        if converged:
            break
    transformed_residual = transformed_values - transformed_design @ beta
    return (
        beta,
        rho,
        n_iteration,
        transformed_values,
        transformed_design,
        transformed_residual,
    )


def _bootstrap_coefficient_deviations(
    beta: np.ndarray,
    rho: float,
    transformed_residual: np.ndarray,
    design: np.ndarray,
    valid: np.ndarray,
    *,
    settings: dict[str, Any],
    random_seed: int,
) -> tuple[np.ndarray, float]:
    """Bootstrap coefficient deviations from a fitted AR(1) innovation model."""

    bootstrap = settings["bootstrap"]
    n_replicate = int(bootstrap["replicates"])
    burn_in = int(bootstrap["burn_in_months"])
    innovations = np.asarray(transformed_residual[valid], dtype="float64")
    innovations = innovations - np.mean(innovations)
    n_valid = int(valid.sum())
    n_parameter = design.shape[1]
    if n_valid > n_parameter:
        innovations *= math.sqrt(n_valid / (n_valid - n_parameter))
    rng = np.random.default_rng(random_seed)
    deviations = np.empty((n_replicate, n_parameter), dtype="float64")
    n_time = design.shape[0]
    adjacent_pairs = int(np.count_nonzero(valid[1:] & valid[:-1]))
    ar1_standard_error = math.sqrt(
        max(0.0, 1.0 - rho * rho) / max(1, adjacent_pairs)
    )
    ar1_quantile = float(
        student_t.ppf(
            (1.0 + float(bootstrap["ar1_confidence_level"])) / 2.0,
            max(1, adjacent_pairs - 1),
        )
    )
    simulation_rho = float(
        np.clip(
            rho + ar1_quantile * ar1_standard_error,
            -float(settings["autocorrelation"]["maximum_absolute_ar1"]),
            float(settings["autocorrelation"]["maximum_absolute_ar1"]),
        )
    )
    for replicate in range(n_replicate):
        draws = rng.choice(innovations, size=n_time + burn_in, replace=True)
        simulated_error = np.empty(draws.size, dtype="float64")
        state = 0.0
        for index, innovation in enumerate(draws):
            state = simulation_rho * state + innovation
            simulated_error[index] = state
        simulated_values = design @ beta + simulated_error[burn_in:]
        simulated_values[~valid] = np.nan
        simulated_beta, _, _, _, _, _ = _fit_prais_winsten(
            simulated_values,
            design,
            valid,
            settings["autocorrelation"],
        )
        deviations[replicate] = simulated_beta - beta
    return deviations, simulation_rho


def _estimate_row(
    contrast: np.ndarray,
    beta: np.ndarray,
    covariance_hac: np.ndarray,
    covariance_iid: np.ndarray,
    bootstrap_deviations: np.ndarray,
    *,
    confidence_level: float,
    degrees_freedom: int,
) -> dict[str, float]:
    estimate = float(contrast @ beta)
    variance_hac = max(0.0, float(contrast @ covariance_hac @ contrast))
    variance_iid = max(0.0, float(contrast @ covariance_iid @ contrast))
    standard_error_hac = math.sqrt(variance_hac)
    standard_error_iid = math.sqrt(variance_iid)
    deviation = bootstrap_deviations @ contrast
    standard_error_bootstrap = float(np.std(deviation, ddof=1))
    alpha = 1.0 - confidence_level
    lower_deviation, upper_deviation = np.quantile(
        deviation, [alpha / 2.0, 1.0 - alpha / 2.0]
    )
    p_value_bootstrap = float(
        (1 + np.count_nonzero(np.abs(deviation) >= abs(estimate)))
        / (deviation.size + 1)
    )
    critical = float(student_t.ppf((1.0 + confidence_level) / 2.0, degrees_freedom))
    conservative_analytic_se = max(standard_error_hac, standard_error_iid)
    if conservative_analytic_se > 0:
        statistic = estimate / conservative_analytic_se
        p_value_analytic = float(
            2.0 * student_t.sf(abs(statistic), degrees_freedom)
        )
    else:
        statistic = math.copysign(math.inf, estimate) if estimate else 0.0
        p_value_analytic = 0.0 if estimate else 1.0
    return {
        "estimate": estimate,
        "standard_error_bootstrap": standard_error_bootstrap,
        "standard_error_hac": standard_error_hac,
        "standard_error_iid": standard_error_iid,
        "ci_low_bootstrap": estimate - upper_deviation,
        "ci_high_bootstrap": estimate - lower_deviation,
        "ci_low_analytic_conservative": estimate - critical * conservative_analytic_se,
        "ci_high_analytic_conservative": estimate + critical * conservative_analytic_se,
        "test_statistic_analytic": statistic,
        "p_value_analytic_conservative": p_value_analytic,
        "p_value": p_value_bootstrap,
    }


def fit_interrupted_series(
    values: Any,
    time: Any,
    *,
    config: dict[str, Any] | None = None,
    config_path: str | Path = DEFAULT_CONFIG,
    registry_path: str | Path = DEFAULT_REGISTRY,
    random_seed: int | None = None,
) -> InterruptedSeriesResult:
    """Fit one P1-P9 segmented monthly series with bootstrap inference."""

    settings = config if config is not None else load_interrupted_config(config_path)
    values = np.asarray(values, dtype="float64")
    time_index = pd.DatetimeIndex(np.asarray(time))
    if values.ndim != 1 or values.size != time_index.size:
        raise ValueError("values and time must be matching one-dimensional arrays")
    design, metadata, fine = build_interrupted_design(
        time_index, config=settings, registry_path=registry_path
    )
    valid = np.isfinite(values)
    n_valid = int(valid.sum())
    valid_fraction = n_valid / values.size
    if n_valid < int(settings["minimum_valid_months"]):
        raise ValueError("Series has fewer than the configured minimum valid months")
    if valid_fraction < float(settings["minimum_valid_fraction"]):
        raise ValueError("Series has less than the configured valid fraction")
    valid_indices = np.flatnonzero(valid)
    span_years = (valid_indices[-1] - valid_indices[0]) / 12.0
    if span_years < float(settings["minimum_span_years"]):
        raise ValueError("Series has less than the configured minimum time span")

    rank = int(np.linalg.matrix_rank(design[valid]))
    n_parameter = design.shape[1]
    if rank != n_parameter:
        raise ValueError(f"Interrupted-series design rank {rank}, expected {n_parameter}")
    autocorrelation = settings["autocorrelation"]
    (
        beta,
        rho,
        n_iteration,
        transformed_values,
        transformed_design,
        transformed_residual,
    ) = _fit_prais_winsten(values, design, valid, autocorrelation)

    x = transformed_design[valid]
    y = transformed_values[valid]
    transformed_rank = int(np.linalg.matrix_rank(x))
    if transformed_rank != n_parameter:
        raise ValueError(
            f"Prais-Winsten design rank {transformed_rank}, expected {n_parameter}"
        )
    singular_values = np.linalg.svd(x, compute_uv=False)
    fitted = np.full(values.shape, np.nan, dtype="float64")
    fitted[valid] = design[valid] @ beta
    residual = values - fitted
    degrees_freedom = n_valid - n_parameter
    if degrees_freedom <= 0:
        raise ValueError("Interrupted-series model has no residual degrees of freedom")

    residual_sum_squares = float(np.sum(residual[valid] ** 2))
    transformed_residual_sum_squares = float(
        np.sum(transformed_residual[valid] ** 2)
    )
    centered_sum_squares = float(
        np.sum((values[valid] - np.mean(values[valid])) ** 2)
    )
    sigma_squared = transformed_residual_sum_squares / degrees_freedom
    bread = np.linalg.inv(x.T @ x)
    covariance_iid = sigma_squared * bread
    covariance_hac = newey_west_covariance(
        transformed_design,
        transformed_residual,
        valid,
        maximum_lag=int(autocorrelation["maximum_lag_months"]),
        finite_sample_correction=bool(autocorrelation["finite_sample_correction"]),
    )
    bootstrap_seed = (
        int(settings["bootstrap"]["random_seed"])
        if random_seed is None
        else int(random_seed)
    )
    bootstrap_deviations, bootstrap_simulation_rho = _bootstrap_coefficient_deviations(
        beta,
        rho,
        transformed_residual,
        design,
        valid,
        settings=settings,
        random_seed=bootstrap_seed,
    )

    coefficient_rows: list[dict[str, Any]] = []
    for index, row in metadata.iterrows():
        contrast = np.zeros(n_parameter, dtype="float64")
        contrast[index] = 1.0
        coefficient_rows.append(
            {
                **row.to_dict(),
                **_estimate_row(
                    contrast,
                    beta,
                    covariance_hac,
                    covariance_iid,
                    bootstrap_deviations,
                    confidence_level=float(settings["confidence_level"]),
                    degrees_freedom=degrees_freedom,
                ),
            }
        )
    coefficients = pd.DataFrame(coefficient_rows)

    coefficient_index = {
        name: index for index, name in enumerate(metadata["coefficient"])
    }
    slope_contrast = np.zeros(n_parameter, dtype="float64")
    slope_contrast[coefficient_index["baseline_slope_P1"]] = 1.0
    period_rows: list[dict[str, Any]] = []
    for row in fine.itertuples():
        slope_name = f"slope_change_{row.period_id}"
        if slope_name in coefficient_index:
            slope_contrast[coefficient_index[slope_name]] += 1.0
        period_rows.append(
            {
                "coefficient": f"period_slope_{row.period_id}",
                "coefficient_type": "period_slope",
                "period_id": row.period_id,
                "boundary": row.start,
                "boundary_rationale": row.boundary_rationale,
                "independently_estimated": bool(row.reliable_for_slope),
                "constraint_note": row.constraint_note,
                **_estimate_row(
                    slope_contrast,
                    beta,
                    covariance_hac,
                    covariance_iid,
                    bootstrap_deviations,
                    confidence_level=float(settings["confidence_level"]),
                    degrees_freedom=degrees_freedom,
                ),
            }
        )
    period_slopes = pd.DataFrame(period_rows)

    adjacent = valid[1:] & valid[:-1]
    if int(adjacent.sum()) >= 3 and np.ptp(residual[:-1][adjacent]) > 0:
        lag1 = float(np.corrcoef(residual[1:][adjacent], residual[:-1][adjacent])[0, 1])
    else:
        lag1 = np.nan
    durbin_watson = (
        float(np.sum((residual[1:][adjacent] - residual[:-1][adjacent]) ** 2))
        / residual_sum_squares
        if residual_sum_squares > 0
        else np.nan
    )
    if int(adjacent.sum()) >= 3 and np.ptp(transformed_residual[:-1][adjacent]) > 0:
        whitened_lag1 = float(
            np.corrcoef(
                transformed_residual[1:][adjacent],
                transformed_residual[:-1][adjacent],
            )[0, 1]
        )
    else:
        whitened_lag1 = np.nan
    summary = {
        "n_time": int(values.size),
        "n_valid": n_valid,
        "valid_fraction": valid_fraction,
        "span_years": span_years,
        "n_parameter": n_parameter,
        "design_rank": rank,
        "transformed_design_rank": transformed_rank,
        "degrees_freedom": degrees_freedom,
        "design_condition_number": float(singular_values[0] / singular_values[-1]),
        "r_squared": 1.0 - residual_sum_squares / centered_sum_squares
        if centered_sum_squares > 0
        else np.nan,
        "residual_lag1_pearson": lag1,
        "estimated_ar1": rho,
        "prais_winsten_iterations": n_iteration,
        "whitened_residual_lag1_pearson": whitened_lag1,
        "durbin_watson": durbin_watson,
        "hac_maximum_lag_months": int(autocorrelation["maximum_lag_months"]),
        "bootstrap_replicates": int(settings["bootstrap"]["replicates"]),
        "bootstrap_random_seed": bootstrap_seed,
        "bootstrap_simulation_ar1": bootstrap_simulation_rho,
        "configuration": json.dumps(settings, sort_keys=True),
    }
    return InterruptedSeriesResult(
        coefficients=coefficients,
        period_slopes=period_slopes,
        summary=summary,
        fitted=fitted,
        residual=residual,
    )
