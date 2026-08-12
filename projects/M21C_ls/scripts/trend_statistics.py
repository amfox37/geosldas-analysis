#!/usr/bin/env python3
"""Robust, autocorrelation-aware trend statistics for M21C monthly series."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import kendalltau, rankdata, theilslopes


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "trend_statistics.json"

STATUS = {
    0: "success",
    1: "insufficient eligible months",
    2: "insufficient valid fraction or valid months",
    3: "insufficient time span",
    4: "insufficient values after seasonal adjustment",
}

METRIC_NAMES = [
    "slope",
    "intercept",
    "slope_ci_low",
    "slope_ci_high",
    "slope_ci_low_nominal",
    "slope_ci_high_nominal",
    "p_value",
    "p_value_original_mk",
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
    "calendar_month_mask",
]


def load_trend_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load and validate the trend-statistics configuration."""

    config = json.loads(Path(path).read_text())
    if config["seasonal_adjustment"] not in {
        "trend_preserving_calendar_month_anomaly",
        "none",
    }:
        raise ValueError(
            "seasonal_adjustment must be "
            "trend_preserving_calendar_month_anomaly or none"
        )
    fraction = float(config["minimum_valid_fraction_of_eligible"])
    if not 0 < fraction <= 1:
        raise ValueError("minimum_valid_fraction_of_eligible must be in (0, 1]")
    if int(config["minimum_valid_months"]) < 3:
        raise ValueError("minimum_valid_months must be at least 3")
    confidence = float(config["theil_sen_confidence_level"])
    if not 0 < confidence < 1:
        raise ValueError("theil_sen_confidence_level must be in (0, 1)")
    fdr_alpha = float(config["fdr"]["alpha"])
    if not 0 < fdr_alpha < 1:
        raise ValueError("FDR alpha must be in (0, 1)")
    if not bool(config["modified_mk"]["positive_significant_lags_only"]):
        raise ValueError("This implementation requires conservative positive-only ACF inflation")
    if float(config["modified_mk"]["minimum_variance_factor"]) < 1:
        raise ValueError("minimum_variance_factor must be at least 1")
    return config


def monthly_time_axis(time: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return exact monthly decimal years and calendar-month numbers."""

    index = pd.DatetimeIndex(np.asarray(time))
    if index.has_duplicates or not index.is_monotonic_increasing:
        raise ValueError("Time coordinate must be unique and increasing")
    expected = pd.date_range(index[0], index[-1], freq="MS")
    if not index.equals(expected):
        raise ValueError("Time coordinate must be a contiguous month-start sequence")
    month_number = index.year * 12 + index.month - 1
    decimal_year = (month_number - month_number[0]).to_numpy(dtype="float64") / 12.0
    return decimal_year, index.month.to_numpy(dtype="int16")


def calendar_month_anomalies(
    values: np.ndarray,
    time_years: np.ndarray,
    months: np.ndarray,
    valid: np.ndarray,
    minimum_samples: int,
) -> tuple[np.ndarray, int]:
    """Remove seasonality while preserving a preliminary exact Sen trend."""

    adjusted = np.full(values.shape, np.nan, dtype="float64")
    calendar_month_mask = 0
    if int(valid.sum()) < 2 or np.ptp(values[valid]) == 0:
        preliminary_slope = 0.0
    else:
        preliminary_slope = float(
            theilslopes(values[valid], x=time_years[valid], method="joint").slope
        )
    detrended = values - preliminary_slope * time_years
    for month in range(1, 13):
        use = valid & (months == month)
        if int(use.sum()) >= minimum_samples:
            seasonal_climatology = float(np.mean(detrended[use]))
            adjusted[use] = values[use] - seasonal_climatology
            calendar_month_mask |= 1 << (month - 1)
    return adjusted, calendar_month_mask


def _tie_counts(values: np.ndarray) -> np.ndarray:
    _, counts = np.unique(values, return_counts=True)
    return counts[counts > 1].astype("float64")


def _mann_kendall_score_and_variance(values: np.ndarray) -> tuple[float, float]:
    """Return Mann-Kendall S and its tie-corrected independent variance."""

    n = values.size
    if n < 2:
        return np.nan, np.nan
    ties = _tie_counts(values)
    n_pairs = n * (n - 1) / 2.0
    tied_pairs = float(np.sum(ties * (ties - 1.0) / 2.0))
    tau = float(kendalltau(np.arange(n), values, nan_policy="omit").statistic)
    if not np.isfinite(tau):
        score = 0.0
    else:
        score = float(np.rint(tau * math.sqrt(n_pairs * (n_pairs - tied_pairs))))
    tie_term = float(np.sum(ties * (ties - 1.0) * (2.0 * ties + 5.0)))
    variance = (n * (n - 1.0) * (2.0 * n + 5.0) - tie_term) / 18.0
    return score, variance


def _normal_two_sided_p(score: float, variance: float) -> float:
    if not np.isfinite(score) or not np.isfinite(variance) or variance <= 0:
        return 1.0
    if score > 0:
        z = (score - 1.0) / math.sqrt(variance)
    elif score < 0:
        z = (score + 1.0) / math.sqrt(variance)
    else:
        z = 0.0
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def _pearson_pair(values_a: np.ndarray, values_b: np.ndarray, minimum_pairs: int) -> float:
    use = np.isfinite(values_a) & np.isfinite(values_b)
    if int(use.sum()) < minimum_pairs:
        return np.nan
    a = values_a[use]
    b = values_b[use]
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def modified_mk_test(
    adjusted: np.ndarray,
    time_years: np.ndarray,
    slope: float,
    intercept: float,
    config: dict[str, Any],
) -> tuple[float, float, float, float, float, int]:
    """Return conservative Hamed-Rao-style modified MK diagnostics.

    The variance correction uses significant positive rank autocorrelation from
    detrended residuals at actual monthly lags. The factor is never allowed
    below one, so autocorrelation handling cannot make a result more significant
    than ordinary Mann-Kendall.
    """

    finite = np.isfinite(adjusted)
    values = adjusted[finite]
    score, variance = _mann_kendall_score_and_variance(values)
    original_p = _normal_two_sided_p(score, variance)
    if values.size < 4 or variance <= 0:
        return 1.0, original_p, 1.0, np.nan, np.nan, 0

    residual = np.full(adjusted.shape, np.nan, dtype="float64")
    residual[finite] = adjusted[finite] - (intercept + slope * time_years[finite])
    residual_values = residual[finite]
    residual_scale = max(1.0, float(np.nanmax(np.abs(adjusted[finite]))))
    if np.ptp(residual_values) <= 1.0e-12 * residual_scale:
        return original_p, original_p, 1.0, np.nan, np.nan, 0
    ranks = np.full(adjusted.shape, np.nan, dtype="float64")
    ranks[finite] = rankdata(residual_values, method="average")

    settings = config["modified_mk"]
    max_lag = min(int(settings["maximum_autocorrelation_lag_months"]), values.size - 3)
    minimum_pairs = int(settings["minimum_pairs_per_lag"])
    acf_alpha = float(settings["autocorrelation_significance_level"])
    critical_z = 1.959963984540054 if acf_alpha == 0.05 else None
    if critical_z is None:
        from scipy.stats import norm

        critical_z = float(norm.ppf(1.0 - acf_alpha / 2.0))

    n = values.size
    weighted_rho = 0.0
    n_positive = 0
    lag1_pearson = _pearson_pair(residual[:-1], residual[1:], minimum_pairs)
    lag1_rank = _pearson_pair(ranks[:-1], ranks[1:], minimum_pairs)
    for lag in range(1, max_lag + 1):
        pair_mask = np.isfinite(ranks[:-lag]) & np.isfinite(ranks[lag:])
        n_pair = int(pair_mask.sum())
        if n_pair < minimum_pairs:
            continue
        rho = _pearson_pair(ranks[:-lag], ranks[lag:], minimum_pairs)
        if not np.isfinite(rho):
            continue
        threshold = critical_z / math.sqrt(n_pair)
        if rho > threshold:
            weight = (n - lag) * (n - lag - 1.0) * (n - lag - 2.0)
            if weight > 0:
                weighted_rho += weight * rho
                n_positive += 1

    denominator = n * (n - 1.0) * (n - 2.0)
    factor = 1.0 + (2.0 * weighted_rho / denominator if denominator > 0 else 0.0)
    factor = max(float(settings["minimum_variance_factor"]), factor)
    corrected_p = _normal_two_sided_p(score, variance * factor)
    return corrected_p, original_p, factor, lag1_pearson, lag1_rank, n_positive


def _empty_result(
    *,
    n_eligible: int,
    n_valid: int,
    n_trend: int,
    valid_fraction: float,
    span_years: float,
    status: int,
    calendar_month_mask: int = 0,
) -> dict[str, float]:
    result = {name: np.nan for name in METRIC_NAMES}
    result.update(
        {
            "n_eligible": float(n_eligible),
            "n_valid": float(n_valid),
            "n_trend": float(n_trend),
            "valid_fraction": float(valid_fraction),
            "span_years": float(span_years),
            "status": float(status),
            "calendar_month_mask": float(calendar_month_mask),
        }
    )
    return result


def analyze_monthly_series(
    values: np.ndarray,
    eligible: np.ndarray,
    time_years: np.ndarray,
    months: np.ndarray,
    config: dict[str, Any],
) -> dict[str, float]:
    """Compute support, Theil-Sen, modified MK, and diagnostic statistics."""

    values = np.asarray(values, dtype="float64")
    eligible = np.asarray(eligible, dtype=bool)
    if values.ndim != 1 or eligible.shape != values.shape:
        raise ValueError("values and eligible must be matching one-dimensional arrays")

    valid = eligible & np.isfinite(values)
    n_eligible = int(eligible.sum())
    n_valid = int(valid.sum())
    valid_fraction = n_valid / n_eligible if n_eligible else np.nan
    if n_eligible < int(config["minimum_valid_months"]):
        return _empty_result(
            n_eligible=n_eligible,
            n_valid=n_valid,
            n_trend=0,
            valid_fraction=valid_fraction,
            span_years=np.nan,
            status=1,
        )
    if (
        n_valid < int(config["minimum_valid_months"])
        or valid_fraction < float(config["minimum_valid_fraction_of_eligible"])
    ):
        return _empty_result(
            n_eligible=n_eligible,
            n_valid=n_valid,
            n_trend=0,
            valid_fraction=valid_fraction,
            span_years=np.nan,
            status=2,
        )

    valid_indices = np.flatnonzero(valid)
    span_years = float(time_years[valid_indices[-1]] - time_years[valid_indices[0]])
    if span_years < float(config["minimum_span_years"]):
        return _empty_result(
            n_eligible=n_eligible,
            n_valid=n_valid,
            n_trend=0,
            valid_fraction=valid_fraction,
            span_years=span_years,
            status=3,
        )

    if config["seasonal_adjustment"] == "trend_preserving_calendar_month_anomaly":
        adjusted, calendar_month_mask = calendar_month_anomalies(
            values,
            time_years,
            months,
            valid,
            int(config["minimum_climatology_samples_per_month"]),
        )
    else:
        adjusted = np.where(valid, values, np.nan)
        calendar_month_mask = sum(
            1 << (month - 1) for month in range(1, 13) if np.any(valid & (months == month))
        )
    trend_mask = np.isfinite(adjusted)
    n_trend = int(trend_mask.sum())
    if n_trend < int(config["minimum_valid_months"]):
        return _empty_result(
            n_eligible=n_eligible,
            n_valid=n_valid,
            n_trend=n_trend,
            valid_fraction=valid_fraction,
            span_years=span_years,
            status=4,
            calendar_month_mask=calendar_month_mask,
        )

    x = time_years[trend_mask]
    y = adjusted[trend_mask]
    if np.ptp(y) == 0:
        slope = 0.0
        intercept = float(y[0])
        low_slope = 0.0
        high_slope = 0.0
    else:
        estimate = theilslopes(
            y,
            x=x,
            alpha=float(config["theil_sen_confidence_level"]),
            method="joint",
        )
        slope = float(estimate.slope)
        intercept = float(estimate.intercept)
        low_slope = float(estimate.low_slope)
        high_slope = float(estimate.high_slope)

    p_value, original_p, factor, lag1_pearson, lag1_rank, n_positive = modified_mk_test(
        adjusted, time_years, slope, intercept, config
    )
    ci_inflation = math.sqrt(factor)
    adjusted_low_slope = slope - (slope - low_slope) * ci_inflation
    adjusted_high_slope = slope + (high_slope - slope) * ci_inflation
    return {
        "slope": slope,
        "intercept": intercept,
        "slope_ci_low": adjusted_low_slope,
        "slope_ci_high": adjusted_high_slope,
        "slope_ci_low_nominal": low_slope,
        "slope_ci_high_nominal": high_slope,
        "p_value": p_value,
        "p_value_original_mk": original_p,
        "mk_variance_factor": factor,
        "lag1_residual_pearson_autocorrelation": lag1_pearson,
        "lag1_rank_autocorrelation": lag1_rank,
        "n_positive_autocorrelation_lags": float(n_positive),
        "n_eligible": float(n_eligible),
        "n_valid": float(n_valid),
        "n_trend": float(n_trend),
        "valid_fraction": float(valid_fraction),
        "span_years": span_years,
        "status": 0.0,
        "calendar_month_mask": float(calendar_month_mask),
    }


def benjamini_hochberg(
    p_values: np.ndarray, alpha: float = 0.05
) -> tuple[np.ndarray, np.ndarray]:
    """Return Benjamini-Hochberg q-values and rejection flags."""

    p_values = np.asarray(p_values, dtype="float64")
    q_values = np.full(p_values.shape, np.nan, dtype="float64")
    significant = np.zeros(p_values.shape, dtype=bool)
    finite_indices = np.flatnonzero(np.isfinite(p_values))
    if finite_indices.size == 0:
        return q_values, significant

    finite_p = p_values[finite_indices]
    order = np.argsort(finite_p)
    sorted_p = finite_p[order]
    m = sorted_p.size
    sorted_q = sorted_p * m / np.arange(1, m + 1, dtype="float64")
    sorted_q = np.minimum.accumulate(sorted_q[::-1])[::-1]
    sorted_q = np.clip(sorted_q, 0.0, 1.0)
    q_for_finite = np.empty(m, dtype="float64")
    q_for_finite[order] = sorted_q
    q_values[finite_indices] = q_for_finite
    significant[finite_indices] = q_for_finite <= alpha
    return q_values, significant


def compute_trend_statistics(
    values: xr.DataArray,
    eligible: xr.DataArray | None = None,
    *,
    config: dict[str, Any] | None = None,
    config_path: str | Path = DEFAULT_CONFIG,
    n_jobs: int = 1,
    parallel_prefer: str = "threads",
) -> xr.Dataset:
    """Compute trend statistics for a `(time, tile)` monthly field."""

    if n_jobs == 0:
        raise ValueError("n_jobs cannot be zero")
    if parallel_prefer not in {"threads", "processes"}:
        raise ValueError("parallel_prefer must be 'threads' or 'processes'")
    if values.dims != ("time", "tile"):
        values = values.transpose("time", "tile")
    if eligible is None:
        eligible = xr.ones_like(values, dtype=bool)
    else:
        eligible = eligible.broadcast_like(values).transpose("time", "tile")
    if values.shape != eligible.shape:
        raise ValueError("values and eligible do not align")

    settings = config if config is not None else load_trend_config(config_path)
    time_years, months = monthly_time_axis(values.time.values)
    value_array = np.asarray(values.values, dtype="float64")
    eligible_array = np.asarray(eligible.values, dtype=bool)

    def analyze_tile(tile_index: int) -> dict[str, float]:
        return analyze_monthly_series(
            value_array[:, tile_index],
            eligible_array[:, tile_index],
            time_years,
            months,
            settings,
        )

    tile_indices = range(value_array.shape[1])
    if n_jobs == 1:
        records = [analyze_tile(index) for index in tile_indices]
    else:
        from joblib import Parallel, delayed

        records = Parallel(n_jobs=n_jobs, prefer=parallel_prefer, batch_size="auto")(
            delayed(analyze_tile)(index) for index in tile_indices
        )

    output_arrays = {
        name: np.asarray([record[name] for record in records], dtype="float64")
        for name in METRIC_NAMES
    }
    q_value, significant = benjamini_hochberg(
        output_arrays["p_value"], alpha=float(settings["fdr"]["alpha"])
    )
    output_arrays["q_value"] = q_value
    calendar_month_mask = output_arrays.pop("calendar_month_mask").astype("int16")
    month_numbers = np.arange(1, 13, dtype="int16")
    calendar_month_used = (
        calendar_month_mask[:, None]
        & (1 << np.arange(12, dtype="int16"))[None, :]
    ) != 0
    ci_excludes_zero = (output_arrays["slope_ci_low"] > 0) | (
        output_arrays["slope_ci_high"] < 0
    )
    inconsistent_significance = significant & ~ci_excludes_zero
    if np.any(inconsistent_significance):
        raise RuntimeError(
            "FDR-significant tiles must exclude zero in the autocorrelation-adjusted "
            f"slope interval; found {int(inconsistent_significance.sum())} inconsistencies"
        )

    coords: dict[str, Any] = {"tile": values.tile.values}
    for name in ("lat", "lon", "tile_area"):
        if name in values.coords:
            coords[name] = values.coords[name].copy(deep=False)
    output = xr.Dataset(
        {name: ("tile", array) for name, array in output_arrays.items()},
        coords=coords,
    )
    output["significant_fdr"] = ("tile", significant)
    output["ci_excludes_zero"] = (
        "tile",
        ci_excludes_zero,
    )
    output["calendar_month_used"] = (
        ("tile", "month"),
        calendar_month_used,
    )
    output = output.assign_coords(month=month_numbers)
    output["n_calendar_months_used"] = output["calendar_month_used"].sum("month").astype(
        "int16"
    )
    output["status"] = output["status"].astype("int16")
    for name in ("n_eligible", "n_valid", "n_trend", "n_positive_autocorrelation_lags"):
        output[name] = output[name].fillna(-1).astype("int32")

    source_units = str(values.attrs.get("units", ""))
    slope_units = f"{source_units} yr-1".strip()
    for name in (
        "slope",
        "slope_ci_low",
        "slope_ci_high",
        "slope_ci_low_nominal",
        "slope_ci_high_nominal",
    ):
        output[name].attrs["units"] = slope_units
    output["slope"].attrs["long_name"] = "exact Theil-Sen slope per year"
    output["slope_ci_low"].attrs[
        "long_name"
    ] = "lower first-order autocorrelation-adjusted Theil-Sen slope confidence limit"
    output["slope_ci_high"].attrs[
        "long_name"
    ] = "upper first-order autocorrelation-adjusted Theil-Sen slope confidence limit"
    output["slope_ci_low_nominal"].attrs[
        "long_name"
    ] = "lower nominal independent-sample Theil-Sen slope confidence limit"
    output["slope_ci_high_nominal"].attrs[
        "long_name"
    ] = "upper nominal independent-sample Theil-Sen slope confidence limit"
    output["intercept"].attrs["units"] = source_units
    for name in ("p_value", "p_value_original_mk", "q_value", "valid_fraction"):
        output[name].attrs["units"] = "1"
    output["mk_variance_factor"].attrs.update(
        units="1",
        long_name="conservative Hamed-Rao-style Mann-Kendall variance inflation",
    )
    output["lag1_residual_pearson_autocorrelation"].attrs.update(
        units="1",
        long_name="lag-1 Pearson autocorrelation of Sen-detrended residuals",
    )
    output["lag1_rank_autocorrelation"].attrs.update(
        units="1",
        long_name="lag-1 rank autocorrelation used by modified Mann-Kendall",
    )
    output["p_value"].attrs["long_name"] = "modified Mann-Kendall two-sided p-value"
    output["p_value_original_mk"].attrs[
        "long_name"
    ] = "ordinary independent Mann-Kendall two-sided p-value for diagnosis only"
    output["q_value"].attrs["long_name"] = "Benjamini-Hochberg adjusted p-value"
    output["span_years"].attrs["units"] = "year"
    output["significant_fdr"].attrs.update(
        long_name="modified Mann-Kendall result significant after BH FDR",
        alpha=float(settings["fdr"]["alpha"]),
    )
    output["ci_excludes_zero"].attrs.update(
        long_name="autocorrelation-adjusted pointwise 95% slope interval excludes zero",
        note="Pointwise interval diagnostic; mapped significance must use significant_fdr",
    )
    output["calendar_month_used"].attrs.update(
        long_name="calendar month retained after seasonal-adjustment sample gate",
        units="1",
    )
    output["n_calendar_months_used"].attrs.update(
        long_name="number of calendar months retained for trend estimation",
        units="1",
    )
    output["status"].attrs.update(
        long_name="trend calculation status",
        flag_values=np.asarray(list(STATUS), dtype="int16"),
        flag_meanings=" ".join(value.replace(" ", "_") for value in STATUS.values()),
    )
    output.attrs = {
        "analysis_variable": str(values.attrs.get("analysis_variable", values.name or "")),
        "source_variable": str(values.attrs.get("source_variable", "")),
        "source_series": str(values.name or ""),
        "source_units": source_units,
        "seasonal_adjustment": settings["seasonal_adjustment"],
        "slope_method": "exact Theil-Sen, scipy.stats.theilslopes, joint intercept",
        "significance_method": settings["modified_mk"]["method"],
        "fdr_method": settings["fdr"]["method"],
        "parallel_execution": (
            "serial" if n_jobs == 1 else f"joblib preferred {parallel_prefer}"
        ),
        "configuration": json.dumps(settings, sort_keys=True),
        "status_codes": json.dumps(STATUS, sort_keys=True),
    }
    return output
