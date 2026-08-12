from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from trend_statistics import (
    _mann_kendall_score_and_variance,
    analyze_monthly_series,
    benjamini_hochberg,
    compute_trend_statistics,
    load_trend_config,
    monthly_time_axis,
)


TIME = pd.date_range("2000-06-01", periods=288, freq="MS")
TIME_YEARS, MONTHS = monthly_time_axis(TIME)


def _config(**overrides: object) -> dict:
    config = copy.deepcopy(load_trend_config())
    config.update(overrides)
    return config


def _ar1_noise(
    rng: np.random.Generator, n_time: int, n_tile: int, phi: float, sigma: float
) -> np.ndarray:
    innovations = rng.normal(scale=sigma, size=(n_time, n_tile))
    noise = np.zeros_like(innovations)
    noise[0] = innovations[0]
    for index in range(1, n_time):
        noise[index] = phi * noise[index - 1] + innovations[index]
    return noise


def _data_array(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(
        values,
        dims=("time", "tile"),
        coords={"time": TIME, "tile": np.arange(values.shape[1])},
        name="delta",
        attrs={
            "units": "m3 m-3",
            "analysis_variable": "SFMC",
            "source_variable": "SFMC",
        },
    )


def test_monthly_time_axis_rejects_non_month_start_or_gap() -> None:
    with pytest.raises(ValueError, match="month-start"):
        monthly_time_axis(TIME + pd.Timedelta(days=14))
    with pytest.raises(ValueError, match="contiguous"):
        monthly_time_axis(TIME.delete(10))


def test_mann_kendall_score_matches_brute_force_with_ties() -> None:
    values = np.array([1.0, 2.0, 2.0, 0.0, 4.0, 3.0])
    brute_score = sum(
        np.sign(values[j] - values[i])
        for i in range(values.size - 1)
        for j in range(i + 1, values.size)
    )
    score, variance = _mann_kendall_score_and_variance(values)
    assert score == brute_score
    assert variance > 0


def test_calendar_anomaly_theil_sen_recovers_known_slope() -> None:
    seasonal_cycle = 0.5 * np.sin(2.0 * np.pi * np.arange(TIME.size) / 12.0)
    values = 0.02 * TIME_YEARS + seasonal_cycle
    result = analyze_monthly_series(
        values,
        np.ones(TIME.size, dtype=bool),
        TIME_YEARS,
        MONTHS,
        _config(),
    )

    assert result["status"] == 0
    assert result["slope"] == pytest.approx(0.02, abs=1.0e-12)
    assert result["slope_ci_low"] <= 0.02 <= result["slope_ci_high"]
    assert result["p_value"] < 1.0e-10
    assert result["mk_variance_factor"] == 1.0
    assert result["calendar_month_mask"] == 4095


def test_valid_fraction_uses_eligible_months_as_denominator() -> None:
    eligible = np.zeros(TIME.size, dtype=bool)
    eligible[TIME.month.isin([12, 1, 2])] = True
    values = 0.01 * TIME_YEARS
    eligible_indices = np.flatnonzero(eligible)
    values[eligible_indices[:10]] = np.nan
    config = _config(minimum_valid_months=50)

    result = analyze_monthly_series(values, eligible, TIME_YEARS, MONTHS, config)

    assert result["status"] == 0
    assert result["n_eligible"] == 72
    assert result["n_valid"] == 62
    assert result["valid_fraction"] == pytest.approx(62.0 / 72.0)


def test_support_thresholds_return_explicit_status() -> None:
    values = 0.01 * TIME_YEARS
    eligible = np.ones(TIME.size, dtype=bool)
    values[:70] = np.nan
    result = analyze_monthly_series(values, eligible, TIME_YEARS, MONTHS, _config())
    assert result["status"] == 2
    assert result["n_valid"] == 218

    short_eligible = np.zeros(TIME.size, dtype=bool)
    short_eligible[:50] = True
    result = analyze_monthly_series(
        0.01 * TIME_YEARS, short_eligible, TIME_YEARS, MONTHS, _config()
    )
    assert result["status"] == 1


def test_positive_autocorrelation_only_inflates_mk_variance() -> None:
    rng = np.random.default_rng(19)
    values = 0.004 * TIME_YEARS + _ar1_noise(rng, TIME.size, 1, 0.9, 0.05)[:, 0]
    result = analyze_monthly_series(
        values, np.ones(TIME.size, dtype=bool), TIME_YEARS, MONTHS, _config()
    )

    assert result["status"] == 0
    assert result["lag1_residual_pearson_autocorrelation"] > 0.5
    assert result["lag1_rank_autocorrelation"] > 0.5
    assert result["mk_variance_factor"] > 1.0
    assert result["n_positive_autocorrelation_lags"] > 0
    assert result["p_value"] >= result["p_value_original_mk"]


def test_benjamini_hochberg_preserves_nan_and_expected_rejections() -> None:
    p_values = np.array([0.001, 0.01, 0.04, 0.2, np.nan])
    q_values, significant = benjamini_hochberg(p_values, alpha=0.05)
    np.testing.assert_allclose(q_values[:4], [0.004, 0.02, 0.0533333333, 0.2])
    np.testing.assert_array_equal(significant, [True, True, False, False, False])
    assert np.isnan(q_values[-1])


def test_known_ar1_trends_are_recovered_and_detected_after_fdr() -> None:
    rng = np.random.default_rng(22)
    n_tile = 24
    slopes = np.linspace(0.015, 0.03, n_tile)
    seasonal = 0.2 * np.sin(2.0 * np.pi * np.arange(TIME.size) / 12.0)
    values = (
        seasonal[:, None]
        + TIME_YEARS[:, None] * slopes[None, :]
        + _ar1_noise(rng, TIME.size, n_tile, 0.6, 0.04)
    )

    result = compute_trend_statistics(_data_array(values), config=_config(), n_jobs=1)

    assert (result["status"] == 0).all()
    np.testing.assert_allclose(result["slope"], slopes, atol=0.003)
    assert int(result["significant_fdr"].sum()) >= 22
    assert bool(result["ci_excludes_zero"].where(result["significant_fdr"], True).all())
    assert (result["p_value"] >= result["p_value_original_mk"]).all()
    assert (
        (result["slope_ci_high"] - result["slope_ci_low"])
        >= (result["slope_ci_high_nominal"] - result["slope_ci_low_nominal"])
    ).all()


def test_autocorrelated_no_trend_field_has_controlled_fdr_false_positives() -> None:
    rng = np.random.default_rng(117)
    n_tile = 48
    seasonal = 0.2 * np.sin(2.0 * np.pi * np.arange(TIME.size) / 12.0)
    values = seasonal[:, None] + _ar1_noise(rng, TIME.size, n_tile, 0.8, 0.08)

    result = compute_trend_statistics(_data_array(values), config=_config(), n_jobs=1)

    assert (result["status"] == 0).all()
    assert int(result["significant_fdr"].sum()) <= 2
    assert float(np.median(result["mk_variance_factor"])) > 1.0


def test_adjusted_ci_widens_under_ar1_and_can_restore_zero() -> None:
    rng = np.random.default_rng(314)
    values = 0.0035 * TIME_YEARS + _ar1_noise(rng, TIME.size, 1, 0.92, 0.06)[:, 0]
    result = analyze_monthly_series(
        values, np.ones(TIME.size, dtype=bool), TIME_YEARS, MONTHS, _config()
    )

    nominal_width = result["slope_ci_high_nominal"] - result["slope_ci_low_nominal"]
    adjusted_width = result["slope_ci_high"] - result["slope_ci_low"]
    assert result["mk_variance_factor"] > 1.0
    assert adjusted_width == pytest.approx(
        nominal_width * np.sqrt(result["mk_variance_factor"])
    )
    if result["p_value"] >= 0.05:
        assert result["slope_ci_low"] <= 0 <= result["slope_ci_high"]


def test_calendar_month_retention_is_exposed_per_tile() -> None:
    eligible = np.ones((TIME.size, 2), dtype=bool)
    eligible[TIME.month == 2, 1] = False
    values = np.column_stack([0.01 * TIME_YEARS, 0.01 * TIME_YEARS])
    result = compute_trend_statistics(
        _data_array(values),
        xr.DataArray(eligible, dims=("time", "tile"), coords={"time": TIME, "tile": [0, 1]}),
        config=_config(),
    )

    assert int(result["n_calendar_months_used"].sel(tile=0)) == 12
    assert int(result["n_calendar_months_used"].sel(tile=1)) == 11
    assert bool(result["calendar_month_used"].sel(tile=0, month=2))
    assert not bool(result["calendar_month_used"].sel(tile=1, month=2))


def test_parallel_and_serial_results_match() -> None:
    rng = np.random.default_rng(9)
    values = 0.01 * TIME_YEARS[:, None] + _ar1_noise(rng, TIME.size, 4, 0.4, 0.03)
    data = _data_array(values)
    serial = compute_trend_statistics(data, config=_config(), n_jobs=1)
    parallel = compute_trend_statistics(data, config=_config(), n_jobs=2)

    xr.testing.assert_allclose(serial, parallel)


def test_invalid_parallel_preference_is_rejected() -> None:
    values = _data_array(np.zeros((TIME.size, 1)))
    with pytest.raises(ValueError, match="parallel_prefer"):
        compute_trend_statistics(values, n_jobs=2, parallel_prefer="invalid")
