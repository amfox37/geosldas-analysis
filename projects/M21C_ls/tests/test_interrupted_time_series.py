from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from interrupted_time_series import (
    build_interrupted_design,
    fit_interrupted_series,
    load_interrupted_config,
)
from trend_statistics import benjamini_hochberg


TIME = pd.date_range("2000-06-01", "2024-05-01", freq="MS")


def _config(**overrides: object) -> dict:
    config = copy.deepcopy(load_interrupted_config())
    config["bootstrap"]["replicates"] = 99
    config.update(overrides)
    return config


def _ar1_noise(
    rng: np.random.Generator, n_time: int, phi: float, sigma: float
) -> np.ndarray:
    innovation = rng.normal(scale=sigma, size=n_time)
    output = np.zeros(n_time, dtype="float64")
    output[0] = innovation[0]
    for index in range(1, n_time):
        output[index] = phi * output[index - 1] + innovation[index]
    return output


def test_design_has_all_levels_and_excludes_p7_slope() -> None:
    design, metadata, fine = build_interrupted_design(TIME)

    assert design.shape == (288, 28)
    assert np.linalg.matrix_rank(design) == design.shape[1]
    levels = metadata.loc[metadata["coefficient_type"] == "level_change", "period_id"]
    slopes = metadata.loc[metadata["coefficient_type"] == "slope_change", "period_id"]
    assert levels.tolist() == [f"P{index}" for index in range(2, 10)]
    assert slopes.tolist() == ["P2", "P3", "P4", "P5", "P6", "P8", "P9"]
    assert not bool(fine.set_index("period_id").loc["P7", "reliable_for_slope"])


def test_noise_free_level_and_slope_changes_are_recovered_exactly() -> None:
    design, metadata, _ = build_interrupted_design(TIME)
    beta = np.zeros(design.shape[1], dtype="float64")
    index = {name: i for i, name in enumerate(metadata["coefficient"])}
    beta[index["intercept"]] = 2.0
    beta[index["baseline_slope_P1"]] = 0.02
    beta[index["calendar_month_07"]] = 0.3
    beta[index["level_change_P6"]] = 0.7
    beta[index["slope_change_P6"]] = -0.05
    beta[index["level_change_P7"]] = 0.2
    beta[index["level_change_P8"]] = -0.4
    beta[index["slope_change_P8"]] = 0.08

    result = fit_interrupted_series(design @ beta, TIME, config=_config())
    estimates = result.coefficients.set_index("coefficient")["estimate"]

    np.testing.assert_allclose(estimates.to_numpy(), beta, atol=1.0e-10)
    period_slopes = result.period_slopes.set_index("period_id")
    assert period_slopes.loc["P6", "estimate"] == pytest.approx(-0.03, abs=1.0e-10)
    assert period_slopes.loc["P7", "estimate"] == pytest.approx(-0.03, abs=1.0e-10)
    assert period_slopes.loc["P8", "estimate"] == pytest.approx(0.05, abs=1.0e-10)
    assert not bool(period_slopes.loc["P7", "independently_estimated"])


def test_known_changes_are_recovered_with_ar1_noise() -> None:
    rng = np.random.default_rng(814)
    design, metadata, _ = build_interrupted_design(TIME)
    beta = np.zeros(design.shape[1], dtype="float64")
    index = {name: i for i, name in enumerate(metadata["coefficient"])}
    beta[index["intercept"]] = 1.0
    beta[index["baseline_slope_P1"]] = 0.01
    beta[index["level_change_P6"]] = 1.0
    beta[index["slope_change_P6"]] = 0.30
    values = design @ beta + _ar1_noise(rng, TIME.size, 0.7, 0.08)

    result = fit_interrupted_series(values, TIME, config=_config())
    estimates = result.coefficients.set_index("coefficient")

    assert estimates.loc["level_change_P6", "estimate"] == pytest.approx(1.0, abs=0.3)
    assert estimates.loc["slope_change_P6", "estimate"] == pytest.approx(0.30, abs=0.08)
    assert estimates.loc["level_change_P6", "p_value"] <= 0.02
    assert estimates.loc["slope_change_P6", "p_value"] <= 0.05
    assert result.summary["residual_lag1_pearson"] > 0.4


def test_gap_aware_hac_and_support_are_explicit() -> None:
    values = 0.01 * np.arange(TIME.size) / 12.0
    values[50:60] = np.nan
    result = fit_interrupted_series(
        values,
        TIME,
        config=_config(minimum_valid_fraction=0.9),
    )
    assert result.summary["n_valid"] == 278
    assert np.isfinite(result.coefficients["standard_error_hac"]).all()

    sparse = values.copy()
    sparse[:80] = np.nan
    with pytest.raises(ValueError, match="minimum valid months"):
        fit_interrupted_series(sparse, TIME)


def test_fixed_seed_ar1_no_transition_has_no_bh_discoveries() -> None:
    rng = np.random.default_rng(20260813)
    design, metadata, _ = build_interrupted_design(TIME)
    beta = np.zeros(design.shape[1], dtype="float64")
    index = {name: i for i, name in enumerate(metadata["coefficient"])}
    beta[index["intercept"]] = 1.0
    beta[index["baseline_slope_P1"]] = 0.01
    for month in range(2, 13):
        beta[index[f"calendar_month_{month:02d}"]] = 0.2 * np.sin(
            2.0 * np.pi * (month - 1) / 12.0
        )

    p_values: list[float] = []
    for _ in range(24):
        values = design @ beta + _ar1_noise(rng, TIME.size, 0.8, 0.08)
        result = fit_interrupted_series(
            values, TIME, config=_config(), random_seed=20260900 + len(p_values)
        )
        transitions = result.coefficients["coefficient_type"].isin(
            ["level_change", "slope_change"]
        )
        p_values.extend(result.coefficients.loc[transitions, "p_value"])

    _, significant = benjamini_hochberg(np.asarray(p_values), alpha=0.05)
    assert int(significant.sum()) == 0
