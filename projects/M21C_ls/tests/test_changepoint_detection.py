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

from changepoint_detection import (
    compare_known_boundaries,
    detect_changepoints,
    known_boundary_table,
    load_changepoint_config,
    seasonal_adjustment,
)


TIME = pd.date_range("2000-06-01", periods=288, freq="MS")


def _fast_config() -> dict:
    config = copy.deepcopy(load_changepoint_config())
    config["penalty"]["multipliers"] = [1.0]
    config["penalty"]["primary_multiplier"] = 1.0
    config["acceptance"]["minimum_penalty_stability_fraction"] = 1.0
    return config


def _ar1(seed: int, phi: float = 0.8) -> np.ndarray:
    rng = np.random.default_rng(seed)
    innovation = rng.normal(size=408)
    values = np.zeros(408)
    for index in range(1, values.size):
        values[index] = phi * values[index - 1] + innovation[index]
    return values[-288:]


def _series(seed: int) -> np.ndarray:
    month = np.arange(288)
    return 10.0 + 0.02 * month / 12.0 + 2.0 * np.sin(2 * np.pi * month / 12) + _ar1(seed)


def test_seasonal_adjustment_preserves_linear_trend() -> None:
    month = np.arange(288)
    values = 4.0 + 0.3 * month / 12.0 + 3.0 * np.sin(2 * np.pi * month / 12)

    adjusted, residual, rho = seasonal_adjustment(values, TIME)

    assert np.polyfit(month / 12.0, adjusted, 1)[0] == pytest.approx(0.3)
    assert np.max(np.abs(residual)) < 1e-10
    assert abs(rho) < 0.98


def test_ar1_null_has_no_accepted_detection() -> None:
    result = detect_changepoints(_series(1), TIME, config=_fast_config())

    assert result.summary["n_accepted_detections"] == 0


def test_strong_level_change_is_recovered() -> None:
    values = _series(2)
    values[120:] += 8.0

    result = detect_changepoints(values, TIME, config=_fast_config())
    accepted = result.detections[result.detections["accepted_detection"]]

    assert np.min(np.abs(accepted["break_index"] - 120)) <= 3


def test_p7_is_detection_exempt_but_p8_is_scored() -> None:
    boundaries = known_boundary_table(TIME)
    p7 = boundaries.set_index("period_id").loc["P7"]
    detection = pd.DataFrame(
        {
            "break_index": [p7["boundary_index"]],
            "break_date": [p7["boundary_date"]],
            "accepted_detection": [True],
        }
    )

    comparison = compare_known_boundaries(detection, TIME, config=_fast_config())
    comparison = comparison.set_index("period_id")

    assert not bool(comparison.loc["P7", "scored"])
    assert bool(comparison.loc["P7", "matched_within_primary_tolerance"])
    assert bool(comparison.loc["P8", "scored"])
