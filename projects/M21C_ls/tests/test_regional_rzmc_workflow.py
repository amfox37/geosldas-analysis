from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_regional_rzmc_transitions import aggregate_adjusted_tile_series  # noqa: E402
from regional_rzmc_common import period_statistics_from_adjusted  # noqa: E402
from trend_statistics import calendar_month_anomalies, monthly_time_axis  # noqa: E402


def test_regional_aggregation_uses_production_tile_adjustment() -> None:
    time = pd.date_range("2000-01-01", periods=60, freq="MS")
    years, months = monthly_time_axis(time)
    seasonal = np.sin(2 * np.pi * (months - 1) / 12)
    values = np.vstack(
        [
            0.2 * years + seasonal,
            -0.1 * years + 2.0 * seasonal + 0.3,
        ]
    )
    weights = np.array([[0.25, 0.75], [1.0, 0.0]])

    raw, adjusted = aggregate_adjusted_tile_series(
        values,
        weights,
        time,
        minimum_samples=3,
        n_jobs=1,
        chunk_size=1,
    )

    expected_tiles = np.vstack(
        [
            calendar_month_anomalies(
                row, years, months, np.ones(row.size, dtype=bool), 3
            )[0]
            for row in values
        ]
    )
    np.testing.assert_allclose(raw, weights @ values)
    np.testing.assert_allclose(adjusted, weights @ expected_tiles)


def test_period_statistics_do_not_adjust_series_again() -> None:
    time = pd.date_range("2000-01-01", periods=48, freq="MS")
    values = np.linspace(-0.02, 0.03, time.size)
    periods = pd.DataFrame(
        {
            "period_id": ["P1", "P2"],
            "start": [time[0], time[24]],
            "end": [time[23], time[-1]],
        }
    )

    result = period_statistics_from_adjusted(values, time, periods)

    assert result["means"]["P1"][0] == np.mean(values[:24])
    assert result["means"]["P2"][0] == np.mean(values[24:])
    np.testing.assert_array_equal(result["monthly"].values, values)
