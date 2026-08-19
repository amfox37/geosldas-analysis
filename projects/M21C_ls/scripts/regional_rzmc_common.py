#!/usr/bin/env python3
"""Shared helpers for the regional RZMC observing-system period analysis.

Self-contained so that the regional analysis does not depend on the changepoint
detection module, which is not part of this analysis. `seasonal_adjustment` is
the same trend-preserving calendar-month adjustment used elsewhere in the
project; `check_matches_reference()` verifies that equivalence numerically.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def _time_years(time: pd.DatetimeIndex) -> np.ndarray:
    month_number = time.year.to_numpy() * 12 + time.month.to_numpy() - 1
    return (month_number - month_number[0]).astype("float64") / 12.0


def seasonal_adjustment(values: Any, time: Any) -> tuple[np.ndarray, np.ndarray, float]:
    """Remove fitted calendar-month effects while retaining the global linear trend."""

    values = np.asarray(values, dtype="float64")
    time = pd.DatetimeIndex(np.asarray(time))
    years = _time_years(time)
    month = time.month.to_numpy()
    design_columns = [np.ones(values.size), years]
    for calendar_month in range(2, 13):
        design_columns.append((month == calendar_month).astype("float64"))
    design = np.column_stack(design_columns)
    beta, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    seasonal_component = design[:, 2:] @ beta[2:]
    adjusted = values - seasonal_component
    global_fit = beta[0] + beta[1] * years
    global_residual = adjusted - global_fit
    pair = global_residual[1:] @ global_residual[:-1]
    denominator = global_residual[:-1] @ global_residual[:-1]
    rho = float(np.clip(pair / denominator, -0.98, 0.98)) if denominator > 0 else 0.0
    return adjusted, global_residual, rho


def check_matches_reference(values: Any, time: Any) -> bool:
    """Confirm this implementation matches the project's reference version."""

    from changepoint_detection import seasonal_adjustment as reference

    a, b, c = seasonal_adjustment(values, time)
    ra, rb, rc = reference(values, time)
    return bool(np.allclose(a, ra) and np.allclose(b, rb) and np.isclose(c, rc))
