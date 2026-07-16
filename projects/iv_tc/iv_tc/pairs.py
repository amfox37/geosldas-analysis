"""Shared records passed between IV/TC readers and statistics code."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

import numpy as np


@dataclass(frozen=True)
class SparseObservation:
    """One day's observation values on sparse grid-cell indices."""

    date: date
    sensor: str
    idx: np.ndarray
    values: np.ndarray
    units: str
    qc_summary: dict[str, int | float | str] = field(default_factory=dict)


@dataclass(frozen=True)
class DailyPair:
    """One day's matched observation/model values."""

    date: date
    sensor: str
    run: str
    idx: np.ndarray
    obs: np.ndarray
    model: np.ndarray
    obs_units: str
    model_units: str

