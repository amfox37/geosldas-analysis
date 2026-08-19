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


# --- shared period-mean statistics (identical settings for every variable) ---

def benjamini_hochberg(p: Any) -> np.ndarray:
    """BH step-up q-values."""

    p = np.asarray(p, dtype="float64"); n = p.size
    order = np.argsort(p); q = np.empty(n); running = 1.0
    for k in range(n - 1, -1, -1):
        running = min(running, p[order[k]] * n / (k + 1))
        q[order[k]] = running
    return q


def period_statistics(values: Any, time: Any, fine) -> dict[str, Any]:
    """Period means and AR(1) effective-sample-size standard errors.

    Seasonality is removed with the trend-preserving calendar-month adjustment.
    A single AR(1) coefficient and residual standard deviation are estimated
    from residuals about the period means, then each period mean is given the
    standard error sd / sqrt(n_eff) with n_eff = n(1-rho)/(1+rho).
    """

    time = pd.DatetimeIndex(np.asarray(time))
    adjusted, _, _ = seasonal_adjustment(values, time)
    series = pd.Series(adjusted, index=time)
    segments = [series[(series.index >= r.start) & (series.index <= r.end)]
                for r in fine.itertuples()]
    residual = np.concatenate([seg.values - seg.mean() for seg in segments])
    centered = residual - residual.mean()
    denominator = max((centered[:-1] ** 2).sum(), 1e-30)
    rho = float(np.clip((centered[1:] * centered[:-1]).sum() / denominator, 0.0, 0.98))
    sd = float(residual.std(ddof=1))
    means = {}
    for r, seg in zip(fine.itertuples(), segments):
        n_eff = max(len(seg) * (1 - rho) / (1 + rho), 1.0)
        means[r.period_id] = (float(seg.mean()), float(sd / np.sqrt(n_eff)), len(seg))
    return {"means": means, "rho": rho, "sd": sd,
            "raw": pd.Series(np.asarray(values, dtype="float64"), index=time)}


def adjacent_differences(stats: dict, order: list, period_ids: list, labels: dict):
    """Adjacent-period differences with BH FDR, one family per transition."""

    from scipy.stats import norm

    rows = []
    for a, b in zip(period_ids[:-1], period_ids[1:]):
        family = []
        for rid in order:
            m1, e1, _ = stats[rid]["means"][a]
            m2, e2, _ = stats[rid]["means"][b]
            d = m2 - m1
            se = float(np.hypot(e1, e2))
            family.append({"region": labels[rid], "region_id": rid,
                           "transition": f"{b}−{a}", "diff": d, "se": se,
                           "p": float(2 * (1 - norm.cdf(abs(d) / se)))})
        for entry, q in zip(family, benjamini_hochberg([f["p"] for f in family])):
            entry["q"] = float(q)
            entry["sig_fdr"] = bool(q < 0.05)
            rows.append(entry)
    return pd.DataFrame(rows)
