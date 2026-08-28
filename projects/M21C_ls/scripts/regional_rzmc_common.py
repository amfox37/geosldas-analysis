#!/usr/bin/env python3
"""Shared statistics for tile-adjusted regional soil-moisture series."""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


# --- shared period-mean statistics (identical settings for every variable) ---

def benjamini_hochberg(p: Any) -> np.ndarray:
    """BH step-up q-values."""

    p = np.asarray(p, dtype="float64"); n = p.size
    order = np.argsort(p); q = np.empty(n); running = 1.0
    for k in range(n - 1, -1, -1):
        running = min(running, p[order[k]] * n / (k + 1))
        q[order[k]] = running
    return q


def period_statistics_from_adjusted(values: Any, time: Any, fine) -> dict[str, Any]:
    """Period means and AR(1) standard errors from an already-adjusted series.

    A single AR(1) coefficient and residual standard deviation are estimated
    from residuals about the period means, then each period mean is given the
    standard error sd / sqrt(n_eff) with n_eff = n(1-rho)/(1+rho).
    """

    time = pd.DatetimeIndex(np.asarray(time))
    values = np.asarray(values, dtype="float64")
    if not np.all(np.isfinite(values)):
        raise ValueError("period statistics require a complete adjusted series")
    series = pd.Series(values, index=time)
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
            "monthly": series}


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
