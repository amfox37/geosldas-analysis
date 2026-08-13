#!/usr/bin/env python3
"""Independent PELT changepoint detection for M21C monthly domain series."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import ruptures as rpt
from packaging.version import Version
from ruptures.base import BaseCost
from ruptures.exceptions import NotEnoughPoints

from m21c_periods import DEFAULT_REGISTRY, load_period_frames


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "changepoint_detection.json"


@dataclass
class ChangepointResult:
    """Detection and preprocessing diagnostics for one monthly series."""

    detections: pd.DataFrame
    penalty_detections: pd.DataFrame
    summary: dict[str, Any]
    seasonal_adjusted: np.ndarray
    global_residual: np.ndarray


class CostARLinear(BaseCost):
    """Segment cost for an AR(1) response with intercept and linear time trend."""

    model = "ar1_linear"

    def __init__(self, min_size: int = 24) -> None:
        self.min_size = int(min_size)
        self.signal: np.ndarray | None = None
        self.covariates: np.ndarray | None = None

    def fit(self, signal: Any) -> "CostARLinear":
        values = np.asarray(signal, dtype="float64").reshape(-1)
        if not np.all(np.isfinite(values)):
            raise ValueError("CostARLinear requires a complete finite series")
        time = np.arange(values.size, dtype="float64") / 12.0
        lagged = np.r_[values[0], values[:-1]]
        self.signal = values
        self.covariates = np.column_stack([np.ones(values.size), time, lagged])
        return self

    def error(self, start: int, end: int) -> float:
        if end - start < self.min_size:
            raise NotEnoughPoints
        if self.signal is None or self.covariates is None:
            raise RuntimeError("CostARLinear must be fitted before error evaluation")
        response = self.signal[start:end]
        design = self.covariates[start:end]
        beta, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
        residual_sum_squares = float(np.sum((response - design @ beta) ** 2))
        variance = max(residual_sum_squares / (end - start), np.finfo(float).tiny)
        return float((end - start) * math.log(variance))


class CostProfileLinear(BaseCost):
    """Gaussian profile-likelihood cost for a response and linear covariates."""

    model = "profile_linear"

    def __init__(self, min_size: int = 24) -> None:
        self.min_size = int(min_size)
        self.signal: np.ndarray | None = None
        self.response: np.ndarray | None = None
        self.covariates: np.ndarray | None = None

    def fit(self, signal: Any) -> "CostProfileLinear":
        signal = np.asarray(signal, dtype="float64")
        if signal.ndim != 2 or signal.shape[1] < 2:
            raise ValueError("CostProfileLinear requires response plus covariates")
        if not np.all(np.isfinite(signal)):
            raise ValueError("CostProfileLinear requires a complete finite signal")
        self.signal = signal
        self.response = signal[:, 0]
        self.covariates = signal[:, 1:]
        return self

    def error(self, start: int, end: int) -> float:
        if end - start < self.min_size:
            raise NotEnoughPoints
        if self.response is None or self.covariates is None:
            raise RuntimeError("CostProfileLinear must be fitted before error evaluation")
        response = self.response[start:end]
        design = self.covariates[start:end]
        beta, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
        residual_sum_squares = float(np.sum((response - design @ beta) ** 2))
        variance = max(residual_sum_squares / (end - start), np.finfo(float).tiny)
        return float((end - start) * math.log(variance))


def load_changepoint_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load and validate the independent changepoint contract."""

    config = json.loads(Path(path).read_text())
    if config["library"] != "ruptures":
        raise ValueError("Only the ruptures changepoint backend is supported")
    installed = Version(rpt.__version__.lstrip("v"))
    required = Version(str(config["minimum_library_version"]))
    if installed < required:
        raise ValueError(f"ruptures {installed} is older than required {required}")
    if int(config["minimum_segment_months"]) < 24:
        raise ValueError("minimum_segment_months must be at least 24")
    if float(config["minimum_valid_fraction"]) != 1.0:
        raise ValueError("This implementation currently requires complete domain series")
    methods = config["methods"]
    if methods["primary"] != "pelt_ar1_linear":
        raise ValueError("Primary method must be pelt_ar1_linear")
    if methods["sensitivity"] != "pelt_prewhitened_linear":
        raise ValueError("Sensitivity method must be pelt_prewhitened_linear")
    if methods["prewhitening_ar1_estimator"] != "median_overlapping_locally_detrended_windows":
        raise ValueError("Unsupported prewhitening AR(1) estimator")
    window = int(methods["ar1_window_months"])
    step = int(methods["ar1_window_step_months"])
    if window < int(config["minimum_segment_months"]) or not 0 < step <= window:
        raise ValueError("Invalid local AR(1) window or step")
    multipliers = [float(value) for value in config["penalty"]["multipliers"]]
    if sorted(set(multipliers)) != multipliers or not multipliers or min(multipliers) <= 0:
        raise ValueError("Penalty multipliers must be unique, increasing, and positive")
    if float(config["penalty"]["primary_multiplier"]) not in multipliers:
        raise ValueError("Primary penalty multiplier must be present in the grid")
    return config


def _time_years(time: pd.DatetimeIndex) -> np.ndarray:
    month_number = time.year.to_numpy() * 12 + time.month.to_numpy() - 1
    return (month_number - month_number[0]).astype("float64") / 12.0


def seasonal_adjustment(
    values: Any, time: Any
) -> tuple[np.ndarray, np.ndarray, float]:
    """Remove fitted month effects while retaining the global linear trend."""

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


def _standardize(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    center = float(np.mean(values))
    scale = float(np.std(values, ddof=1))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Series has no finite variation for changepoint detection")
    return (values - center) / scale, center, scale


def robust_local_ar1(values: np.ndarray, *, window: int, step: int) -> float:
    """Estimate persistence without letting long structural changes dominate."""

    estimates: list[float] = []
    local_time = np.arange(window, dtype="float64")
    for start in range(0, values.size - window + 1, step):
        segment = values[start : start + window]
        beta = np.polyfit(local_time, segment, 1)
        residual = segment - np.polyval(beta, local_time)
        denominator = float(residual[:-1] @ residual[:-1])
        if denominator > 0:
            estimates.append(float(residual[1:] @ residual[:-1] / denominator))
    if not estimates:
        raise ValueError("No local windows were available for AR(1) estimation")
    return float(np.clip(np.median(estimates), -0.98, 0.98))


def _prais_winsten_signal(values: np.ndarray, rho: float) -> np.ndarray:
    n = values.size
    years = np.arange(n, dtype="float64") / 12.0
    design = np.column_stack([np.ones(n), years])
    transformed_values = np.empty(n, dtype="float64")
    transformed_design = np.empty_like(design)
    first_scale = math.sqrt(max(0.0, 1.0 - rho * rho))
    transformed_values[0] = first_scale * values[0]
    transformed_design[0] = first_scale * design[0]
    transformed_values[1:] = values[1:] - rho * values[:-1]
    transformed_design[1:] = design[1:] - rho * design[:-1]
    return np.column_stack([transformed_values, transformed_design])


def _method_breaks(
    standardized: np.ndarray,
    *,
    method: str,
    rho: float,
    minimum_segment: int,
    multiplier: float,
) -> tuple[list[int], float]:
    n = standardized.size
    if method == "pelt_ar1_linear":
        parameter_count = 3
        signal: np.ndarray = standardized
        cost: BaseCost = CostARLinear(min_size=minimum_segment)
    elif method == "pelt_prewhitened_linear":
        parameter_count = 2
        signal = _prais_winsten_signal(standardized, rho)
        cost = CostProfileLinear(min_size=minimum_segment)
    else:
        raise ValueError(f"Unknown changepoint method: {method}")
    penalty = float(multiplier * parameter_count * math.log(n))
    endpoints = rpt.Pelt(
        custom_cost=cost,
        min_size=minimum_segment,
        jump=1,
    ).fit(signal).predict(pen=penalty)
    return [int(endpoint) for endpoint in endpoints if endpoint < n], penalty


def _nearest_distance(index: int, candidates: list[int]) -> float:
    if not candidates:
        return np.inf
    return float(min(abs(index - candidate) for candidate in candidates))


def detect_changepoints(
    values: Any,
    time: Any,
    *,
    config: dict[str, Any] | None = None,
    config_path: str | Path = DEFAULT_CONFIG,
) -> ChangepointResult:
    """Detect stable, cross-method changepoints in one complete monthly series."""

    settings = config if config is not None else load_changepoint_config(config_path)
    values = np.asarray(values, dtype="float64")
    time_index = pd.DatetimeIndex(np.asarray(time))
    if values.ndim != 1 or values.size != time_index.size:
        raise ValueError("values and time must be matching one-dimensional arrays")
    if values.size != int(settings["record_length_months"]):
        raise ValueError("Series length differs from the changepoint contract")
    if not np.all(np.isfinite(values)):
        raise ValueError("Independent changepoint detection requires complete finite values")
    adjusted, global_residual, global_rho = seasonal_adjustment(values, time_index)
    rho = robust_local_ar1(
        adjusted,
        window=int(settings["methods"]["ar1_window_months"]),
        step=int(settings["methods"]["ar1_window_step_months"]),
    )
    standardized, center, scale = _standardize(adjusted)
    methods = [settings["methods"]["primary"], settings["methods"]["sensitivity"]]
    multipliers = [float(value) for value in settings["penalty"]["multipliers"]]
    minimum_segment = int(settings["minimum_segment_months"])
    break_grid: dict[tuple[str, float], list[int]] = {}
    penalties: dict[tuple[str, float], float] = {}
    for method in methods:
        for multiplier in multipliers:
            breaks, penalty = _method_breaks(
                standardized,
                method=method,
                rho=rho,
                minimum_segment=minimum_segment,
                multiplier=multiplier,
            )
            break_grid[(method, multiplier)] = breaks
            penalties[(method, multiplier)] = penalty

    penalty_rows = [
        {
            "method": method,
            "penalty_multiplier": multiplier,
            "penalty": penalties[(method, multiplier)],
            "break_index": position,
            "break_date": time_index[position],
        }
        for method in methods
        for multiplier in multipliers
        for position in break_grid[(method, multiplier)]
    ]
    penalty_detections = pd.DataFrame(
        penalty_rows,
        columns=[
            "method",
            "penalty_multiplier",
            "penalty",
            "break_index",
            "break_date",
        ],
    )

    primary_multiplier = float(settings["penalty"]["primary_multiplier"])
    match_tolerance = int(
        settings["acceptance"]["penalty_match_tolerance_months"]
    )
    consensus_tolerance = int(
        settings["acceptance"]["method_consensus_tolerance_months"]
    )
    stability_floor = float(
        settings["acceptance"]["minimum_penalty_stability_fraction"]
    )
    rows: list[dict[str, Any]] = []
    for method in methods:
        other_method = next(candidate for candidate in methods if candidate != method)
        primary_breaks = break_grid[(method, primary_multiplier)]
        other_breaks = break_grid[(other_method, primary_multiplier)]
        for position in primary_breaks:
            stable_count = sum(
                _nearest_distance(position, break_grid[(method, multiplier)])
                <= match_tolerance
                for multiplier in multipliers
            )
            stability = stable_count / len(multipliers)
            other_distance = _nearest_distance(position, other_breaks)
            consensus = other_distance <= consensus_tolerance
            rows.append(
                {
                    "method": method,
                    "break_index": position,
                    "break_date": time_index[position],
                    "primary_penalty_multiplier": primary_multiplier,
                    "primary_penalty": penalties[(method, primary_multiplier)],
                    "penalty_stability_fraction": stability,
                    "nearest_other_method_distance_months": other_distance,
                    "method_consensus": consensus,
                    "accepted_detection": (
                        method == settings["methods"]["primary"]
                        and stability >= stability_floor
                        and consensus
                    ),
                }
            )
    detections = pd.DataFrame(rows)
    if detections.empty:
        detections = pd.DataFrame(
            columns=[
                "method",
                "break_index",
                "break_date",
                "primary_penalty_multiplier",
                "primary_penalty",
                "penalty_stability_fraction",
                "nearest_other_method_distance_months",
                "method_consensus",
                "accepted_detection",
            ]
        )
    summary = {
        "n_time": values.size,
        "n_valid": int(np.isfinite(values).sum()),
        "estimated_global_ar1": global_rho,
        "estimated_local_detrended_ar1": rho,
        "standardization_center": center,
        "standardization_scale": scale,
        "n_primary_detections": len(
            break_grid[(settings["methods"]["primary"], primary_multiplier)]
        ),
        "n_sensitivity_detections": len(
            break_grid[(settings["methods"]["sensitivity"], primary_multiplier)]
        ),
        "n_accepted_detections": int(detections["accepted_detection"].sum()),
        "configuration": json.dumps(settings, sort_keys=True),
    }
    return ChangepointResult(
        detections=detections,
        penalty_detections=penalty_detections,
        summary=summary,
        seasonal_adjusted=adjusted,
        global_residual=global_residual,
    )


def known_boundary_table(
    time: Any,
    *,
    registry_path: str | Path = DEFAULT_REGISTRY,
) -> pd.DataFrame:
    """Return P2-P9 boundary indices and detection-scoring exemptions."""

    time_index = pd.DatetimeIndex(np.asarray(time))
    _, fine, _, _ = load_period_frames(registry_path)
    rows: list[dict[str, Any]] = []
    for row in fine.iloc[1:].itertuples():
        boundary = pd.Timestamp(row.start)
        matches = np.flatnonzero(time_index == boundary)
        if matches.size != 1:
            raise ValueError(f"Boundary {boundary:%Y-%m} is absent from the time axis")
        rows.append(
            {
                "period_id": row.period_id,
                "boundary_date": boundary,
                "boundary_index": int(matches[0]),
                "boundary_rationale": row.boundary_rationale,
                "detection_exempt": bool(row.changepoint_detection_exempt),
            }
        )
    return pd.DataFrame(rows)


def compare_known_boundaries(
    detections: pd.DataFrame,
    time: Any,
    *,
    config: dict[str, Any] | None = None,
    config_path: str | Path = DEFAULT_CONFIG,
    registry_path: str | Path = DEFAULT_REGISTRY,
) -> pd.DataFrame:
    """Compare accepted detections one-to-one with known P2-P9 boundaries."""

    settings = config if config is not None else load_changepoint_config(config_path)
    boundaries = known_boundary_table(time, registry_path=registry_path)
    accepted = detections.loc[detections["accepted_detection"]].copy()
    accepted = accepted.sort_values("break_index")
    unmatched = set(accepted.index)
    rows: list[dict[str, Any]] = []
    primary_tolerance = int(
        settings["known_boundary_comparison"]["primary_tolerance_months"]
    )
    sensitivity_tolerance = int(
        settings["known_boundary_comparison"]["sensitivity_tolerance_months"]
    )
    for boundary in boundaries.itertuples():
        candidates = sorted(
            (
                abs(int(accepted.loc[index, "break_index"]) - boundary.boundary_index),
                index,
            )
            for index in unmatched
        )
        if candidates and candidates[0][0] <= sensitivity_tolerance:
            distance, matched_index = candidates[0]
            matched = accepted.loc[matched_index]
            unmatched.remove(matched_index)
            detected_date: pd.Timestamp | pd.NaT = pd.Timestamp(matched["break_date"])
            signed_offset = int(matched["break_index"]) - boundary.boundary_index
        else:
            distance = np.nan
            detected_date = pd.NaT
            signed_offset = np.nan
        rows.append(
            {
                **boundary._asdict(),
                "detected_date": detected_date,
                "signed_offset_months": signed_offset,
                "absolute_offset_months": distance,
                "matched_within_primary_tolerance": (
                    bool(distance <= primary_tolerance) if np.isfinite(distance) else False
                ),
                "matched_within_sensitivity_tolerance": np.isfinite(distance),
                "scored": not bool(boundary.detection_exempt),
            }
        )
    return pd.DataFrame(rows)
