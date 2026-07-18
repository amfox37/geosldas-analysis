"""Sparse triple collocation following the existing MATLAB TC scripts."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
import re

import numpy as np

from .climatology import (
    M36_GRID_SIZE,
    PENTAD_CENTER_DOYS,
    PentadClimatology,
    pentad_index,
)
from .generation import load_daily_pair_npz
from .pairs import DailyPair


TC_PAIR_COLUMNS = ((0, 1), (0, 2), (1, 2))
TC_SEASONS = ("annual", "summer", "winter")


@dataclass(frozen=True)
class TripleCollocation:
    """Sparse three-series covariance and TC statistics."""

    idx: np.ndarray
    count: np.ndarray
    variance: np.ndarray
    covariance: np.ndarray
    correlation: np.ndarray
    r2: np.ndarray
    error_variance: np.ndarray
    source_names: tuple[str, str, str]
    source_units: tuple[str, str, str]
    source_scales: tuple[float, float, float]
    run: str
    start_date: date
    end_date: date
    season: str
    min_count: int
    r2_floor: float
    primary_sensor: str
    secondary_sensor: str
    model_values_source: str
    model_climatology_source: str
    primary_climatology_start_date: date
    primary_climatology_end_date: date
    primary_climatology_window_days: int
    primary_climatology_min_count: int
    secondary_climatology_start_date: date
    secondary_climatology_end_date: date
    secondary_climatology_window_days: int
    secondary_climatology_min_count: int
    primary_days: int
    secondary_days: int
    candidate_days: int
    paired_days: int
    contributing_days: int
    missing_days: int
    model_compared_values: int
    model_max_abs_difference: float
    grid_size: int = M36_GRID_SIZE


def tc_output_path(
    output_root: Path | str,
    source_names: Sequence[str],
    run_name: str,
    start_date: date,
    end_date: date,
    *,
    season: str = "annual",
    min_count: int = 20,
) -> Path:
    """Return the standard sparse TC output path."""

    if len(source_names) != 3:
        raise ValueError("Triple collocation requires exactly three source names")
    triplet = "__".join(_path_token(name) for name in source_names)
    return (
        Path(output_root)
        / "step4_tc"
        / triplet
        / run_name
        / f"{start_date:%Y%m%d}_{end_date:%Y%m%d}_{season}_n{min_count}.npz"
    )


def compute_triple_collocation(
    primary_pairs: Iterable[DailyPair],
    primary_climatology: PentadClimatology,
    secondary_pairs: Iterable[DailyPair],
    secondary_climatology: PentadClimatology,
    *,
    primary_name: str,
    secondary_name: str,
    model_name: str = "model",
    primary_scale: float = 1.0,
    secondary_scale: float = 1.0,
    model_scale: float = 1.0,
    model_values_source: str = "primary",
    model_climatology_source: str = "secondary",
    start_date: date | None = None,
    end_date: date | None = None,
    season: str = "annual",
    min_count: int = 20,
    r2_floor: float = 0.001,
    strict_missing: bool = True,
    model_agreement_tolerance: float | None = None,
) -> TripleCollocation:
    """Compute MATLAB two-pair TC statistics from in-memory daily pairs."""

    primary = _pairs_by_date(primary_pairs, "primary")
    secondary = _pairs_by_date(secondary_pairs, "secondary")
    return _compute_from_stores(
        tuple(primary),
        primary.__getitem__,
        primary_climatology,
        tuple(secondary),
        secondary.__getitem__,
        secondary_climatology,
        primary_name=primary_name,
        secondary_name=secondary_name,
        model_name=model_name,
        primary_scale=primary_scale,
        secondary_scale=secondary_scale,
        model_scale=model_scale,
        model_values_source=model_values_source,
        model_climatology_source=model_climatology_source,
        start_date=start_date,
        end_date=end_date,
        season=season,
        min_count=min_count,
        r2_floor=r2_floor,
        strict_missing=strict_missing,
        model_agreement_tolerance=model_agreement_tolerance,
    )


def compute_triple_collocation_from_paths(
    primary_paths: Sequence[Path | str],
    primary_climatology: PentadClimatology,
    secondary_paths: Sequence[Path | str],
    secondary_climatology: PentadClimatology,
    **kwargs,
) -> TripleCollocation:
    """Compute MATLAB two-pair TC statistics from daily pair paths."""

    primary = _paths_by_date(primary_paths, "primary")
    secondary = _paths_by_date(secondary_paths, "secondary")

    def load_primary(day: date) -> DailyPair:
        return _load_pair_for_date(primary[day], day)

    def load_secondary(day: date) -> DailyPair:
        return _load_pair_for_date(secondary[day], day)

    return _compute_from_stores(
        tuple(primary),
        load_primary,
        primary_climatology,
        tuple(secondary),
        load_secondary,
        secondary_climatology,
        **kwargs,
    )


def save_triple_collocation_npz(result: TripleCollocation, path: Path | str) -> None:
    """Write a sparse generic TC result."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        idx0=np.asarray(result.idx, dtype=np.int32),
        N_sm=np.asarray(result.count, dtype=np.int32),
        variance=np.asarray(result.variance, dtype=np.float64),
        covariance=np.asarray(result.covariance, dtype=np.float64),
        correlation=np.asarray(result.correlation, dtype=np.float64),
        R2_TC=np.asarray(result.r2, dtype=np.float64),
        sigma2=np.asarray(result.error_variance, dtype=np.float64),
        source_names=np.asarray(result.source_names),
        source_units=np.asarray(result.source_units),
        source_scales=np.asarray(result.source_scales, dtype=np.float64),
        covariance_pairs=np.asarray(TC_PAIR_COLUMNS, dtype=np.int8),
        run=np.asarray(result.run),
        start_date=np.asarray(result.start_date.strftime("%Y%m%d")),
        end_date=np.asarray(result.end_date.strftime("%Y%m%d")),
        season=np.asarray(result.season),
        Nmin=np.asarray(result.min_count, dtype=np.int32),
        r2_floor=np.asarray(result.r2_floor, dtype=np.float64),
        primary_sensor=np.asarray(result.primary_sensor),
        secondary_sensor=np.asarray(result.secondary_sensor),
        model_values_source=np.asarray(result.model_values_source),
        model_climatology_source=np.asarray(result.model_climatology_source),
        primary_climatology_start_date=np.asarray(
            result.primary_climatology_start_date.strftime("%Y%m%d")
        ),
        primary_climatology_end_date=np.asarray(
            result.primary_climatology_end_date.strftime("%Y%m%d")
        ),
        primary_climatology_window_days=np.asarray(
            result.primary_climatology_window_days, dtype=np.int16
        ),
        primary_climatology_min_count=np.asarray(
            result.primary_climatology_min_count, dtype=np.int32
        ),
        secondary_climatology_start_date=np.asarray(
            result.secondary_climatology_start_date.strftime("%Y%m%d")
        ),
        secondary_climatology_end_date=np.asarray(
            result.secondary_climatology_end_date.strftime("%Y%m%d")
        ),
        secondary_climatology_window_days=np.asarray(
            result.secondary_climatology_window_days, dtype=np.int16
        ),
        secondary_climatology_min_count=np.asarray(
            result.secondary_climatology_min_count, dtype=np.int32
        ),
        primary_days=np.asarray(result.primary_days, dtype=np.int32),
        secondary_days=np.asarray(result.secondary_days, dtype=np.int32),
        candidate_days=np.asarray(result.candidate_days, dtype=np.int32),
        paired_days=np.asarray(result.paired_days, dtype=np.int32),
        contributing_days=np.asarray(result.contributing_days, dtype=np.int32),
        missing_days=np.asarray(result.missing_days, dtype=np.int32),
        model_compared_values=np.asarray(result.model_compared_values, dtype=np.int64),
        model_max_abs_difference=np.asarray(
            result.model_max_abs_difference, dtype=np.float64
        ),
        grid_size=np.asarray(result.grid_size, dtype=np.int32),
    )


def load_triple_collocation_npz(path: Path | str) -> TripleCollocation:
    """Load a sparse generic TC result."""

    with np.load(Path(path)) as z:
        return TripleCollocation(
            idx=np.asarray(z["idx0"], dtype=np.int64),
            count=np.asarray(z["N_sm"], dtype=np.int64),
            variance=np.asarray(z["variance"], dtype=np.float64),
            covariance=np.asarray(z["covariance"], dtype=np.float64),
            correlation=np.asarray(z["correlation"], dtype=np.float64),
            r2=np.asarray(z["R2_TC"], dtype=np.float64),
            error_variance=np.asarray(z["sigma2"], dtype=np.float64),
            source_names=tuple(str(value) for value in z["source_names"]),
            source_units=tuple(str(value) for value in z["source_units"]),
            source_scales=tuple(float(value) for value in z["source_scales"]),
            run=str(z["run"]),
            start_date=_parse_yyyymmdd(str(z["start_date"])),
            end_date=_parse_yyyymmdd(str(z["end_date"])),
            season=str(z["season"]),
            min_count=int(z["Nmin"]),
            r2_floor=float(z["r2_floor"]),
            primary_sensor=str(z["primary_sensor"]),
            secondary_sensor=str(z["secondary_sensor"]),
            model_values_source=str(z["model_values_source"]),
            model_climatology_source=str(z["model_climatology_source"]),
            primary_climatology_start_date=_parse_yyyymmdd(
                str(z["primary_climatology_start_date"])
            ),
            primary_climatology_end_date=_parse_yyyymmdd(
                str(z["primary_climatology_end_date"])
            ),
            primary_climatology_window_days=int(
                z["primary_climatology_window_days"]
            ),
            primary_climatology_min_count=int(z["primary_climatology_min_count"]),
            secondary_climatology_start_date=_parse_yyyymmdd(
                str(z["secondary_climatology_start_date"])
            ),
            secondary_climatology_end_date=_parse_yyyymmdd(
                str(z["secondary_climatology_end_date"])
            ),
            secondary_climatology_window_days=int(
                z["secondary_climatology_window_days"]
            ),
            secondary_climatology_min_count=int(
                z["secondary_climatology_min_count"]
            ),
            primary_days=int(z["primary_days"]),
            secondary_days=int(z["secondary_days"]),
            candidate_days=int(z["candidate_days"]),
            paired_days=int(z["paired_days"]),
            contributing_days=int(z["contributing_days"]),
            missing_days=int(z["missing_days"]),
            model_compared_values=int(z["model_compared_values"]),
            model_max_abs_difference=float(z["model_max_abs_difference"]),
            grid_size=int(z["grid_size"]),
        )


def triple_collocation_npz_is_valid(path: Path | str) -> bool:
    """Return whether a file contains shape-consistent sparse TC fields."""

    path = Path(path)
    if not path.exists():
        return False
    required = (
        "idx0",
        "N_sm",
        "variance",
        "covariance",
        "correlation",
        "R2_TC",
        "sigma2",
        "source_names",
        "Nmin",
    )
    try:
        with np.load(path) as z:
            if not all(field in z for field in required):
                return False
            n = np.asarray(z["idx0"]).size
            return bool(
                np.asarray(z["idx0"]).shape == (n,)
                and np.asarray(z["N_sm"]).shape == (n,)
                and np.asarray(z["variance"]).shape == (n, 3)
                and np.asarray(z["covariance"]).shape == (n, 3)
                and np.asarray(z["correlation"]).shape == (n, 3)
                and np.asarray(z["R2_TC"]).shape == (n, 3)
                and np.asarray(z["sigma2"]).shape == (n, 3)
                and np.asarray(z["source_names"]).shape == (3,)
            )
    except Exception:
        return False


def _compute_from_stores(
    primary_dates: Sequence[date],
    load_primary: Callable[[date], DailyPair],
    primary_climatology: PentadClimatology,
    secondary_dates: Sequence[date],
    load_secondary: Callable[[date], DailyPair],
    secondary_climatology: PentadClimatology,
    *,
    primary_name: str,
    secondary_name: str,
    model_name: str = "model",
    primary_scale: float = 1.0,
    secondary_scale: float = 1.0,
    model_scale: float = 1.0,
    model_values_source: str = "primary",
    model_climatology_source: str = "secondary",
    start_date: date | None = None,
    end_date: date | None = None,
    season: str = "annual",
    min_count: int = 20,
    r2_floor: float = 0.001,
    strict_missing: bool = True,
    model_agreement_tolerance: float | None = None,
) -> TripleCollocation:
    _validate_options(
        primary_name,
        secondary_name,
        model_name,
        primary_scale,
        secondary_scale,
        model_scale,
        model_values_source,
        model_climatology_source,
        season,
        min_count,
        r2_floor,
        model_agreement_tolerance,
    )
    _validate_climatologies(primary_climatology, secondary_climatology)
    primary_date_set = set(primary_dates)
    secondary_date_set = set(secondary_dates)
    all_dates = primary_date_set | secondary_date_set
    if not all_dates:
        raise ValueError("At least one daily pair is required")
    analysis_start = min(all_dates) if start_date is None else start_date
    analysis_end = max(all_dates) if end_date is None else end_date
    if analysis_end < analysis_start:
        raise ValueError(
            f"end date {analysis_end} is before start date {analysis_start}"
        )

    grid_size = primary_climatology.grid_size
    value_sum = np.zeros((3, grid_size), dtype=np.float64)
    value2_sum = np.zeros((3, grid_size), dtype=np.float64)
    product_sum = np.zeros((3, grid_size), dtype=np.float64)
    count = np.zeros(grid_size, dtype=np.int32)
    candidate_days = (analysis_end - analysis_start).days + 1
    paired_days = 0
    contributing_days = 0
    missing_days = 0
    model_compared_values = 0
    model_max_abs_difference = 0.0

    day = analysis_start
    while day <= analysis_end:
        if day not in primary_date_set or day not in secondary_date_set:
            missing_days += 1
            if strict_missing:
                missing = []
                if day not in primary_date_set:
                    missing.append("primary")
                if day not in secondary_date_set:
                    missing.append("secondary")
                raise FileNotFoundError(
                    f"Missing {' and '.join(missing)} daily pair on {day}"
                )
            day += timedelta(days=1)
            continue

        primary = load_primary(day)
        secondary = load_secondary(day)
        _validate_pair(primary, primary_climatology, "primary")
        _validate_pair(secondary, secondary_climatology, "secondary")
        if primary.run != secondary.run:
            raise ValueError(
                f"Run mismatch on {day}: {primary.run!r} != {secondary.run!r}"
            )
        paired_days += 1
        if not _season_includes(day, season):
            day += timedelta(days=1)
            continue

        common, primary_pos, secondary_pos = np.intersect1d(
            np.asarray(primary.idx, dtype=np.int64),
            np.asarray(secondary.idx, dtype=np.int64),
            assume_unique=True,
            return_indices=True,
        )
        if common.size == 0:
            day += timedelta(days=1)
            continue

        primary_model = np.asarray(primary.model, dtype=np.float64)[primary_pos]
        secondary_model = np.asarray(secondary.model, dtype=np.float64)[secondary_pos]
        compare = np.isfinite(primary_model) & np.isfinite(secondary_model)
        if np.any(compare):
            differences = np.abs(primary_model[compare] - secondary_model[compare])
            model_compared_values += int(differences.size)
            model_max_abs_difference = max(
                model_max_abs_difference, float(np.max(differences))
            )
            if (
                model_agreement_tolerance is not None
                and np.any(differences > model_agreement_tolerance)
            ):
                raise ValueError(
                    f"Model values differ by up to {np.max(differences):.6g} "
                    f"on {day}, exceeding tolerance {model_agreement_tolerance}"
                )

        column = pentad_index(day)
        primary_obs = (
            np.asarray(primary.obs, dtype=np.float64)[primary_pos]
            * primary_scale
        )
        secondary_obs = (
            np.asarray(secondary.obs, dtype=np.float64)[secondary_pos]
            * secondary_scale
        )
        model = primary_model if model_values_source == "primary" else secondary_model
        model = model * model_scale
        primary_obs -= _climatology_values(
            primary_climatology, common, column, "obs"
        ) * primary_scale
        secondary_obs -= _climatology_values(
            secondary_climatology, common, column, "obs"
        ) * secondary_scale
        model_climatology = (
            primary_climatology
            if model_climatology_source == "primary"
            else secondary_climatology
        )
        model -= _climatology_values(
            model_climatology, common, column, "model"
        ) * model_scale

        values = np.vstack((primary_obs, model, secondary_obs))
        valid = np.all(np.isfinite(values), axis=0)
        if not np.any(valid):
            day += timedelta(days=1)
            continue
        contributing_days += 1
        idx = common[valid]
        values = values[:, valid]
        value_sum[:, idx] += values
        value2_sum[:, idx] += values**2
        product_sum[0, idx] += values[0] * values[1]
        product_sum[1, idx] += values[0] * values[2]
        product_sum[2, idx] += values[1] * values[2]
        count[idx] += 1
        day += timedelta(days=1)

    variance, covariance, correlation, r2, error_variance = _statistics_from_sums(
        value_sum, value2_sum, product_sum, count, min_count, r2_floor
    )
    occupied = count > 0
    source_units = (
        primary_climatology.obs_units,
        primary_climatology.model_units
        if model_values_source == "primary"
        else secondary_climatology.model_units,
        secondary_climatology.obs_units,
    )
    return TripleCollocation(
        idx=np.flatnonzero(occupied),
        count=count[occupied],
        variance=variance[:, occupied].T,
        covariance=covariance[:, occupied].T,
        correlation=correlation[:, occupied].T,
        r2=r2[:, occupied].T,
        error_variance=error_variance[:, occupied].T,
        source_names=(primary_name, model_name, secondary_name),
        source_units=source_units,
        source_scales=(primary_scale, model_scale, secondary_scale),
        run=primary_climatology.run,
        start_date=analysis_start,
        end_date=analysis_end,
        season=season,
        min_count=min_count,
        r2_floor=r2_floor,
        primary_sensor=primary_climatology.sensor,
        secondary_sensor=secondary_climatology.sensor,
        model_values_source=model_values_source,
        model_climatology_source=model_climatology_source,
        primary_climatology_start_date=primary_climatology.start_date,
        primary_climatology_end_date=primary_climatology.end_date,
        primary_climatology_window_days=primary_climatology.window_days,
        primary_climatology_min_count=primary_climatology.min_count,
        secondary_climatology_start_date=secondary_climatology.start_date,
        secondary_climatology_end_date=secondary_climatology.end_date,
        secondary_climatology_window_days=secondary_climatology.window_days,
        secondary_climatology_min_count=secondary_climatology.min_count,
        primary_days=sum(
            analysis_start <= value <= analysis_end for value in primary_date_set
        ),
        secondary_days=sum(
            analysis_start <= value <= analysis_end for value in secondary_date_set
        ),
        candidate_days=candidate_days,
        paired_days=paired_days,
        contributing_days=contributing_days,
        missing_days=missing_days,
        model_compared_values=model_compared_values,
        model_max_abs_difference=(
            model_max_abs_difference if model_compared_values else np.nan
        ),
        grid_size=grid_size,
    )


def _statistics_from_sums(
    value_sum: np.ndarray,
    value2_sum: np.ndarray,
    product_sum: np.ndarray,
    count: np.ndarray,
    min_count: int,
    r2_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        nn = np.where(count >= min_count, count, np.nan).astype(np.float64)
        mean = value_sum / nn
        mean[np.isinf(mean)] = np.nan
        variance = value2_sum / nn - mean**2
        variance[variance < 0.0] = np.nan
        covariance = np.empty_like(product_sum)
        covariance[0] = product_sum[0] / nn - mean[0] * mean[1]
        covariance[1] = product_sum[1] / nn - mean[0] * mean[2]
        covariance[2] = product_sum[2] / nn - mean[1] * mean[2]
        correlation = np.empty_like(covariance)
        correlation[0] = covariance[0] / np.sqrt(variance[0] * variance[1])
        correlation[1] = covariance[1] / np.sqrt(variance[0] * variance[2])
        correlation[2] = covariance[2] / np.sqrt(variance[1] * variance[2])

        r2 = np.empty_like(variance)
        r2[0] = covariance[0] * covariance[1] / covariance[2] / variance[0]
        r2[1] = covariance[0] * covariance[2] / covariance[1] / variance[1]
        r2[2] = covariance[2] * covariance[1] / covariance[0] / variance[2]
        r2[r2 < r2_floor] = np.nan
        r2[r2 > 1.0] = 1.0

        error_variance = np.empty_like(variance)
        error_variance[0] = (
            variance[0] - covariance[1] * covariance[0] / covariance[2]
        )
        error_variance[1] = (
            variance[1] - covariance[2] * covariance[0] / covariance[1]
        )
        error_variance[2] = (
            variance[2] - covariance[1] * covariance[2] / covariance[0]
        )
    return variance, covariance, correlation, r2, error_variance


def _climatology_values(
    climatology: PentadClimatology, idx: np.ndarray, column: int, field: str
) -> np.ndarray:
    rows = np.searchsorted(climatology.idx, idx)
    found = rows < climatology.idx.size
    found[found] &= climatology.idx[rows[found]] == idx[found]
    values = np.full(idx.size, np.nan, dtype=np.float64)
    source = climatology.model if field == "model" else climatology.obs
    values[found] = source[rows[found], column]
    return values


def _validate_climatologies(
    primary: PentadClimatology, secondary: PentadClimatology
) -> None:
    for label, climatology in (("primary", primary), ("secondary", secondary)):
        idx = np.asarray(climatology.idx)
        expected = (idx.size, PENTAD_CENTER_DOYS.size)
        if idx.ndim != 1 or np.any(np.diff(idx) <= 0):
            raise ValueError(f"{label} climatology indices must be sorted unique")
        if np.any(idx < 0) or np.any(idx >= climatology.grid_size):
            raise ValueError(f"{label} climatology index is outside the grid")
        if climatology.obs.shape != expected or climatology.model.shape != expected:
            raise ValueError(f"{label} climatology arrays must have shape {expected}")
        if not np.array_equal(climatology.center_doy, PENTAD_CENTER_DOYS):
            raise ValueError(f"{label} climatology has non-MATLAB pentad centers")
    if primary.grid_size != secondary.grid_size:
        raise ValueError("Climatology grid sizes differ")
    if primary.run != secondary.run:
        raise ValueError(f"Climatology runs differ: {primary.run!r} != {secondary.run!r}")


def _validate_pair(pair: DailyPair, climatology: PentadClimatology, label: str) -> None:
    idx = np.asarray(pair.idx)
    obs = np.asarray(pair.obs)
    model = np.asarray(pair.model)
    if idx.ndim != 1 or obs.ndim != 1 or model.ndim != 1:
        raise ValueError(f"{label} pair arrays must be one-dimensional on {pair.date}")
    if not (idx.size == obs.size == model.size):
        raise ValueError(f"{label} pair arrays have unequal lengths on {pair.date}")
    if np.unique(idx).size != idx.size:
        raise ValueError(f"{label} pair has duplicate grid indices on {pair.date}")
    if np.any(idx < 0) or np.any(idx >= climatology.grid_size):
        raise ValueError(f"{label} pair index is outside the grid on {pair.date}")
    metadata = (pair.sensor, pair.run, pair.obs_units, pair.model_units)
    expected = (
        climatology.sensor,
        climatology.run,
        climatology.obs_units,
        climatology.model_units,
    )
    if metadata != expected:
        raise ValueError(
            f"{label} pair metadata does not match climatology on {pair.date}; "
            f"expected {expected}, got {metadata}"
        )


def _validate_options(
    primary_name: str,
    secondary_name: str,
    model_name: str,
    primary_scale: float,
    secondary_scale: float,
    model_scale: float,
    model_values_source: str,
    model_climatology_source: str,
    season: str,
    min_count: int,
    r2_floor: float,
    model_agreement_tolerance: float | None,
) -> None:
    names = (primary_name.strip(), model_name.strip(), secondary_name.strip())
    if any(not name for name in names) or len(set(names)) != 3:
        raise ValueError("TC source names must be non-empty and unique")
    scales = (primary_scale, model_scale, secondary_scale)
    if any(not np.isfinite(scale) or scale == 0 for scale in scales):
        raise ValueError("TC source scales must be finite and nonzero")
    if model_values_source not in ("primary", "secondary"):
        raise ValueError("model_values_source must be primary or secondary")
    if model_climatology_source not in ("primary", "secondary"):
        raise ValueError("model_climatology_source must be primary or secondary")
    if season not in TC_SEASONS:
        raise ValueError(f"season must be one of {TC_SEASONS}")
    if min_count < 0:
        raise ValueError("min_count must be non-negative")
    if not np.isfinite(r2_floor):
        raise ValueError("r2_floor must be finite")
    if model_agreement_tolerance is not None and model_agreement_tolerance < 0:
        raise ValueError("model_agreement_tolerance must be non-negative")


def _season_includes(day: date, season: str) -> bool:
    if season == "annual":
        return True
    if season == "summer":
        return 5 <= day.month <= 9
    return day.month >= 10 or day.month <= 4


def _pairs_by_date(pairs: Iterable[DailyPair], label: str) -> dict[date, DailyPair]:
    result: dict[date, DailyPair] = {}
    for pair in pairs:
        if pair.date in result:
            raise ValueError(f"Duplicate {label} daily pair date: {pair.date}")
        result[pair.date] = pair
    return result


def _paths_by_date(paths: Sequence[Path | str], label: str) -> dict[date, Path]:
    result: dict[date, Path] = {}
    for value in paths:
        path = Path(value)
        day = _parse_yyyymmdd(path.stem)
        if day in result:
            raise ValueError(f"Duplicate {label} daily pair date: {day}")
        result[day] = path
    return result


def _load_pair_for_date(path: Path, day: date) -> DailyPair:
    pair = load_daily_pair_npz(path)
    if pair.date != day:
        raise ValueError(
            f"Pair metadata date {pair.date} does not match filename date {day}"
        )
    return pair


def _path_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    if not token:
        raise ValueError(f"Cannot make output path token from {value!r}")
    return token


def _parse_yyyymmdd(value: str) -> date:
    value = value.strip()
    if len(value) != 8 or not value.isdigit():
        raise ValueError(f"Expected YYYYMMDD date, got {value!r}")
    return date(int(value[:4]), int(value[4:6]), int(value[6:8]))
