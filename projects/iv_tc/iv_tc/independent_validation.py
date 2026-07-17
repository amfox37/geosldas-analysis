"""MATLAB-compatible independent validation for sparse daily pairs."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np

from .climatology import (
    M36_GRID_SIZE,
    PENTAD_CENTER_DOYS,
    PentadClimatology,
    pentad_index,
)
from .generation import canonical_sensor, load_daily_pair_npz
from .pairs import DailyPair


MATLAB_UNDEFINED_CORRELATION = -9999.0


@dataclass(frozen=True)
class IndependentValidation:
    """Sparse MATLAB step-4 IVd/IVs statistics."""

    idx: np.ndarray
    count: np.ndarray
    r2_ivd_model: np.ndarray
    r2_ivd_obs: np.ndarray
    r2_ivs_model: np.ndarray
    r2_ivs_obs: np.ndarray
    r_model_obs: np.ndarray
    sensor: str
    run: str
    obs_units: str
    model_units: str
    start_date: date
    end_date: date
    lag_days: int
    min_count: int
    source_days: int
    candidate_days: int
    paired_days: int
    contributing_days: int
    climatology_start_date: date
    climatology_end_date: date
    climatology_window_days: int
    climatology_min_count: int
    grid_size: int = M36_GRID_SIZE


def iv_output_path(
    output_root: Path | str,
    sensor: str,
    run_name: str,
    start_date: date,
    end_date: date,
    lag_days: int = 2,
    min_count: int = 100,
) -> Path:
    """Return the standard sparse step-4 output path."""

    return (
        Path(output_root)
        / "step4_iv"
        / canonical_sensor(sensor)
        / run_name
        / (
            f"{start_date:%Y%m%d}_{end_date:%Y%m%d}"
            f"_lag{lag_days}_n{min_count}.npz"
        )
    )


def compute_independent_validation(
    pairs: Iterable[DailyPair],
    climatology: PentadClimatology,
    *,
    lag_days: int = 2,
    min_count: int = 100,
    start_date: date | None = None,
    end_date: date | None = None,
) -> IndependentValidation:
    """Compute step-4 statistics from in-memory daily pairs."""

    pair_by_date: dict[date, DailyPair] = {}
    for pair in pairs:
        if pair.date in pair_by_date:
            raise ValueError(f"Duplicate daily pair date: {pair.date}")
        pair_by_date[pair.date] = pair
    return _compute_from_store(
        tuple(pair_by_date),
        pair_by_date.__getitem__,
        climatology,
        lag_days=lag_days,
        min_count=min_count,
        start_date=start_date,
        end_date=end_date,
    )


def compute_independent_validation_from_paths(
    paths: Sequence[Path | str],
    climatology: PentadClimatology,
    *,
    lag_days: int = 2,
    min_count: int = 100,
    start_date: date | None = None,
    end_date: date | None = None,
) -> IndependentValidation:
    """Compute step-4 statistics while loading only each required pair."""

    path_by_date: dict[date, Path] = {}
    for value in paths:
        path = Path(value)
        day = _parse_yyyymmdd(path.stem)
        if day in path_by_date:
            raise ValueError(f"Duplicate daily pair date: {day}")
        path_by_date[day] = path

    def load(day: date) -> DailyPair:
        pair = load_daily_pair_npz(path_by_date[day])
        if pair.date != day:
            raise ValueError(
                f"Pair metadata date {pair.date} does not match filename date {day}"
            )
        return pair

    return _compute_from_store(
        tuple(path_by_date),
        load,
        climatology,
        lag_days=lag_days,
        min_count=min_count,
        start_date=start_date,
        end_date=end_date,
    )


def save_independent_validation_npz(
    result: IndependentValidation, path: Path | str
) -> None:
    """Write sparse IV fields using MATLAB variable names where applicable."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        idx0=np.asarray(result.idx, dtype=np.int32),
        N_sm=np.asarray(result.count, dtype=np.int32),
        R2_ivd_mod=np.asarray(result.r2_ivd_model, dtype=np.float64),
        R2_ivd_obs=np.asarray(result.r2_ivd_obs, dtype=np.float64),
        R2_ivs_mod=np.asarray(result.r2_ivs_model, dtype=np.float64),
        R2_ivs_obs=np.asarray(result.r2_ivs_obs, dtype=np.float64),
        R_mod_obs=np.asarray(result.r_model_obs, dtype=np.float64),
        Nmin=np.asarray(result.min_count, dtype=np.int32),
        Nlag=np.asarray(result.lag_days, dtype=np.int16),
        sensor=np.asarray(result.sensor),
        run=np.asarray(result.run),
        obs_units=np.asarray(result.obs_units),
        model_units=np.asarray(result.model_units),
        start_date=np.asarray(result.start_date.strftime("%Y%m%d")),
        end_date=np.asarray(result.end_date.strftime("%Y%m%d")),
        source_days=np.asarray(result.source_days, dtype=np.int32),
        candidate_days=np.asarray(result.candidate_days, dtype=np.int32),
        paired_days=np.asarray(result.paired_days, dtype=np.int32),
        contributing_days=np.asarray(result.contributing_days, dtype=np.int32),
        climatology_start_date=np.asarray(
            result.climatology_start_date.strftime("%Y%m%d")
        ),
        climatology_end_date=np.asarray(
            result.climatology_end_date.strftime("%Y%m%d")
        ),
        climatology_window_days=np.asarray(
            result.climatology_window_days, dtype=np.int16
        ),
        climatology_min_count=np.asarray(
            result.climatology_min_count, dtype=np.int32
        ),
        grid_size=np.asarray(result.grid_size, dtype=np.int32),
        undefined_correlation=np.asarray(
            MATLAB_UNDEFINED_CORRELATION, dtype=np.float64
        ),
    )


def load_independent_validation_npz(path: Path | str) -> IndependentValidation:
    """Load a sparse step-4 result."""

    with np.load(Path(path)) as z:
        return IndependentValidation(
            idx=np.asarray(z["idx0"], dtype=np.int64),
            count=np.asarray(z["N_sm"], dtype=np.int64),
            r2_ivd_model=np.asarray(z["R2_ivd_mod"], dtype=np.float64),
            r2_ivd_obs=np.asarray(z["R2_ivd_obs"], dtype=np.float64),
            r2_ivs_model=np.asarray(z["R2_ivs_mod"], dtype=np.float64),
            r2_ivs_obs=np.asarray(z["R2_ivs_obs"], dtype=np.float64),
            r_model_obs=np.asarray(z["R_mod_obs"], dtype=np.float64),
            sensor=str(z["sensor"]),
            run=str(z["run"]),
            obs_units=str(z["obs_units"]),
            model_units=str(z["model_units"]),
            start_date=_parse_yyyymmdd(str(z["start_date"])),
            end_date=_parse_yyyymmdd(str(z["end_date"])),
            lag_days=int(z["Nlag"]),
            min_count=int(z["Nmin"]),
            source_days=int(z["source_days"]),
            candidate_days=int(z["candidate_days"]),
            paired_days=int(z["paired_days"]),
            contributing_days=int(z["contributing_days"]),
            climatology_start_date=_parse_yyyymmdd(
                str(z["climatology_start_date"])
            ),
            climatology_end_date=_parse_yyyymmdd(
                str(z["climatology_end_date"])
            ),
            climatology_window_days=int(z["climatology_window_days"]),
            climatology_min_count=int(z["climatology_min_count"]),
            grid_size=int(z["grid_size"]),
        )


def independent_validation_npz_is_valid(path: Path | str) -> bool:
    """Return whether a file contains shape-consistent sparse IV fields."""

    path = Path(path)
    if not path.exists():
        return False
    fields = (
        "idx0",
        "N_sm",
        "R2_ivd_mod",
        "R2_ivd_obs",
        "R2_ivs_mod",
        "R2_ivs_obs",
        "R_mod_obs",
        "Nmin",
        "Nlag",
    )
    try:
        with np.load(path) as z:
            if not all(field in z for field in fields):
                return False
            arrays = [np.asarray(z[field]) for field in fields[:7]]
    except Exception:
        return False
    return bool(
        arrays[0].ndim == 1
        and all(array.ndim == 1 for array in arrays[1:])
        and all(array.size == arrays[0].size for array in arrays[1:])
    )


def _compute_from_store(
    available_dates: Sequence[date],
    load_pair: Callable[[date], DailyPair],
    climatology: PentadClimatology,
    *,
    lag_days: int,
    min_count: int,
    start_date: date | None,
    end_date: date | None,
) -> IndependentValidation:
    _validate_options(lag_days, min_count)
    _validate_climatology(climatology)
    dates = set(available_dates)
    if not dates:
        raise ValueError("At least one daily pair is required")

    analysis_start = min(dates) if start_date is None else start_date
    analysis_end = max(dates) if end_date is None else end_date
    if analysis_end < analysis_start:
        raise ValueError(
            f"end date {analysis_end} is before start date {analysis_start}"
        )

    grid_size = climatology.grid_size
    sums = [np.zeros(grid_size, dtype=np.float64) for _ in range(9)]
    (
        model_sum,
        model2_sum,
        obs_sum,
        obs2_sum,
        model_obs_sum,
        model_model_lag_sum,
        obs_obs_lag_sum,
        model_obs_lag_sum,
        model_lag_obs_sum,
    ) = sums
    count = np.zeros(grid_size, dtype=np.int32)

    candidate_days = max((analysis_end - analysis_start).days - lag_days + 1, 0)
    paired_days = 0
    contributing_days = 0
    pair_cache: dict[date, DailyPair] = {}

    def cached_load(value: date) -> DailyPair:
        if value not in pair_cache:
            pair_cache[value] = load_pair(value)
        return pair_cache[value]

    day = analysis_start + timedelta(days=lag_days)
    while day <= analysis_end:
        previous_day = day - timedelta(days=lag_days)
        if day not in dates or previous_day not in dates:
            day += timedelta(days=1)
            continue

        current = cached_load(day)
        previous = cached_load(previous_day)
        _validate_pair(current, climatology)
        _validate_pair(previous, climatology)
        paired_days += 1

        common, current_pos, previous_pos = np.intersect1d(
            np.asarray(current.idx, dtype=np.int64),
            np.asarray(previous.idx, dtype=np.int64),
            assume_unique=True,
            return_indices=True,
        )
        if common.size == 0:
            day += timedelta(days=1)
            continue

        current_col = pentad_index(day)
        previous_col = pentad_index(previous_day)
        model_now = np.asarray(current.model, dtype=np.float64)[current_pos]
        obs_now = np.asarray(current.obs, dtype=np.float64)[current_pos]
        model_previous = np.asarray(previous.model, dtype=np.float64)[previous_pos]
        obs_previous = np.asarray(previous.obs, dtype=np.float64)[previous_pos]
        model_now -= _climatology_values(climatology, common, current_col, "model")
        obs_now -= _climatology_values(climatology, common, current_col, "obs")
        model_previous -= _climatology_values(
            climatology, common, previous_col, "model"
        )
        obs_previous -= _climatology_values(
            climatology, common, previous_col, "obs"
        )

        # MATLAB selects only on the current model anomaly. Any NaN in the
        # other three terms then propagates through that cell's running sums.
        valid = ~np.isnan(model_now)
        if not np.any(valid):
            day += timedelta(days=1)
            continue
        contributing_days += 1
        idx = common[valid]
        model_now = model_now[valid]
        obs_now = obs_now[valid]
        model_previous = model_previous[valid]
        obs_previous = obs_previous[valid]

        model_sum[idx] += model_now
        model2_sum[idx] += model_now**2
        obs_sum[idx] += obs_now
        obs2_sum[idx] += obs_now**2
        model_obs_sum[idx] += model_now * obs_now
        model_model_lag_sum[idx] += model_now * model_previous
        obs_obs_lag_sum[idx] += obs_now * obs_previous
        model_obs_lag_sum[idx] += model_now * obs_previous
        model_lag_obs_sum[idx] += model_previous * obs_now
        count[idx] += 1
        for cached_day in tuple(pair_cache):
            if cached_day <= previous_day:
                del pair_cache[cached_day]
        day += timedelta(days=1)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        nn = np.where(count >= min_count, count, np.nan).astype(np.float64)
        model_mean = model_sum / nn
        obs_mean = obs_sum / nn
        model_mean[np.isinf(model_mean)] = np.nan
        obs_mean[np.isinf(obs_mean)] = np.nan

        c_model_model = model2_sum / nn - model_mean**2
        c_obs_obs = obs2_sum / nn - obs_mean**2
        c_model_model[c_model_model < 0.0] = np.nan
        c_obs_obs[c_obs_obs < 0.0] = np.nan

        c_model_obs = model_obs_sum / nn - model_mean * obs_mean
        r_model_obs = c_model_obs / np.sqrt(c_model_model * c_obs_obs)
        # A zero-variance model or obs sample (denominator exactly 0) gives a
        # meaningless +/-inf "correlation" when the covariance isn't also
        # exactly 0 -- not caught by the isnan check below, since inf is not
        # NaN. Treat it the same as the undefined 0/0 case.
        r_model_obs[~np.isfinite(r_model_obs)] = MATLAB_UNDEFINED_CORRELATION
        c_model_obs[c_model_obs < 0.0] = np.nan

        c_model_model_lag = model_model_lag_sum / nn - model_mean**2
        c_obs_obs_lag = obs_obs_lag_sum / nn - obs_mean**2
        c_model_model_lag[c_model_model_lag < 0.0] = np.nan
        c_obs_obs_lag[c_obs_obs_lag < 0.0] = np.nan
        c_model_obs_lag = model_obs_lag_sum / nn - model_mean * obs_mean
        c_model_lag_obs = model_lag_obs_sum / nn - model_mean * obs_mean

        scale_ivd = np.sqrt(c_model_model_lag / c_obs_obs_lag)
        scale_ivs_obs = c_model_obs_lag / c_obs_obs_lag

        r2_ivd_model = c_model_obs * scale_ivd / c_model_model
        r2_ivd_obs = c_model_obs / c_obs_obs / scale_ivd
        r2_ivs_model = c_model_obs * scale_ivs_obs / c_model_model
        r2_ivs_obs = c_model_obs / c_obs_obs / scale_ivs_obs

    for field in (r2_ivd_model, r2_ivd_obs, r2_ivs_model, r2_ivs_obs):
        field[field < 0.0] = np.nan
        field[field > 1.0] = 1.0

    occupied = count > 0
    idx = np.flatnonzero(occupied)
    return IndependentValidation(
        idx=idx,
        count=count[occupied],
        r2_ivd_model=r2_ivd_model[occupied],
        r2_ivd_obs=r2_ivd_obs[occupied],
        r2_ivs_model=r2_ivs_model[occupied],
        r2_ivs_obs=r2_ivs_obs[occupied],
        r_model_obs=r_model_obs[occupied],
        sensor=climatology.sensor,
        run=climatology.run,
        obs_units=climatology.obs_units,
        model_units=climatology.model_units,
        start_date=analysis_start,
        end_date=analysis_end,
        lag_days=lag_days,
        min_count=min_count,
        source_days=sum(analysis_start <= value <= analysis_end for value in dates),
        candidate_days=candidate_days,
        paired_days=paired_days,
        contributing_days=contributing_days,
        climatology_start_date=climatology.start_date,
        climatology_end_date=climatology.end_date,
        climatology_window_days=climatology.window_days,
        climatology_min_count=climatology.min_count,
        grid_size=grid_size,
    )


def _climatology_values(
    climatology: PentadClimatology,
    idx: np.ndarray,
    column: int,
    field: str,
) -> np.ndarray:
    rows = np.searchsorted(climatology.idx, idx)
    found = rows < climatology.idx.size
    found[found] &= climatology.idx[rows[found]] == idx[found]
    values = np.full(idx.size, np.nan, dtype=np.float64)
    source = climatology.model if field == "model" else climatology.obs
    values[found] = source[rows[found], column]
    return values


def _validate_options(lag_days: int, min_count: int) -> None:
    if lag_days <= 0:
        raise ValueError("lag_days must be positive")
    if min_count < 0:
        raise ValueError("min_count must be non-negative")


def _validate_climatology(climatology: PentadClimatology) -> None:
    idx = np.asarray(climatology.idx)
    expected = (idx.size, PENTAD_CENTER_DOYS.size)
    if idx.ndim != 1 or np.any(np.diff(idx) <= 0):
        raise ValueError("Climatology indices must be one-dimensional and sorted unique")
    if climatology.obs.shape != expected or climatology.model.shape != expected:
        raise ValueError(f"Climatology arrays must have shape {expected}")
    if not np.array_equal(climatology.center_doy, PENTAD_CENTER_DOYS):
        raise ValueError("Climatology pentad centers do not match MATLAB days 3:5:365")
    if np.any(idx < 0) or np.any(idx >= climatology.grid_size):
        raise ValueError("Climatology index is outside its grid")


def _validate_pair(pair: DailyPair, climatology: PentadClimatology) -> None:
    idx = np.asarray(pair.idx)
    obs = np.asarray(pair.obs)
    model = np.asarray(pair.model)
    if idx.ndim != 1 or obs.ndim != 1 or model.ndim != 1:
        raise ValueError(f"Daily pair arrays must be one-dimensional on {pair.date}")
    if not (idx.size == obs.size == model.size):
        raise ValueError(f"Daily pair arrays have unequal lengths on {pair.date}")
    if np.any(idx < 0) or np.any(idx >= climatology.grid_size):
        raise ValueError(f"Daily pair index is outside the grid on {pair.date}")
    if np.unique(idx).size != idx.size:
        raise ValueError(f"Daily pair contains duplicate grid indices on {pair.date}")
    metadata = (pair.sensor, pair.run, pair.obs_units, pair.model_units)
    expected = (
        climatology.sensor,
        climatology.run,
        climatology.obs_units,
        climatology.model_units,
    )
    if metadata != expected:
        raise ValueError(
            "Daily pair metadata does not match climatology; "
            f"expected {expected}, got {metadata} on {pair.date}"
        )


def _parse_yyyymmdd(value: str) -> date:
    value = value.strip()
    if len(value) != 8 or not value.isdigit():
        raise ValueError(f"Expected YYYYMMDD date, got {value!r}")
    return date(int(value[:4]), int(value[4:6]), int(value[6:8]))
