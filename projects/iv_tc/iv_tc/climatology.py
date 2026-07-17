"""Sparse pentad climatologies for IV/TC daily pairs."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np

from .generation import canonical_sensor, load_daily_pair_npz, pair_output_path
from .pairs import DailyPair


M36_GRID_SIZE = 964 * 406
PENTAD_CENTER_DOYS = np.arange(3, 366, 5, dtype=np.int16)


@dataclass(frozen=True)
class PentadClimatology:
    """Matched observation/model climatologies on sparse M36 cells."""

    idx: np.ndarray
    obs: np.ndarray
    model: np.ndarray
    count: np.ndarray
    center_doy: np.ndarray
    sensor: str
    run: str
    obs_units: str
    model_units: str
    start_date: date
    end_date: date
    window_days: int
    min_count: int
    source_days: int
    grid_size: int = M36_GRID_SIZE


def normal_year_doy(day: date) -> int:
    """Map a date to a 365-day calendar, assigning Feb 29 to Feb 28."""

    month = day.month
    month_day = 28 if month == 2 and day.day == 29 else day.day
    return date(2017, month, month_day).timetuple().tm_yday


def pentad_index(day: date) -> int:
    """Return the zero-based pentad used by the MATLAB downstream scripts."""

    return (normal_year_doy(day) - 1) // 5


def matlab_default_min_count(start_date: date, end_date: date) -> int:
    """Return MATLAB's ``4 * year_span`` threshold for an inclusive range."""

    if end_date < start_date:
        raise ValueError(f"end date {end_date} is before start date {start_date}")
    end_exclusive = end_date + timedelta(days=1)
    return 4 * (end_exclusive.year - start_date.year)


def climatology_output_path(
    output_root: Path | str,
    sensor: str,
    run_name: str,
    start_date: date,
    end_date: date,
    window_days: int = 31,
) -> Path:
    """Return the standard step-3 climatology output path."""

    return (
        Path(output_root)
        / "step3_climatology"
        / canonical_sensor(sensor)
        / run_name
        / f"{start_date:%Y%m%d}_{end_date:%Y%m%d}_w{window_days}.npz"
    )


def daily_pair_paths(
    pair_root: Path | str,
    sensor: str,
    run_name: str,
    start_date: date,
    end_date: date,
) -> tuple[list[Path], list[date]]:
    """Collect existing pair paths and dates missing from an inclusive range."""

    if end_date < start_date:
        raise ValueError(f"end date {end_date} is before start date {start_date}")

    paths: list[Path] = []
    missing: list[date] = []
    day = start_date
    while day <= end_date:
        path = pair_output_path(pair_root, sensor, run_name, day)
        if path.exists():
            paths.append(path)
        else:
            missing.append(day)
        day += timedelta(days=1)
    return paths, missing


def compute_pentad_climatology(
    pairs: Iterable[DailyPair],
    *,
    min_count: int | None = None,
    window_days: int = 31,
    grid_size: int = M36_GRID_SIZE,
) -> PentadClimatology:
    """Compute a MATLAB-step3-compatible climatology from daily pairs."""

    pair_list = list(pairs)
    return _compute_from_pair_source(
        lambda: iter(pair_list),
        min_count=min_count,
        window_days=window_days,
        grid_size=grid_size,
    )


def compute_pentad_climatology_from_paths(
    paths: Sequence[Path | str],
    *,
    min_count: int | None = None,
    window_days: int = 31,
    grid_size: int = M36_GRID_SIZE,
) -> PentadClimatology:
    """Compute a climatology while streaming daily pair files twice."""

    pair_paths = tuple(Path(path) for path in paths)
    return _compute_from_pair_source(
        lambda: (load_daily_pair_npz(path) for path in pair_paths),
        min_count=min_count,
        window_days=window_days,
        grid_size=grid_size,
    )


def save_pentad_climatology_npz(climatology: PentadClimatology, path: Path | str) -> None:
    """Write a sparse climatology with MATLAB-compatible array names."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        idx0=np.asarray(climatology.idx, dtype=np.int32),
        obs_sm_clim=np.asarray(climatology.obs, dtype=np.float32),
        mod_sm_clim=np.asarray(climatology.model, dtype=np.float32),
        N_sm_clim=np.asarray(climatology.count, dtype=np.int32),
        pentad_center_doy=np.asarray(climatology.center_doy, dtype=np.int16),
        sensor=np.asarray(climatology.sensor),
        run=np.asarray(climatology.run),
        obs_units=np.asarray(climatology.obs_units),
        model_units=np.asarray(climatology.model_units),
        start_date=np.asarray(climatology.start_date.strftime("%Y%m%d")),
        end_date=np.asarray(climatology.end_date.strftime("%Y%m%d")),
        window_days=np.asarray(climatology.window_days, dtype=np.int16),
        Nday_min=np.asarray(climatology.min_count, dtype=np.int32),
        source_days=np.asarray(climatology.source_days, dtype=np.int32),
        grid_size=np.asarray(climatology.grid_size, dtype=np.int32),
    )


def load_pentad_climatology_npz(path: Path | str) -> PentadClimatology:
    """Load a sparse step-3 climatology file."""

    with np.load(Path(path)) as z:
        return PentadClimatology(
            idx=np.asarray(z["idx0"], dtype=np.int64),
            obs=np.asarray(z["obs_sm_clim"], dtype=np.float64),
            model=np.asarray(z["mod_sm_clim"], dtype=np.float64),
            count=np.asarray(z["N_sm_clim"], dtype=np.int64),
            center_doy=np.asarray(z["pentad_center_doy"], dtype=np.int64),
            sensor=str(z["sensor"]),
            run=str(z["run"]),
            obs_units=str(z["obs_units"]),
            model_units=str(z["model_units"]),
            start_date=_parse_yyyymmdd(str(z["start_date"])),
            end_date=_parse_yyyymmdd(str(z["end_date"])),
            window_days=int(z["window_days"]),
            min_count=int(z["Nday_min"]),
            source_days=int(z["source_days"]),
            grid_size=int(z["grid_size"]),
        )


def climatology_npz_is_valid(path: Path | str) -> bool:
    """Return True when a file contains a shape-consistent climatology."""

    path = Path(path)
    if not path.exists():
        return False
    try:
        with np.load(path) as z:
            required = ("idx0", "obs_sm_clim", "mod_sm_clim", "N_sm_clim", "pentad_center_doy")
            if not all(key in z for key in required):
                return False
            idx = np.asarray(z["idx0"])
            obs = np.asarray(z["obs_sm_clim"])
            model = np.asarray(z["mod_sm_clim"])
            count = np.asarray(z["N_sm_clim"])
            centers = np.asarray(z["pentad_center_doy"])
    except Exception:
        return False
    return bool(
        idx.ndim == 1
        and centers.ndim == 1
        and obs.shape == model.shape == count.shape == (idx.size, centers.size)
        and centers.size == 73
    )


def _compute_from_pair_source(
    pair_source: Callable[[], Iterator[DailyPair]],
    *,
    min_count: int | None,
    window_days: int,
    grid_size: int,
) -> PentadClimatology:
    _validate_options(min_count, window_days, grid_size)

    occupied = np.zeros(grid_size, dtype=bool)
    metadata: tuple[str, str, str, str] | None = None
    seen_dates: set[date] = set()
    start_date: date | None = None
    end_date: date | None = None

    for pair in pair_source():
        _validate_pair(pair, grid_size)
        current_metadata = (pair.sensor, pair.run, pair.obs_units, pair.model_units)
        if metadata is None:
            metadata = current_metadata
        elif current_metadata != metadata:
            raise ValueError(
                "All daily pairs must have identical sensor, run, and units; "
                f"expected {metadata}, got {current_metadata} on {pair.date}"
            )
        if pair.date in seen_dates:
            raise ValueError(f"Duplicate daily pair date: {pair.date}")
        seen_dates.add(pair.date)
        start_date = pair.date if start_date is None else min(start_date, pair.date)
        end_date = pair.date if end_date is None else max(end_date, pair.date)
        occupied[np.asarray(pair.idx, dtype=np.int64)] = True

    if metadata is None or start_date is None or end_date is None:
        raise ValueError("At least one daily pair is required")

    idx = np.flatnonzero(occupied)
    shape = (idx.size, PENTAD_CENTER_DOYS.size)
    obs_sum = np.zeros(shape, dtype=np.float64)
    model_sum = np.zeros(shape, dtype=np.float64)
    count = np.zeros(shape, dtype=np.int32)
    half_window = window_days // 2

    for pair in pair_source():
        pair_idx = np.asarray(pair.idx, dtype=np.int64)
        obs = np.asarray(pair.obs, dtype=np.float64)
        model = np.asarray(pair.model, dtype=np.float64)
        good = np.isfinite(obs) & np.isfinite(model)
        if not np.any(good):
            continue
        rows = np.searchsorted(idx, pair_idx[good])
        obs = obs[good]
        model = model[good]
        source_doy = normal_year_doy(pair.date)
        distance = np.abs(PENTAD_CENTER_DOYS.astype(np.int32) - source_doy)
        distance = np.minimum(distance, 365 - distance)
        for column in np.flatnonzero(distance <= half_window):
            obs_sum[rows, column] += obs
            model_sum[rows, column] += model
            count[rows, column] += 1

    threshold = matlab_default_min_count(start_date, end_date) if min_count is None else min_count
    valid = (count >= threshold) & (count > 0)
    obs_clim = np.full(shape, np.nan, dtype=np.float32)
    model_clim = np.full(shape, np.nan, dtype=np.float32)
    obs_clim[valid] = (obs_sum[valid] / count[valid]).astype(np.float32)
    model_clim[valid] = (model_sum[valid] / count[valid]).astype(np.float32)

    sensor, run, obs_units, model_units = metadata
    return PentadClimatology(
        idx=idx,
        obs=obs_clim,
        model=model_clim,
        count=count,
        center_doy=PENTAD_CENTER_DOYS.copy(),
        sensor=sensor,
        run=run,
        obs_units=obs_units,
        model_units=model_units,
        start_date=start_date,
        end_date=end_date,
        window_days=window_days,
        min_count=threshold,
        source_days=len(seen_dates),
        grid_size=grid_size,
    )


def _validate_options(min_count: int | None, window_days: int, grid_size: int) -> None:
    if window_days <= 0 or window_days > 365 or window_days % 2 == 0:
        raise ValueError("window_days must be an odd integer from 1 through 365")
    if min_count is not None and min_count < 0:
        raise ValueError("min_count must be non-negative")
    if grid_size <= 0:
        raise ValueError("grid_size must be positive")


def _validate_pair(pair: DailyPair, grid_size: int) -> None:
    idx = np.asarray(pair.idx)
    obs = np.asarray(pair.obs)
    model = np.asarray(pair.model)
    if idx.ndim != 1 or obs.ndim != 1 or model.ndim != 1:
        raise ValueError(f"Daily pair arrays must be one-dimensional on {pair.date}")
    if not (idx.size == obs.size == model.size):
        raise ValueError(f"Daily pair arrays have unequal lengths on {pair.date}")
    if np.any(idx < 0) or np.any(idx >= grid_size):
        raise ValueError(f"Daily pair index is outside grid_size={grid_size} on {pair.date}")
    if np.unique(idx).size != idx.size:
        raise ValueError(f"Daily pair contains duplicate grid indices on {pair.date}")


def _parse_yyyymmdd(value: str) -> date:
    value = value.strip()
    return date(int(value[:4]), int(value[4:6]), int(value[6:8]))
