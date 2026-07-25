"""ISMN archive readers for in-situ soil-moisture validation.

The `ismn` Python package is not installed on Discover, so this module reads
the local ISMN "header_values" (.stm) archive directly. The per-sensor
inventory comes from the archive's own `python_metadata/ISMN_data.csv`, which
carries each sensor's network, station, coordinates, depth range, time range,
and the archive-relative path to its .stm file -- enough to work across every
network without hardcoding a network list.

Surface and root-zone construction follow the M21C_ls ISMN workflow
(`projects/M21C_ls/ISMN_methods_readme.md`), except that root zone uses the
generic profile-weighted rule rather than per-network strict layer rules, so
it generalizes to networks that workflow never covered.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import csv

import numpy as np
import pandas as pd


SM_VARIABLE = "soil_moisture"
DEFAULT_REQUIRED_FLAG = "G"

# The metadata CSV has a two-row header; every column is a (name, sub) pair.
_REQUIRED_COLUMNS = (
    ("network", "val"),
    ("station", "val"),
    ("latitude", "val"),
    ("longitude", "val"),
    ("variable", "val"),
    ("variable", "depth_from"),
    ("variable", "depth_to"),
    ("timerange_from", "val"),
    ("timerange_to", "val"),
    ("file_path", "val"),
)


def _parse_float(value) -> float:
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "nat"}:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def _parse_datetime(value):
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "nat"}:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def station_key(network: str, station: str) -> str:
    """Stable station identifier; station names are not unique across networks."""

    return f"{str(network).strip()}::{str(station).strip()}"


def load_sm_sensor_inventory(metadata_csv, archive_root=None) -> pd.DataFrame:
    """Return one row per soil-moisture sensor across every network."""

    metadata_csv = Path(metadata_csv)
    if not metadata_csv.exists():
        raise FileNotFoundError(f"ISMN metadata CSV not found: {metadata_csv}")

    archive_root = Path(archive_root) if archive_root is not None else metadata_csv.parent.parent

    rows = []
    with metadata_csv.open(newline="") as handle:
        reader = csv.reader(handle)
        header_top = next(reader)
        header_sub = next(reader)
        col_idx = {col: i for i, col in enumerate(zip(header_top, header_sub))}

        missing = [col for col in _REQUIRED_COLUMNS if col not in col_idx]
        if missing:
            raise KeyError(f"Missing required metadata columns: {missing}")

        for record in reader:
            if record[col_idx[("variable", "val")]].strip() != SM_VARIABLE:
                continue

            network = record[col_idx[("network", "val")]].strip()
            station = record[col_idx[("station", "val")]].strip()
            depth_from = _parse_float(record[col_idx[("variable", "depth_from")]])
            depth_to = _parse_float(record[col_idx[("variable", "depth_to")]])
            rel_path = record[col_idx[("file_path", "val")]].strip()

            rows.append(
                {
                    "station_key": station_key(network, station),
                    "network": network,
                    "station": station,
                    "lat": _parse_float(record[col_idx[("latitude", "val")]]),
                    "lon": _parse_float(record[col_idx[("longitude", "val")]]),
                    "depth_from": depth_from,
                    "depth_to": depth_to,
                    "depth_center": 0.5 * (depth_from + depth_to),
                    "timerange_from": _parse_datetime(record[col_idx[("timerange_from", "val")]]),
                    "timerange_to": _parse_datetime(record[col_idx[("timerange_to", "val")]]),
                    "path": archive_root / rel_path,
                }
            )

    if not rows:
        raise RuntimeError(f"No {SM_VARIABLE} sensors found in {metadata_csv}")

    return pd.DataFrame(rows)


def filter_inventory_by_window(inventory: pd.DataFrame, start, end) -> pd.DataFrame:
    """Keep sensors whose reported time range overlaps [start, end]."""

    start = pd.Timestamp(start)
    end = pd.Timestamp(end)

    t_from = pd.to_datetime(inventory["timerange_from"], errors="coerce")
    t_to = pd.to_datetime(inventory["timerange_to"], errors="coerce")

    # Missing bounds are kept: the .stm file is the authority, not the metadata.
    overlaps = (t_from.isna() | (t_from <= end)) & (t_to.isna() | (t_to >= start))
    return inventory.loc[overlaps].copy().reset_index(drop=True)


def read_stm_series(path, required_flag: str = DEFAULT_REQUIRED_FLAG) -> pd.Series:
    """Read one .stm sensor file into a QC-filtered, duplicate-averaged series.

    Data lines are ``YYYY/MM/DD HH:MM value ismn_flag provider_flag``; only
    records carrying the requested ISMN quality flag are retained.
    """

    frame = pd.read_csv(
        path,
        sep=r"\s+",
        skiprows=1,
        header=None,
        names=["date", "time", "value", "flag", "provider_flag"],
        dtype=str,
        on_bad_lines="skip",
        engine="c",
    )

    if frame.empty:
        return pd.Series(index=pd.DatetimeIndex([], dtype="datetime64[ns]"), dtype=float)

    stamps = pd.to_datetime(
        frame["date"].str.strip() + " " + frame["time"].str.strip(),
        format="%Y/%m/%d %H:%M",
        errors="coerce",
    )
    values = pd.to_numeric(frame["value"], errors="coerce")
    flags = frame["flag"].astype(str).str.strip().str.upper()

    keep = flags.eq(str(required_flag).strip().upper()) & stamps.notna() & values.notna()
    if not keep.any():
        return pd.Series(index=pd.DatetimeIndex([], dtype="datetime64[ns]"), dtype=float)

    series = pd.Series(values[keep].to_numpy(dtype=float), index=pd.DatetimeIndex(stamps[keep]))
    return series.groupby(series.index).mean().sort_index()


def depth_layer_weights(depth_centers, z_top: float = 0.0, z_bottom: float = 1.0) -> dict:
    """Layer thickness each sensor depth represents within [z_top, z_bottom]."""

    depths = np.array(depth_centers, dtype=float)
    depths = np.sort(np.unique(depths[np.isfinite(depths)]))
    if depths.size == 0:
        return {}

    weights = {}
    for i, depth in enumerate(depths):
        upper = z_top if i == 0 else 0.5 * (depths[i - 1] + depth)
        lower = z_bottom if i == depths.size - 1 else 0.5 * (depth + depths[i + 1])
        thickness = max(0.0, float(min(lower, z_bottom) - max(upper, z_top)))
        if thickness > 0:
            weights[float(depth)] = thickness
    return weights


def weighted_rootzone(
    depth_frame: pd.DataFrame,
    z_top: float = 0.0,
    z_bottom: float = 1.0,
    min_sensors: int = 3,
    shallow_max: float = 0.20,
    deep_min: float = 0.50,
) -> pd.Series:
    """Profile-weighted root-zone series over [z_top, z_bottom].

    A timestamp yields a value only when it has at least ``min_sensors`` finite
    layers spanning both a shallow (<= ``shallow_max``) and a deep
    (>= ``deep_min``) sensor, so partial profiles never masquerade as root zone.
    """

    if depth_frame is None or depth_frame.shape[1] == 0:
        return pd.Series(index=pd.DatetimeIndex([], dtype="datetime64[ns]"), dtype=float)

    weight_map = depth_layer_weights([float(c) for c in depth_frame.columns], z_top, z_bottom)
    use_depths = [float(c) for c in depth_frame.columns if float(c) in weight_map]
    if not use_depths:
        return pd.Series(index=depth_frame.index, dtype=float)

    values = depth_frame[use_depths].to_numpy(dtype=float)
    weights = np.array([weight_map[d] for d in use_depths], dtype=float)
    depths = np.array(use_depths, dtype=float)

    finite = np.isfinite(values)
    weight_sum = np.where(finite, weights[None, :], 0.0).sum(axis=1)
    weighted_sum = np.nansum(np.where(finite, values * weights[None, :], np.nan), axis=1)

    good = (
        (finite.sum(axis=1) >= int(min_sensors))
        & np.any(finite & (depths[None, :] <= float(shallow_max)), axis=1)
        & np.any(finite & (depths[None, :] >= float(deep_min)), axis=1)
        & (weight_sum > 0)
    )

    rootzone = np.full(values.shape[0], np.nan, dtype=float)
    rootzone[good] = weighted_sum[good] / weight_sum[good]

    series = pd.Series(rootzone, index=depth_frame.index)
    return series.groupby(series.index).mean().sort_index()


def build_station_series(
    sensor_rows: pd.DataFrame,
    required_flag: str = DEFAULT_REQUIRED_FLAG,
    target_surface_depth_m: float = 0.05,
    surface_max_depth_m: float = 0.10,
    daily_shift_hours: int = 12,
    rz_top_m: float = 0.0,
    rz_bottom_m: float = 1.0,
    rz_min_sensors: int = 3,
    rz_shallow_max_m: float = 0.20,
    rz_deep_min_m: float = 0.50,
    start=None,
    end=None,
) -> dict:
    """Build daily surface and root-zone series for one ISMN station.

    Surface uses the available sensor depth nearest ``target_surface_depth_m``
    among sensors no deeper than ``surface_max_depth_m``; a station with no
    such sensor gets an empty surface series rather than a mislabeled deep one.
    """

    # Sub-daily records are trimmed to the analysis window as they are read;
    # ISMN histories run decades longer than any one experiment.
    t0 = pd.Timestamp(start) if start is not None else None
    t1 = pd.Timestamp(end) + pd.Timedelta(days=1) if end is not None else None

    depth_series: dict[float, list[pd.Series]] = {}
    for row in sensor_rows.itertuples():
        depth = float(row.depth_center)
        if not np.isfinite(depth):
            continue
        try:
            series = read_stm_series(row.path, required_flag=required_flag)
        except (OSError, ValueError, pd.errors.ParserError):
            continue
        if t0 is not None:
            series = series[series.index >= t0]
        if t1 is not None:
            series = series[series.index <= t1]
        if series.empty:
            continue
        depth_series.setdefault(depth, []).append(series)

    if not depth_series:
        return {"ok": False, "reason": "no usable QC-filtered sensors"}

    # Collapse replicate sensors at the same nominal depth.
    merged: dict[float, pd.Series] = {}
    for depth, series_list in depth_series.items():
        if len(series_list) == 1:
            combined = series_list[0]
        else:
            combined = pd.concat(series_list, axis=1).mean(axis=1)
        merged[depth] = combined.groupby(combined.index).mean().sort_index()

    depths = sorted(merged)

    surface_candidates = [d for d in depths if d <= float(surface_max_depth_m)]
    if surface_candidates:
        surface_depth = min(surface_candidates, key=lambda d: abs(d - float(target_surface_depth_m)))
        surface_native = merged[surface_depth]
    else:
        surface_depth = np.nan
        surface_native = pd.Series(index=pd.DatetimeIndex([], dtype="datetime64[ns]"), dtype=float)

    depth_frame = pd.concat([merged[d].rename(d) for d in depths], axis=1, join="outer").sort_index()
    rootzone_native = weighted_rootzone(
        depth_frame,
        z_top=rz_top_m,
        z_bottom=rz_bottom_m,
        min_sensors=rz_min_sensors,
        shallow_max=rz_shallow_max_m,
        deep_min=rz_deep_min_m,
    )

    def to_daily(series: pd.Series) -> pd.Series:
        if series.empty:
            return series
        daily = series.resample("D").mean()
        if daily_shift_hours:
            daily.index = daily.index + pd.Timedelta(hours=int(daily_shift_hours))
        return daily.dropna()

    return {
        "ok": True,
        "surface_daily": to_daily(surface_native),
        "rz_daily": to_daily(rootzone_native),
        "surface_depth_m": float(surface_depth) if np.isfinite(surface_depth) else np.nan,
        "depths_m": depths,
        "n_depths": len(depths),
    }
