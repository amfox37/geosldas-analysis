"""Readers for sparse observations and GEOSldas model fixtures."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import struct

import numpy as np

from .pairs import DailyPair, SparseObservation


M36_NX = 964
SMOSIC_DATE_RE = re.compile(r"smos_ic_sm_m36_(\d{8})\.nc$")


@dataclass(frozen=True)
class ModelTileValues:
    """One GEOSldas daily model variable on tile space."""

    values: np.ndarray
    units: str


@dataclass(frozen=True)
class RepresentativeTiles:
    """Representative tile chosen for each M36 EASE grid cell."""

    m36_linear: np.ndarray
    tile_index: np.ndarray


def read_smosic_sparse(path: Path | str) -> SparseObservation:
    """Read one preprocessed SMOS-IC sparse daily NetCDF file."""

    path = Path(path)
    Dataset = _dataset_cls()

    with Dataset(path, "r") as ds:
        if "idx_EASEv2_lonxlat" not in ds.variables or "sm_obs" not in ds.variables:
            raise KeyError(f"{path} is missing SMOS-IC sparse variables")

        idx = _filled(ds.variables["idx_EASEv2_lonxlat"][:]).astype(np.int64)
        values = _filled(ds.variables["sm_obs"][:]).astype(np.float64)
        units = getattr(ds.variables["sm_obs"], "units", "")
        date_value = getattr(ds, "date", None)

        qc_summary: dict[str, int | float | str] = {"n_points": int(idx.size)}
        if "coverage_frac" in ds.variables:
            coverage = _filled(ds.variables["coverage_frac"][:]).astype(np.float64)
            finite_cov = coverage[np.isfinite(coverage)]
            if finite_cov.size:
                qc_summary["coverage_min"] = float(finite_cov.min())
                qc_summary["coverage_mean"] = float(finite_cov.mean())
                qc_summary["coverage_max"] = float(finite_cov.max())

        for attr in ("qc_scene_flag_max", "qc_tb_rmse_max", "qc_sm_min", "qc_sm_max", "min_coverage_frac"):
            if hasattr(ds, attr):
                value = getattr(ds, attr)
                if isinstance(value, (int, float, str)):
                    qc_summary[attr] = value

    return SparseObservation(
        date=_parse_smosic_date(path, date_value),
        sensor="SMOS-IC",
        idx=idx,
        values=values,
        units=units,
        qc_summary=qc_summary,
    )


def read_tilecoord(path: Path | str) -> dict[str, np.ndarray]:
    """Read a GEOSldas little-endian binary tilecoord file."""

    path = Path(path)
    tile_coord: dict[str, np.ndarray] = {}
    int_fields = {"tile_id", "typ", "pfaf", "i_indg", "j_indg"}
    fields = [
        "tile_id",
        "typ",
        "pfaf",
        "com_lon",
        "com_lat",
        "min_lon",
        "max_lon",
        "min_lat",
        "max_lat",
        "i_indg",
        "j_indg",
        "frac_cell",
        "frac_pfaf",
        "area",
        "elev",
    ]

    with path.open("rb") as fp:
        tag = _read_record_tag(fp)
        if tag != 4:
            raise ValueError(f"{path} N_tile record has tag {tag}, expected 4")
        n_tile = struct.unpack("<i", _read_exact(fp, 4))[0]
        end_tag = _read_record_tag(fp)
        if end_tag != tag:
            raise ValueError(f"{path} N_tile record has mismatched record tags {tag} and {end_tag}")
        tile_coord["N_tile"] = np.asarray(n_tile, dtype=np.int32)

        for field in fields:
            dtype = np.dtype("<i4") if field in int_fields else np.dtype("<f4")
            expected_bytes = n_tile * dtype.itemsize
            tag = _read_record_tag(fp)
            if tag != expected_bytes:
                raise ValueError(f"{path} field {field} has record tag {tag}, expected {expected_bytes}")
            tile_coord[field] = np.frombuffer(_read_exact(fp, expected_bytes), dtype=dtype).copy()
            end_tag = _read_record_tag(fp)
            if end_tag != tag:
                raise ValueError(f"{path} field {field} has mismatched record tags {tag} and {end_tag}")

    return tile_coord


def representative_tiles_from_tilecoord(tilecoord: dict[str, np.ndarray], nx: int = M36_NX) -> RepresentativeTiles:
    """Choose the deterministic representative tile for each M36 cell."""

    i_indg = np.asarray(tilecoord["i_indg"], dtype=np.int64)
    j_indg = np.asarray(tilecoord["j_indg"], dtype=np.int64)
    tile_id = np.asarray(tilecoord["tile_id"], dtype=np.int64)
    frac_cell = np.asarray(tilecoord["frac_cell"], dtype=np.float64)

    frac_sort = np.where(np.isfinite(frac_cell), frac_cell, -np.inf)
    m36_linear = j_indg * int(nx) + i_indg

    order = np.lexsort((tile_id, -frac_sort, m36_linear))
    sorted_linear = m36_linear[order]
    if sorted_linear.size == 0:
        return RepresentativeTiles(
            m36_linear=np.array([], dtype=np.int64),
            tile_index=np.array([], dtype=np.int64),
        )

    first = np.empty(sorted_linear.size, dtype=bool)
    first[0] = True
    first[1:] = sorted_linear[1:] != sorted_linear[:-1]
    rep_idx = order[first]

    return RepresentativeTiles(
        m36_linear=m36_linear[rep_idx].astype(np.int64),
        tile_index=rep_idx.astype(np.int64),
    )


def read_model_tile_values(path: Path | str, variable: str = "SFMC") -> ModelTileValues:
    """Read one GEOSldas daily model variable from a tile-space NetCDF file."""

    path = Path(path)
    Dataset = _dataset_cls()

    with Dataset(path, "r") as ds:
        if variable not in ds.variables:
            raise KeyError(f"{path} is missing model variable {variable!r}")
        var = ds.variables[variable]
        data = _filled(var[:])
        if data.ndim == 2:
            data = data[0, :]
        elif data.ndim != 1:
            raise ValueError(f"{path}:{variable} has unsupported shape {data.shape}")
        units = getattr(var, "units", "")

    return ModelTileValues(values=np.asarray(data, dtype=np.float64), units=units)


def pair_sparse_observation_with_model(
    observation: SparseObservation,
    model_path: Path | str,
    tilecoord_path: Path | str,
    run: str,
    model_variable: str = "SFMC",
    nx: int = M36_NX,
) -> DailyPair:
    """Match a sparse observation to representative-tile model values."""

    tilecoord = read_tilecoord(tilecoord_path)
    reps = representative_tiles_from_tilecoord(tilecoord, nx=nx)
    model = read_model_tile_values(model_path, variable=model_variable)

    order = np.argsort(reps.m36_linear)
    sorted_linear = reps.m36_linear[order]
    sorted_tile_index = reps.tile_index[order]

    pos = np.searchsorted(sorted_linear, observation.idx)
    in_range = pos < sorted_linear.size
    matched = np.zeros(observation.idx.size, dtype=bool)
    matched[in_range] = sorted_linear[pos[in_range]] == observation.idx[in_range]

    tile_index = sorted_tile_index[pos[matched]]
    obs = observation.values[matched]
    mod = model.values[tile_index]
    good = np.isfinite(obs) & np.isfinite(mod)

    return DailyPair(
        date=observation.date,
        sensor=observation.sensor,
        run=run,
        idx=observation.idx[matched][good],
        obs=obs[good],
        model=mod[good],
        obs_units=observation.units,
        model_units=model.units,
    )


def read_smosic_model_pair(
    smosic_path: Path | str,
    model_path: Path | str,
    tilecoord_path: Path | str,
    run: str,
    model_variable: str = "SFMC",
    nx: int = M36_NX,
) -> DailyPair:
    """Read SMOS-IC and GEOSldas files and return their matched daily pair."""

    obs = read_smosic_sparse(smosic_path)
    return pair_sparse_observation_with_model(
        observation=obs,
        model_path=model_path,
        tilecoord_path=tilecoord_path,
        run=run,
        model_variable=model_variable,
        nx=nx,
    )


def _parse_smosic_date(path: Path, date_value: object | None) -> datetime.date:
    if date_value is not None:
        text = str(date_value)
    else:
        match = SMOSIC_DATE_RE.search(path.name)
        if match is None:
            raise ValueError(f"Cannot infer SMOS-IC date from {path}")
        text = match.group(1)
    return datetime.strptime(text, "%Y%m%d").date()


def _filled(arr: np.ndarray | np.ma.MaskedArray) -> np.ndarray:
    if np.ma.isMaskedArray(arr):
        return np.asarray(arr.filled(np.nan))
    return np.asarray(arr)


def _read_exact(fp, n: int) -> bytes:
    data = fp.read(n)
    if len(data) != n:
        raise EOFError(f"Expected {n} bytes, got {len(data)}")
    return data


def _read_record_tag(fp) -> int:
    return struct.unpack("<i", _read_exact(fp, 4))[0]


def _dataset_cls():
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise ImportError("netCDF4 is required for IV/TC NetCDF readers") from exc
    return Dataset
