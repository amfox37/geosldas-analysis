"""GEOS-LDAS Tb scaling-parameter sequential-binary reader and writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np

from .io import DateTime


@dataclass(frozen=True)
class TbScalingFile:
    asc_flag: int
    ndata_min: int
    start_time: DateTime
    end_time: DateTime
    angles: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    tile_id: np.ndarray
    data: np.ndarray


def _write_record(fp: BinaryIO, values: np.ndarray, dtype: str) -> None:
    payload = np.asarray(values, dtype=dtype).tobytes(order="C")
    tag = np.asarray([len(payload)], dtype=">i4").tobytes()
    fp.write(tag)
    fp.write(payload)
    fp.write(tag)


def _read_record(fp: BinaryIO, dtype: str, count: int | None = None) -> np.ndarray:
    tag_data = fp.read(4)
    if len(tag_data) != 4:
        raise EOFError("Unexpected EOF before Fortran record")
    tag = int(np.frombuffer(tag_data, dtype=">i4")[0])
    payload = fp.read(tag)
    if len(payload) != tag:
        raise EOFError("Unexpected EOF in Fortran record")
    end_tag_data = fp.read(4)
    if len(end_tag_data) != 4 or int(np.frombuffer(end_tag_data, dtype=">i4")[0]) != tag:
        raise ValueError("Mismatched Fortran record tags")
    values = np.frombuffer(payload, dtype=dtype).copy()
    if count is not None and values.size != count:
        raise ValueError(f"Record contains {values.size} values; expected {count}")
    return values


def _time_values(value: DateTime) -> np.ndarray:
    return np.array([value.year, value.month, value.day, value.hour, value.minute])


def write_tb_scaling_file(
    path: str | Path,
    *,
    asc_flag: int,
    ndata_min: int,
    start_time: DateTime,
    end_time: DateTime,
    angles: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    tile_id: np.ndarray,
    data: np.ndarray,
    overwrite: bool = True,
) -> None:
    """Write the 14-field format consumed by ``scale_obs_Tb_zscore``."""

    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    angles = np.atleast_1d(np.asarray(angles, dtype=np.float32))
    lon = np.asarray(lon, dtype=np.float32)
    lat = np.asarray(lat, dtype=np.float32)
    tile_id = np.asarray(tile_id)
    if tile_id.ndim == 1:
        tile_id = tile_id[:, None]
    data = np.asarray(data, dtype=np.float32)
    n_grid = lon.size
    expected_shape = (14, n_grid, angles.size)
    if data.shape != expected_shape:
        raise ValueError(f"data shape is {data.shape}; expected {expected_shape}")
    if lat.size != n_grid or tile_id.shape[0] != n_grid:
        raise ValueError("lon, lat, tile_id, and data grid dimensions do not match")
    if asc_flag not in (0, 1):
        raise ValueError("asc_flag must be 0 (descending) or 1 (ascending)")

    with path.open("wb") as fp:
        _write_record(fp, np.array([asc_flag, ndata_min]), ">i4")
        _write_record(fp, _time_values(start_time), ">i4")
        _write_record(fp, _time_values(end_time), ">i4")
        _write_record(fp, np.array([n_grid, angles.size, tile_id.shape[1]]), ">i4")
        _write_record(fp, angles, ">f4")
        _write_record(fp, lon, ">f4")
        _write_record(fp, lat, ">f4")
        for column in range(tile_id.shape[1]):
            _write_record(fp, tile_id[:, column], ">i4")
        for field in range(data.shape[0]):
            for angle in range(data.shape[2]):
                _write_record(fp, data[field, :, angle], ">f4")


def read_tb_scaling_file(path: str | Path) -> TbScalingFile:
    """Read and validate a GEOS-LDAS Tb scaling-parameter file."""

    path = Path(path)
    with path.open("rb") as fp:
        header = _read_record(fp, ">i4", 2)
        start = _read_record(fp, ">i4", 5)
        end = _read_record(fp, ">i4", 5)
        dimensions = _read_record(fp, ">i4", 3)
        n_grid, n_angle, n_tile_columns = (int(value) for value in dimensions)
        angles = _read_record(fp, ">f4", n_angle)
        lon = _read_record(fp, ">f4", n_grid)
        lat = _read_record(fp, ">f4", n_grid)
        tile_columns = [_read_record(fp, ">i4", n_grid) for _ in range(n_tile_columns)]
        tile_id = np.column_stack(tile_columns)
        data = np.empty((14, n_grid, n_angle), dtype=np.float32)
        for field in range(14):
            for angle in range(n_angle):
                data[field, :, angle] = _read_record(fp, ">f4", n_grid)
        if fp.read(1):
            raise ValueError(f"Unexpected trailing data in {path}")
    if tile_id.shape[1] == 1:
        tile_id = tile_id[:, 0]
    return TbScalingFile(
        asc_flag=int(header[0]),
        ndata_min=int(header[1]),
        start_time=DateTime(*[int(value) for value in start], second=0),
        end_time=DateTime(*[int(value) for value in end], second=0),
        angles=angles,
        lon=lon,
        lat=lat,
        tile_id=tile_id,
        data=data,
    )
