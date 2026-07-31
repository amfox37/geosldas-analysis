"""Readers for GEOS-LDAS little-endian tile coordinate/grid files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np


@dataclass(frozen=True)
class TileCoordinates:
    tile_id: np.ndarray
    typ: np.ndarray
    pfaf: np.ndarray
    com_lon: np.ndarray
    com_lat: np.ndarray
    min_lon: np.ndarray
    max_lon: np.ndarray
    min_lat: np.ndarray
    max_lat: np.ndarray
    i_indg: np.ndarray
    j_indg: np.ndarray
    frac_cell: np.ndarray
    frac_pfaf: np.ndarray
    area: np.ndarray
    elev: np.ndarray

    @property
    def n_tile(self) -> int:
        return int(self.tile_id.size)


@dataclass(frozen=True)
class TileGrid:
    gridtype: str
    ind_base: int
    i_dir: int
    j_dir: int
    n_lon: int
    n_lat: int
    i_offg: int
    j_offg: int
    ll_lon: float
    ll_lat: float
    ur_lon: float
    ur_lat: float
    dlon: float
    dlat: float


def _read_exact(fp: BinaryIO, nbytes: int) -> bytes:
    value = fp.read(nbytes)
    if len(value) != nbytes:
        raise EOFError(f"Expected {nbytes} bytes, found {len(value)}")
    return value


def _read_record(fp: BinaryIO, dtype: str, count: int) -> np.ndarray:
    tag = int(np.frombuffer(_read_exact(fp, 4), dtype="<i4")[0])
    expected = np.dtype(dtype).itemsize * count
    if tag != expected:
        raise ValueError(f"Fortran record is {tag} bytes; expected {expected}")
    values = np.frombuffer(_read_exact(fp, expected), dtype=dtype).copy()
    end_tag = int(np.frombuffer(_read_exact(fp, 4), dtype="<i4")[0])
    if end_tag != tag:
        raise ValueError(f"Mismatched Fortran record tags: {tag} and {end_tag}")
    return values


def read_tilecoord(path: str | Path) -> TileCoordinates:
    """Read a GEOS-LDAS ``ldas_tilecoord.bin`` file."""

    path = Path(path)
    with path.open("rb") as fp:
        n_tile = int(_read_record(fp, "<i4", 1)[0])
        int_fields = {"tile_id", "typ", "pfaf", "i_indg", "j_indg"}
        fields = [
            "tile_id", "typ", "pfaf", "com_lon", "com_lat", "min_lon", "max_lon",
            "min_lat", "max_lat", "i_indg", "j_indg", "frac_cell", "frac_pfaf",
            "area", "elev",
        ]
        values = {
            field: _read_record(fp, "<i4" if field in int_fields else "<f4", n_tile)
            for field in fields
        }
        if fp.read(1):
            raise ValueError(f"Unexpected trailing data in {path}")
    return TileCoordinates(**values)


def _read_tilegrid_record(fp: BinaryIO) -> TileGrid:
    tag = int(np.frombuffer(_read_exact(fp, 4), dtype="<i4")[0])
    expected = 40 + 7 * 4 + 6 * 4
    if tag != expected:
        raise ValueError(f"Tile-grid record is {tag} bytes; expected {expected}")
    gridtype = _read_exact(fp, 40).decode("ascii").rstrip("\x00 ")
    integers = np.frombuffer(_read_exact(fp, 7 * 4), dtype="<i4")
    floats = np.frombuffer(_read_exact(fp, 6 * 4), dtype="<f4")
    end_tag = int(np.frombuffer(_read_exact(fp, 4), dtype="<i4")[0])
    if end_tag != tag:
        raise ValueError(f"Mismatched tile-grid record tags: {tag} and {end_tag}")
    return TileGrid(
        gridtype=gridtype,
        ind_base=int(integers[0]),
        i_dir=int(integers[1]),
        j_dir=int(integers[2]),
        n_lon=int(integers[3]),
        n_lat=int(integers[4]),
        i_offg=int(integers[5]),
        j_offg=int(integers[6]),
        ll_lon=float(floats[0]),
        ll_lat=float(floats[1]),
        ur_lon=float(floats[2]),
        ur_lat=float(floats[3]),
        dlon=float(floats[4]),
        dlat=float(floats[5]),
    )


def read_tilegrids(path: str | Path) -> tuple[TileGrid, TileGrid]:
    """Read global and domain grids from ``ldas_tilegrids.bin``."""

    path = Path(path)
    with path.open("rb") as fp:
        global_grid = _read_tilegrid_record(fp)
        domain_grid = _read_tilegrid_record(fp)
        if fp.read(1):
            raise ValueError(f"Unexpected trailing data in {path}")
    return global_grid, domain_grid
