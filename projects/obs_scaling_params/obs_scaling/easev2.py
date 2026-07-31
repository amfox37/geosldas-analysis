"""EASE-Grid 2.0 coordinate transforms used by GEOS-LDAS tile scaling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EaseGrid:
    scale_m: float
    n_cols: int
    n_rows: int

    @property
    def col_origin(self) -> float:
        return (self.n_cols - 1) / 2.0

    @property
    def row_origin(self) -> float:
        return (self.n_rows - 1) / 2.0


_GRIDS = {
    "M36": EaseGrid(36032.220840584, 964, 406),
    "M09": EaseGrid(9008.055210146, 3856, 1624),
    "M03": EaseGrid(3002.6850700487, 11568, 4872),
    "M01": EaseGrid(1000.89502334956, 34704, 14616),
}

_RADIUS_M = 6378137.0
_ECCENTRICITY = 0.081819190843
_E2 = _ECCENTRICITY**2
_PHI1 = np.deg2rad(30.0)
_KZ = np.cos(_PHI1) / np.sqrt(1.0 - _E2 * np.sin(_PHI1) ** 2)


def _grid(grid_id: str) -> EaseGrid:
    try:
        return _GRIDS[grid_id]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported EASEv2 grid {grid_id!r}; expected one of {sorted(_GRIDS)}"
        ) from exc


def latlon_to_ind(
    lat: np.ndarray | float,
    lon: np.ndarray | float,
    grid_id: str,
    *,
    rounded: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return zero-based EASEv2 row and column indices for latitude/longitude."""

    grid = _grid(grid_id)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    dlon = lon.copy()
    dlon = np.where(dlon < -180.0, dlon + 360.0, dlon)
    dlon = np.where(dlon > 180.0, dlon - 360.0, dlon)
    phi = np.deg2rad(lat)
    lam = np.deg2rad(dlon)
    sin_phi = np.sin(phi)
    q = (1.0 - _E2) * (
        sin_phi / (1.0 - _E2 * sin_phi**2)
        - np.log((1.0 - _ECCENTRICITY * sin_phi) / (1.0 + _ECCENTRICITY * sin_phi))
        / (2.0 * _ECCENTRICITY)
    )
    x = _RADIUS_M * _KZ * lam
    y = _RADIUS_M * q / (2.0 * _KZ)
    row = grid.row_origin - y / grid.scale_m
    col = grid.col_origin + x / grid.scale_m
    if rounded:
        # MATLAB round() is half-away-from-zero. EASE indices are non-negative.
        row = np.floor(row + 0.5).astype(np.int64)
        col = np.floor(col + 0.5).astype(np.int64)
    return row, col


def ind_to_latlon(
    row: np.ndarray | float,
    col: np.ndarray | float,
    grid_id: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return latitude/longitude at zero-based EASEv2 row and column indices."""

    grid = _grid(grid_id)
    row = np.asarray(row, dtype=np.float64)
    col = np.asarray(col, dtype=np.float64)
    x = (col - grid.col_origin) * grid.scale_m
    y = (grid.row_origin - row) * grid.scale_m
    e4 = _ECCENTRICITY**4
    e6 = _ECCENTRICITY**6
    qp = (1.0 - _E2) * (
        1.0 / (1.0 - _E2)
        - np.log((1.0 - _ECCENTRICITY) / (1.0 + _ECCENTRICITY))
        / (2.0 * _ECCENTRICITY)
    )
    beta = np.arcsin(2.0 * y * _KZ / (_RADIUS_M * qp))
    lam = x / (_RADIUS_M * _KZ)
    phi = (
        beta
        + (_E2 / 3.0 + 31.0 * e4 / 180.0 + 517.0 * e6 / 5040.0) * np.sin(2.0 * beta)
        + (23.0 * e4 / 360.0 + 251.0 * e6 / 3780.0) * np.sin(4.0 * beta)
        + 761.0 * e6 / 45360.0 * np.sin(6.0 * beta)
    )
    lat = np.rad2deg(phi)
    lon = np.rad2deg(lam)
    lon = np.where(lon < -180.0, lon + 360.0, lon)
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    return lat, lon


def grid_id_from_gridtype(gridtype: str) -> str:
    """Translate a GEOS-LDAS EASEv2 gridtype string to an EASE grid id."""

    normalized = gridtype.strip().replace("-", "_")
    for grid_id in _GRIDS:
        if normalized.endswith(grid_id):
            return grid_id
    raise ValueError(f"Unsupported GEOS-LDAS tile grid type {gridtype!r}")
