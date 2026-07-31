"""NetCDF output for owner-tile scaling climatologies."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Sequence

import netCDF4 as nc
import numpy as np

from .io import DateTime
from .tile_io import TileGrid


FILL_VALUE = -9999.0
FIELD_NAMES = ("o_mean", "o_std", "m_mean", "m_std", "n_data", "m_min", "m_max")


def write_netcdf_owner_grid(
    fname: str | Path,
    *,
    data: np.ndarray,
    tile_i: np.ndarray,
    tile_j: np.ndarray,
    grid: TileGrid,
    pentads: Sequence[int],
    start_time: Sequence[DateTime],
    end_time: Sequence[DateTime],
    window_days: int,
    ndata_min: int,
    std_epsilon: float,
    species_id: int,
    species_description: str,
    overwrite: bool = True,
) -> None:
    """Write local-tile statistics onto a dense global EASE owner grid."""

    path = Path(fname)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists and overwrite is False")
    path.parent.mkdir(parents=True, exist_ok=True)

    values = np.asarray(data, dtype=np.float64)
    pentad_values = np.asarray(pentads, dtype=np.int32)
    tile_i = np.asarray(tile_i, dtype=np.int64)
    tile_j = np.asarray(tile_j, dtype=np.int64)
    expected = (len(FIELD_NAMES), tile_i.size, pentad_values.size)
    if values.shape != expected:
        raise ValueError(f"data shape must be {expected}; got {values.shape}")
    if tile_j.size != tile_i.size:
        raise ValueError("tile_i and tile_j must have equal length")
    if len(start_time) != pentad_values.size or len(end_time) != pentad_values.size:
        raise ValueError("start_time/end_time must contain one value per pentad")
    if np.any(tile_i < 0) or np.any(tile_i >= grid.n_lon):
        raise ValueError("Owner-tile column indices fall outside the global grid")
    if np.any(tile_j < 0) or np.any(tile_j >= grid.n_lat):
        raise ValueError("Owner-tile row indices fall outside the global grid")
    linear = tile_j * grid.n_lon + tile_i
    if np.unique(linear).size != linear.size:
        raise ValueError("More than one model tile maps to an owner-grid cell")

    epoch = datetime(1950, 1, 1)

    def days_since_1950(items: Sequence[DateTime]) -> np.ndarray:
        return np.asarray(
            [
                (
                    datetime(item.year, item.month, item.day, item.hour, item.minute, item.second)
                    - epoch
                ).total_seconds()
                / 86400.0
                for item in items
            ],
            dtype=np.float64,
        )

    with nc.Dataset(path, "w", format="NETCDF4") as dataset:
        dataset.setncatts(
            {
                "title": "CYGNSS L1 owner-tile z-score scaling climatology",
                "schema_version": "1.0",
                "gridtype": grid.gridtype,
                "source_ind_base": grid.ind_base,
                "source_i_offg": grid.i_offg,
                "source_j_offg": grid.j_offg,
                "window_days": window_days,
                "ndata_min": ndata_min,
                "std_epsilon": std_epsilon,
                "standard_deviation_convention": "population (ddof=0)",
                "species_id": species_id,
                "species_description": species_description,
                "spatial_key": "global owner tile i_indg/j_indg",
                "sampling_contract": "at most one selected observation per owner tile and cycle",
            }
        )
        dataset.createDimension("pentad", pentad_values.size)
        dataset.createDimension("y", grid.n_lat)
        dataset.createDimension("x", grid.n_lon)

        pentad_var = dataset.createVariable("pentad", "i4", ("pentad",))
        pentad_var[:] = pentad_values
        pentad_var.setncatts({"long_name": "pentad of year", "units": "1", "axis": "T"})

        x_var = dataset.createVariable("x", "i4", ("x",))
        x_var[:] = np.arange(grid.n_lon, dtype=np.int32)
        x_var.setncatts(
            {"long_name": "zero-based global EASE-grid column index", "units": "1", "axis": "X"}
        )
        y_var = dataset.createVariable("y", "i4", ("y",))
        y_var[:] = np.arange(grid.n_lat, dtype=np.int32)
        y_var.setncatts(
            {"long_name": "zero-based global EASE-grid row index", "units": "1", "axis": "Y"}
        )

        time_attrs = {"units": "days since 1950-01-01 00:00:00", "calendar": "standard"}
        start_var = dataset.createVariable("start_time", "f8", ("pentad",))
        end_var = dataset.createVariable("end_time", "f8", ("pentad",))
        start_var[:] = days_since_1950(start_time)
        end_var[:] = days_since_1950(end_time)
        start_var.setncatts(time_attrs)
        end_var.setncatts(time_attrs)

        descriptions = {
            "o_mean": "Mean unscaled CYGNSS L1 observation",
            "o_std": "Population standard deviation of unscaled CYGNSS L1 observations",
            "m_mean": "Mean CYGNSS L1 model forecast equivalent",
            "m_std": "Population standard deviation of CYGNSS L1 model forecast equivalents",
            "n_data": "Number of paired observation and forecast samples",
            "m_min": "Minimum CYGNSS L1 model forecast equivalent over all pentads",
            "m_max": "Maximum CYGNSS L1 model forecast equivalent over all pentads",
        }
        chunks = (1, min(64, grid.n_lat), min(128, grid.n_lon))
        variables: dict[str, nc.Variable] = {}
        for name in FIELD_NAMES[:5]:
            variable = dataset.createVariable(
                name,
                "f8",
                ("pentad", "y", "x"),
                fill_value=FILL_VALUE,
                zlib=True,
                complevel=5,
                chunksizes=chunks,
            )
            variable.setncatts(
                {"long_name": descriptions[name], "units": "1" if name == "n_data" else "dB"}
            )
            variables[name] = variable

        for field, name in enumerate(FIELD_NAMES[:5]):
            for p_index in range(pentad_values.size):
                dense = np.full((grid.n_lat, grid.n_lon), FILL_VALUE, dtype=np.float64)
                tile_values = values[field, :, p_index]
                valid = np.isfinite(tile_values)
                dense[tile_j[valid], tile_i[valid]] = tile_values[valid]
                variables[name][p_index, :, :] = dense

        for field, name in ((5, "m_min"), (6, "m_max")):
            variable = dataset.createVariable(
                name,
                "f8",
                ("y", "x"),
                fill_value=FILL_VALUE,
                zlib=True,
                complevel=5,
                chunksizes=chunks[1:],
            )
            variable.setncatts({"long_name": descriptions[name], "units": "dB"})
            finite = np.isfinite(values[field])
            has_data = finite.any(axis=1)
            if name == "m_min":
                tile_values = np.where(finite, values[field], np.inf).min(axis=1)
            else:
                tile_values = np.where(finite, values[field], -np.inf).max(axis=1)
            tile_values[~has_data] = np.nan
            dense = np.full((grid.n_lat, grid.n_lon), FILL_VALUE, dtype=np.float64)
            valid = np.isfinite(tile_values)
            dense[tile_j[valid], tile_i[valid]] = tile_values[valid]
            variable[:, :] = dense


__all__ = ["FIELD_NAMES", "FILL_VALUE", "write_netcdf_owner_grid"]
