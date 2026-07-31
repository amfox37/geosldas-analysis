"""CYGNSS L1 scaling climatology on the GEOS-LDAS owner-tile grid."""

from __future__ import annotations

from dataclasses import dataclass
from math import floor
from pathlib import Path
from typing import Sequence

import numpy as np

from .io import (
    DateTime,
    ObsFcstAnaRecord,
    ObsParam,
    augment_date_time,
    days_in_month,
    read_obs_fcst_ana,
)
from .owner_grid_io import FIELD_NAMES, write_netcdf_owner_grid
from .tb_tile_stats import _date_range_tags, _obsfcstana_path
from .tile_io import TileCoordinates, TileGrid, read_tilecoord, read_tilegrids


NODATA = -9999.0
N_PENTADS = 73
CYGNSS_L1_DESCRIPTION = "CYGNSS_L1_DDM3X5_CROP_SCALAR"
CYGNSS_L1_VARNAME = "cygl1scal"
CYGNSS_L1_UNITS = "dB"


@dataclass(frozen=True)
class OwnerGridMapping:
    """Map experiment tile positions to zero-based global owner-grid cells."""

    i: np.ndarray
    j: np.ndarray
    linear: np.ndarray
    grid: TileGrid


def build_owner_grid_mapping(
    tile_coord: TileCoordinates,
    global_grid: TileGrid,
) -> OwnerGridMapping:
    """Convert global tile indices to the dense file's zero-based coordinates."""

    if global_grid.gridtype != "EASEv2_M36":
        raise ValueError(
            f"CYGNSS L1 owner scaling requires EASEv2_M36; got {global_grid.gridtype!r}"
        )
    i = tile_coord.i_indg.astype(np.int64) - global_grid.i_offg - global_grid.ind_base
    j = tile_coord.j_indg.astype(np.int64) - global_grid.j_offg - global_grid.ind_base
    if np.any(i < 0) or np.any(i >= global_grid.n_lon):
        raise ValueError("tilecoord i_indg values fall outside the global tile grid")
    if np.any(j < 0) or np.any(j >= global_grid.n_lat):
        raise ValueError("tilecoord j_indg values fall outside the global tile grid")
    linear = j * global_grid.n_lon + i
    if np.unique(linear).size != linear.size:
        raise ValueError("More than one model tile occupies an owner-grid cell")
    return OwnerGridMapping(i=i, j=j, linear=linear, grid=global_grid)


def select_cygnss_l1_species(
    obs_params: Sequence[ObsParam],
    description: str = CYGNSS_L1_DESCRIPTION,
) -> ObsParam:
    """Resolve and validate exactly one CYGNSS L1 species from ldas_obsparam."""

    matches = [param for param in obs_params if param.descr == description]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one ldas_obsparam entry named {description!r}; "
            f"found {[param.species for param in matches]}"
        )
    selected = matches[0]
    expected = (CYGNSS_L1_VARNAME, CYGNSS_L1_UNITS, CYGNSS_L1_VARNAME, CYGNSS_L1_UNITS)
    actual = (selected.varname, selected.units, selected.fcstvarname, selected.fcstunits)
    if actual != expected:
        raise ValueError(f"Unexpected CYGNSS L1 variable metadata: {actual}; expected {expected}")
    if _logical_value(selected.scale, "scale"):
        raise ValueError("CYGNSS L1 scaling climatology requires an unscaled monitoring run")
    return selected


def _logical_value(value: str, name: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"t", "true", ".true.", "1"}:
        return True
    if normalized in {"f", "false", ".false.", "0"}:
        return False
    raise ValueError(f"Unrecognized {name} logical value: {value!r}")


def _validate_ofa_metadata(record: ObsFcstAnaRecord, selected: ObsParam, path: Path) -> None:
    if record.species_metadata is None:
        raise ValueError(f"ObsFcstAna file lacks per-species scaling metadata: {path}")
    matches = [item for item in record.species_metadata if item.species == selected.species]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one metadata entry for species {selected.species} in {path}; "
            f"found {len(matches)}"
        )
    metadata = matches[0]
    actual = (
        metadata.descr,
        metadata.varname,
        metadata.units,
        metadata.fcstvarname,
        metadata.fcstunits,
    )
    expected = (
        selected.descr,
        CYGNSS_L1_VARNAME,
        CYGNSS_L1_UNITS,
        CYGNSS_L1_VARNAME,
        CYGNSS_L1_UNITS,
    )
    if actual != expected:
        raise ValueError(f"CYGNSS L1 metadata mismatch in {path}: {actual}; expected {expected}")
    if metadata.scale:
        raise ValueError(f"CYGNSS L1 observations are already scaled in {path}")


def _compute_window_statistics(
    ring: np.ndarray,
    *,
    ndata_min: int,
    std_epsilon: float,
) -> np.ndarray:
    count = ring[4].sum(axis=-1)
    obs_sum = ring[0].sum(axis=-1)
    obs_sum2 = ring[1].sum(axis=-1)
    mod_sum = ring[2].sum(axis=-1)
    mod_sum2 = ring[3].sum(axis=-1)
    obs_mean = np.divide(obs_sum, count, out=np.full_like(obs_sum, np.nan), where=count > 0)
    mod_mean = np.divide(mod_sum, count, out=np.full_like(mod_sum, np.nan), where=count > 0)
    obs_second = np.divide(
        obs_sum2, count, out=np.full_like(obs_sum2, np.nan), where=count > 0
    )
    mod_second = np.divide(
        mod_sum2, count, out=np.full_like(mod_sum2, np.nan), where=count > 0
    )
    obs_std = np.sqrt(np.maximum(obs_second - obs_mean**2, 0.0))
    mod_std = np.sqrt(np.maximum(mod_second - mod_mean**2, 0.0))
    model_min = ring[5].min(axis=-1)
    model_max = ring[6].max(axis=-1)
    model_min[~np.isfinite(model_min)] = np.nan
    model_max[~np.isfinite(model_max)] = np.nan

    output = np.stack((obs_mean, obs_std, mod_mean, mod_std, count, model_min, model_max))
    valid = (count >= ndata_min) & (obs_std > std_epsilon) & (mod_std > std_epsilon)
    output[:, ~valid] = np.nan
    return output


def generate_cygnss_l1_scaling_params(
    *,
    run_months: Sequence[int],
    exp_path: str | Path,
    exp_run: str,
    domain: str,
    start_year: Sequence[int],
    end_year: Sequence[int],
    dt_assim: int,
    t0_assim: int,
    obs_params: Sequence[ObsParam],
    window_days: int,
    ndata_min: int,
    prefix: str,
    description: str = CYGNSS_L1_DESCRIPTION,
    out_dir: str = "cygnss_l1_z_score_clim",
    std_epsilon: float = 1e-6,
    tilecoord_path: str | Path | None = None,
    tilegrids_path: str | Path | None = None,
    overwrite: bool = True,
) -> Path:
    """Generate one 73-pentad CYGNSS L1 owner-tile scaling file."""

    if window_days % 5 or window_days % 10 == 0:
        raise ValueError("window_days must be an odd number of pentads (for example, 75)")
    if ndata_min < 1:
        raise ValueError("ndata_min must be at least 1")
    if std_epsilon < 0:
        raise ValueError("std_epsilon cannot be negative")
    if dt_assim <= 0:
        raise ValueError("dt_assim must be positive")
    if len(run_months) != len(start_year) or len(run_months) != len(end_year):
        raise ValueError("run_months, start_year, and end_year must have equal length")

    selected = select_cygnss_l1_species(obs_params, description)
    inpath = Path(exp_path) / exp_run / "output" / domain
    rc_out = inpath / "rc_out"
    tilecoord = Path(tilecoord_path) if tilecoord_path else rc_out / f"{exp_run}.ldas_tilecoord.bin"
    tilegrids = Path(tilegrids_path) if tilegrids_path else rc_out / f"{exp_run}.ldas_tilegrids.bin"
    tile_coord = read_tilecoord(tilecoord)
    global_grid, _ = read_tilegrids(tilegrids)
    mapping = build_owner_grid_mapping(tile_coord, global_grid)

    # sum(obs), sum(obs^2), sum(model), sum(model^2), count, model min, model max
    ring = np.zeros((len(FIELD_NAMES), tile_coord.n_tile, window_days), dtype=np.float64)
    ring[5] = np.inf
    ring[6] = -np.inf
    data_out = np.full((len(FIELD_NAMES), tile_coord.n_tile, N_PENTADS), np.nan)
    default_time = DateTime(2014, 1, 1, 0, 0, 0)
    start_times = [default_time for _ in range(N_PENTADS)]
    end_times = [default_time for _ in range(N_PENTADS)]

    count_days = 0
    t0_assim %= dt_assim
    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month(2014, month) + 1):
            if count_days < window_days:
                count_days += 1
            day_index = count_days - 1
            last_time = DateTime(2014, month, day, 0, 0, 0)
            for seconds in range(t0_assim, 86400, dt_assim):
                hour = seconds // 3600
                minute = (seconds % 3600) // 60
                last_time = DateTime(2014, month, day, hour, minute, 0)
                for year in range(start_year[imonth], end_year[imonth] + 1):
                    actual_time = DateTime(year, month, day, hour, minute, 0)
                    path = _obsfcstana_path(inpath, exp_run, actual_time, "nc4")
                    if path is None or not path.exists():
                        continue
                    record = read_obs_fcst_ana(path)
                    if record is None:
                        continue
                    _validate_ofa_metadata(record, selected, path)
                    selected_mask = record.obs_species == selected.species
                    if not np.any(selected_mask):
                        continue
                    model_indices = record.obs_tilenum[selected_mask].astype(np.int64) - 1
                    if np.any(model_indices < 0) or np.any(model_indices >= tile_coord.n_tile):
                        raise ValueError(
                            f"Species {selected.species} contains out-of-range tile "
                            f"numbers in {path}"
                        )
                    unique, multiplicity = np.unique(model_indices, return_counts=True)
                    duplicate = multiplicity > 1
                    if np.any(duplicate):
                        details = ", ".join(
                            f"tile {int(tile + 1)} x{int(count)}"
                            for tile, count in zip(unique[duplicate], multiplicity[duplicate])
                        )
                        raise ValueError(f"Duplicate CYGNSS L1 owner tiles in {path}: {details}")

                    obs = record.obs_obs[selected_mask].astype(np.float64)
                    model = record.obs_fcst[selected_mask].astype(np.float64)
                    valid = (
                        np.isfinite(obs)
                        & np.isfinite(model)
                        & (np.abs(obs - NODATA) > 1e-4)
                        & (np.abs(model - NODATA) > 1e-4)
                    )
                    indices = model_indices[valid]
                    obs = obs[valid]
                    model = model[valid]
                    np.add.at(ring[0, :, day_index], indices, obs)
                    np.add.at(ring[1, :, day_index], indices, obs**2)
                    np.add.at(ring[2, :, day_index], indices, model)
                    np.add.at(ring[3, :, day_index], indices, model**2)
                    np.add.at(ring[4, :, day_index], indices, 1.0)
                    np.minimum.at(ring[5, :, day_index], indices, model)
                    np.maximum.at(ring[6, :, day_index], indices, model)

            if count_days >= window_days:
                center_time = augment_date_time(-floor(window_days * 86400 / 2), last_time)
                doy = int(center_time.dofyr)
                if (doy + 2) % 5 == 0:
                    pentad = (doy + 2) // 5
                    data_out[:, :, pentad - 1] = _compute_window_statistics(
                        ring, ndata_min=ndata_min, std_epsilon=std_epsilon
                    )
                    start_times[pentad - 1] = augment_date_time(-window_days * 86400, last_time)
                    end_times[pentad - 1] = last_time

                ring[..., :-1] = ring[..., 1:]
                ring[:5, :, -1] = 0.0
                ring[5, :, -1] = np.inf
                ring[6, :, -1] = -np.inf

    doy_start, doy_end, _, _ = _date_range_tags(run_months, start_year, end_year)
    output_path = inpath / "stats" / out_dir / (
        f"{prefix}{min(start_year)}_doy{doy_start}_{max(end_year)}_doy{doy_end}"
        f"_W_{window_days}d_Nmin_{ndata_min}_sp_{description}_all_pentads.nc4"
    )
    write_netcdf_owner_grid(
        output_path,
        data=data_out,
        tile_i=mapping.i,
        tile_j=mapping.j,
        grid=mapping.grid,
        pentads=np.arange(1, N_PENTADS + 1),
        start_time=start_times,
        end_time=end_times,
        window_days=window_days,
        ndata_min=ndata_min,
        std_epsilon=std_epsilon,
        species_id=selected.species,
        species_description=description,
        overwrite=overwrite,
    )
    return output_path


__all__ = [
    "CYGNSS_L1_DESCRIPTION",
    "OwnerGridMapping",
    "build_owner_grid_mapping",
    "generate_cygnss_l1_scaling_params",
    "select_cygnss_l1_species",
]
