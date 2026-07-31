"""Tile-space SMAP Tb scaling climatology ported from the MATLAB workflow."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
import shutil
from typing import Sequence

import numpy as np

from .easev2 import grid_id_from_gridtype, ind_to_latlon, latlon_to_ind
from .io import DateTime, ObsParam, augment_date_time, days_in_month, read_obs_fcst_ana
from .seqbin import write_tb_scaling_file
from .tile_io import TileCoordinates, TileGrid, read_tilecoord, read_tilegrids


NODATA = -9999.0
N_FIELDS = 14


@dataclass(frozen=True)
class AdministeringTiles:
    """Mapping from model-tile positions to output scaling-tile positions."""

    model_to_output: np.ndarray
    model_indices: np.ndarray
    tile_id: np.ndarray
    lon: np.ndarray
    lat: np.ndarray

    @property
    def n_tiles(self) -> int:
        return int(self.model_indices.size)


@dataclass(frozen=True)
class TbSpecies:
    species: int
    description: str
    pol_index: int
    angle_index: int


def _distance_km_to_degrees(distance_km: float, lat: np.ndarray) -> tuple[np.ndarray, float]:
    distance_y = distance_km * 180.0 / np.pi / 6371.0
    distance_x = distance_y / np.cos(np.deg2rad(lat))
    return distance_x, distance_y


def build_administering_tiles(
    tile_coord: TileCoordinates,
    tile_grid: TileGrid,
    *,
    convert_grid: str | None = "EASEv2_M36",
    shift_lon: float = 0.01,
    shift_lat: float = 0.005,
    fov_km: float = 20.0,
) -> AdministeringTiles:
    """Reproduce MATLAB's M36 administering-tile selection."""

    if convert_grid is None:
        model_indices = np.arange(tile_coord.n_tile, dtype=np.int64)
    else:
        if convert_grid != "EASEv2_M36":
            raise ValueError("The initial Tb port supports convert_grid='EASEv2_M36' only")
        m36_row, m36_col = latlon_to_ind(
            tile_coord.com_lat, tile_coord.com_lon, "M36", rounded=True
        )
        m36_cells = np.column_stack((m36_row, m36_col))
        unique_cells = np.unique(m36_cells, axis=0)
        target_lat, target_lon = ind_to_latlon(unique_cells[:, 0], unique_cells[:, 1], "M36")
        target_lat = target_lat + shift_lat
        target_lon = target_lon + shift_lon

        local_i = tile_coord.i_indg.astype(np.int64) - tile_grid.i_offg - tile_grid.ind_base
        local_j = tile_coord.j_indg.astype(np.int64) - tile_grid.j_offg - tile_grid.ind_base
        if np.any(local_i < 0) or np.any(local_i >= tile_grid.n_lon):
            raise ValueError("tilecoord i_indg values fall outside the domain tile grid")
        if np.any(local_j < 0) or np.any(local_j >= tile_grid.n_lat):
            raise ValueError("tilecoord j_indg values fall outside the domain tile grid")
        linear = local_j * tile_grid.n_lon + local_i
        if np.unique(linear).size != linear.size:
            raise ValueError("More than one tile occupies an EASE grid cell")
        cell_to_tile = {int(cell): index for index, cell in enumerate(linear)}

        model_grid_id = grid_id_from_gridtype(tile_grid.gridtype)
        target_row, target_col = latlon_to_ind(
            target_lat, target_lon, model_grid_id, rounded=True
        )
        target_i = target_col - tile_grid.i_offg - tile_grid.ind_base
        target_j = target_row - tile_grid.j_offg - tile_grid.ind_base
        distance_x, distance_y = _distance_km_to_degrees(fov_km, target_lat)
        selected: list[int] = []
        for n in range(target_lat.size):
            if target_i[n] < 0 or target_i[n] >= tile_grid.n_lon:
                continue
            if target_j[n] < 0 or target_j[n] >= tile_grid.n_lat:
                continue
            model_index = cell_to_tile.get(int(target_j[n] * tile_grid.n_lon + target_i[n]))
            if model_index is None:
                continue
            outside_bbox = (
                target_lon[n] < tile_coord.min_lon[model_index]
                or target_lon[n] > tile_coord.max_lon[model_index]
                or target_lat[n] < tile_coord.min_lat[model_index]
                or target_lat[n] > tile_coord.max_lat[model_index]
            )
            too_far = (
                abs(target_lon[n] - tile_coord.com_lon[model_index]) > distance_x[n]
                or abs(target_lat[n] - tile_coord.com_lat[model_index]) > distance_y
            )
            if not (outside_bbox and too_far):
                selected.append(model_index)
        model_indices = np.asarray(selected, dtype=np.int64)

    if np.unique(model_indices).size != model_indices.size:
        raise ValueError("Administering-tile selection produced duplicate model tiles")
    model_to_output = np.full(tile_coord.n_tile, -1, dtype=np.int64)
    model_to_output[model_indices] = np.arange(model_indices.size, dtype=np.int64)
    return AdministeringTiles(
        model_to_output=model_to_output,
        model_indices=model_indices,
        tile_id=tile_coord.tile_id[model_indices].astype(np.int32),
        lon=tile_coord.com_lon[model_indices].astype(np.float32),
        lat=tile_coord.com_lat[model_indices].astype(np.float32),
    )


def select_tb_species(
    obs_params: Sequence[ObsParam],
    *,
    description_contains: str,
    orbit: int,
    angles: Sequence[float],
) -> list[TbSpecies]:
    """Select exactly one H and V Tb species for every requested angle."""

    selected: list[TbSpecies] = []
    for pol in (1, 2):
        for angle_index, angle in enumerate(angles):
            matches = [
                param
                for param in obs_params
                if param.varname == "Tb"
                and description_contains in param.descr
                and param.orbit == orbit
                and param.pol == pol
                and np.any(np.isclose(param.ang, angle, atol=1e-6, rtol=0.0))
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one Tb species for orbit={orbit}, pol={pol}, angle={angle}; "
                    f"found {[param.descr for param in matches]}"
                )
            selected.append(
                TbSpecies(matches[0].species, matches[0].descr, pol - 1, angle_index)
            )
    return selected


def _obsfcstana_path(base: Path, exp_run: str, dt: DateTime, file_format: str) -> Path | None:
    stem = (
        base / "ana" / "ens_avg" / f"Y{dt.year:04d}" / f"M{dt.month:02d}"
        / f"{exp_run}.ens_avg.ldas_ObsFcstAna.{dt.year:04d}{dt.month:02d}{dt.day:02d}_"
        f"{dt.hour:02d}{dt.minute:02d}z"
    )
    if file_format == "auto":
        for suffix in (".bin", ".nc4"):
            candidate = Path(f"{stem}{suffix}")
            if candidate.exists():
                return candidate
        return None
    return Path(f"{stem}.{file_format}")


def _date_range_tags(run_months: Sequence[int], start_year: Sequence[int], end_year: Sequence[int]):
    start_year_array = np.asarray(start_year)
    end_year_array = np.asarray(end_year)
    months = np.asarray(run_months)
    first_month = int(months[start_year_array == start_year_array.min()].min())
    last_month = int(months[end_year_array == end_year_array.max()].max())
    cumulative = np.concatenate(([0], np.cumsum([days_in_month(2014, m) for m in range(1, 13)])))
    doy_start = int(cumulative[first_month - 1] + 1)
    doy_end = int(cumulative[last_month])
    pentad_start = int(ceil(doy_start / 5))
    pentad_end = int(floor(doy_end / 5))
    return doy_start, doy_end, pentad_start, pentad_end


def _compute_output(
    ring: np.ndarray,
    ndata_min: int,
    middle_index: int,
) -> np.ndarray:
    data_out = np.full((N_FIELDS, ring.shape[2], ring.shape[3]), NODATA, dtype=np.float32)
    for pol in range(2):
        offset = pol * 5
        count = ring[4, pol].sum(axis=-1)
        valid = count >= ndata_min
        obs_mean = np.full_like(count, np.nan)
        obs_second = np.full_like(count, np.nan)
        mod_mean = np.full_like(count, np.nan)
        mod_second = np.full_like(count, np.nan)
        np.divide(ring[0, pol].sum(axis=-1), count, out=obs_mean, where=count > 0)
        np.divide(ring[1, pol].sum(axis=-1), count, out=obs_second, where=count > 0)
        np.divide(ring[2, pol].sum(axis=-1), count, out=mod_mean, where=count > 0)
        np.divide(ring[3, pol].sum(axis=-1), count, out=mod_second, where=count > 0)
        obs_std = np.sqrt(np.maximum(obs_second - obs_mean**2, 0.0))
        mod_std = np.sqrt(np.maximum(mod_second - mod_mean**2, 0.0))
        values = (obs_mean, obs_std, mod_mean, mod_std, count)
        for field, value in enumerate(values):
            data_out[offset + field] = np.where(valid, value, NODATA)
        center_count = ring[4, pol, :, :, middle_index]
        center_obs = np.full_like(center_count, np.nan)
        center_mod = np.full_like(center_count, np.nan)
        np.divide(
            ring[0, pol, :, :, middle_index], center_count,
            out=center_obs, where=center_count > 0,
        )
        np.divide(
            ring[2, pol, :, :, middle_index], center_count,
            out=center_mod, where=center_count > 0,
        )
        data_out[10 + 2 * pol] = np.where(center_count > 0, center_obs, NODATA)
        data_out[11 + 2 * pol] = np.where(center_count > 0, center_mod, NODATA)
    return data_out


def generate_tb_scaling_params(
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
    description_contains: str,
    orbit: int,
    angles: Sequence[float],
    window_days: int,
    ndata_min: int,
    prefix: str,
    convert_grid: str | None = "EASEv2_M36",
    obsfcstana_format: str = "bin",
    print_each_doy: bool = False,
) -> list[Path]:
    """Generate SMAP Tb tile-space scaling files for one orbit."""

    if orbit not in (1, 2):
        raise ValueError("orbit must be 1 (ascending) or 2 (descending)")
    if window_days % 5 or window_days % 10 == 0:
        raise ValueError("window_days must be an odd number of pentads (for example, 75)")
    if ndata_min < 1:
        raise ValueError("ndata_min must be at least 1")
    if len(run_months) != len(start_year) or len(run_months) != len(end_year):
        raise ValueError("run_months, start_year, and end_year must have equal length")
    if obsfcstana_format not in {"auto", "bin", "nc4"}:
        raise ValueError("obsfcstana_format must be auto, bin, or nc4")

    inpath = Path(exp_path) / exp_run / "output" / domain
    tilecoord_path = inpath / "rc_out" / f"{exp_run}.ldas_tilecoord.bin"
    tilegrids_path = inpath / "rc_out" / f"{exp_run}.ldas_tilegrids.bin"
    tile_coord = read_tilecoord(tilecoord_path)
    _, tile_grid = read_tilegrids(tilegrids_path)
    administering = build_administering_tiles(
        tile_coord, tile_grid, convert_grid=convert_grid
    )
    species = select_tb_species(
        obs_params,
        description_contains=description_contains,
        orbit=orbit,
        angles=angles,
    )
    species_lookup = {item.species: item for item in species}
    print(
        f"Selected {administering.n_tiles:,} administering tiles "
        f"from {tile_coord.n_tile:,} model tiles"
    )
    print("Tb species:", ", ".join(f"{item.species}:{item.description}" for item in species))

    # sum(obs), sum(obs^2), sum(model), sum(model^2), count
    ring = np.zeros(
        (5, 2, administering.n_tiles, len(angles), window_days), dtype=np.float64
    )
    count_days = 0
    t0_assim %= dt_assim
    is_m36_native = "M36" in domain and all(
        "_E" not in item.description for item in species
    )
    tolerance = 1e-3 if is_m36_native else 2.0
    output_dir = inpath / "stats" / "z_score_clim"
    output_dir.mkdir(parents=True, exist_ok=True)
    doy_start, doy_end, pentad_start, pentad_end = _date_range_tags(
        run_months, start_year, end_year
    )
    base = (
        f"{prefix}{min(start_year)}_doy{doy_start}_{max(end_year)}_doy{doy_end}"
        f"_hscale_0.00_W_{window_days}d_Nmin_{ndata_min}_{'A' if orbit == 1 else 'D'}"
    )
    base_p = (
        f"{prefix}{min(start_year)}_p{pentad_start}_{max(end_year)}_p{pentad_end}"
        f"_hscale_0.00_W_{round(window_days / 5)}p_Nmin_{ndata_min}_{'A' if orbit == 1 else 'D'}"
    )
    written: list[Path] = []

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
                    path = _obsfcstana_path(inpath, exp_run, actual_time, obsfcstana_format)
                    if path is None:
                        continue
                    record = read_obs_fcst_ana(path)
                    if record is None:
                        continue
                    valid_fcst = np.isfinite(record.obs_fcst)
                    for species_id, config in species_lookup.items():
                        mask = valid_fcst & (record.obs_species == species_id)
                        if not np.any(mask):
                            continue
                        model_indices = record.obs_tilenum[mask].astype(np.int64) - 1
                        if np.any(model_indices < 0) or np.any(model_indices >= tile_coord.n_tile):
                            raise ValueError(
                                f"Species {species_id} contains out-of-range tile numbers in {path}"
                            )
                        if np.unique(model_indices).size != model_indices.size:
                            raise ValueError(
                                f"Species {species_id} has multiple observations on one tile in {path}"
                            )
                        output_indices = administering.model_to_output[model_indices]
                        if np.any(output_indices < 0):
                            raise ValueError(
                                f"Species {species_id} uses {np.count_nonzero(output_indices < 0)} "
                                f"non-administering tiles in {path}"
                            )
                        obs_lon = record.obs_lon[mask]
                        obs_lat = record.obs_lat[mask]
                        if (
                            np.any(np.abs(tile_coord.com_lon[model_indices] - obs_lon) > tolerance)
                            or np.any(
                                np.abs(tile_coord.com_lat[model_indices] - obs_lat) > tolerance
                            )
                        ):
                            raise ValueError(
                                f"Observation and administering-tile coordinates disagree in {path}"
                            )
                        obs = record.obs_obs[mask].astype(np.float64)
                        model = record.obs_fcst[mask].astype(np.float64)
                        valid_obs = np.isfinite(obs) & (np.abs(obs - NODATA) > 1e-4)
                        pol = config.pol_index
                        angle = config.angle_index
                        np.add.at(
                            ring[0, pol, :, angle, day_index],
                            output_indices[valid_obs], obs[valid_obs],
                        )
                        np.add.at(
                            ring[1, pol, :, angle, day_index],
                            output_indices[valid_obs], obs[valid_obs] ** 2,
                        )
                        np.add.at(ring[2, pol, :, angle, day_index], output_indices, model)
                        np.add.at(ring[3, pol, :, angle, day_index], output_indices, model**2)
                        np.add.at(ring[4, pol, :, angle, day_index], output_indices[valid_obs], 1.0)

            if count_days < window_days:
                continue
            start_time = augment_date_time(-window_days * 86400, last_time)
            center_time = augment_date_time(-floor(window_days * 86400 / 2), last_time)
            doy = int(center_time.dofyr)
            data_out = _compute_output(
                ring, ndata_min, window_days - floor(window_days / 2) - 1
            )
            doy_path = output_dir / f"{base}_DOY{doy:03d}.bin"
            is_middle_of_pentad = (doy + 2) % 5 == 0
            if print_each_doy or is_middle_of_pentad:
                write_tb_scaling_file(
                    doy_path,
                    asc_flag=1 if orbit == 1 else 0,
                    # MATLAB passes zero in this legacy header/version slot.
                    ndata_min=0,
                    start_time=start_time,
                    end_time=last_time,
                    angles=np.asarray(angles),
                    lon=administering.lon,
                    lat=administering.lat,
                    tile_id=administering.tile_id,
                    data=data_out,
                )
                written.append(doy_path)
            if is_middle_of_pentad:
                pentad = (doy + 2) // 5
                pentad_path = output_dir / f"{base_p}_p{pentad:02d}.bin"
                shutil.copyfile(doy_path, pentad_path)
                written.append(pentad_path)
            ring[..., :-1] = ring[..., 1:]
            ring[..., -1] = 0.0

    return written
