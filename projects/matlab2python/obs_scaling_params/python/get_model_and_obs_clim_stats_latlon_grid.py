"""
Python port of `get_model_and_obs_clim_stats_latlon_grid.m`.

This follows the MATLAB implementation closely so the output NetCDF files can
be used interchangeably for observation scaling.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
from typing import Sequence

import numpy as np

try:
    from .obs_scaling_utils import (
        DateTime,
        augment_date_time,
        read_obs_fcst_ana,
        write_netcdf_latlon_grid,
    )
except ImportError:
    import sys
    from pathlib import Path

    THIS_DIR = Path(__file__).resolve().parent
    if str(THIS_DIR) not in sys.path:
        sys.path.append(str(THIS_DIR))
    from obs_scaling_utils import (  # type: ignore
        DateTime,
        augment_date_time,
        read_obs_fcst_ana,
        write_netcdf_latlon_grid,
    )


@dataclass
class GridDefinition:
    ll_lon: float
    ll_lat: float
    d_lon: float
    d_lat: float
    n_lon: int
    n_lat: int
    ll_lons: np.ndarray
    ll_lats: np.ndarray
    obsnum: np.ndarray
    i_out: np.ndarray
    j_out: np.ndarray


@dataclass
class _DedupEntry:
    scnt: int
    cell: int
    idx: int
    obs: float
    obs2: float
    mod: float
    mod2: float
    count: float
    cycle_id: int


def _build_grid(resol: float) -> GridDefinition:
    ll_lon = -180.0
    ll_lat = -90.0
    d_lon = resol
    d_lat = resol
    n_lon = int(round(360.0 / d_lon))
    n_lat = int(round(180.0 / d_lat))
    ll_lons = np.linspace(ll_lon, ll_lon + (n_lon - 1) * d_lon, n_lon)
    ll_lats = np.linspace(ll_lat, ll_lat + (n_lat - 1) * d_lat, n_lat)
    obsnum = np.arange(n_lon * n_lat)
    j_out, i_out = np.divmod(obsnum, n_lon)
    return GridDefinition(
        ll_lon=ll_lon,
        ll_lat=ll_lat,
        d_lon=d_lon,
        d_lat=d_lat,
        n_lon=n_lon,
        n_lat=n_lat,
        ll_lons=ll_lons,
        ll_lats=ll_lats,
        obsnum=obsnum,
        i_out=i_out,
        j_out=j_out,
    )


def get_model_and_obs_clim_stats_latlon_grid(
    species_names: Sequence[str],
    run_months: Sequence[int],
    exp_path: str,
    exp_run: str,
    domain: str,
    start_year: Sequence[int],
    end_year: Sequence[int],
    dt_assim: int,
    t0_assim: int,
    species: Sequence[int],
    combine_species_stats: bool,
    resol: float,
    w_days: int,
    ndata_min: int,
    prefix: str,
    print_each_DOY: bool,
    print_each_pentad: bool,
    print_all_pentads: bool,
    out_dir: str,
    enable_dedup: bool = True,
) -> None:
    nodata = -9999.0
    overwrite = True
    n_fields = 7
    n_pentads = 73
    lonlat_q = 1e3
    obs_q = 1e4
    dedup_keep_cycles = 3

    if combine_species_stats:
        species_groups = [np.array(species, dtype=int)]
    else:
        species_groups = [np.array([s], dtype=int) for s in species]

    inpath = Path(exp_path) / exp_run / "output" / domain
    outpath = inpath / "stats" / out_dir
    outpath.mkdir(parents=True, exist_ok=True)

    start_year = np.array(start_year, dtype=int)
    end_year = np.array(end_year, dtype=int)
    run_months = np.array(run_months, dtype=int)

    grid = _build_grid(resol)
    n_species = len(species_groups)
    n_gridcells = grid.n_lon * grid.n_lat

    shape = (n_species, n_gridcells, w_days)
    o_data_sum = np.full(shape, np.nan)
    m_data_sum = np.full(shape, np.nan)
    o_data_sum2 = np.full(shape, np.nan)
    m_data_sum2 = np.full(shape, np.nan)
    m_data_min = np.full(shape, np.nan)
    m_data_max = np.full(shape, np.nan)
    n_data = np.full(shape, np.nan)

    data_out = np.full((n_species, n_fields, n_gridcells, n_pentads), np.nan)
    start_time_p = [DateTime(2014, 1, 1, 0, 0, 0) for _ in range(n_pentads)]
    end_time_p = [DateTime(2014, 1, 1, 0, 0, 0) for _ in range(n_pentads)]

    data2d = np.full((n_fields, n_gridcells), np.nan)

    count = 0

    start_mask = start_year == start_year.min()
    end_mask = end_year == end_year.max()
    mi_m = int(run_months[start_mask].min())
    ma_m = int(run_months[end_mask].max())

    days_non_leap = np.array([_days_in_month_non_leap(m) for m in range(1, 13)])
    cumulative_days = np.concatenate(([0], np.cumsum(days_non_leap)))
    if mi_m > 1:
        d1 = int(cumulative_days[mi_m - 1] + 1)
        p1 = int(ceil(d1 / 5))
    else:
        d1 = 1
        p1 = 1
    d2 = int(cumulative_days[ma_m])
    p2 = int(floor(d2 / 5))

    fname_out_base_d = (
        f"{outpath}/{prefix}{start_year.min()}_doy{d1}_{end_year.max()}_doy{d2}"
        f"_W_{w_days}d_Nmin_{ndata_min}"
    )
    fname_out_base_p = (
        f"{outpath}/{prefix}{start_year.min()}_p{p1}_{end_year.max()}_p{p2}"
        f"_W_{int(round(w_days / 5))}p_Nmin_{ndata_min}"
    )

    cycle_id = 0
    dedup_tracker: dict[tuple[int, ...], _DedupEntry] = {}

    for imonth, month in enumerate(run_months.tolist()):
        month_days = _days_in_month_non_leap(month)
        for day in range(1, month_days + 1):
            print(f"Processing {month:02d}-{day:02d}")
            if count < w_days:
                count += 1
            idx = min(count, w_days) - 1

            for seconds_in_day in range(t0_assim, 86400, dt_assim):
                hour = seconds_in_day // 3600
                minute = (seconds_in_day % 3600) // 60
                second = seconds_in_day % 60
                if enable_dedup:
                    cycle_id += 1
                    _prune_dedup(dedup_tracker, cycle_id, dedup_keep_cycles)
                for year in range(start_year[imonth], end_year[imonth] + 1):
                    fname = (
                        inpath
                        / "ana"
                        / "ens_avg"
                        / f"Y{year:04d}"
                        / f"M{month:02d}"
                        / f"{exp_run}.ens_avg.ldas_ObsFcstAna.{year:04d}{month:02d}{day:02d}_{hour:02d}{minute:02d}z.bin"
                    )
                    record = read_obs_fcst_ana(fname)
                    if record is None:
                        continue
                    valid = ~np.isnan(record.obs_fcst)
                    if not np.any(valid):
                        continue

                    for scnt, species_group in enumerate(species_groups):
                        if combine_species_stats:
                            mask = np.isin(record.obs_species, species_group)
                        else:
                            mask = record.obs_species == species_group[0]
                        mask &= valid
                        if not np.any(mask):
                            continue
                        obs_lon = record.obs_lon[mask]
                        obs_lat = record.obs_lat[mask]
                        obs_obs = record.obs_obs[mask]
                        obs_fcst = record.obs_fcst[mask]

                        i_idx = np.floor((obs_lon - grid.ll_lon) / grid.d_lon).astype(int)
                        j_idx = np.floor((obs_lat - grid.ll_lat) / grid.d_lat).astype(int)
                        inside = (
                            (i_idx >= 0)
                            & (i_idx < grid.n_lon)
                            & (j_idx >= 0)
                            & (j_idx < grid.n_lat)
                        )
                        if not np.any(inside):
                            continue
                        i_idx = i_idx[inside]
                        j_idx = j_idx[inside]
                        obs_obs = obs_obs[inside]
                        obs_fcst = obs_fcst[inside]
                        obs_cells = j_idx * grid.n_lon + i_idx

                        obs_species = record.obs_species[mask][inside]
                        obs_tile = record.obs_tilenum[mask][inside]
                        obs_lon_sel = obs_lon[inside]
                        obs_lat_sel = obs_lat[inside]

                        for cell, obs_val, fcst_val, tile, species_id, lon, lat in zip(
                            obs_cells, obs_obs, obs_fcst, obs_tile, obs_species, obs_lon_sel, obs_lat_sel
                        ):
                            obs_contrib = np.nan if np.isnan(obs_val) else float(obs_val)
                            obs_sq = np.nan if np.isnan(obs_val) else obs_contrib * obs_contrib
                            mod_contrib = np.nan if np.isnan(fcst_val) else float(fcst_val)
                            mod_sq = np.nan if np.isnan(fcst_val) else mod_contrib * mod_contrib

                            key = None
                            obs_valid = not np.isnan(obs_contrib)
                            mod_valid = not np.isnan(mod_contrib)
                            obs_delta = obs_contrib if obs_valid else 0.0
                            obs2_delta = obs_sq if obs_valid else 0.0
                            mod_delta = mod_contrib if mod_valid else 0.0
                            mod2_delta = mod_sq if mod_valid else 0.0
                            n_delta = 1.0 if obs_valid else 0.0

                            if enable_dedup and obs_valid:
                                key = (
                                    int(tile),
                                    int(species_id),
                                    int(round(lon * lonlat_q)),
                                    int(round(lat * lonlat_q)),
                                    int(round(obs_contrib * obs_q)),
                                    int(year),
                                )
                                prev = dedup_tracker.get(key)
                                if prev is not None:
                                    _accumulate(o_data_sum, prev.scnt, prev.cell, prev.idx, -prev.obs)
                                    _accumulate(m_data_sum, prev.scnt, prev.cell, prev.idx, -prev.mod)
                                    _accumulate(o_data_sum2, prev.scnt, prev.cell, prev.idx, -prev.obs2)
                                    _accumulate(m_data_sum2, prev.scnt, prev.cell, prev.idx, -prev.mod2)
                                    _accumulate(n_data, prev.scnt, prev.cell, prev.idx, -prev.count)

                            if obs_valid:
                                _accumulate(o_data_sum, scnt, cell, idx, obs_delta)
                                _accumulate(o_data_sum2, scnt, cell, idx, obs2_delta)
                                _accumulate(n_data, scnt, cell, idx, n_delta)
                            if mod_valid:
                                _accumulate(m_data_sum, scnt, cell, idx, mod_delta)
                                _accumulate(m_data_sum2, scnt, cell, idx, mod2_delta)
                                m_data_min[scnt, cell, idx] = _nan_min(m_data_min[scnt, cell, idx], mod_contrib)
                                m_data_max[scnt, cell, idx] = _nan_max(m_data_max[scnt, cell, idx], mod_contrib)

                            if key is not None:
                                dedup_tracker[key] = _DedupEntry(
                                    scnt=scnt,
                                    cell=cell,
                                    idx=idx,
                                    obs=obs_delta,
                                    obs2=obs2_delta,
                                    mod=mod_delta,
                                    mod2=mod2_delta,
                                    count=n_delta,
                                    cycle_id=cycle_id,
                                )

            if count >= w_days:
                end_time = DateTime(2014, month, day, hour, minute, second)
                start_time = augment_date_time(-int(w_days * 86400), end_time)
                for scnt in range(n_species):
                    window_n = np.nansum(n_data[scnt, :, :], axis=1)
                    obs_sum = np.nansum(o_data_sum[scnt, :, :], axis=1)
                    obs_sum2 = np.nansum(o_data_sum2[scnt, :, :], axis=1)
                    mod_sum = np.nansum(m_data_sum[scnt, :, :], axis=1)
                    mod_sum2 = np.nansum(m_data_sum2[scnt, :, :], axis=1)

                    obs_mean = np.divide(obs_sum, window_n, out=np.full_like(obs_sum, np.nan), where=window_n > 0)
                    obs_std = np.sqrt(
                        np.maximum(
                            np.divide(obs_sum2, window_n, out=np.full_like(obs_sum2, np.nan), where=window_n > 0)
                            - obs_mean**2,
                            0.0,
                        )
                    )
                    mod_mean = np.divide(mod_sum, window_n, out=np.full_like(mod_sum, np.nan), where=window_n > 0)
                    mod_std = np.sqrt(
                        np.maximum(
                            np.divide(mod_sum2, window_n, out=np.full_like(mod_sum2, np.nan), where=window_n > 0)
                            - mod_mean**2,
                            0.0,
                        )
                    )
                    data2d[0, :] = obs_mean
                    data2d[1, :] = obs_std
                    data2d[2, :] = mod_mean
                    data2d[3, :] = mod_std
                    data2d[4, :] = window_n
                    with np.errstate(all="ignore"):
                        data2d[5, :] = np.nanmin(m_data_min[scnt, :, :], axis=1)
                        data2d[6, :] = np.nanmax(m_data_max[scnt, :, :], axis=1)

                    invalid = window_n < ndata_min
                    data2d[:, invalid] = np.nan

                    doy = augment_date_time(-int(w_days * 86400 / 2.0), end_time).dofyr
                    pentad = floor((doy + 2) / 5)

                    if print_each_DOY:
                        species_label = "ALL" if combine_species_stats else species_names[scnt]
                        fname_out = f"{fname_out_base_d}_sp_{species_label}_DOY{doy:03d}.nc4"
                        write_netcdf_latlon_grid(
                            fname_out,
                            grid.i_out,
                            grid.j_out,
                            grid.ll_lons,
                            grid.ll_lats,
                            data2d,
                            pentad,
                            [start_time],
                            [end_time],
                            overwrite,
                            n_fields,
                            grid.ll_lon,
                            grid.ll_lat,
                            grid.d_lon,
                            grid.d_lat,
                        )

                    if (doy + 2) % 5 == 0:
                        data_out[scnt, :, :, pentad - 1] = data2d
                        start_time_p[pentad - 1] = start_time
                        end_time_p[pentad - 1] = end_time
                        if print_each_pentad:
                            species_label = "ALL" if combine_species_stats else species_names[scnt]
                            fname_out = f"{fname_out_base_p}_sp_{species_label}_p{pentad:02d}.nc4"
                            write_netcdf_latlon_grid(
                                fname_out,
                                grid.i_out,
                                grid.j_out,
                                grid.ll_lons,
                                grid.ll_lats,
                                data2d,
                                pentad,
                                [start_time],
                                [end_time],
                                overwrite,
                                n_fields,
                                grid.ll_lon,
                                grid.ll_lat,
                                grid.d_lon,
                                grid.d_lat,
                            )

                o_data_sum[:, :, :-1] = o_data_sum[:, :, 1:]
                m_data_sum[:, :, :-1] = m_data_sum[:, :, 1:]
                o_data_sum2[:, :, :-1] = o_data_sum2[:, :, 1:]
                m_data_sum2[:, :, :-1] = m_data_sum2[:, :, 1:]
                m_data_min[:, :, :-1] = m_data_min[:, :, 1:]
                m_data_max[:, :, :-1] = m_data_max[:, :, 1:]
                n_data[:, :, :-1] = n_data[:, :, 1:]
                o_data_sum[:, :, -1] = np.nan
                m_data_sum[:, :, -1] = np.nan
                o_data_sum2[:, :, -1] = np.nan
                m_data_sum2[:, :, -1] = np.nan
                m_data_min[:, :, -1] = np.nan
                m_data_max[:, :, -1] = np.nan
                n_data[:, :, -1] = np.nan
                data2d[:] = np.nan

    if print_all_pentads:
        pentads = np.arange(1, n_pentads + 1)
        for scnt in range(n_species):
            species_label = "ALL" if combine_species_stats else species_names[scnt]
            fname_out = f"{fname_out_base_d}_sp_{species_label}_all_pentads.nc4"
            write_netcdf_latlon_grid(
                fname_out,
                grid.i_out,
                grid.j_out,
                grid.ll_lons,
                grid.ll_lats,
                data_out[scnt],
                pentads,
                start_time_p,
                end_time_p,
                overwrite,
                n_fields,
                grid.ll_lon,
                grid.ll_lat,
                grid.d_lon,
                grid.d_lat,
            )


def _day_of_year_bounds(month: int, start: bool) -> int:
    fake_year = 2014
    if start:
        start_day = 1
    else:
        start_day = _days_in_month_non_leap(month)
    date = DateTime(fake_year, month, start_day, 0, 0, 0)
    return augment_date_time(0, date).dofyr or 1


def _pentad_bounds(month: int, start: bool) -> int:
    doy = _day_of_year_bounds(month, start)
    return floor((doy + 2) / 5)


def _days_in_month_non_leap(month: int) -> int:
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return days[month - 1]


def _accumulate(arr: np.ndarray, scnt: int, cell: int, idx: int, delta: float) -> None:
    if delta == 0:
        return
    current = arr[scnt, cell, idx]
    if np.isnan(current):
        current = 0.0
    arr[scnt, cell, idx] = current + delta


def _prune_dedup(dedup_tracker: dict[tuple[int, ...], _DedupEntry], cycle_id: int, keep_cycles: int) -> None:
    if not dedup_tracker:
        return
    to_remove = [key for key, entry in dedup_tracker.items() if cycle_id - entry.cycle_id > keep_cycles]
    for key in to_remove:
        dedup_tracker.pop(key, None)


def _nan_min(value: float, candidate: float) -> float:
    if np.isnan(candidate):
        return value
    if np.isnan(value):
        return candidate
    return min(value, candidate)


def _nan_max(value: float, candidate: float) -> float:
    if np.isnan(candidate):
        return candidate
    if np.isnan(value):
        return candidate
    return max(value, candidate)
