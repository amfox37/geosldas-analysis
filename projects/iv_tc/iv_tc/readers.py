"""Readers for sparse observations and GEOSldas model fixtures."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
import importlib
from pathlib import Path
import re
import struct
import sys
import types

import numpy as np

from .pairs import DailyPair, SparseObservation


M36_NX = 964
SMOSIC_DATE_RE = re.compile(r"smos_ic_sm_m36_(\d{8})\.nc$")
SMAPL3_DATE_RE = re.compile(r"SMAP_L3_SM_P_(\d{8})_R\d+_.*\.h5$")
ASCAT_H119_H120_DATE_RE = re.compile(r"ASCAT_HSAF_H119_SM_(\d{8})_AD\.mat$")
SMAP_L3_ORBITS = (
    ("AM", "", ""),
    ("PM", "_dca_pm", "_pm"),
)


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


def read_smap_l3_sparse(path: Path | str, qc: bool = True, nx: int = M36_NX) -> SparseObservation:
    """Read one SMAP L3 SPL3SMP v009/R19240 file as sparse M36 observations."""

    path = Path(path)
    Dataset = _dataset_cls()

    grids: list[np.ndarray] = []
    units = ""
    qc_summary: dict[str, int | float | str] = {
        "qc": int(qc),
        "qc_source": "Save_SMPL3_LDAS_tavg24_nc4_daily.m",
    }

    with Dataset(path, "r") as ds:
        for orbit, retrieval_suffix, ancillary_suffix in SMAP_L3_ORBITS:
            group_name = f"Soil_Moisture_Retrieval_Data_{orbit}"
            if group_name not in ds.groups:
                continue

            group = ds.groups[group_name]
            sm_name = _first_variable_name(
                group,
                f"soil_moisture{retrieval_suffix}",
                f"soil_moisture{ancillary_suffix}",
                "soil_moisture",
            )
            sm_var = group.variables[sm_name]
            sm = _filled(sm_var[:]).astype(np.float64)
            sm[sm < 0] = np.nan
            units = units or _soil_moisture_units(getattr(sm_var, "units", ""))

            raw_finite = int(np.count_nonzero(np.isfinite(sm)))
            if qc:
                qf = _int_variable(group, f"retrieval_qual_flag{retrieval_suffix}", f"retrieval_qual_flag{ancillary_suffix}")
                surface_status = _int_variable(group, f"grid_surface_status{ancillary_suffix}")
                surface_flag = _int_variable(group, f"surface_flag{ancillary_suffix}")
                rfi_h = _int_variable(group, f"tb_qual_flag_h{ancillary_suffix}")
                rfi_v = _int_variable(group, f"tb_qual_flag_v{ancillary_suffix}")

                bad = (
                    ((qf & 4) != 0)
                    | ((surface_flag & 256) != 0)
                    | ((surface_flag & 32) != 0)
                    | ((surface_flag & 64) != 0)
                    | (surface_status != 0)
                    | ((rfi_h & 1) != 0)
                    | ((rfi_v & 1) != 0)
                )
                sm[bad] = np.nan

            kept = int(np.count_nonzero(np.isfinite(sm)))
            qc_summary[f"{orbit.lower()}_raw_finite"] = raw_finite
            qc_summary[f"{orbit.lower()}_kept"] = kept
            grids.append(sm)

    if not grids:
        raise KeyError(f"{path} has no SMAP L3 AM/PM retrieval groups")

    stack = np.stack(grids, axis=0)
    finite = np.isfinite(stack)
    counts = finite.sum(axis=0)
    summed = np.where(finite, stack, 0.0).sum(axis=0)
    daily = np.full(counts.shape, np.nan, dtype=np.float64)
    daily[counts > 0] = summed[counts > 0] / counts[counts > 0]

    row, col = np.nonzero(np.isfinite(daily))
    idx = row.astype(np.int64) * int(nx) + col.astype(np.int64)
    values = daily[row, col].astype(np.float64)
    qc_summary["n_points"] = int(idx.size)

    return SparseObservation(
        date=_parse_smap_l3_date(path),
        sensor="SMAP L3",
        idx=idx,
        values=values,
        units=units,
        qc_summary=qc_summary,
    )


def read_ascat_h121_sparse(
    h121_paths: Path | str | list[Path | str],
    day: date | datetime,
    tilecoord_path: Path | str,
    qc: dict | None = None,
    nx: int = M36_NX,
) -> SparseObservation:
    """Read H121/H139 NetCDF swaths and aggregate QC'd obs to daily M36 cells."""

    read_h121, form_super_obs, qc_default = _ascat_h121_tools()
    day_dt = _as_datetime(day)
    dirs = _h121_day_directories(h121_paths, day_dt.date())
    tilecoord = read_tilecoord(tilecoord_path)

    lat_parts = []
    lon_parts = []
    ssm_parts = []
    n_raw_qc = 0
    qc_used = qc_default if qc is None else qc

    for hdir in dirs:
        obs = read_h121(str(hdir), day_dt, qc=qc_used)
        n_raw_qc += int(obs["ssm"].size)
        if obs["ssm"].size:
            lat_parts.append(obs["lat"])
            lon_parts.append(obs["lon"])
            ssm_parts.append(obs["ssm"])

    if ssm_parts:
        lat = np.concatenate(lat_parts)
        lon = np.concatenate(lon_parts)
        ssm = np.concatenate(ssm_parts)
    else:
        lat = np.array([], dtype=np.float64)
        lon = np.array([], dtype=np.float64)
        ssm = np.array([], dtype=np.float64)

    super_obs = form_super_obs(lat, lon, ssm, tile_coord=tilecoord)
    ij = np.asarray(super_obs["ij"], dtype=np.int64)
    values = np.asarray(super_obs["ssm"], dtype=np.float64)
    if ij.size:
        idx = ij[:, 1] * int(nx) + ij[:, 0]
    else:
        idx = np.array([], dtype=np.int64)

    counts = np.asarray(super_obs["count"], dtype=np.int64)
    qc_summary: dict[str, int | float | str] = {
        "source_dirs": len(dirs),
        "qc_obs": n_raw_qc,
        "n_points": int(idx.size),
        "qc_source": "projects/ascat_da/lib/qc.py:QC_DEFAULT_H121",
        "superob_source": "projects/ascat_da/lib/superob.py:form_super_obs",
    }
    if counts.size:
        qc_summary["superob_count_mean"] = float(counts.mean())
        qc_summary["superob_count_max"] = int(counts.max())

    return SparseObservation(
        date=day_dt.date(),
        sensor="ASCAT H121",
        idx=idx.astype(np.int64),
        values=values,
        units="percent saturation",
        qc_summary=qc_summary,
    )


def read_ascat_h119_h120_sparse(
    mat_path: Path | str,
    auxiliary_path: Path | str,
    tilecoord_path: Path | str,
    method: str = "linear",
    fill_nearest: bool = True,
    nx: int = M36_NX,
) -> SparseObservation:
    """Read one processed ASCAT H119/H120 daily .mat file as sparse M36 obs."""

    mat_path = Path(mat_path)
    target_i, target_j, target_lon, target_lat = _m36_model_land_targets(tilecoord_path)
    lon_gpi, lat_gpi = _read_ascat_h119_h120_land_info(auxiliary_path)
    mat = _load_mat_fields(mat_path, ("sm_tile", "conf_flag_tile"))

    sm = np.asarray(mat["sm_tile"], dtype=np.float64).reshape(-1)
    conf = np.asarray(mat["conf_flag_tile"]).reshape(-1)
    if sm.shape[0] != lon_gpi.shape[0]:
        raise RuntimeError(
            f"ASCAT H119/H120 sm_tile len {sm.shape[0]} != land-grid len {lon_gpi.shape[0]} ({mat_path})"
        )
    if conf.shape[0] != sm.shape[0]:
        raise RuntimeError(f"ASCAT H119/H120 conf_flag_tile len {conf.shape[0]} != sm_tile len {sm.shape[0]}")

    raw_finite = int(np.count_nonzero(np.isfinite(sm)))
    sm[sm > 100] = np.nan
    sm[conf >= 1] = np.nan
    kept = int(np.count_nonzero(np.isfinite(sm)))

    interp = _interpolate_to_targets(
        lon_gpi=lon_gpi,
        lat_gpi=lat_gpi,
        values=sm,
        target_lon=target_lon,
        target_lat=target_lat,
        method=method,
        fill_nearest=fill_nearest,
    )
    good = np.isfinite(interp)
    idx = target_j[good].astype(np.int64) * int(nx) + target_i[good].astype(np.int64)

    return SparseObservation(
        date=_parse_ascat_h119_h120_date(mat_path),
        sensor="ASCAT H119/H120",
        idx=idx.astype(np.int64),
        values=interp[good].astype(np.float64),
        units="percent saturation",
        qc_summary={
            "raw_finite": raw_finite,
            "kept": kept,
            "n_targets": int(target_i.size),
            "n_points": int(idx.size),
            "qc_source": "Save_ASCAT_LDAS_tavg24_nc4_daily.m",
            "interp_method": method,
            "fill_nearest": int(fill_nearest),
        },
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


def read_smap_l3_model_pair(
    smap_path: Path | str,
    model_path: Path | str,
    tilecoord_path: Path | str,
    run: str,
    model_variable: str = "SFMC",
    qc: bool = True,
    nx: int = M36_NX,
) -> DailyPair:
    """Read SMAP L3 and GEOSldas files and return their matched daily pair."""

    obs = read_smap_l3_sparse(smap_path, qc=qc, nx=nx)
    return pair_sparse_observation_with_model(
        observation=obs,
        model_path=model_path,
        tilecoord_path=tilecoord_path,
        run=run,
        model_variable=model_variable,
        nx=nx,
    )


def read_ascat_h121_model_pair(
    h121_paths: Path | str | list[Path | str],
    day: date | datetime,
    model_path: Path | str,
    tilecoord_path: Path | str,
    run: str,
    qc: dict | None = None,
    model_variable: str = "SFMC",
    nx: int = M36_NX,
) -> DailyPair:
    """Read ASCAT H121/H139 and GEOSldas files and return their matched daily pair."""

    obs = read_ascat_h121_sparse(
        h121_paths=h121_paths,
        day=day,
        tilecoord_path=tilecoord_path,
        qc=qc,
        nx=nx,
    )
    return pair_sparse_observation_with_model(
        observation=obs,
        model_path=model_path,
        tilecoord_path=tilecoord_path,
        run=run,
        model_variable=model_variable,
        nx=nx,
    )


def read_ascat_h119_h120_model_pair(
    mat_path: Path | str,
    auxiliary_path: Path | str,
    model_path: Path | str,
    tilecoord_path: Path | str,
    run: str,
    method: str = "linear",
    fill_nearest: bool = True,
    model_variable: str = "SFMC",
    nx: int = M36_NX,
) -> DailyPair:
    """Read ASCAT H119/H120 and GEOSldas files and return their matched pair."""

    obs = read_ascat_h119_h120_sparse(
        mat_path=mat_path,
        auxiliary_path=auxiliary_path,
        tilecoord_path=tilecoord_path,
        method=method,
        fill_nearest=fill_nearest,
        nx=nx,
    )
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


def _parse_smap_l3_date(path: Path) -> datetime.date:
    match = SMAPL3_DATE_RE.search(path.name)
    if match is None:
        raise ValueError(f"Cannot infer SMAP L3 date from {path}")
    return datetime.strptime(match.group(1), "%Y%m%d").date()


def _parse_ascat_h119_h120_date(path: Path) -> datetime.date:
    match = ASCAT_H119_H120_DATE_RE.search(path.name)
    if match is None:
        raise ValueError(f"Cannot infer ASCAT H119/H120 date from {path}")
    return datetime.strptime(match.group(1), "%Y%m%d").date()


def _m36_model_land_targets(tilecoord_path: Path | str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tilecoord = read_tilecoord(tilecoord_path)
    ij = np.column_stack(
        (
            np.asarray(tilecoord["i_indg"], dtype=np.int64),
            np.asarray(tilecoord["j_indg"], dtype=np.int64),
        )
    )
    ij = np.unique(ij, axis=0)
    target_i = ij[:, 0]
    target_j = ij[:, 1]
    target_lon, target_lat = _ease2_m36_lon_lat_for_ij(target_i, target_j)
    return target_i, target_j, target_lon, target_lat


def _read_ascat_h119_h120_land_info(auxiliary_path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    Dataset = _dataset_cls()
    with Dataset(auxiliary_path, "r") as ds:
        land_flag = np.asarray(ds.variables["land_flag"][:]).reshape(-1)
        lon = _filled(ds.variables["lon"][:]).astype(np.float64).reshape(-1)
        lat = _filled(ds.variables["lat"][:]).astype(np.float64).reshape(-1)
    keep = land_flag == 1
    return lon[keep], lat[keep]


def _load_mat_fields(path: Path, fields: tuple[str, ...]) -> dict[str, np.ndarray]:
    try:
        from scipy import io as scipy_io
    except ImportError as exc:
        raise ImportError("scipy is required for ASCAT H119/H120 .mat readers") from exc

    try:
        data = scipy_io.loadmat(path, squeeze_me=True, struct_as_record=False)
        return {field: np.asarray(data[field]) for field in fields}
    except NotImplementedError:
        try:
            import h5py
        except ImportError as exc:
            raise ImportError("h5py is required for MATLAB v7.3 ASCAT H119/H120 .mat files") from exc

        out: dict[str, np.ndarray] = {}
        with h5py.File(path, "r") as h5:
            for field in fields:
                out[field] = np.asarray(h5[field]).squeeze()
        return out


def _interpolate_to_targets(
    lon_gpi: np.ndarray,
    lat_gpi: np.ndarray,
    values: np.ndarray,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
    method: str,
    fill_nearest: bool,
) -> np.ndarray:
    try:
        from scipy.interpolate import griddata
    except ImportError as exc:
        raise ImportError("scipy is required for ASCAT H119/H120 interpolation") from exc

    method = method.lower()
    min_points = 3 if method in {"linear", "cubic"} else 1
    finite = np.isfinite(values) & np.isfinite(lon_gpi) & np.isfinite(lat_gpi)
    if int(np.count_nonzero(finite)) < min_points:
        return np.full(target_lon.shape, np.nan, dtype=np.float64)

    points = np.column_stack((lon_gpi[finite], lat_gpi[finite]))
    out = np.asarray(
        griddata(points, values[finite], (target_lon, target_lat), method=method),
        dtype=np.float64,
    )
    if method == "linear" and fill_nearest:
        missing = ~np.isfinite(out)
        if np.any(missing):
            out[missing] = griddata(points, values[finite], (target_lon[missing], target_lat[missing]), method="nearest")
    return out


def _ease2_m36_lon_lat_for_ij(i: np.ndarray, j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return EASE2 M36 longitude/latitude for zero-based M36 cell indices."""

    i = np.asarray(i, dtype=np.float64)
    j = np.asarray(j, dtype=np.float64)
    nlon = 964
    nlat = 406
    r0 = (nlon - 1) / 2.0
    s0 = (nlat - 1) / 2.0
    map_scale_m = 36032.220840584

    x = (i - r0) * map_scale_m
    y = (s0 - j) * map_scale_m

    a = 6378137.0
    e = 0.081819190843
    e2 = e**2
    e4 = e**4
    e6 = e**6
    phi1 = np.deg2rad(30.0)
    kz = np.cos(phi1) / np.sqrt(1.0 - e2 * np.sin(phi1) ** 2)
    qp = (1.0 - e2) * ((1.0 / (1.0 - e2)) - (1.0 / (2.0 * e)) * np.log((1.0 - e) / (1.0 + e)))

    beta = np.arcsin(2.0 * y * kz / (a * qp))
    lam = x / (a * kz)
    phi = (
        beta
        + (((e2 / 3.0) + ((31.0 / 180.0) * e4) + ((517.0 / 5040.0) * e6)) * np.sin(2.0 * beta))
        + ((((23.0 / 360.0) * e4) + ((251.0 / 3780.0) * e6)) * np.sin(4.0 * beta))
        + (((761.0 / 45360.0) * e6) * np.sin(6.0 * beta))
    )

    lat = np.rad2deg(phi)
    lon = np.rad2deg(lam)
    lon = np.where(lon < -180.0, lon + 360.0, lon)
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    return lon.astype(np.float64), lat.astype(np.float64)


def _filled(arr: np.ndarray | np.ma.MaskedArray) -> np.ndarray:
    if np.ma.isMaskedArray(arr):
        return np.asarray(arr.filled(np.nan))
    return np.asarray(arr)


def _first_variable_name(group, *names: str) -> str:
    for name in names:
        if name and name in group.variables:
            return name
    available = ", ".join(sorted(group.variables))
    raise KeyError(f"None of {names!r} found in group {group.path}; available: {available}")


def _int_variable(group, *names: str) -> np.ndarray:
    name = _first_variable_name(group, *names)
    var = group.variables[name]
    fill = getattr(var, "_FillValue", 0)
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        arr = arr.filled(fill)
    return np.asarray(arr, dtype=np.int64)


def _soil_moisture_units(units: str) -> str:
    if units in {"cm**3/cm**3", "m3/m3", "m3 m-3"}:
        return "m3 m-3"
    return units


def _as_datetime(value: date | datetime) -> datetime:
    if isinstance(value, datetime):
        return value
    return datetime(value.year, value.month, value.day)


def _as_path_list(paths: Path | str | list[Path | str]) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(path) for path in paths]


def _h121_day_directories(paths: Path | str | list[Path | str], day: date) -> list[Path]:
    yyyymmdd = day.strftime("%Y%m%d")
    dirs: list[Path] = []
    for path in _as_path_list(paths):
        if path.is_file():
            dirs.append(path.parent)
        elif any(path.glob(f"*_{yyyymmdd}*.nc")):
            dirs.append(path)
        else:
            dirs.extend(sorted({candidate.parent for candidate in path.glob(f"**/*_{yyyymmdd}*.nc")}))

    unique: dict[str, Path] = {}
    for path in dirs:
        unique[str(path.resolve())] = path
    return sorted(unique.values())


def _ascat_h121_tools():
    package = "_iv_tc_ascat_da_lib"
    repo_root = Path(__file__).resolve().parents[3]
    lib_dir = repo_root / "projects" / "ascat_da" / "lib"
    if not lib_dir.exists():
        raise ImportError(f"ASCAT DA support library not found at {lib_dir}")

    if package not in sys.modules:
        module = types.ModuleType(package)
        module.__path__ = [str(lib_dir)]
        sys.modules[package] = module

    qc_mod = importlib.import_module(f"{package}.qc")
    readers_mod = importlib.import_module(f"{package}.readers")
    superob_mod = importlib.import_module(f"{package}.superob")
    return readers_mod.read_h121, superob_mod.form_super_obs, qc_mod.QC_DEFAULT_H121


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
