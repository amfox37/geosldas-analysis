#!/usr/bin/env python3
"""
Preprocess SMOS-IC daily soil moisture to GEOSldas M36 land cells.

Design:
1) Read daily SMOS-IC ASC/DES files.
2) Apply QC on native 25 km grid:
   - Scene_Flags <= threshold
   - TB-RMSE <= threshold
   - Soil_Moisture in [sm_min, sm_max]
3) Merge ASC+DES per day.
4) Remap to M36 using precomputed conservative area-overlap weights
   between rectilinear EASE2 cells (source) and M36 EASE2 cells (target).
5) Write sparse daily outputs (idx_EASEv2_lonxlat, sm_obs) for reuse.

Outputs:
- mapping cache (.npz)
- daily sparse .nc files
- manifest CSV
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


DATE_RE = re.compile(r"_(20\d{6})T")


@dataclass(frozen=True)
class M36Grid:
    nx: int = 964
    ny: int = 406
    dx_m: float = 36032.220840584

    @property
    def r0(self) -> float:
        return (self.nx - 1) / 2.0

    @property
    def s0(self) -> float:
        return (self.ny - 1) / 2.0

    def x_edges(self) -> np.ndarray:
        # x_center(i) = (i - r0) * dx
        # x_edge(i)   = (i - 0.5 - r0) * dx
        return (np.arange(self.nx + 1, dtype=np.float64) - 0.5 - self.r0) * self.dx_m

    def y_edges_desc(self) -> np.ndarray:
        # y_center(j) = (s0 - j) * dx
        # y_edge(j)   = (s0 - (j - 0.5)) * dx
        return (self.s0 - (np.arange(self.ny + 1, dtype=np.float64) - 0.5)) * self.dx_m


@dataclass(frozen=True)
class SourceLayout:
    raw_shape: tuple[int, int]
    nlat: int
    nlon: int
    lat_axis: int
    lon_axis: int
    transpose: bool
    flip_lat: bool
    flip_lon: bool


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Preprocess SMOS-IC daily files to M36 sparse cache.")
    ap.add_argument("--smos-root", required=True, help="Root directory containing SMOS-IC folders")
    ap.add_argument("--tilecoord", required=True, help="GEOSldas tilecoord binary file")
    ap.add_argument("--out-dir", required=True, help="Output directory for cache + daily files")
    ap.add_argument("--start-date", default=None, help="Optional YYYY-MM-DD")
    ap.add_argument("--end-date", default=None, help="Optional YYYY-MM-DD")
    ap.add_argument("--scene-flag-max", type=int, default=1, help="QC threshold: Scene_Flags <= value")
    ap.add_argument("--tb-rmse-max", type=float, default=8.0, help="QC threshold: RMSE <= value (K)")
    ap.add_argument("--sm-min", type=float, default=0.0, help="QC threshold: Soil_Moisture >= value")
    ap.add_argument("--sm-max", type=float, default=1.0, help="QC threshold: Soil_Moisture <= value")
    ap.add_argument(
        "--min-coverage-frac",
        type=float,
        default=0.05,
        help="Minimum target-cell area coverage fraction required to keep remapped value",
    )
    ap.add_argument(
        "--mapping-cache",
        default=None,
        help="Optional mapping cache .npz path (default: <out-dir>/smos25km_to_m36_conservative_mapping.npz)",
    )
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing daily outputs")
    ap.add_argument("--rebuild-mapping", action="store_true", help="Recompute mapping even if cache exists")
    return ap.parse_args()


def find_repo_root(start: Path) -> Path:
    p = start.resolve()
    for c in [p] + list(p.parents):
        if (c / "common/python/io/read_GEOSldas.py").exists():
            return c
    raise FileNotFoundError("Could not find repo root containing common/python/io/read_GEOSldas.py")


def parse_date_from_name(name: str) -> date | None:
    m = DATE_RE.search(name)
    if not m:
        return None
    return datetime.strptime(m.group(1), "%Y%m%d").date()


def build_file_index(root: Path) -> tuple[dict[date, Path], dict[date, Path]]:
    """Build ASC/DES daily file indexes from expected SMOS-IC folders."""
    asc_dirs = [
        root / "SMOS_IC_V2_ASC_2015_2017",
        root / "SMOS_IC_V2_ASC_2018_2021",
        root / "SMOS_IC_V2_ASC",
    ]
    des_dirs = [
        root / "SMOS_IC_V2_DES_2015_2017",
        root / "SMOS_IC_V2_DES_2018_2021",
        root / "SMOS_IC_OPER_DES",
    ]

    asc: dict[date, Path] = {}
    des: dict[date, Path] = {}

    for d in asc_dirs:
        if not d.exists():
            continue
        for f in sorted(d.glob("*.nc")):
            dt = parse_date_from_name(f.name)
            if dt is None:
                continue
            asc[dt] = f

    for d in des_dirs:
        if not d.exists():
            continue
        for f in sorted(d.glob("*.nc")):
            dt = parse_date_from_name(f.name)
            if dt is None:
                continue
            des[dt] = f

    if not asc and not des:
        raise FileNotFoundError(f"No SMOS-IC .nc files found under {root}")

    return asc, des


def choose_representative_tile_per_cell(tc: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Deterministically choose one tile per M36 cell.

    Rule:
    - highest frac_cell
    - tie-break: smallest tile_id
    """
    i_indg = tc["i_indg"].astype(np.int64)
    j_indg = tc["j_indg"].astype(np.int64)
    tile_id = tc["tile_id"].astype(np.int64)
    frac_cell = np.asarray(tc["frac_cell"], dtype=np.float64)

    frac_sort = np.where(np.isfinite(frac_cell), frac_cell, -np.inf)

    # Use tilecoord index space only to identify unique cells. Do not infer
    # target-grid geometry from these maxima (some domains are partial subsets).
    nx_local = int(i_indg.max()) + 1
    ny_local = int(j_indg.max()) + 1
    cell_code = j_indg * nx_local + i_indg

    order = np.lexsort((tile_id, -frac_sort, cell_code))
    code_sorted = cell_code[order]

    first = np.empty(code_sorted.size, dtype=bool)
    first[0] = True
    first[1:] = code_sorted[1:] != code_sorted[:-1]
    rep_idx = order[first]

    rep_i = i_indg[rep_idx].astype(np.int32)
    rep_j = j_indg[rep_idx].astype(np.int32)

    return {
        "rep_idx": rep_idx,
        "rep_i": rep_i,
        "rep_j": rep_j,
        "rep_tile_id": tile_id[rep_idx].astype(np.int32),
        "rep_frac_cell": frac_cell[rep_idx].astype(np.float32),
        "rep_lat": np.asarray(tc["com_lat"][rep_idx], dtype=np.float64),
        "rep_lon": np.asarray(tc["com_lon"][rep_idx], dtype=np.float64),
        "nx_local": np.int32(nx_local),
        "ny_local": np.int32(ny_local),
    }


def read_source_georef(sample_nc: Path) -> tuple[int, int, np.ndarray]:
    """Read source grid shape and GeoTransform from one SMOS-IC file."""
    with Dataset(sample_nc, "r") as ds:
        if "lon" not in ds.dimensions or "lat" not in ds.dimensions:
            raise KeyError(f"{sample_nc} missing lon/lat dimensions")
        nlon = int(len(ds.dimensions["lon"]))
        nlat = int(len(ds.dimensions["lat"]))
        if "crs" not in ds.variables or not hasattr(ds.variables["crs"], "GeoTransform"):
            raise KeyError(f"{sample_nc} missing crs GeoTransform")
        gt = np.asarray([float(v) for v in str(ds.variables["crs"].GeoTransform).split()], dtype=np.float64)
    return nlat, nlon, gt


def _interval_overlaps(edges_asc: np.ndarray, low: float, high: float) -> tuple[np.ndarray, np.ndarray]:
    """Return overlapping source interval indices and overlap lengths."""
    if not (high > low):
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.float64)

    n = edges_asc.size - 1
    i0 = int(np.searchsorted(edges_asc, low, side="right") - 1)
    i1 = int(np.searchsorted(edges_asc, high, side="left"))
    i0 = max(i0, 0)
    i1 = min(i1, n)
    if i1 <= i0:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.float64)

    idx = np.arange(i0, i1, dtype=np.int32)
    lo = np.maximum(low, edges_asc[idx])
    hi = np.minimum(high, edges_asc[idx + 1])
    ol = hi - lo
    ok = ol > 0.0
    return idx[ok], ol[ok]


def _precompute_1d_overlap_maps(
    src_x_edges: np.ndarray,
    src_y_edges_desc: np.ndarray,
    tgt: M36Grid,
) -> tuple[list[tuple[np.ndarray, np.ndarray]], list[tuple[np.ndarray, np.ndarray]]]:
    """
    Precompute source-overlap maps for each M36 column and row.

    Returns:
    - x_map[i] -> (src_cols, overlap_x_m)
    - y_map[j] -> (src_rows_orig, overlap_y_m)
    """
    x_map: list[tuple[np.ndarray, np.ndarray]] = []
    y_map: list[tuple[np.ndarray, np.ndarray]] = []

    tgt_x_edges = tgt.x_edges()  # ascending
    tgt_y_edges_desc = tgt.y_edges_desc()  # descending

    src_y_edges_asc = src_y_edges_desc[::-1]
    src_ny = src_y_edges_desc.size - 1

    for i in range(tgt.nx):
        low = min(tgt_x_edges[i], tgt_x_edges[i + 1])
        high = max(tgt_x_edges[i], tgt_x_edges[i + 1])
        src_idx, ol = _interval_overlaps(src_x_edges, low, high)
        x_map.append((src_idx, ol))

    for j in range(tgt.ny):
        low = min(tgt_y_edges_desc[j], tgt_y_edges_desc[j + 1])
        high = max(tgt_y_edges_desc[j], tgt_y_edges_desc[j + 1])
        src_idx_asc, ol = _interval_overlaps(src_y_edges_asc, low, high)
        # Convert asc-row indexing back to native source row indexing.
        src_rows_orig = (src_ny - 1 - src_idx_asc).astype(np.int32)
        y_map.append((src_rows_orig, ol))

    return x_map, y_map


def build_conservative_mapping(
    sample_nc: Path,
    rep: dict[str, np.ndarray],
    tgt: M36Grid,
    cache_path: Path,
) -> dict[str, np.ndarray]:
    """Build source->target conservative overlap mapping and write cache."""
    with Dataset(sample_nc, "r") as ds:
        if "lon" not in ds.variables or "lat" not in ds.variables:
            raise KeyError(f"{sample_nc} missing lon/lat variables")
        nlon = len(ds.dimensions["lon"])
        nlat = len(ds.dimensions["lat"])
        if "crs" not in ds.variables or not hasattr(ds.variables["crs"], "GeoTransform"):
            raise KeyError(f"{sample_nc} missing crs GeoTransform")
        gt = [float(v) for v in str(ds.variables["crs"].GeoTransform).split()]

    # GDAL GeoTransform: x0 dx rx y0 ry dy, where x0/y0 are top-left corner edges.
    x0, dx, rx, y0, ry, dy = gt
    if abs(rx) > 1e-12 or abs(ry) > 1e-12:
        raise RuntimeError("SMOS-IC GeoTransform includes rotation terms; expected 0")

    src_x_edges = x0 + np.arange(nlon + 1, dtype=np.float64) * dx
    src_y_edges_desc = y0 + np.arange(nlat + 1, dtype=np.float64) * dy

    if not np.all(np.diff(src_x_edges) > 0):
        raise RuntimeError("Source x edges are not strictly increasing")
    # dy is expected negative for north->south rows.
    if not np.all(np.diff(src_y_edges_desc) < 0):
        raise RuntimeError("Source y edges are not strictly decreasing")

    x_map, y_map = _precompute_1d_overlap_maps(src_x_edges, src_y_edges_desc, tgt)

    rep_i = rep["rep_i"].astype(np.int32)
    rep_j = rep["rep_j"].astype(np.int32)
    if np.any(rep_i < 0) or np.any(rep_i >= tgt.nx) or np.any(rep_j < 0) or np.any(rep_j >= tgt.ny):
        raise RuntimeError(
            "Representative tilecoord indices are outside expected M36 global bounds: "
            f"i=[{int(rep_i.min())},{int(rep_i.max())}] vs [0,{tgt.nx-1}], "
            f"j=[{int(rep_j.min())},{int(rep_j.max())}] vs [0,{tgt.ny-1}]"
        )
    m36_linear = rep_j.astype(np.int64) * int(tgt.nx) + rep_i.astype(np.int64)

    tgt_ids: list[np.ndarray] = []
    src_ids: list[np.ndarray] = []
    areas: list[np.ndarray] = []

    for t, (i, j) in enumerate(zip(rep_i.tolist(), rep_j.tolist())):
        sx, ox = x_map[i]
        sy, oy = y_map[j]
        if sx.size == 0 or sy.size == 0:
            continue

        src_grid = sy[:, None].astype(np.int64) * nlon + sx[None, :].astype(np.int64)
        area_grid = oy[:, None] * ox[None, :]

        src_flat = src_grid.ravel(order="C")
        area_flat = area_grid.ravel(order="C")
        ok = area_flat > 0.0
        if not np.any(ok):
            continue

        src_flat = src_flat[ok]
        area_flat = area_flat[ok]
        tgt_flat = np.full(src_flat.size, t, dtype=np.int32)

        tgt_ids.append(tgt_flat)
        src_ids.append(src_flat.astype(np.int64))
        areas.append(area_flat.astype(np.float64))

    if not tgt_ids:
        raise RuntimeError("No overlap mapping entries were built")

    map_tgt = np.concatenate(tgt_ids)
    map_src = np.concatenate(src_ids)
    map_area = np.concatenate(areas)

    target_cell_area = tgt.dx_m * tgt.dx_m

    cache = {
        "map_tgt": map_tgt,
        "map_src": map_src,
        "map_area": map_area,
        "n_tgt": np.int32(rep_i.size),
        "src_nlat": np.int32(nlat),
        "src_nlon": np.int32(nlon),
        "target_cell_area": np.float64(target_cell_area),
        "rep_i": rep_i,
        "rep_j": rep_j,
        "m36_linear": m36_linear.astype(np.int64),
        "rep_tile_id": rep["rep_tile_id"].astype(np.int32),
        "rep_frac_cell": rep["rep_frac_cell"].astype(np.float32),
        "rep_lat": rep["rep_lat"].astype(np.float32),
        "rep_lon": rep["rep_lon"].astype(np.float32),
        "m36_nx": np.int32(int(tgt.nx)),
        "m36_ny": np.int32(int(tgt.ny)),
        "source_gt": np.asarray(gt, dtype=np.float64),
    }

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, **cache)
    return cache


def load_mapping(cache_path: Path) -> dict[str, np.ndarray]:
    with np.load(cache_path) as z:
        return {k: z[k] for k in z.files}


def _filled(arr: np.ndarray | np.ma.MaskedArray) -> np.ndarray:
    if np.ma.isMaskedArray(arr):
        return np.asarray(arr.filled(np.nan))
    return np.asarray(arr)


def _infer_axis_by_name_or_size(var: np.ndarray, dims: tuple[str, ...], name: str, expected_size: int) -> int:
    if name in dims:
        return int(dims.index(name))
    candidates = [i for i, n in enumerate(var.shape) if int(n) == int(expected_size)]
    if len(candidates) == 1:
        return int(candidates[0])
    raise RuntimeError(f"Could not infer axis for '{name}' from dims={dims}, shape={var.shape}")


def _mono_sign(v: np.ndarray) -> int:
    dv = np.diff(np.asarray(v, dtype=np.float64))
    dv = dv[np.isfinite(dv)]
    if dv.size == 0:
        return 0
    med = float(np.nanmedian(dv))
    if med > 0.0:
        return 1
    if med < 0.0:
        return -1
    return 0


def _coord_line_1d_or_2d(
    arr: np.ndarray,
    canonical_shape: tuple[int, int],
    axis: int,
    transpose_2d: bool = False,
) -> np.ndarray | None:
    a = _filled(arr).astype(np.float64, copy=False)
    if a.ndim == 1:
        if a.size == canonical_shape[axis]:
            return a
        return None
    if a.ndim == 2:
        if a.shape == canonical_shape:
            pass
        elif transpose_2d and a.shape == (canonical_shape[1], canonical_shape[0]):
            a = a.T
        else:
            return None
        other_axis = 1 - axis
        return np.nanmean(a, axis=other_axis)
    return None


def infer_source_layout(sample_nc: Path) -> SourceLayout:
    with Dataset(sample_nc, "r") as ds:
        sm_var = ds.variables["Soil_Moisture"]
        dims = tuple(str(d) for d in sm_var.dimensions)
        sm_shape = tuple(int(n) for n in sm_var.shape)
        if len(sm_shape) != 2:
            raise RuntimeError(f"Expected 2D Soil_Moisture, got shape={sm_shape}")
        if "lat" not in ds.dimensions or "lon" not in ds.dimensions:
            raise KeyError(f"{sample_nc} missing lat/lon dimensions")

        nlat = int(len(ds.dimensions["lat"]))
        nlon = int(len(ds.dimensions["lon"]))
        lat_axis = _infer_axis_by_name_or_size(sm_var, dims, "lat", nlat)
        lon_axis = _infer_axis_by_name_or_size(sm_var, dims, "lon", nlon)
        if lat_axis == lon_axis:
            raise RuntimeError(f"Invalid Soil_Moisture axis mapping: lat_axis={lat_axis}, lon_axis={lon_axis}")

        if (lat_axis, lon_axis) not in ((0, 1), (1, 0)):
            raise RuntimeError(
                f"Unsupported axis order for Soil_Moisture: dims={dims}, lat_axis={lat_axis}, lon_axis={lon_axis}"
            )
        transpose = (lat_axis, lon_axis) == (1, 0)
        canonical_shape = (nlat, nlon)

        flip_lat = False
        if "lat" in ds.variables:
            lat_line = _coord_line_1d_or_2d(ds.variables["lat"][:], canonical_shape, axis=0, transpose_2d=transpose)
            if lat_line is not None:
                # Canonical source rows must be north->south (lat descending).
                s = _mono_sign(lat_line)
                flip_lat = s > 0

        flip_lon = False
        if "lon" in ds.variables:
            lon_line = _coord_line_1d_or_2d(ds.variables["lon"][:], canonical_shape, axis=1, transpose_2d=transpose)
            if lon_line is not None:
                # Canonical source cols must be west->east (lon ascending).
                s = _mono_sign(lon_line)
                flip_lon = s < 0

    return SourceLayout(
        raw_shape=sm_shape,
        nlat=nlat,
        nlon=nlon,
        lat_axis=lat_axis,
        lon_axis=lon_axis,
        transpose=transpose,
        flip_lat=flip_lat,
        flip_lon=flip_lon,
    )


def canonicalize_source_array(arr: np.ndarray, layout: SourceLayout) -> np.ndarray:
    a = _filled(arr).astype(np.float64, copy=False)
    if tuple(int(n) for n in a.shape) != tuple(int(n) for n in layout.raw_shape):
        raise RuntimeError(f"Unexpected source shape {a.shape}; expected {layout.raw_shape}")
    if layout.transpose:
        a = a.T
    if layout.flip_lat:
        a = a[::-1, :]
    if layout.flip_lon:
        a = a[:, ::-1]
    return np.asarray(a, dtype=np.float64)


def _read_smos_fields(path: Path, layout: SourceLayout) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(path, "r") as ds:
        sm = ds.variables["Soil_Moisture"][:]
        sf = ds.variables["Scene_Flags"][:]
        rmse = ds.variables["RMSE"][:]

    return (
        canonicalize_source_array(sm, layout),
        canonicalize_source_array(sf, layout),
        canonicalize_source_array(rmse, layout),
    )


def apply_qc(
    sm: np.ndarray,
    sf: np.ndarray,
    rmse: np.ndarray,
    scene_flag_max: int,
    tb_rmse_max: float,
    sm_min: float,
    sm_max: float,
) -> np.ndarray:
    good = (
        np.isfinite(sm)
        & np.isfinite(sf)
        & np.isfinite(rmse)
        & (sf <= float(scene_flag_max))
        & (rmse <= float(tb_rmse_max))
        & (sm >= float(sm_min))
        & (sm <= float(sm_max))
    )
    out = np.full(sm.shape, np.nan, dtype=np.float64)
    out[good] = sm[good]
    return out


def merge_asc_des(sm_asc: np.ndarray | None, sm_des: np.ndarray | None) -> np.ndarray | None:
    if sm_asc is None and sm_des is None:
        return None
    if sm_asc is None:
        return sm_des
    if sm_des is None:
        return sm_asc

    fa = np.isfinite(sm_asc)
    fd = np.isfinite(sm_des)

    out = np.full(sm_asc.shape, np.nan, dtype=np.float64)
    both = fa & fd
    out[both] = 0.5 * (sm_asc[both] + sm_des[both])
    out[fa & ~fd] = sm_asc[fa & ~fd]
    out[~fa & fd] = sm_des[~fa & fd]
    return out


def remap_source_to_m36(
    src_sm: np.ndarray,
    mapping: dict[str, np.ndarray],
    min_coverage_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Area-weighted conservative remap from source grid to representative M36 cells."""
    src_flat = src_sm.ravel(order="C")

    map_tgt = mapping["map_tgt"].astype(np.int64)
    map_src = mapping["map_src"].astype(np.int64)
    map_area = mapping["map_area"].astype(np.float64)
    n_tgt = int(mapping["n_tgt"])
    target_cell_area = float(mapping["target_cell_area"])

    vals = src_flat[map_src]
    good = np.isfinite(vals)

    num = np.bincount(map_tgt[good], weights=map_area[good] * vals[good], minlength=n_tgt)
    den = np.bincount(map_tgt[good], weights=map_area[good], minlength=n_tgt)

    out = np.full(n_tgt, np.nan, dtype=np.float64)
    coverage = den / target_cell_area
    ok = den > (float(min_coverage_frac) * target_cell_area)
    out[ok] = num[ok] / den[ok]
    return out, coverage


def write_daily_sparse_netcdf(
    out_nc: Path,
    idx: np.ndarray,
    vals: np.ndarray,
    cov: np.ndarray,
    tag: str,
    args: argparse.Namespace,
) -> None:
    """Write one day of sparse M36 values in NetCDF format."""
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(out_nc, "w", format="NETCDF4") as ds:
        n = int(idx.size)
        ds.createDimension("n_points", n)

        v_idx = ds.createVariable("idx_EASEv2_lonxlat", "i4", ("n_points",), zlib=True, complevel=4)
        v_idx.long_name = "EASEv2 M36 linear index, 0-based (column-major flattening)"
        v_idx[:] = idx.astype(np.int32, copy=False)

        v_sm = ds.createVariable("sm_obs", "f4", ("n_points",), zlib=True, complevel=4)
        v_sm.long_name = "SMOS-IC daily soil moisture on M36"
        v_sm.units = "m3 m-3"
        v_sm[:] = vals.astype(np.float32, copy=False)

        v_cov = ds.createVariable("coverage_frac", "f4", ("n_points",), zlib=True, complevel=4)
        v_cov.long_name = "Conservative overlap coverage fraction for each retained M36 cell"
        v_cov.units = "1"
        v_cov[:] = cov.astype(np.float32, copy=False)

        ds.date = tag
        ds.qc_scene_flag_max = int(args.scene_flag_max)
        ds.qc_tb_rmse_max = float(args.tb_rmse_max)
        ds.qc_sm_min = float(args.sm_min)
        ds.qc_sm_max = float(args.sm_max)
        ds.min_coverage_frac = float(args.min_coverage_frac)


def date_range(start: date, end: date) -> list[date]:
    out: list[date] = []
    d = start
    while d <= end:
        out.append(d)
        d += timedelta(days=1)
    return out


def format_day(d: date) -> str:
    return d.strftime("%Y%m%d")


def main() -> None:
    args = parse_args()

    smos_root = Path(args.smos_root).expanduser().resolve()
    tilecoord = Path(args.tilecoord).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()

    out_dir.mkdir(parents=True, exist_ok=True)

    mapping_cache = (
        Path(args.mapping_cache).expanduser().resolve()
        if args.mapping_cache
        else out_dir / "smos25km_to_m36_conservative_mapping.npz"
    )

    repo_root = find_repo_root(Path(__file__).resolve())
    sys.path.insert(0, str(repo_root / "common/python/io"))
    from read_GEOSldas import read_tilecoord  # type: ignore

    asc_idx, des_idx = build_file_index(smos_root)
    all_dates = sorted(set(asc_idx.keys()) | set(des_idx.keys()))

    if not all_dates:
        raise RuntimeError("No dates found in SMOS-IC indexes")

    if args.start_date:
        d0 = datetime.strptime(args.start_date, "%Y-%m-%d").date()
    else:
        d0 = all_dates[0]

    if args.end_date:
        d1 = datetime.strptime(args.end_date, "%Y-%m-%d").date()
    else:
        d1 = all_dates[-1]

    if d1 < d0:
        raise ValueError(f"Invalid date range: {d0}..{d1}")

    run_dates = [d for d in date_range(d0, d1) if d in asc_idx or d in des_idx]

    print(f"SMOS root: {smos_root}")
    print(f"Date range requested: {d0} .. {d1}")
    print(f"Dates with ASC/DES files: {len(run_dates)}")
    print(f"QC: Scene_Flags <= {args.scene_flag_max}, RMSE <= {args.tb_rmse_max} K, {args.sm_min} <= SM <= {args.sm_max}")

    sample_path = asc_idx.get(run_dates[0], None) or des_idx.get(run_dates[0], None)
    if sample_path is None:
        raise RuntimeError("Could not find sample SMOS file")
    sample_nlat, sample_nlon, sample_gt = read_source_georef(sample_path)
    source_layout = infer_source_layout(sample_path)
    print(
        "Source layout normalization: "
        f"shape={source_layout.raw_shape}, axes(lat,lon)=({source_layout.lat_axis},{source_layout.lon_axis}), "
        f"transpose={source_layout.transpose}, flip_lat={source_layout.flip_lat}, flip_lon={source_layout.flip_lon}"
    )

    tc = read_tilecoord(str(tilecoord))
    rep = choose_representative_tile_per_cell(tc)
    tgt = M36Grid()  # always use full global M36 geometry (964x406)
    print(
        "Tilecoord representative index range: "
        f"i=[{int(rep['rep_i'].min())},{int(rep['rep_i'].max())}] "
        f"j=[{int(rep['rep_j'].min())},{int(rep['rep_j'].max())}] "
        f"(local max extents nx={int(rep['nx_local'])}, ny={int(rep['ny_local'])})"
    )

    mapping: dict[str, np.ndarray] | None = None
    cache_reasons: list[str] = []
    if mapping_cache.exists() and not args.rebuild_mapping:
        mapping = load_mapping(mapping_cache)
        map_nlat = int(mapping["src_nlat"])
        map_nlon = int(mapping["src_nlon"])
        if map_nlat != int(sample_nlat) or map_nlon != int(sample_nlon):
            cache_reasons.append(
                f"source dims cache=({map_nlat},{map_nlon}) current=({sample_nlat},{sample_nlon})"
            )
        map_m36_nx = int(mapping.get("m36_nx", np.int32(-1)))
        map_m36_ny = int(mapping.get("m36_ny", np.int32(-1)))
        if map_m36_nx != int(tgt.nx) or map_m36_ny != int(tgt.ny):
            cache_reasons.append(
                f"target dims cache=({map_m36_nx},{map_m36_ny}) expected=({tgt.nx},{tgt.ny})"
            )
        if "source_gt" in mapping:
            gt_cache = np.asarray(mapping["source_gt"], dtype=np.float64).reshape(-1)
            if gt_cache.size != sample_gt.size or not np.allclose(gt_cache, sample_gt, rtol=0.0, atol=1e-6):
                cache_reasons.append("source GeoTransform differs from current sample file")
        else:
            cache_reasons.append("cache missing source_gt metadata")

        if cache_reasons:
            print("Mapping cache is incompatible and will be rebuilt:")
            for r in cache_reasons:
                print(f"  - {r}")
            mapping = None
        else:
            print(f"Loaded mapping cache: {mapping_cache}")

    if mapping is None:
        mapping = build_conservative_mapping(
            sample_nc=sample_path,
            rep=rep,
            tgt=tgt,
            cache_path=mapping_cache,
        )
        print(f"Built mapping cache: {mapping_cache}")

    map_nlat = int(mapping["src_nlat"])
    map_nlon = int(mapping["src_nlon"])
    if map_nlat != int(source_layout.nlat) or map_nlon != int(source_layout.nlon):
        raise RuntimeError(
            f"Mapping source dims ({map_nlat},{map_nlon}) do not match source layout ({source_layout.nlat},{source_layout.nlon}). "
            "Rebuild mapping cache with --rebuild-mapping."
        )
    if int(mapping["m36_nx"]) != int(tgt.nx) or int(mapping["m36_ny"]) != int(tgt.ny):
        raise RuntimeError(
            f"Mapping target dims ({int(mapping['m36_nx'])},{int(mapping['m36_ny'])}) do not match expected ({tgt.nx},{tgt.ny}). "
            "Rebuild mapping cache with --rebuild-mapping."
        )

    m36_linear = mapping["m36_linear"].astype(np.int64)

    manifest_path = out_dir / "smos_ic_m36_daily_manifest.csv"
    fieldnames = [
        "date",
        "asc_file",
        "des_file",
        "n_valid_asc_native",
        "n_valid_des_native",
        "n_valid_combined_native",
        "n_valid_m36",
        "mean_coverage_m36",
        "output_nc",
    ]

    rows: list[dict[str, str]] = []

    for n, d in enumerate(run_dates, start=1):
        tag = format_day(d)
        out_nc = out_dir / f"smos_ic_sm_m36_{tag}.nc"

        if out_nc.exists() and not args.overwrite:
            rows.append(
                {
                    "date": tag,
                    "asc_file": str(asc_idx.get(d, "")),
                    "des_file": str(des_idx.get(d, "")),
                    "n_valid_asc_native": "",
                    "n_valid_des_native": "",
                    "n_valid_combined_native": "",
                    "n_valid_m36": "",
                    "mean_coverage_m36": "",
                    "output_nc": str(out_nc),
                }
            )
            continue

        asc_file = asc_idx.get(d, None)
        des_file = des_idx.get(d, None)

        sm_asc = None
        sm_des = None
        n_asc = 0
        n_des = 0

        if asc_file is not None:
            sm, sf, rmse = _read_smos_fields(asc_file, source_layout)
            sm_asc = apply_qc(
                sm=sm,
                sf=sf,
                rmse=rmse,
                scene_flag_max=args.scene_flag_max,
                tb_rmse_max=args.tb_rmse_max,
                sm_min=args.sm_min,
                sm_max=args.sm_max,
            )
            n_asc = int(np.isfinite(sm_asc).sum())

        if des_file is not None:
            sm, sf, rmse = _read_smos_fields(des_file, source_layout)
            sm_des = apply_qc(
                sm=sm,
                sf=sf,
                rmse=rmse,
                scene_flag_max=args.scene_flag_max,
                tb_rmse_max=args.tb_rmse_max,
                sm_min=args.sm_min,
                sm_max=args.sm_max,
            )
            n_des = int(np.isfinite(sm_des).sum())

        sm_comb = merge_asc_des(sm_asc, sm_des)
        if sm_comb is None:
            continue

        n_comb = int(np.isfinite(sm_comb).sum())

        sm_m36, coverage = remap_source_to_m36(
            src_sm=sm_comb,
            mapping=mapping,
            min_coverage_frac=float(args.min_coverage_frac),
        )

        valid_tgt = np.isfinite(sm_m36)
        idx = m36_linear[valid_tgt].astype(np.int32)
        vals = sm_m36[valid_tgt].astype(np.float32)
        cov = coverage[valid_tgt].astype(np.float32)

        write_daily_sparse_netcdf(
            out_nc=out_nc,
            idx=idx,
            vals=vals,
            cov=cov,
            tag=tag,
            args=args,
        )

        rows.append(
            {
                "date": tag,
                "asc_file": str(asc_file) if asc_file else "",
                "des_file": str(des_file) if des_file else "",
                "n_valid_asc_native": str(n_asc),
                "n_valid_des_native": str(n_des),
                "n_valid_combined_native": str(n_comb),
                "n_valid_m36": str(int(valid_tgt.sum())),
                "mean_coverage_m36": f"{float(np.nanmean(coverage[valid_tgt])):.6f}" if np.any(valid_tgt) else "",
                "output_nc": str(out_nc),
            }
        )

        if (n % 50) == 0 or n == len(run_dates):
            print(f"Processed {n}/{len(run_dates)} days ({tag})")

    with open(manifest_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote manifest: {manifest_path}")
    print(f"Done. Days listed: {len(rows)}")


if __name__ == "__main__":
    main()
