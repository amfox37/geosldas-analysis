#!/usr/bin/env python3
"""
Nearest-neighbor regridding of IMS categorical data onto GEOSldas M36 grid.

Design goals (explicitly enforced):
1) Categorical-safe: nearest-neighbor only (no interpolation, no averaging).
2) One category per EASE cell: choose one deterministic representative tile
   per EASE (y, x) cell before mapping IMS.
3) Guard rails:
   - max-distance threshold (km) to reject suspicious source matches
   - mapping diagnostics and metadata in output

Typical usage:
  python regrid_ims_to_m36_nearest.py \
    --ims-nc ../output/ims_snowcover_24km_2024.nc4 \
    --tilecoord /discover/.../LS_OLv8_M36.ldas_tilecoord.bin \
    --out-nc ../output/ims_snowcover_24km_2024_on_m36_nearest.nc4

Batch yearly usage:
  python regrid_ims_to_m36_nearest.py \
    --ims-dir /gpfsm/dnb06/projects/p163/IMS \
    --out-dir /gpfsm/dnb06/projects/p163/IMS \
    --year-start 2000 \
    --year-end 2023 \
    --tilecoord /discover/.../LS_OLv8_M36.ldas_tilecoord.bin
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

try:
    # cKDTree is fast and robust for nearest-neighbor search.
    from scipy.spatial import cKDTree
except Exception as exc:  # pragma: no cover
    raise ImportError(
        "scipy is required for nearest-neighbor mapping (cKDTree). "
        "Install scipy in your environment."
    ) from exc


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Nearest-neighbor regrid of IMS categorical grids to GEOSldas M36 cells."
    )
    ap.add_argument("--ims-nc", help="Input IMS NetCDF with categorical grid (time,y,x)")
    ap.add_argument("--ims-dir", help="Directory containing yearly IMS NetCDF files for batch mode")
    ap.add_argument("--ims-var", default="ims_snowcover", help="Input IMS variable name (default: ims_snowcover)")
    ap.add_argument(
        "--ims-template",
        default="ims_snowcover_24km_{year}.nc4",
        help="Batch mode input filename template (must include {year})",
    )
    ap.add_argument("--tilecoord", required=True, help="GEOSldas tilecoord binary file (.ldas_tilecoord.bin)")
    ap.add_argument("--out-nc", help="Output NetCDF path on M36 grid (single-file mode)")
    ap.add_argument("--out-dir", help="Batch mode output directory (default: --ims-dir)")
    ap.add_argument(
        "--out-template",
        default="ims_snowcover_24km_{year}_on_m36_nearest.nc4",
        help="Batch mode output filename template (must include {year})",
    )
    ap.add_argument("--year-start", type=int, default=2000, help="Batch mode first year (default: 2000)")
    ap.add_argument("--year-end", type=int, default=2023, help="Batch mode last year (default: 2023)")
    ap.add_argument(
        "--max-distance-km",
        type=float,
        default=60.0,
        help="Maximum allowed nearest-source distance in km (default: 60)",
    )
    ap.add_argument(
        "--fill-value",
        type=int,
        default=-32768,
        help="Fill value for unmapped/missing categorical cells (default: -32768)",
    )
    ap.add_argument(
        "--overwrite-output",
        action="store_true",
        help="Overwrite output NetCDF if it already exists",
    )
    ap.add_argument(
        "--skip-existing-output",
        action="store_true",
        help="Batch mode: skip files whose outputs already exist",
    )
    ap.add_argument(
        "--skip-missing-input",
        action="store_true",
        help="Batch mode: skip missing yearly IMS inputs instead of failing",
    )
    ap.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Batch mode: continue to next file if one file fails",
    )
    ap.add_argument(
        "--copy-time-values",
        action="store_true",
        help=(
            "If possible, copy 1D time variable values from input. "
            "Checks [time_dim, time, doy, day_of_year], else writes 0..N-1."
        ),
    )
    return ap.parse_args()


def find_repo_root(start: Path) -> Path:
    """Find repo root by locating common/python/io/read_GEOSldas.py upward."""
    p = start.resolve()
    candidates = [p] + list(p.parents)
    for c in candidates:
        if (c / "common/python/io/read_GEOSldas.py").exists():
            return c
    raise FileNotFoundError(
        "Could not locate repo root containing common/python/io/read_GEOSldas.py"
    )


def wrap_lon_180(lon: np.ndarray) -> np.ndarray:
    """Wrap longitude to [-180, 180)."""
    return ((lon + 180.0) % 360.0) - 180.0


def latlon_to_unit_xyz(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """
    Convert lat/lon (degrees) to unit-sphere xyz.

    Using xyz avoids dateline issues in nearest-neighbor search.
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    clat = np.cos(lat)
    x = clat * np.cos(lon)
    y = clat * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack([x, y, z])


def chord_to_km(chord_dist: np.ndarray, radius_km: float = 6371.0) -> np.ndarray:
    """
    Convert unit-sphere chord distance to great-circle distance in km.
    """
    # chord = 2*sin(theta/2)  =>  theta = 2*asin(chord/2)
    x = np.clip(chord_dist / 2.0, 0.0, 1.0)
    theta = 2.0 * np.arcsin(x)
    return radius_km * theta


def choose_representative_tile_per_cell(tc: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Deterministically pick one tile per EASE cell.

    Rule:
    - Primary: highest frac_cell in the cell
    - Tie-break: smallest tile_id

    This guarantees exactly one target category cell per (y, x), with no averaging.
    """
    i_indg = tc["i_indg"].astype(np.int64)
    j_indg = tc["j_indg"].astype(np.int64)
    tile_id = tc["tile_id"].astype(np.int64)
    frac_cell = np.asarray(tc["frac_cell"], dtype=np.float64)

    # Replace NaN frac with -inf so valid fractions always win.
    frac_sort = np.where(np.isfinite(frac_cell), frac_cell, -np.inf)

    nx = int(i_indg.max()) + 1
    ny = int(j_indg.max()) + 1
    cell_code = j_indg * nx + i_indg

    # Sort by:
    # 1) cell_code ascending
    # 2) frac_cell descending (use -frac for lexsort)
    # 3) tile_id ascending
    order = np.lexsort((tile_id, -frac_sort, cell_code))
    code_sorted = cell_code[order]

    # Keep first row for each cell in sorted order.
    first = np.empty(code_sorted.size, dtype=bool)
    first[0] = True
    first[1:] = code_sorted[1:] != code_sorted[:-1]
    rep_idx = order[first]

    out = {
        "rep_idx": rep_idx,
        "rep_i": i_indg[rep_idx].astype(np.int32),
        "rep_j": j_indg[rep_idx].astype(np.int32),
        "rep_tile_id": tile_id[rep_idx].astype(np.int32),
        "rep_frac_cell": frac_cell[rep_idx].astype(np.float32),
        "rep_lat": np.asarray(tc["com_lat"][rep_idx], dtype=np.float64),
        "rep_lon": np.asarray(tc["com_lon"][rep_idx], dtype=np.float64),
        "rep_elev": np.asarray(tc["elev"][rep_idx], dtype=np.float64),
        "nx": np.int32(nx),
        "ny": np.int32(ny),
    }
    return out


def build_jobs(args: argparse.Namespace) -> list[tuple[int | None, Path, Path]]:
    """Build single-file or batch job list."""
    if bool(args.ims_nc) != bool(args.out_nc):
        raise ValueError("Single-file mode requires both --ims-nc and --out-nc together.")

    if args.ims_nc and args.out_nc:
        if args.ims_dir or args.out_dir:
            raise ValueError("Use either single-file args (--ims-nc/--out-nc) or batch args (--ims-dir/--out-dir).")
        return [(None, Path(args.ims_nc), Path(args.out_nc))]

    if not args.ims_dir:
        raise ValueError(
            "No inputs specified. Provide --ims-nc/--out-nc for one file, "
            "or --ims-dir with --year-start/--year-end for batch mode."
        )
    if args.year_end < args.year_start:
        raise ValueError(f"Invalid year range: {args.year_start}..{args.year_end}")

    if "{year}" not in args.ims_template or "{year}" not in args.out_template:
        raise ValueError("Both --ims-template and --out-template must include '{year}'.")

    ims_dir = Path(args.ims_dir)
    out_dir = Path(args.out_dir) if args.out_dir else ims_dir

    jobs: list[tuple[int | None, Path, Path]] = []
    for year in range(args.year_start, args.year_end + 1):
        ims_name = args.ims_template.format(year=year)
        out_name = args.out_template.format(year=year)
        jobs.append((year, ims_dir / ims_name, out_dir / out_name))
    return jobs


def select_time_values(
    ds_in: Dataset,
    time_dim: str,
    n_time: int,
    copy_time_values: bool,
) -> tuple[np.ndarray, dict[str, str], str | None]:
    """
    Pick 1D time coordinate values robustly across input variants.

    Priority: [time_dim variable, time, doy, day_of_year].
    """
    if not copy_time_values:
        return np.arange(n_time, dtype=np.float64), {"long_name": "time_index"}, None

    for name in (time_dim, "time", "doy", "day_of_year"):
        if not name or name not in ds_in.variables:
            continue
        v = ds_in.variables[name]
        if getattr(v, "ndim", 0) != 1 or len(v) != n_time:
            continue
        values = np.asarray(v[:], dtype=np.float64)
        attrs: dict[str, str] = {}
        for a in ("units", "calendar", "long_name", "standard_name", "axis"):
            if hasattr(v, a):
                attrs[a] = str(getattr(v, a))
        if name == "doy" and "long_name" not in attrs:
            attrs["long_name"] = "day_of_year"
        return values, attrs, name

    return np.arange(n_time, dtype=np.float64), {"long_name": "time_index"}, None


def regrid_one_file(
    ims_nc: Path,
    out_nc: Path,
    args: argparse.Namespace,
    tilecoord_path: Path,
    rep: dict[str, np.ndarray],
) -> None:
    """Regrid one IMS file to the M36 grid using nearest-neighbor only."""
    nx = int(rep["nx"])
    ny = int(rep["ny"])
    rep_i = rep["rep_i"]
    rep_j = rep["rep_j"]
    rep_lat = np.asarray(rep["rep_lat"], dtype=np.float64)
    rep_lon = wrap_lon_180(np.asarray(rep["rep_lon"], dtype=np.float64))

    print(f"IMS input: {ims_nc}")
    print(f"IMS var: {args.ims_var}")
    print(f"Output: {out_nc}")

    with Dataset(ims_nc, "r") as ds_in:
        if args.ims_var not in ds_in.variables:
            raise KeyError(f"IMS var '{args.ims_var}' not found in {ims_nc}")
        if "lat" not in ds_in.variables or "lon" not in ds_in.variables:
            raise KeyError("IMS input must contain 2D 'lat' and 'lon' variables")

        ims_var = ds_in.variables[args.ims_var]
        ims_dims = ims_var.dimensions
        if len(ims_dims) != 3:
            raise ValueError(
                f"IMS var {args.ims_var} must be 3D (time,y,x), found dims={ims_dims}"
            )
        time_dim = str(ims_dims[0])

        src_lat = np.asarray(ds_in.variables["lat"][:], dtype=np.float64)
        src_lon = wrap_lon_180(np.asarray(ds_in.variables["lon"][:], dtype=np.float64))
        if src_lat.ndim != 2 or src_lon.ndim != 2:
            raise ValueError("IMS lat/lon must be 2D")
        if src_lat.shape != src_lon.shape:
            raise ValueError("IMS lat/lon shape mismatch")

        n_time = int(ims_var.shape[0])
        src_ny, src_nx = src_lat.shape

        transpose_input_spatial = False
        if tuple(ims_var.shape[1:]) == tuple(src_lat.shape):
            transpose_input_spatial = False
        elif tuple(ims_var.shape[1:]) == tuple(src_lat.T.shape):
            transpose_input_spatial = True
        else:
            raise ValueError(
                f"IMS var spatial shape {ims_var.shape[1:]} does not match lat/lon shape "
                f"{src_lat.shape} (or transpose {src_lat.T.shape})"
            )

        print(f"IMS var dims: {ims_dims}, shape={ims_var.shape}")
        print(f"IMS source grid: ny={src_ny}, nx={src_nx}")
        if transpose_input_spatial:
            print("Detected input spatial order (x, y); transposing each slice to (y, x).")

        src_valid = np.isfinite(src_lat) & np.isfinite(src_lon)
        src_flat_idx = np.flatnonzero(src_valid.ravel())
        if src_flat_idx.size == 0:
            raise RuntimeError("No finite IMS source lat/lon points found")

        src_lat_v = src_lat.ravel()[src_flat_idx]
        src_lon_v = src_lon.ravel()[src_flat_idx]
        src_xyz = latlon_to_unit_xyz(src_lat_v, src_lon_v)
        tree = cKDTree(src_xyz)

        tgt_xyz = latlon_to_unit_xyz(rep_lat, rep_lon)
        chord, nn_idx = tree.query(tgt_xyz, k=1)
        dist_km = chord_to_km(np.asarray(chord, dtype=np.float64))

        within = dist_km <= float(args.max_distance_km)
        matched_src_flat = np.full(rep_i.size, -1, dtype=np.int64)
        matched_src_flat[within] = src_flat_idx[np.asarray(nn_idx[within], dtype=np.int64)]

        matched_src_y = np.full(rep_i.size, -1, dtype=np.int32)
        matched_src_x = np.full(rep_i.size, -1, dtype=np.int32)
        if np.any(within):
            yy, xx = np.unravel_index(matched_src_flat[within], (src_ny, src_nx))
            matched_src_y[within] = yy.astype(np.int32)
            matched_src_x[within] = xx.astype(np.int32)

        print(
            "Nearest mapping stats:",
            f"matched={int(np.sum(within))}/{within.size}",
            f"min_km={float(np.nanmin(dist_km)):.3f}",
            f"p50_km={float(np.nanmedian(dist_km)):.3f}",
            f"max_km={float(np.nanmax(dist_km)):.3f}",
        )

        fill_i32 = np.int32(-1)
        map_src_y_grid = np.full((ny, nx), fill_i32, dtype=np.int32)
        map_src_x_grid = np.full((ny, nx), fill_i32, dtype=np.int32)
        map_dist_grid = np.full((ny, nx), np.nan, dtype=np.float32)
        map_tile_id_grid = np.full((ny, nx), fill_i32, dtype=np.int32)
        map_frac_grid = np.full((ny, nx), np.nan, dtype=np.float32)
        map_elev_grid = np.full((ny, nx), np.nan, dtype=np.float32)
        m36_lat_grid = np.full((ny, nx), np.nan, dtype=np.float32)
        m36_lon_grid = np.full((ny, nx), np.nan, dtype=np.float32)
        land_mask_grid = np.zeros((ny, nx), dtype=np.int8)
        within_grid = np.zeros((ny, nx), dtype=np.int8)

        map_tile_id_grid[rep_j, rep_i] = rep["rep_tile_id"].astype(np.int32)
        map_frac_grid[rep_j, rep_i] = rep["rep_frac_cell"].astype(np.float32)
        map_elev_grid[rep_j, rep_i] = rep["rep_elev"].astype(np.float32)
        m36_lat_grid[rep_j, rep_i] = rep_lat.astype(np.float32)
        m36_lon_grid[rep_j, rep_i] = rep_lon.astype(np.float32)
        land_mask_grid[rep_j, rep_i] = 1

        map_src_y_grid[rep_j, rep_i] = matched_src_y
        map_src_x_grid[rep_j, rep_i] = matched_src_x
        map_dist_grid[rep_j, rep_i] = dist_km.astype(np.float32)
        within_grid[rep_j, rep_i] = within.astype(np.int8)

        time_values, time_attrs, time_source_name = select_time_values(
            ds_in=ds_in,
            time_dim=time_dim,
            n_time=n_time,
            copy_time_values=bool(args.copy_time_values),
        )
        if time_source_name:
            print(f"Time coordinate copied from input variable: {time_source_name}")
        else:
            print("Time coordinate source variable not found; writing 0..N-1 index.")

        ds_out = Dataset(out_nc, "w", format="NETCDF4")
        try:
            ds_out.createDimension("time", n_time)
            ds_out.createDimension("y", ny)
            ds_out.createDimension("x", nx)

            tvar_out = ds_out.createVariable("time", "f8", ("time",))
            tvar_out[:] = time_values
            for k, v in time_attrs.items():
                setattr(tvar_out, k, v)

            lat_out = ds_out.createVariable("lat", "f4", ("y", "x"), zlib=True, complevel=4)
            lon_out = ds_out.createVariable("lon", "f4", ("y", "x"), zlib=True, complevel=4)
            lat_out[:] = m36_lat_grid
            lon_out[:] = m36_lon_grid
            lat_out.units = "degrees_north"
            lon_out.units = "degrees_east"

            land_out = ds_out.createVariable("land_mask", "i1", ("y", "x"), zlib=True, complevel=4)
            land_out[:] = land_mask_grid
            land_out.long_name = "M36 land cell mask from tilecoord representative tiles (1=land)"

            within_out = ds_out.createVariable(
                "within_max_distance", "i1", ("y", "x"), zlib=True, complevel=4
            )
            within_out[:] = within_grid
            within_out.long_name = "Nearest IMS source point is within max distance threshold (1=yes)"

            dist_out = ds_out.createVariable(
                "nearest_distance_km", "f4", ("y", "x"), zlib=True, complevel=4
            )
            dist_out[:] = map_dist_grid
            dist_out.units = "km"

            sy_out = ds_out.createVariable("source_y", "i4", ("y", "x"), zlib=True, complevel=4, fill_value=-1)
            sx_out = ds_out.createVariable("source_x", "i4", ("y", "x"), zlib=True, complevel=4, fill_value=-1)
            sy_out[:] = map_src_y_grid
            sx_out[:] = map_src_x_grid
            sy_out.long_name = "IMS source y index used by nearest map (-1 if no valid map)"
            sx_out.long_name = "IMS source x index used by nearest map (-1 if no valid map)"

            tid_out = ds_out.createVariable(
                "rep_tile_id", "i4", ("y", "x"), zlib=True, complevel=4, fill_value=-1
            )
            fcf_out = ds_out.createVariable("rep_frac_cell", "f4", ("y", "x"), zlib=True, complevel=4)
            elv_out = ds_out.createVariable("rep_elev_m", "f4", ("y", "x"), zlib=True, complevel=4)
            tid_out[:] = map_tile_id_grid
            fcf_out[:] = map_frac_grid
            elv_out[:] = map_elev_grid
            fcf_out.long_name = "Representative tile frac_cell used to select one tile per M36 cell"
            elv_out.units = "m"

            cat_out = ds_out.createVariable(
                "ims_category",
                "i2",
                ("time", "y", "x"),
                zlib=True,
                complevel=4,
                fill_value=np.int16(args.fill_value),
            )
            cat_out.long_name = "IMS categorical code mapped to M36 by nearest-neighbor (no averaging)"
            cat_out.comment = (
                "One category per M36 cell. Representative tile per cell selected by "
                "max frac_cell then min tile_id. Source match rejected if distance > max_distance_km."
            )

            valid_map = within & (matched_src_y >= 0) & (matched_src_x >= 0)
            tgt_j = rep_j[valid_map]
            tgt_i = rep_i[valid_map]
            src_j = matched_src_y[valid_map]
            src_i = matched_src_x[valid_map]

            print(f"Time steps to process: {n_time}")
            for t in range(n_time):
                out_slice = np.full((ny, nx), np.int16(args.fill_value), dtype=np.int16)

                ims_slice = np.asarray(ims_var[t, :, :], dtype=np.float64)
                if transpose_input_spatial:
                    ims_slice = ims_slice.T
                src_vals = ims_slice[src_j, src_i]

                ok = np.isfinite(src_vals)
                if np.any(ok):
                    out_slice[tgt_j[ok], tgt_i[ok]] = np.rint(src_vals[ok]).astype(np.int16)

                cat_out[t, :, :] = out_slice

                if (t + 1) % 30 == 0 or (t + 1) == n_time:
                    print(f"  wrote {t + 1}/{n_time} time slices")

            ds_out.title = "IMS categorical data on GEOSldas M36 grid (nearest-neighbor)"
            ds_out.history = f"{datetime.utcnow().isoformat()}Z created by regrid_ims_to_m36_nearest.py"
            ds_out.source_ims_file = str(ims_nc)
            ds_out.source_ims_var = str(args.ims_var)
            if time_source_name:
                ds_out.source_time_var = str(time_source_name)
            ds_out.source_tilecoord = str(tilecoord_path)
            ds_out.regrid_method = "nearest"
            ds_out.no_averaging = "true"
            ds_out.max_distance_km = float(args.max_distance_km)
            ds_out.fill_value_category = int(args.fill_value)
            ds_out.representative_tile_rule = "max frac_cell, tie-break min tile_id"
            ds_out.input_ims_dims = ",".join(str(d) for d in ims_dims)
            ds_out.input_ims_spatial_transposed = "true" if transpose_input_spatial else "false"

        finally:
            ds_out.close()

    print(f"Wrote output: {out_nc}")


def main() -> None:
    args = parse_args()
    tilecoord_path = Path(args.tilecoord)

    if not tilecoord_path.exists():
        raise FileNotFoundError(f"Tilecoord file not found: {tilecoord_path}")

    jobs = build_jobs(args)
    print(f"Prepared {len(jobs)} job(s).")
    print(f"Tilecoord: {tilecoord_path}")
    print(f"Max distance km: {args.max_distance_km}")
    if len(jobs) > 1:
        print(f"Batch years: {args.year_start}..{args.year_end}")

    repo_root = find_repo_root(Path(__file__).resolve())
    sys.path.insert(0, str(repo_root / "common/python/io"))
    from read_GEOSldas import read_tilecoord  # type: ignore

    tc = read_tilecoord(str(tilecoord_path))
    rep = choose_representative_tile_per_cell(tc)
    print(f"Tile count in tilecoord: {len(tc['tile_id'])}")
    print(f"Unique M36 cells with representative tile: {rep['rep_i'].size}")
    print(f"M36 grid shape from tilecoord: ny={int(rep['ny'])}, nx={int(rep['nx'])}")

    done = 0
    skipped = 0
    failed = 0

    for idx, (year, ims_nc, out_nc) in enumerate(jobs, start=1):
        tag = f"[{idx}/{len(jobs)}]"
        year_txt = f" year={year}" if year is not None else ""
        print(f"{tag}{year_txt} {ims_nc.name} -> {out_nc.name}")

        if not ims_nc.exists():
            msg = f"Input file not found: {ims_nc}"
            if args.skip_missing_input:
                print(f"  skip: {msg}")
                skipped += 1
                continue
            raise FileNotFoundError(msg)

        if out_nc.exists():
            if args.skip_existing_output:
                print(f"  skip: output already exists: {out_nc}")
                skipped += 1
                continue
            if not args.overwrite_output:
                raise FileExistsError(f"Output exists (use --overwrite-output or --skip-existing-output): {out_nc}")

        out_nc.parent.mkdir(parents=True, exist_ok=True)

        try:
            regrid_one_file(
                ims_nc=ims_nc,
                out_nc=out_nc,
                args=args,
                tilecoord_path=tilecoord_path,
                rep=rep,
            )
            done += 1
        except Exception as exc:
            failed += 1
            print(f"  failed: {exc}")
            if not args.continue_on_error:
                raise

    print("Done.")
    print(f"Summary: completed={done}, skipped={skipped}, failed={failed}, total={len(jobs)}")


if __name__ == "__main__":
    main()
