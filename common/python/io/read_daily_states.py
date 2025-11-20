#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Daily OL/DA pipeline (speed-optimized):
- Recursively discovers daily *_1200z.nc4 files for a prefix
- Opens lazily with Dask using a fast open_mfdataset configuration
- Applies frozen/snow masking (TSOIL1 < TEMP or FRLANDSNO > EPS)
- Enforces per-tile minimum valid-day fraction in BOTH experiments
- Builds smoothed DOY climatology from OL (CNTL)
- Computes daily anomalies for OL and DA with the SAME climatology
- Writes outputs with compression, chunking, and visible progress bars
- Logs sizes, chunks, kept tiles/days, and timing of open_mfdataset

Usage (Jupyter):
    %run /path/to/daily_da_pipeline_fast.py \
      --ol-root /discover/nobackup/.../LS_OLv8_M36/.../ens_avg \
      --ol-prefix LS_OLv8_M36 \
      --da-root /discover/nobackup/.../LS_DAv8_M36/.../ens_avg \
      --da-prefix LS_DAv8_M36 \
      --vars SFMC RZMC PRECTOTCORRLAND FRLANDSNO TSOIL1 \
      --anom-vars SFMC RZMC \
      --start-date 2000-01-01 --end-date 2024-12-31 \
      --chunks time:31,tile:16384 \
      --write-engine h5netcdf \
      --verbose 1
"""

import os
import sys
import re
import glob
import time
import argparse
import logging
import numpy as np
import xarray as xr

# --- Dask setup (threads typically best for NetCDF/HDF5 I/O) ---
import dask
from dask.diagnostics import ProgressBar
dask.config.set(scheduler="threads")

# Progress helper for discovery loops (tolerate missing widgets)
try:
    # Force text progress even in notebooks to avoid ipywidgets warnings
    os.environ.setdefault("TQDM_NOTEBOOK", "0")
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else []

# ------------------------
# Tunables / defaults
# ------------------------
FILL_MASK_THRESHOLD = 1e10
DEFAULT_READ_ENGINE  = "netcdf4"   # robust and fast for many small HDF5 files
DEFAULT_WRITE_ENGINE = "h5netcdf"

# Default chunking: ~1 month per chunk in time; 16k tiles per block
_DEFAULT_CHUNKS = {"time": 31, "tile": 16384}

# Mask thresholds (replicating user's choices)
TEMP_THRESH_K   = 275.15   # 2 °C
SNOW_EPS        = 1e-2     # 1% snow cover
MIN_VALID_FRAC  = 0.7      # require >=70% valid days per gridcell

# File name pattern: ...YYYYMMDD_1200z.nc4
_TS_RE = re.compile(r"\.(\d{8})_1200z\.nc4$")

# ------------------------
# Logging helpers
# ------------------------
def setup_logger(verbosity: int = 1) -> logging.Logger:
    level = logging.INFO if verbosity == 1 else (logging.DEBUG if verbosity >= 2 else logging.WARNING)
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
        force=True,
    )
    return logging.getLogger("daily-da")

def stamp(log: logging.Logger, msg: str):
    log.info(msg)

# ------------------------
# Discovery / parsing
# ------------------------
def _parse_ts(bname: str) -> np.datetime64 | None:
    m = _TS_RE.search(bname)
    if not m:
        return None
    ymd = m.group(1)
    return np.datetime64(f"{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}T12:00")

def collect_daily_files(root_dir: str, file_prefix: str,
                        start_date: np.datetime64, end_date: np.datetime64,
                        log: logging.Logger | None = None):
    """
    Recursively glob all *_1200z.nc4 files, parse dates, filter to [start_date, end_date], sort.
    Shows a tqdm bar during the scan.
    """
    pattern = os.path.join(root_dir, "**", f"{file_prefix}.tavg24_1d_lnd_Nt.*_1200z.nc4")
    hits = glob.glob(pattern, recursive=True)
    if not hits:
        raise FileNotFoundError(f"No daily *_1200z.nc4 under {root_dir} for {file_prefix}")

    start = np.datetime64(str(start_date), "ns")
    end   = np.datetime64(str(end_date), "ns")

    files, times = [], []
    it = tqdm(hits, desc=f"[{file_prefix}] scanning files", unit="file")
    for p in it:
        ts = _parse_ts(os.path.basename(p))
        if ts is not None and start <= ts <= end:
            files.append(p)
            times.append(ts)

    if not files:
        raise FileNotFoundError(f"No daily files for {file_prefix} within [{start_date} .. {end_date}]")

    order = np.argsort(np.asarray(times))
    files = [files[i] for i in order]
    times = np.asarray(times, dtype="datetime64[ns]")[order]

    if log:
        stamp(log, f"[{file_prefix}] in-range files: {len(files)} / hits: {len(hits)}")
        stamp(log, f"[{file_prefix}] time span: {str(times[0])} … {str(times[-1])}")

    return files, times

# ------------------------
# Fast open_mfdataset
# ------------------------
def speed_open_mfdataset(files, varnames, engine="netcdf4", chunks={"time":31, "tile":16384}, log: logging.Logger | None = None):
    """
    Fast open for thousands of GEOS daily files:
    - minimal metadata coordination
    - no CF decoding (we assign 'time' from filenames later)
    - thread-safe HDF5 access (lock=True)
    """
    want = set(varnames)
    def _pre(ds):
        keep = [v for v in ds.variables if v in want]
        if not keep:
            raise KeyError(f"None of requested variables present. Requested: {varnames}")
        return ds[keep]

    if log:
        stamp(log, f"[open] calling open_mfdataset on {len(files)} files …")
    t0 = time.perf_counter()
    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        preprocess=_pre,
        engine=engine,
        # -------- performance-critical knobs --------
        parallel=False,            # avoid HDF5 handle churn; Dask still parallelizes tasks
        lock=True,                 # safe on shared FS
        chunks=chunks,             # lazy chunking
        data_vars="minimal",       # avoid global metadata merging
        coords="minimal",
        compat="override",
        # Skip CF decoding to reduce per-file header work; we reassign 'time'
        mask_and_scale=False,
        decode_times=False,
        decode_coords=False,
        use_cftime=False,
    )
    t1 = time.perf_counter()
    if log:
        stamp(log, f"[open] done in {t1 - t0:.1f}s | chunks={chunks} | engine={engine}")
    return ds

# ------------------------
# Opening helpers
# ------------------------
def _open_first_for_latlon(path: str, engine: str):
    with xr.open_dataset(path, engine=engine, chunks={}) as ds0:
        if "lat" not in ds0 or "lon" not in ds0:
            raise KeyError("Expected 'lat' and 'lon' in input files.")
        lat = ds0["lat"].values
        lon = ds0["lon"].values
    return lat, lon

def build_lsm_dataset_daily(
    root_dir: str,
    file_prefix: str,
    varnames: list[str],
    start_date: str = "2000-01-01",
    end_date: str   = "2024-12-31",
    read_engine: str = DEFAULT_READ_ENGINE,
    chunks: dict | None = None,
    log: logging.Logger | None = None,
) -> xr.Dataset:
    files, times = collect_daily_files(root_dir, file_prefix,
                                       np.datetime64(start_date),
                                       np.datetime64(end_date),
                                       log=log)
    if chunks is None:
        chunks = _DEFAULT_CHUNKS

    # lat/lon from the first file
    try:
        lat, lon = _open_first_for_latlon(files[0], read_engine)
    except Exception as e:
        if read_engine != "netcdf4":
            if log: stamp(log, f"[warn] lat/lon read failed with {read_engine}: {e} → retry netcdf4")
            lat, lon = _open_first_for_latlon(files[0], "netcdf4")
        else:
            raise

    # Fast open
    ds = speed_open_mfdataset(files, varnames, engine=read_engine, chunks=chunks, log=log)

    # Attach coordinates
    ds = ds.assign_coords({
        "time": ("time", times),
        "lat":  ("tile", lat),
        "lon":  ("tile", lon),
    })

    # Mask giant fill values lazily
    for v in varnames:
        if v in ds:
            ds[v] = ds[v].where(ds[v] < FILL_MASK_THRESHOLD)

    if log:
        stamp(log, f"[{file_prefix}] dims={dict(ds.dims)}")
        try:
            stamp(log, f"[{file_prefix}] chunks={ds.chunks}")
        except Exception:
            pass
    return ds

# ------------------------
# Masking utilities
# ------------------------
def apply_frozen_snow_mask(sm_da: xr.DataArray, tsoil: xr.DataArray, frsnow: xr.DataArray,
                           temp_thresh: float = TEMP_THRESH_K, snow_eps: float = SNOW_EPS) -> xr.DataArray:
    mask = (tsoil < temp_thresh) | (frsnow > snow_eps)
    return sm_da.where(~mask)

def per_tile_min_valid_mask(da: xr.DataArray, min_valid_frac: float = MIN_VALID_FRAC) -> xr.DataArray:
    """Return boolean over 'tile' where the fraction of valid (non-NaN) days ≥ threshold."""
    frac = xr.where(xr.apply_ufunc(np.isfinite, da), 1.0, 0.0).mean("time")
    return (frac >= min_valid_frac)

# ------------------------
# Climatology + anomalies
# ------------------------
def _cyclic_rolling_mean(da: xr.DataArray, window: int) -> xr.DataArray:
    half = window // 2
    pad = xr.concat([da.isel(dayofyear=slice(-half, None)), da, da.isel(dayofyear=slice(0, half))],
                    dim="dayofyear")
    sm  = pad.rolling(dayofyear=window, center=True, min_periods=1).mean()
    return sm.isel(dayofyear=slice(half, -half))

def build_cntl_climatology_doy_smooth(ds_cntl: xr.Dataset, clim_vars: list[str], window: int = 31) -> xr.Dataset:
    doy = ds_cntl.time.dt.dayofyear
    clim = xr.Dataset()
    for v in clim_vars:
        cr = ds_cntl[v].groupby(doy).mean("time", skipna=True)
        # Ensure 1..366; map 366 -> 365 before smoothing
        if 366 not in cr["dayofyear"]:
            full = xr.DataArray(np.arange(1, 367), dims="dayofyear", name="dayofyear")
            cr = cr.reindex(dayofyear=full)
        cr = cr.where(cr["dayofyear"] != 366, cr.sel(dayofyear=365))
        clim[f"{v}_clim"] = _cyclic_rolling_mean(cr, window=window)
    clim = clim.assign_coords(dayofyear=clim["dayofyear"], tile=ds_cntl["tile"])
    clim.attrs.update(baseline="CNTL masked daily DOY mean", smoothing=f"cyclic {window}-day")
    return clim

def compute_daily_anomalies_from_clim(ds: xr.Dataset, clim: xr.Dataset, vars_: list[str]) -> xr.Dataset:
    doy = xr.where(ds.time.dt.dayofyear == 366, 365, ds.time.dt.dayofyear)
    out = xr.Dataset()
    for v in vars_:
        base = clim[f"{v}_clim"].sel(dayofyear=doy).transpose("time", "tile")
        out[v] = ds[v] - base
    return out.assign_coords(time=ds["time"], tile=ds["tile"], lat=ds["lat"], lon=ds["lon"])

# ------------------------
# Writing + logging helpers
# ------------------------
def _estimate_nbytes(ds: xr.Dataset) -> int:
    n = 0
    for v in ds.data_vars:
        da = ds[v]
        if da.dtype.kind in "fcui":
            try:
                n += da.data.nbytes
            except Exception:
                n += int(da.size * da.dtype.itemsize)
    return n

def write_compressed(ds: xr.Dataset, out_path: str, write_engine: str = DEFAULT_WRITE_ENGINE,
                     chunks: dict | None = None, log: logging.Logger | None = None):
    if chunks is None:
        chunks = _DEFAULT_CHUNKS
    tlen = int(ds.sizes.get("time", 1))
    ntiles = int(ds.sizes.get("tile", 1))
    time_chunk = min(tlen, int(chunks.get("time", 31)))
    tile_chunk = min(ntiles, int(chunks.get("tile", 16384)))

    comp = dict(zlib=True, complevel=4)
    encoding = {}
    for v in ds.data_vars:
        if ds[v].dims == ("time", "tile"):
            encoding[v] = {**comp, "chunksizes": (time_chunk, tile_chunk)}

    if log:
        est = _estimate_nbytes(ds)
        stamp(log, f"→ Writing {out_path} ({write_engine}); chunks (time,tile)=({time_chunk},{tile_chunk}); est raw ~ {est/1e9:.2f} GB")

    delayed = ds.to_netcdf(out_path, engine=write_engine, encoding=encoding, compute=False)
    with ProgressBar():
        dask.compute(delayed)
    if log:
        stamp(log, f"✓ Wrote {out_path}")

# ------------------------
# CLI parsing
# ------------------------
def parse_chunk_flag(s: str) -> dict:
    out = {}
    if not s:
        return out
    for kv in s.split(","):
        if not kv.strip():
            continue
        k, v = kv.split(":")
        out[k.strip()] = int(v.strip())
    return out

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Daily OL/DA with Dask (fast open); frozen/snow masks; CNTL DOY climatology; daily anomalies; min-valid filtering; progress/logging."
    )
    # Roots/prefixes
    p.add_argument("--ol-root", default="/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg")
    p.add_argument("--ol-prefix", default="LS_OLv8_M36")
    p.add_argument("--da-root", default="/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg")
    p.add_argument("--da-prefix", default="LS_DAv8_M36")

    # Variables: which to read; which to use for climatology/anomalies
    p.add_argument("--vars", nargs="+", default=["SFMC", "RZMC", "PRECTOTCORRLAND", "FRLANDSNO", "TSOIL1"])
    p.add_argument("--anom-vars", nargs="+", default=["SFMC", "RZMC"], help="Variables to build climatology/anomalies for")

    # Time, I/O, chunking
    p.add_argument("--start-date", default="2000-01-01")
    p.add_argument("--end-date",   default="2024-12-31")
    p.add_argument("--read-engine", choices=["h5netcdf", "netcdf4"], default=DEFAULT_READ_ENGINE)
    p.add_argument("--write-engine", choices=["h5netcdf", "netcdf4"], default=DEFAULT_WRITE_ENGINE)
    p.add_argument("--chunks", default="time:31,tile:16384", help='Chunk spec, e.g., "time:31,tile:16384"')

    # Climatology smoothing + masking thresholds
    p.add_argument("--clim-window", type=int, default=31, help="Cyclic smoothing window (days)")
    p.add_argument("--temp-K", type=float, default=TEMP_THRESH_K, help="Frozen mask: TSOIL1 < TEMP_K")
    p.add_argument("--snow-eps", type=float, default=SNOW_EPS, help="Snow mask: FRLANDSNO > SNOW_EPS")
    p.add_argument("--min-valid-frac", type=float, default=MIN_VALID_FRAC,
                   help="Require ≥ fraction of valid days per tile in BOTH expts")

    # Outputs
    p.add_argument("--out-clim", default="OLv8_daily_climatology_DOY_smooth_masked.nc")
    p.add_argument("--out-ol-anom", default="OLv8_daily_anomalies_masked.nc")
    p.add_argument("--out-da-anom", default="DAv8_daily_anomalies_masked.nc")

    # Verbosity
    p.add_argument("--verbose", type=int, default=1, help="0=warn, 1=info, 2=debug")

    return p.parse_args(argv)

# ------------------------
# Main pipeline
# ------------------------
def main(argv=None):
    args = parse_args(argv)
    log = setup_logger(args.verbose)
    chunks = parse_chunk_flag(args.chunks)

    # Environment niceties (safe defaults on many HPC FS)
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    stamp(log, "Starting daily DA pipeline")
    stamp(log, f"OL prefix/root: {args.ol_prefix} | {args.ol_root}")
    stamp(log, f"DA prefix/root: {args.da_prefix} | {args.da_root}")
    stamp(log, f"Time range: {args.start_date} .. {args.end_date}")
    stamp(log, f"Vars (read): {args.vars}")
    stamp(log, f"Vars (anom): {args.anom_vars}")
    stamp(log, f"Chunks: {chunks}")

    # === Open OL & DA daily stacks ===
    stamp(log, "Opening OL dataset …")
    ds_ol = build_lsm_dataset_daily(
        root_dir=args.ol_root,
        file_prefix=args.ol_prefix,
        varnames=args.vars,
        start_date=args.start_date,
        end_date=args.end_date,
        read_engine=args.read_engine,
        chunks=chunks,
        log=log,
    )
    stamp(log, "Opening DA dataset …")
    ds_da = build_lsm_dataset_daily(
        root_dir=args.da_root,
        file_prefix=args.da_prefix,
        varnames=args.vars,
        start_date=args.start_date,
        end_date=args.end_date,
        read_engine=args.read_engine,
        chunks=chunks,
        log=log,
    )

    # Log dims/chunks
    stamp(log, f"OL dims: {dict(ds_ol.dims)}")
    stamp(log, f"DA dims: {dict(ds_da.dims)}")
    try:
        stamp(log, f"OL chunks: {ds_ol.chunks}")
        stamp(log, f"DA chunks: {ds_da.chunks}")
    except Exception:
        pass

    # --- Apply frozen/snow masks to anomaly vars (OL & DA) ---
    stamp(log, f"Mask thresholds: TSOIL1<{args.temp_K}K OR FRLANDSNO>{args.snow_eps}; per-tile min valid frac={args.min_valid_frac:.0%}")

    def mask_vars(ds: xr.Dataset) -> xr.Dataset:
        out = xr.Dataset(coords=ds.coords)
        for v in args.anom_vars:
            sm = ds[v]
            tsoil = ds["TSOIL1"] if "TSOIL1" in ds else xr.full_like(sm, np.nan)
            frsnow = ds["FRLANDSNO"] if "FRLANDSNO" in ds else xr.full_like(sm, np.nan)
            out[v] = apply_frozen_snow_mask(sm, tsoil, frsnow,
                                            temp_thresh=args.temp_K, snow_eps=args.snow_eps)
        return out

    stamp(log, "Applying frozen/snow mask to OL anomaly vars …")
    ol_masked = mask_vars(ds_ol)
    stamp(log, "Applying frozen/snow mask to DA anomaly vars …")
    da_masked = mask_vars(ds_da)

    if args.verbose >= 2:
        stamp(log, "Persist masked OL/DA vars (optional, speeds up later steps)…")
        with ProgressBar():
            ol_masked = ol_masked.persist()
            da_masked = da_masked.persist()

    # --- Enforce per-tile minimum valid-day fraction (separately), then intersect ---
    stamp(log, "Computing per-tile valid-day fractions …")
    keep_tile_ol = xr.ones_like(ol_masked[args.anom_vars[0]].isel(time=0), dtype=bool)
    keep_tile_da = xr.ones_like(da_masked[args.anom_vars[0]].isel(time=0), dtype=bool)
    for v in args.anom_vars:
        keep_tile_ol = keep_tile_ol & per_tile_min_valid_mask(ol_masked[v], args.min_valid_frac)
        keep_tile_da = keep_tile_da & per_tile_min_valid_mask(da_masked[v], args.min_valid_frac)
    keep_tile = keep_tile_ol & keep_tile_da

    n_tiles_total = int(keep_tile.size)
    n_tiles_kept  = int(keep_tile.sum().compute())
    stamp(log, f"Tiles kept (≥{args.min_valid_frac:.0%} valid days in BOTH expts): {n_tiles_kept}/{n_tiles_total}")

    ol_masked = ol_masked.where(keep_tile, np.nan)
    da_masked = da_masked.where(keep_tile, np.nan)

    # --- Build OL climatology (from masked OL) ---
    stamp(log, "Building CNTL (OL) DOY climatology with cyclic smoothing …")
    clim = build_cntl_climatology_doy_smooth(ol_masked, args.anom_vars, window=args.clim_window)
    write_compressed(clim, args.out_clim, write_engine=args.write_engine,
                     chunks={"dayofyear": 366, "tile": _DEFAULT_CHUNKS["tile"]}, log=log)

    # --- Compute anomalies for OL & DA using OL climatology ---
    stamp(log, "Computing daily anomalies (OL) …")
    anom_ol = compute_daily_anomalies_from_clim(
        ol_masked.assign_coords(lat=ds_ol["lat"], lon=ds_ol["lon"]),
        clim, args.anom_vars
    )
    stamp(log, "Computing daily anomalies (DA) …")
    anom_da = compute_daily_anomalies_from_clim(
        da_masked.assign_coords(lat=ds_da["lat"], lon=ds_da["lon"]),
        clim, args.anom_vars
    )

    # --- Align time and enforce strict common day validity across vars & expts ---
    stamp(log, "Aligning OL/DA anomaly times & enforcing common-day validity …")
    anom_ol, anom_da = xr.align(anom_ol, anom_da, join="inner")
    valid_day = xr.full_like(anom_ol["time"], True, dtype=bool)
    for v in args.anom_vars:
        valid_day = valid_day & xr.apply_ufunc(np.isfinite, anom_ol[v]).all("tile")
        valid_day = valid_day & xr.apply_ufunc(np.isfinite, anom_da[v]).all("tile")
    n_days_total = anom_ol.sizes["time"]
    n_days_kept  = int(valid_day.sum().compute())
    stamp(log, f"Days kept after strict validity: {n_days_kept}/{n_days_total}")
    anom_ol = anom_ol.where(valid_day)
    anom_da = anom_da.where(valid_day)

    # --- Write anomalies ---
    stamp(log, "Writing daily anomaly files …")
    write_compressed(anom_ol, args.out_ol_anom, write_engine=args.write_engine, chunks=chunks, log=log)
    write_compressed(anom_da, args.out_da_anom, write_engine=args.write_engine, chunks=chunks, log=log)

    stamp(log, "Done.")
    stamp(log, f"Anomaly vars: {list(anom_da.data_vars)}")
    stamp(log, f"Dims (DA anomalies): {dict(anom_da.dims)}")

if __name__ == "__main__":
    sys.exit(main())
