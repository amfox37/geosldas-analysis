# === User parameters ===

# Directory or glob for input LDAS lnd files (edit this)
INPUT_GLOB = "/discover/nobackup/projects/land_da/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/Y????/M??/OLv8_M36_cd.tavg24_1d_lnd*.*"  # e.g., "/discover/nobackup/.../lnd/Y*/M*/*.nc4"

# Output folder (created if it doesn't exist)
OUT_DIR = "/discover/nobackup/projects/land_da/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg"

# Time subset (inclusive). Set to None to skip filtering.
START_DATE = "2018-08-01"   # or None
END_DATE   = "2024-06-30"   # or None

# Variables to extract and process
VARS = [
    "GRN", "LAI",
    "GWETPROF", "GWETROOT", "GWETTOP",
    "PRMC", "RZMC", "SFMC",
    "PRECTOTCORRLAND", "QINFILLAND",
    "SHLAND", "LHLAND", "EVLAND",
]

# Chunking for xarray/dask (None = let xarray decide; or provide dict like {"time": 64, "tile": 4096})
CHUNKS = {"time": 64, "tile": 32768}

# Compression settings for NetCDF (set zlib=True for smaller files)
ENCODING_COMP = {"zlib": True, "complevel": 4, "shuffle": True}

# =======================

import os
from glob import glob
import numpy as np
import xarray as xr
import warnings
import re
from pathlib import Path
import pandas as pd

# Helpers for filename-based time filtering
DATE_RE = re.compile(r"\.(\d{8})_")   # captures YYYYMMDD in filename pieces like .20221222_
FV = 1e15

def clean(ds, name):
    """Mask out large fill values common in GEOS-LDAS outputs."""
    v = ds[name]
    return v.where(v < FV)

def collect_files(base: Path, start: pd.Timestamp, end: pd.Timestamp, file_glob: str):
    files = []
    for ydir in sorted(base.glob("Y*/")):
        for mdir in sorted(ydir.glob("M*/")):
            for f in mdir.glob(file_glob):
                m = DATE_RE.search(f.name)
                if not m:
                    continue
                dt = pd.to_datetime(m.group(1), format="%Y%m%d")
                if start <= dt <= end:
                    files.append(str(f))
    files.sort()
    return files

os.makedirs(OUT_DIR, exist_ok=True)

# prepare file list using filename dates if START_DATE/END_DATE provided
start_ts = pd.to_datetime(START_DATE) if START_DATE is not None else None
end_ts = pd.to_datetime(END_DATE) if END_DATE is not None else None

if start_ts is not None and end_ts is not None:
    # derive base dir and per-month filename glob from INPUT_GLOB
    s = INPUT_GLOB
    if "/Y" in s:
        base_str, rest = s.split("/Y", 1)
        base = Path(base_str)
        file_glob = Path(rest).name
    else:
        base = Path(s).parent
        file_glob = Path(s).name
    files = collect_files(base, start_ts, end_ts, file_glob)
else:
    files = sorted(glob(INPUT_GLOB))

if not files:
    raise FileNotFoundError(f"No input files found for pattern/time range: {INPUT_GLOB}")

# Open multi-file dataset (combine='nested' required). Avoid parallel=True which can trigger libnetcdf assertions.
try:
    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        parallel=False,        # avoid netCDF4 parallel internals that crash in some builds
        decode_times=True,
        chunks=CHUNKS,
        engine="netcdf4",
    )
except Exception as e_open:
    warnings.warn(f"open_mfdataset with netcdf4 failed: {e_open!s}. Falling back to manual open/concat.")
    # try manual safe open/concat, first with netcdf4, then with h5netcdf if needed
    def _try_manual(engine_name):
        try:
            ds_list = [xr.open_dataset(f, engine=engine_name, chunks=CHUNKS, decode_times=True) for f in files]
            return xr.concat(ds_list, dim="time", data_vars="minimal", coords="minimal")
        except Exception as e:
            raise
    try:
        ds = _try_manual("netcdf4")
    except Exception as e_manual_netcdf4:
        warnings.warn(f"Manual open with netcdf4 failed: {e_manual_netcdf4!s}. Retrying with h5netcdf.")
        try:
            ds = _try_manual("h5netcdf")
        except Exception as e_manual_h5:
            # re-raise a helpful error including original failure
            raise RuntimeError(
                "Failed to open files with nested combine via netcdf4 and h5netcdf. "
                f"Original open_mfdataset error: {e_open!s}. Manual netcdf4 error: {e_manual_netcdf4!s}. "
                f"h5netcdf error: {e_manual_h5!s}."
            )

# Ensure time is sorted before slicing/groupby
ds = ds.sortby("time")

# Keep only the variables of interest that are actually present
vars_present = [v for v in VARS if v in ds.data_vars]
if not vars_present:
    raise ValueError("None of the requested variables are present in the input files.")
ds = ds[vars_present]

# Mask using native _FillValue if present; fallback to 1e14 guard
for v in vars_present:
    da = ds[v]
    # preserve attributes by using DataArray.where
    fv = da.attrs.get("_FillValue", None)
    if fv is not None:
        da = da.where((da != fv) & np.isfinite(da))
    else:
        da = da.where(np.isfinite(da) & (np.abs(da) < 1e14))
    ds[v] = da

# Annual means (calendar-year)
annual_means = ds.groupby("time.year").mean("time", skipna=True, keep_attrs=True)
annual_means.attrs.update(ds.attrs)

# Monthly climatology (across all years)
monthly_climo = ds.groupby("time.month").mean("time", skipna=True, keep_attrs=True)
monthly_climo.attrs.update(ds.attrs)

# Grand mean & std across years (from annual_means)
grand_mean = annual_means.mean(dim="year", skipna=True, keep_attrs=True)
grand_std  = annual_means.std(dim="year", skipna=True, ddof=1, keep_attrs=True)

# Combine into one Dataset with clear suffixes
gm = xr.Dataset()
for v in annual_means.data_vars:
    gm[f"{v}_mean"] = grand_mean[v]
    gm[f"{v}_std"]  = grand_std[v]
gm.attrs.update(annual_means.attrs)

def build_encoding(ds_like):
    enc = {}
    for v in ds_like.data_vars:
        enc[v] = dict(ENCODING_COMP)
        enc[v]["dtype"] = "float32"
        enc[v]["_FillValue"] = np.float32(1.0e15)
    return enc

annual_path = os.path.join(OUT_DIR, "annual_means.nc")
climo_path  = os.path.join(OUT_DIR, "monthly_climatology.nc")
grand_path  = os.path.join(OUT_DIR, "grand_stats.nc")  # renamed

annual_means.astype("float32").to_netcdf(
    annual_path, format="NETCDF4", encoding=build_encoding(annual_means)
)
monthly_climo.astype("float32").to_netcdf(
    climo_path, format="NETCDF4", encoding=build_encoding(monthly_climo)
)
gm.astype("float32").to_netcdf(
    grand_path, format="NETCDF4", encoding=build_encoding(gm)  # <-- fixed
)

# close dataset to free resources (especially important when using dask/xarray)
ds.close()