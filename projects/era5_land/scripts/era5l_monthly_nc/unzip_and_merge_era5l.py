#!/usr/bin/env python3
import os, re, sys, glob, shutil, zipfile, tempfile, datetime as dt
from pathlib import Path

import numpy as np
import xarray as xr

# ---------- settings ----------
IN_DIR  = Path(".")                        # folder with era5l_monthly_YYYY.nc (ZIPs)
OUT_ALL = Path("ERA5L_monthly_merged.nc") # final merged NetCDF
KEEP_UNZIPPED = True                       # keep era5l_monthly_YYYY_unzipped.nc files

# ---------- helpers ----------
def magic_bytes(path, n=8):
    with open(path, "rb") as f:
        return f.read(n)

def looks_like_zip(path):
    m = magic_bytes(path, 4)
    return m == b"PK\x03\x04"  # ZIP local file header

def looks_like_netcdf(path):
    m = magic_bytes(path)
    return m.startswith(b"CDF") or m.startswith(b"\x89HDF")  # classic / netCDF4

def looks_like_grib(path):
    return magic_bytes(path, 4) == b"GRIB"

def year_from_name(p: Path):
    m = re.search(r"(\d{4})", p.name)
    return m.group(1) if m else None

def unzip_inner_nc(zip_path: Path, out_path: Path) -> Path:
    """Extract inner NetCDF from CDS ZIP (named data_stream-moda.nc) and write to out_path."""
    with zipfile.ZipFile(zip_path, "r") as z:
        # If you ever get multiple files, pick the first .nc
        nc_members = [n for n in z.namelist() if n.lower().endswith(".nc")]
        if not nc_members:
            raise RuntimeError(f"No .nc member inside ZIP: {zip_path}")
        inner = nc_members[0]
        with z.open(inner) as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return out_path

def ensure_netcdf(nc_path: Path) -> Path:
    """If file is GRIB, convert to NetCDF; if NetCDF, return as-is; else raise."""
    if looks_like_netcdf(nc_path):
        return nc_path
    if looks_like_grib(nc_path):
        # convert with cfgrib → netcdf via xarray
        print(f"[info] {nc_path.name}: GRIB detected, converting to NetCDF via cfgrib ...")
        ds = xr.open_dataset(nc_path, engine="cfgrib")
        tmp = nc_path.with_suffix(".converted.nc")
        ds.to_netcdf(tmp)
        ds.close()
        nc_path.unlink()
        tmp.rename(nc_path)
        return nc_path
    # else: maybe an HTML/JSON error slipped into the ZIP (rare)
    raise RuntimeError(f"{nc_path} is not NetCDF/GRIB (first bytes: {magic_bytes(nc_path,16)!r})")

# ---------- main ----------
def main():
    zips = sorted(IN_DIR.glob("era5l_monthly_*.nc"))
    if not zips:
        print("No era5l_monthly_*.nc files found (remember these are ZIPs).")
        sys.exit(1)

    unzipped_paths = []
    for z in zips:
        y = year_from_name(z) or "unknown"
        out = IN_DIR / f"era5l_monthly_{y}_unzipped.nc"
        if looks_like_zip(z):
            print(f"[unzip] {z.name} -> {out.name}")
            unzip_inner_nc(z, out)
        else:
            # If someone already replaced the ZIP with the inner file, just copy/rename
            print(f"[warn] {z.name} does not look like ZIP; copying to {out.name}")
            shutil.copy2(z, out)
        # Validate/convert to proper NetCDF
        ensure_netcdf(out)
        unzipped_paths.append(out)

    # Merge all years
    paths = sorted(map(str, unzipped_paths))
    print("[merge] Concatenating by coordinates:", paths)
    # compat='override' tolerates minor attr differences across yearly files
    ds = xr.open_mfdataset(
        paths,
        combine="by_coords",
        coords="minimal",      # <-- add this
        compat="override",     # keep this
        join="exact"           # optional: stricter alignment
    )


    # Ensure time monotonic & unique; drop accidental dups
    ds = ds.sortby("valid_time")
    _, uniq_idx = np.unique(ds["valid_time"].values, return_index=True)
    if len(uniq_idx) < ds.sizes["valid_time"]:
        print(f"[info] Dropping {ds.sizes['valid_time'] - len(uniq_idx)} duplicate timesteps")
        ds = ds.isel(time=np.sort(uniq_idx))

    # Reasonable compression
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    for c in ("latitude", "longitude", "time"):
        if c in ds.coords:
            enc[c] = {"zlib": True, "complevel": 4}

    # Write final file
    ds.to_netcdf(OUT_ALL, encoding=enc)
    ds.close()
    print(f"[done] Wrote {OUT_ALL} with {OUT_ALL.stat().st_size/1e6:.1f} MB")

    if not KEEP_UNZIPPED:
        for p in unzipped_paths:
            try: p.unlink()
            except Exception: pass

if __name__ == "__main__":
    main()
