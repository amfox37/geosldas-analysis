#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import xarray as xr

FILL_MASK_THRESHOLD = 1e10
DEFAULT_READ_ENGINE  = "h5netcdf"   # we'll fallback to netcdf4 if needed
DEFAULT_WRITE_ENGINE = "h5netcdf"   # you can switch to netcdf4 if you prefer

def collect_monthly_files(root_dir: str, file_prefix: str, start_year: int, end_year: int):
    files, dates = [], []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            fn = f"{file_prefix}.tavg24_1d_lnd_Nt.monthly.{year:04d}{month:02d}.nc4"
            fpath = os.path.join(root_dir, f"Y{year:04d}", f"M{month:02d}", fn)
            if os.path.exists(fpath):
                files.append(fpath)
                dates.append(np.datetime64(f"{year:04d}-{month:02d}-01"))
    return files, np.asarray(dates, dtype="datetime64[ns]")

def _open_first_for_latlon(path, engine):
    # read lat/lon as plain arrays
    with xr.open_dataset(path, engine=engine) as ds0:
        if "lat" not in ds0 or "lon" not in ds0:
            raise KeyError("Expected 'lat' and 'lon' coordinates in input files.")
        lat = ds0["lat"].values
        lon = ds0["lon"].values
    return lat, lon

def _safe_open_mfdataset(files, preprocess, engine):
    """
    Open many files robustly. Try the requested engine first (parallel=False),
    then fall back to netcdf4 if the first attempt hits an HDF5/h5netcdf quirk.
    """
    try:
        return xr.open_mfdataset(
            files,
            combine="nested",
            concat_dim="time",
            preprocess=preprocess,
            engine=engine,
            parallel=False,   # <- safer for h5netcdf; avoids invalid file identifier issues
        )
    except Exception as e:
        if engine != "netcdf4":
            print(f"[warn] open_mfdataset failed with engine={engine}: {e}")
            print("[info] Falling back to engine='netcdf4' for reading …")
            return xr.open_mfdataset(
                files,
                combine="nested",
                concat_dim="time",
                preprocess=preprocess,
                engine="netcdf4",
                parallel=False,
            )
        raise

def build_lsm_dataset(
    root_dir: str,
    file_prefix: str,
    varnames: list[str],
    start_year: int = 2000,
    end_year: int = 2024,
    read_engine: str = DEFAULT_READ_ENGINE,
) -> xr.Dataset:
    """Open monthly tavg24_1d_lnd_Nt files and return a Dataset with requested variables."""
    files, dates = collect_monthly_files(root_dir, file_prefix, start_year, end_year)
    if not files:
        raise FileNotFoundError(
            f"No monthly files found under {root_dir} for prefix {file_prefix}"
        )
    print(f"[{file_prefix}] Found {len(files)} monthly files.")

    # lat/lon from the first file using the same engine (fall back if needed)
    try:
        lat, lon = _open_first_for_latlon(files[0], read_engine)
    except Exception as e:
        if read_engine != "netcdf4":
            print(f"[warn] Failed reading lat/lon with engine={read_engine}: {e}")
            print("[info] Retrying first file with engine='netcdf4' …")
            lat, lon = _open_first_for_latlon(files[0], "netcdf4")
        else:
            raise

    def _pre(ds):
        keep = [v for v in varnames if v in ds.variables]
        if not keep:
            raise KeyError(f"None of requested variables present. Requested: {varnames}")
        return ds[keep]

    ds = _safe_open_mfdataset(files, preprocess=_pre, engine=read_engine)

    ds = ds.assign_coords({
        "time": ("time", dates),
        "lat":  ("tile", lat),
        "lon":  ("tile", lon),
    })

    for v in varnames:
        if v in ds:
            ds[v] = ds[v].where(ds[v] < FILL_MASK_THRESHOLD)

    ds.attrs.update(
        source_root=root_dir,
        file_prefix=file_prefix,
        note="Monthly means from tavg24_1d_lnd_Nt.*; lat/lon attached on tile."
    )
    return ds

def write_compressed(ds: xr.Dataset, out_path: str, write_engine: str = DEFAULT_WRITE_ENGINE):
    # pick reasonable chunksizes aligned to dims (only for 2D time/tile vars)
    tlen = int(ds.sizes.get("time", 1))
    ntiles = int(ds.sizes.get("tile", 1))
    time_chunk = min(tlen, 12)
    tile_chunk = min(ntiles, 100000)

    comp = dict(zlib=True, complevel=4)
    encoding = {}
    for v in ds.data_vars:
        dims = ds[v].dims
        if dims == ("time", "tile"):
            encoding[v] = {**comp, "chunksizes": (time_chunk, tile_chunk)}
        # leave 1D vars alone (no chunksizes)

    print(f"Writing {out_path} with {write_engine}, chunks (time,tile)=({time_chunk},{tile_chunk}) …")
    ds.to_netcdf(out_path, engine=write_engine, encoding=encoding)
    print(f"Wrote: {out_path}")

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Extract monthly state variables from OL/DA tavg24_1d_lnd_Nt.* files and save compressed NetCDFs."
    )
    # Hardwired defaults: your paths/prefixes/vars
    p.add_argument("--ol-root", default="/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg")
    p.add_argument("--ol-prefix", default="LS_OLv8_M36")
    p.add_argument("--da-root", default="/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg")
    p.add_argument("--da-prefix", default="LS_DAv8_M36")
    p.add_argument("--vars", nargs="+",
                   default=["SFMC", "RZMC", "PRECTOTCORRLAND", "FRLANDSNO", "TSOIL1", "SNOMASLAND"])
    p.add_argument("--start-year", type=int, default=2000)
    p.add_argument("--end-year", type=int, default=2024)
    p.add_argument("--read-engine", choices=["h5netcdf", "netcdf4"], default=DEFAULT_READ_ENGINE)
    p.add_argument("--write-engine", choices=["h5netcdf", "netcdf4"], default=DEFAULT_WRITE_ENGINE)
    p.add_argument("--out-ol", default="OLv8_land_variables_2000_2024_compressed.nc")
    p.add_argument("--out-da", default="DAv8_land_variables_2000_2024_compressed.nc")
    return p.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)

    print("=== Building OL dataset ===")
    ds_ol = build_lsm_dataset(
        root_dir=args.ol_root,
        file_prefix=args.ol_prefix,
        varnames=args.vars,
        start_year=args.start_year,
        end_year=args.end_year,
        read_engine=args.read_engine,
    )
    write_compressed(ds_ol, args.out_ol, write_engine=args.write_engine)

    print("=== Building DA dataset ===")
    ds_da = build_lsm_dataset(
        root_dir=args.da_root,
        file_prefix=args.da_prefix,
        varnames=args.vars,
        start_year=args.start_year,
        end_year=args.end_year,
        read_engine=args.read_engine,
    )
    write_compressed(ds_da, args.out_da, write_engine=args.write_engine)

    print("Done.")
    print("OL vars:", list(ds_ol.data_vars))
    print("DA vars:", list(ds_da.data_vars))
    print("Dims:", dict(ds_da.dims))

if __name__ == "__main__":
    sys.exit(main())

