#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import xarray as xr
import dask

# ---------------------
# Tunables / defaults
# ---------------------
FILL = 1e10  # mask huge fill values like 1e15
DEFAULT_TIME_CHUNK = 12
DEFAULT_TILE_CHUNK = 100_000

def _collect_monthly_files(root_dir, file_template, start_year=2000, end_year=2024):
    files, dates = [], []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            fn = file_template.format(YYYY=year, MM=month)
            fpath = os.path.join(root_dir, f"Y{year:04d}", f"M{month:02d}", fn)
            if os.path.exists(fpath):
                files.append(fpath)
                dates.append(np.datetime64(f"{year:04d}-{month:02d}-01"))
    return files, np.array(dates, dtype="datetime64[ns]")

def _load_lat_lon_from(first_file, engine="h5netcdf"):
    # Only to fetch lat/lon on tile
    with xr.open_dataset(first_file, engine=engine) as ds0:
        # Return plain 1-D numpy arrays to avoid coord/Index baggage
        lat = np.asarray(ds0["lat"].values)
        lon = np.asarray(ds0["lon"].values)
    return lat, lon

def build_fcst_ana_increments(
    root_dir,
    start_year=2000,
    end_year=2024,
    file_prefix="LS_DAv8_M36",
    engine="h5netcdf",
    chunks=None,
):
    """
    From LS_*inst3_1d_lndfcstana_Nt.monthly.YYYYMM.nc4:
      SFMC_INC = SFMC_ANA - SFMC_FCST (m3 m-3)
      RZMC_INC = RZMC_ANA - RZMC_FCST (m3 m-3)
    """
    tmpl = f"{file_prefix}.inst3_1d_lndfcstana_Nt.monthly.{{YYYY:04d}}{{MM:02d}}.nc4"
    files, dates = _collect_monthly_files(root_dir, tmpl, start_year, end_year)
    if not files:
        raise FileNotFoundError("No fcst/ana monthly files found.")
    lat, lon = _load_lat_lon_from(files[0], engine=engine)

    vars_needed = ["SFMC_FCST", "RZMC_FCST", "SFMC_ANA", "RZMC_ANA"]

    def _pre(ds):
        keep = [v for v in vars_needed if v in ds.variables]
        return ds[keep]

    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        preprocess=_pre,
        engine=engine,
        parallel=True,
        chunks=chunks,
    ).assign_coords({
        "time": ("time", dates),
        "lat": ("tile", lat),
        "lon": ("tile", lon),
    })

    # Mask fills lazily
    for v in vars_needed:
        if v in ds:
            ds[v] = ds[v].where(ds[v] < FILL)

    # increments as DataArrays, keep laziness
    sfmc_inc = (ds["SFMC_ANA"] - ds["SFMC_FCST"]).rename("SFMC_INC")
    rzmc_inc = (ds["RZMC_ANA"] - ds["RZMC_FCST"]).rename("RZMC_INC")

    # build a clean Dataset with numpy coords
    out = xr.Dataset(
        coords={
            "time": ("time", ds["time"].values),
            "tile": ("tile", ds["tile"].values),
            "lat":  ("tile", lat),   # numpy arrays returned by _load_lat_lon_from
            "lon":  ("tile", lon),
        },
    )

    # attach variables using raw array (.data keeps dask)
    out["SFMC_INC"] = (("time", "tile"), sfmc_inc.data)
    out["RZMC_INC"] = (("time", "tile"), rzmc_inc.data)

    # attrs
    out["SFMC_INC"].attrs.update(long_name="soil_moisture_surface_increment (ANA−FCST)", units="m3 m-3")
    out["RZMC_INC"].attrs.update(long_name="soil_moisture_rootzone_increment (ANA−FCST)", units="m3 m-3")

    return out

def build_snowmass_increment(
    root_dir,
    start_year=2000,
    end_year=2024,
    file_prefix="LS_DAv8_M36",
    engine="h5netcdf",
    chunks=None,
):
    """
    From LS_*catch_progn_incr.monthly.YYYYMM.nc4:
      SNOWMASS_INCR = WESNN1_INCR + WESNN2_INCR + WESNN3_INCR (kg m-2)
    """
    tmpl = f"{file_prefix}.catch_progn_incr.monthly.{{YYYY:04d}}{{MM:02d}}.nc4"
    files, dates = _collect_monthly_files(root_dir, tmpl, start_year, end_year)
    if not files:
        raise FileNotFoundError("No catch_progn_incr monthly files found.")
    lat, lon = _load_lat_lon_from(files[0], engine=engine)

    vars_needed = ["WESNN1_INCR", "WESNN2_INCR", "WESNN3_INCR"]

    def _pre(ds):
        keep = [v for v in vars_needed if v in ds.variables]
        return ds[keep]

    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        preprocess=_pre,
        engine=engine,
        parallel=True,
        chunks=chunks,
    ).assign_coords({
        "time": ("time", dates),
        "lat": ("tile", lat),
        "lon": ("tile", lon),
    })

    # Mask fills lazily
    for v in vars_needed:
        if v in ds:
            ds[v] = ds[v].where(ds[v] < FILL)

    snow_incr = (
        ds["WESNN1_INCR"].fillna(0.0)
    + ds["WESNN2_INCR"].fillna(0.0)
    + ds["WESNN3_INCR"].fillna(0.0)
    ).rename("SNOWMASS_INCR")

    out = xr.Dataset(
        coords={
            "time": ("time", ds["time"].values),
            "tile": ("tile", ds["tile"].values),
            "lat":  ("tile", lat),
            "lon":  ("tile", lon),
        },
    )
    out["SNOWMASS_INCR"] = (("time", "tile"), snow_incr.data)
    out["SNOWMASS_INCR"].attrs.update(long_name="increment_snow_mass_total (sum of layers)", units="kg m-2")

    return out

def build_monthly_increments(
    fcst_ana_root,
    snow_incr_root,
    start_year=2000,
    end_year=2024,
    file_prefix="LS_DAv8_M36",
    out_nc="LS_monthly_increments_2000_2024.nc",
    engine="h5netcdf",
    time_chunk=DEFAULT_TIME_CHUNK,
    tile_chunk=DEFAULT_TILE_CHUNK,
    compress=True,
    complevel=2,
    scheduler="processes",
    write_scheduler="threads",
):
    # Configure Dask
    dask.config.set(scheduler=scheduler)

    chunks = {"time": time_chunk, "tile": tile_chunk}

    # Build lazily
    ds_sm = build_fcst_ana_increments(
        fcst_ana_root, start_year, end_year, file_prefix=file_prefix, engine=engine, chunks=chunks
    )
    ds_sn = build_snowmass_increment(
        snow_incr_root, start_year, end_year, file_prefix=file_prefix, engine=engine, chunks=chunks
    )

    # Align and merge on (time, tile); keep lazy
    ds_sm, ds_sn = xr.align(ds_sm, ds_sn, join="outer")
    out = xr.merge([ds_sm, ds_sn])

    # Rechunk once for write
    out = out.chunk({"time": time_chunk, "tile": tile_chunk})

    # Encoding: align NetCDF chunksizes with Dask chunks; light/optional compression
    comp = {"zlib": bool(compress), "complevel": int(complevel)} if compress else {"zlib": False}
    encoding = {v: {**comp, "chunksizes": (time_chunk, tile_chunk)} for v in out.data_vars}

    # Write chunk-aligned NetCDF
    with dask.config.set(scheduler=write_scheduler):
        out.to_netcdf(out_nc, engine=engine, encoding=encoding)
    return out

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Build monthly ANA−FCST increments and snow mass increments (chunked NetCDF).")
    p.add_argument("--fcst-ana-root", default="/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
                   help="Root with inst3_1d_lndfcstana_Nt.monthly.YYYYMM.nc4 in YYYYY/MM folders")
    p.add_argument("--snow-incr-root", default="/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
                   help="Root with catch_progn_incr.monthly.YYYYMM.nc4 in YYYYY/MM folders")
    p.add_argument("--file-prefix", default="LS_DAv8_M36", help="Common file prefix")
    p.add_argument("--start-year", type=int, default=2000)
    p.add_argument("--end-year", type=int, default=2024)
    p.add_argument("--out-nc", default="LS_monthly_increments_2000_2024.nc")
    p.add_argument("--engine", choices=["h5netcdf", "netcdf4"], default="h5netcdf", help="NetCDF engine")
    p.add_argument("--time-chunk", type=int, default=DEFAULT_TIME_CHUNK, help="Dask/NetCDF chunk size for time")
    p.add_argument("--tile-chunk", type=int, default=DEFAULT_TILE_CHUNK, help="Dask/NetCDF chunk size for tile")
    p.add_argument("--no-compress", action="store_true", help="Disable zlib compression for faster writes")
    p.add_argument("--complevel", type=int, default=2, help="Compression level (1–9) if compression is enabled")
    p.add_argument("--scheduler", choices=["threads", "processes", "single-threaded"], default="processes", help="Dask scheduler")
    return p.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)
    print("Building monthly increments with args:")
    for k, v in vars(args).items():
        print(f"  {k}: {v}")

    out = build_monthly_increments(
        fcst_ana_root=args.fcst_ana_root,
        snow_incr_root=args.snow_incr_root,
        start_year=args.start_year,
        end_year=args.end_year,
        file_prefix=args.file_prefix,
        out_nc=args.out_nc,
        engine=args.engine,
        time_chunk=args.time_chunk,
        tile_chunk=args.tile_chunk,
        compress=not args.no_compress,
        complevel=args.complevel,
        scheduler=args.scheduler,
    )

    print(f"Done. Wrote: {args.out_nc}")
    print("Variables:", list(out.data_vars))
    print("Dims:", dict(out.dims))

if __name__ == "__main__":
    sys.exit(main())
