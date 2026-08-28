"""
Build the two delivery files requested in the Discover brief's packaging
specification (Section 2a of docs/notes/discover_peatclsm_budget_brief.md):

    OLv8_peat_fsw_2000_2024_compressed.nc
    DAv8_peat_fsw_2000_2024_compressed.nc

PEATCLSM_FSWCHANGE, monthly mean, (time=288, tile=112573), float32, zlib
complevel 4, chunks (12, 100000). Non-peat tiles are written as 0.0 (not
NaN/fill), matching the Fortran source's FSW_CHANGE=0. initialization.
lat/lon/time coordinates are copied verbatim from the existing audited
OLv8_water_budget_2000_2024_compressed.nc so tile ordering and the monthly
time axis are guaranteed identical, per the brief's requirement not to
re-sort.

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import glob
import hashlib
import os
from datetime import datetime, timezone

import numpy as np
from netCDF4 import Dataset

SOIL_PARAM = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/soil_param.dat"
POROS_COL = 6  # 0-indexed column 7
POROS_THRESHOLD = 0.90

CATCHMENT_DEF = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def"
AREA_COL = 6  # 0-indexed last column, km2

REFERENCE_PRODUCT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_water_budget_2000_2024_compressed.nc"

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states"

EXPERIMENTS = {
    "OL": {
        "cat_ens_avg": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
        "file_prefix": "LS_OLv8_M36",
        "out_name": "OLv8_peat_fsw_2000_2024_compressed.nc",
    },
    "DA": {
        "cat_ens_avg": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
        "file_prefix": "LS_DAv8_M36",
        "out_name": "DAv8_peat_fsw_2000_2024_compressed.nc",
    },
}


def load_poros():
    return np.loadtxt(SOIL_PARAM, usecols=(POROS_COL,))


def load_area_m2():
    return np.loadtxt(CATCHMENT_DEF, usecols=(AREA_COL,), skiprows=1) * 1.0e6


def find_monthly_files(base):
    pattern = os.path.join(base, "Y*", "M*", "*tavg24_1d_lnd_Nt.monthly.*.nc4")
    return sorted(glob.glob(pattern))


def load_fsw(path, n_tile, peat_mask):
    with Dataset(path) as ds:
        raw = ds.variables["PEATCLSM_FSWCHANGE"][0, :]
        val = np.ma.filled(raw.astype(np.float64), np.nan)

    undefined_on_peat = peat_mask & np.isnan(val)
    if undefined_on_peat.any():
        raise RuntimeError(
            f"{path}: PEATCLSM_FSWCHANGE undefined on {undefined_on_peat.sum()} "
            "peat tiles — unexpected, aborting rather than silently zeroing."
        )

    out = np.zeros(n_tile, dtype=np.float32)
    out[peat_mask] = val[peat_mask].astype(np.float32)
    return out


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    poros = load_poros()
    n_tile = poros.size
    peat_mask = poros >= POROS_THRESHOLD
    area_m2 = load_area_m2()
    assert area_m2.size == n_tile

    with Dataset(REFERENCE_PRODUCT) as ref:
        ref_time = np.array(ref.variables["time"][:])
        ref_lat = np.array(ref.variables["lat"][:])
        ref_lon = np.array(ref.variables["lon"][:])
        time_units = ref.variables["time"].units
        time_calendar = ref.variables["time"].calendar
    n_time = ref_time.size
    assert ref_lat.size == n_tile and ref_lon.size == n_tile

    print(f"tiles: {n_tile}, peat tiles (POROS>={POROS_THRESHOLD}): {peat_mask.sum()} "
          f"({100*peat_mask.mean():.3f}%)")
    print(f"time axis: {n_time} months, units='{time_units}'")

    report = {}

    for label, cfg in EXPERIMENTS.items():
        files = find_monthly_files(cfg["cat_ens_avg"])
        print(f"\n=== {label}: {len(files)} monthly files ===")
        if len(files) != n_time:
            raise RuntimeError(
                f"{label}: found {len(files)} files, expected {n_time} to match "
                "the reference product's time axis"
            )

        data = np.empty((n_time, n_tile), dtype=np.float32)
        for i, f in enumerate(files):
            data[i, :] = load_fsw(f, n_tile, peat_mask)

        out_path = os.path.join(OUT_DIR, cfg["out_name"])
        with Dataset(out_path, "w", format="NETCDF4") as ds:
            ds.createDimension("time", n_time)
            ds.createDimension("tile", n_tile)

            v_time = ds.createVariable("time", "i8", ("time",))
            v_time[:] = ref_time
            v_time.units = time_units
            v_time.calendar = time_calendar

            v_lat = ds.createVariable("lat", "f4", ("tile",), fill_value=np.nan)
            v_lat[:] = ref_lat
            v_lon = ds.createVariable("lon", "f4", ("tile",), fill_value=np.nan)
            v_lon[:] = ref_lon

            v_fsw = ds.createVariable(
                "PEATCLSM_FSWCHANGE", "f4", ("time", "tile"),
                fill_value=np.nan, zlib=True, complevel=4,
                chunksizes=(12, 100000),
            )
            v_fsw[:, :] = data
            v_fsw.long_name = "free_surface_water_on_peat_flux"
            v_fsw.units = "kg m-2 s-1"
            v_fsw.cell_methods = "time: mean"
            v_fsw.coordinates = "lat lon"

            ds.CreatedBy = ""
            ds.Date = datetime.now(timezone.utc).strftime("%Y-%m-%d at %H:%M:%S UTC")
            ds.source_root = cfg["cat_ens_avg"]
            ds.file_prefix = cfg["file_prefix"]
            ds.note = (
                "Monthly means from tavg24_1d_lnd_Nt.*; lat/lon/time copied from "
                "the audited *_water_budget_2000_2024_compressed.nc product "
                "(guarantees identical tile order and time axis). "
                "PEATCLSM_FSWCHANGE is 0.0 (not NaN) on non-peat tiles "
                f"(POROS<{POROS_THRESHOLD}), matching the model's own "
                "FSW_CHANGE=0. initialization in catchment.F90."
            )

        checksum = sha256sum(out_path)
        size_mb = os.path.getsize(out_path) / 1e6
        print(f"  wrote {out_path} ({size_mb:.1f} MB)")
        print(f"  sha256: {checksum}")

        # --- report stats ---
        nonzero_ever = (data != 0.0).any(axis=0)
        n_nonzero_tiles = int(nonzero_ever.sum())

        # area-weighted 24-yr mean, kg m-2 yr-1: mean flux (kg m-2 s-1) * seconds/yr
        seconds_per_year = 365.25 * 86400.0
        mean_flux_per_tile = data.mean(axis=0)  # kg m-2 s-1, time-mean per tile
        mean_rate_per_tile = mean_flux_per_tile * seconds_per_year  # kg m-2 yr-1

        aw_all = np.average(mean_rate_per_tile, weights=area_m2)
        aw_peat = np.average(mean_rate_per_tile[peat_mask], weights=area_m2[peat_mask])

        print(f"  tiles with any non-zero PEATCLSM_FSWCHANGE: {n_nonzero_tiles}")
        print(f"  area-weighted 24-yr mean, all land : {aw_all:.4f} kg m-2 yr-1")
        print(f"  area-weighted 24-yr mean, peat only : {aw_peat:.4f} kg m-2 yr-1")

        report[label] = {
            "path": out_path,
            "sha256": checksum,
            "n_nonzero_tiles": n_nonzero_tiles,
            "aw_all": aw_all,
            "aw_peat": aw_peat,
        }

    print("\n=== summary ===")
    for label, r in report.items():
        print(f"{label}: {r['path']}")
        print(f"  sha256: {r['sha256']}")
        print(f"  nonzero peat tiles: {r['n_nonzero_tiles']}")
        print(f"  area-weighted mean (all land / peat only): "
              f"{r['aw_all']:.4f} / {r['aw_peat']:.4f} kg m-2 yr-1")


if __name__ == "__main__":
    main()
