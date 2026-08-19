"""
Package the restart-reconstructed TWLAND water-year boundary values (Task 2
of the Discover brief, docs/notes/discover_peatclsm_budget_brief.md) as two
delivery files, one per experiment:

    OLv8_twland_wy_endpoints_2000_2006.nc
    DAv8_twland_wy_endpoints_2000_2006.nc

Each holds TWLAND at Oct 1 of every year 2000-2006 (the 7 boundaries bounding
water years WY2001-WY2006), reconstructed exactly from catch_internal_rst
restarts via

    TWLAND = cdcr2/(1-wpwet) - catdef + rzexc + srfexc + capac + sum(wesnn)

ensemble-mean over all 24 members, per tile. This replaces the Sep/Oct
monthly-mean approximation for this window since the true restart-based
values are available at zero extra cost (see the companion doc,
docs/notes/twland_endpoint_reconstruction.md, for the full explanation and
the aggregate-vs-annual accuracy comparison against the approximation).

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import hashlib
import os
from datetime import datetime, timezone

import numpy as np
from netCDF4 import Dataset

N_ENS = 24
N_TILE = 112573
YEARS = list(range(2000, 2007))  # Oct 1 of each year, 2000..2006

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states"
REFERENCE_PRODUCT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_water_budget_2000_2024_compressed.nc"

EXPERIMENTS = {
    "OL": {
        "restart_root": "/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs",
        "prefix": "LS_OLv8_M36",
        "out_name": "OLv8_twland_wy_endpoints_2000_2006.nc",
        "source_experiment": "M21C_land_sweeper_OLv8_M36 (predecessor leg of LS_OLv8_M36_v2, via its RESTART_PATH)",
    },
    "DA": {
        "restart_root": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs",
        "prefix": "LS_DAv8_M36",
        "out_name": "DAv8_twland_wy_endpoints_2000_2006.nc",
        "source_experiment": "LS_DAv8_M36_v2 (predecessor leg of LS_DAv8_M36_v3)",
    },
}


def twland_one_member(path):
    with Dataset(path) as ds:
        catdef = ds.variables["CATDEF"][:].astype(np.float64)
        rzexc = ds.variables["RZEXC"][:].astype(np.float64)
        srfexc = ds.variables["SRFEXC"][:].astype(np.float64)
        capac = ds.variables["CAPAC"][:].astype(np.float64)
        wesnn = (ds.variables["WESNN1"][:].astype(np.float64)
                 + ds.variables["WESNN2"][:].astype(np.float64)
                 + ds.variables["WESNN3"][:].astype(np.float64))
        cdcr2 = ds.variables["CDCR2"][:].astype(np.float64)
        wpwet = ds.variables["WPWET"][:].astype(np.float64)
    return cdcr2 / (1.0 - wpwet) - catdef + rzexc + srfexc + capac + wesnn


def twland_ens_mean(restart_root, prefix, year):
    acc = np.zeros(N_TILE, dtype=np.float64)
    for e in range(N_ENS):
        path = (f"{restart_root}/ens{e:04d}/Y{year}/M10/"
                f"{prefix}.catch_internal_rst.{year}1001_0000")
        acc += twland_one_member(path)
    return acc / N_ENS


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    with Dataset(REFERENCE_PRODUCT) as ref:
        ref_lat = np.array(ref.variables["lat"][:])
        ref_lon = np.array(ref.variables["lon"][:])
    assert ref_lat.size == N_TILE and ref_lon.size == N_TILE

    # CF time: days since 2000-06-01, matching the existing product's convention
    boundary_days = [
        int((np.datetime64(f"{y}-10-01") - np.datetime64("2000-06-01"))
            / np.timedelta64(1, "D"))
        for y in YEARS
    ]

    for label, cfg in EXPERIMENTS.items():
        print(f"\n=== {label}: reconstructing TWLAND at Oct 1, {YEARS[0]}-{YEARS[-1]} ===")
        data = np.empty((len(YEARS), N_TILE), dtype=np.float32)
        for i, y in enumerate(YEARS):
            tw = twland_ens_mean(cfg["restart_root"], cfg["prefix"], y)
            data[i, :] = tw.astype(np.float32)
            print(f"  Y{y}-10-01: domain-mean TWLAND = {tw.mean():.4f} kg m-2")

        out_path = os.path.join(OUT_DIR, cfg["out_name"])
        with Dataset(out_path, "w", format="NETCDF4") as ds:
            ds.createDimension("wy_boundary", len(YEARS))
            ds.createDimension("tile", N_TILE)

            v_time = ds.createVariable("time", "i8", ("wy_boundary",))
            v_time[:] = np.array(boundary_days, dtype=np.int64)
            v_time.units = "days since 2000-06-01 00:00:00"
            v_time.calendar = "proleptic_gregorian"
            v_time.long_name = "water_year_boundary_date_Oct1_00Z"

            v_lat = ds.createVariable("lat", "f4", ("tile",), fill_value=np.nan)
            v_lat[:] = ref_lat
            v_lon = ds.createVariable("lon", "f4", ("tile",), fill_value=np.nan)
            v_lon[:] = ref_lon

            v_tw = ds.createVariable(
                "TWLAND", "f4", ("wy_boundary", "tile"),
                fill_value=np.nan, zlib=True, complevel=4,
            )
            v_tw[:, :] = data
            v_tw.long_name = "total_water_storage_land"
            v_tw.units = "kg m-2"
            v_tw.cell_methods = "restart snapshot at 00Z Oct 1, ensemble-mean of 24 members"
            v_tw.coordinates = "lat lon"

            ds.CreatedBy = ""
            ds.Date = datetime.now(timezone.utc).strftime("%Y-%m-%d at %H:%M:%S UTC")
            ds.source_root = cfg["restart_root"]
            ds.source_experiment = cfg["source_experiment"]
            ds.ensemble_size = N_ENS
            ds.formula = ("TWLAND = CDCR2/(1-WPWET) - CATDEF + RZEXC + SRFEXC "
                           "+ CAPAC + WESNN1 + WESNN2 + WESNN3, all read from "
                           "catch_internal_rst.<YYYY>1001_0000 restarts, "
                           "ensemble-mean over 24 members.")
            ds.note = (
                "Exact restart-based TWLAND at each water-year boundary "
                "(Oct 1, 00Z), replacing the Sep/Oct monthly-mean "
                "approximation for WY2001-WY2006. See "
                "docs/notes/twland_endpoint_reconstruction.md for the full "
                "explanation, provenance, and accuracy comparison against "
                "the approximation (aggregate 6-yr match is close, ~2%, "
                "but individual water years differ by up to ~20%)."
            )

        checksum = sha256sum(out_path)
        size_kb = os.path.getsize(out_path) / 1e3
        print(f"  wrote {out_path} ({size_kb:.1f} KB)")
        print(f"  sha256: {checksum}")


if __name__ == "__main__":
    main()
