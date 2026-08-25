"""
Sum the DA analysis-increment water mass (catch_progn_incr) over WY2001-WY2006
to test whether it explains the gap between integrated WCHANGELAND and the
true restart-based storage change found in close_da_water_budget.py.

catch_progn_incr files carry the direct state increments applied by the
analysis step (CATDEF_INCR, RZEXC_INCR, SRFEXC_INCR, CAPAC_INCR, WESNN1-3_INCR),
each in kg m-2, at up to 8 three-hourly analysis windows per day. This mass
injection happens at the analysis step, outside the physics timestep, so it
is invisible to WCHANGELAND (a pure physics-flux tendency) but IS captured by
the true restart-to-restart TWLAND difference. The mirror term to WTOT's
formula (CDCR2, WPWET are static, never incremented):

    dWTOT_from_increments = -CATDEF_INCR + RZEXC_INCR + SRFEXC_INCR
                              + CAPAC_INCR + WESNN1_INCR + WESNN2_INCR + WESNN3_INCR

Extracts one calendar year of catch_progn_incr files at a time from the tar
archive, sums, deletes, moves to next year.

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import os
import shutil
import tarfile

import numpy as np
from netCDF4 import Dataset

N_TILE = 112573
SCRATCH = "/gpfsm/dnb34/amfox/scratch_da_incr"
TAR_ROOT = "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/SMAP_EASEv2_M36_GLOBAL_cat_tars"

area_m2 = np.loadtxt(
    "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def",
    usecols=(6,), skiprows=1) * 1.0e6


def aw_mean(field):
    valid = np.isfinite(field)
    return float(np.average(field[valid], weights=area_m2[valid]))


YEAR_MONTHS = {
    2000: [10, 11, 12],
    2001: list(range(1, 13)),
    2002: list(range(1, 13)),
    2003: list(range(1, 13)),
    2004: list(range(1, 13)),
    2005: list(range(1, 13)),
    2006: list(range(1, 10)),
}

INCR_VARS = ["CATDEF_INCR", "RZEXC_INCR", "SRFEXC_INCR", "CAPAC_INCR",
             "WESNN1_INCR", "WESNN2_INCR", "WESNN3_INCR"]


def main():
    cum_incr = np.zeros(N_TILE, dtype=np.float64)
    n_files_total = 0
    n_nonzero_windows = 0

    for year, months in YEAR_MONTHS.items():
        tarpath = os.path.join(TAR_ROOT, f"Y{year}.tar")
        year_dir = os.path.join(SCRATCH, f"Y{year}")
        os.makedirs(year_dir, exist_ok=True)

        with tarfile.open(tarpath, "r") as tf:
            members = [
                m for m in tf.getmembers()
                if "catch_progn_incr" in m.name
                and any(f"/M{mo:02d}/" in m.name for mo in months)
            ]
            print(f"Y{year}: extracting {len(members)} daily catch_progn_incr files...", flush=True)
            tf.extractall(path=SCRATCH, members=members)

        day_files = []
        for mo in months:
            month_dir = os.path.join(SCRATCH, f"Y{year}", f"M{mo:02d}")
            if os.path.isdir(month_dir):
                for fn in sorted(os.listdir(month_dir)):
                    if "catch_progn_incr" in fn:
                        day_files.append(os.path.join(month_dir, fn))

        for path in day_files:
            with Dataset(path) as ds:
                n_t = ds.dimensions["time"].size
                for v in INCR_VARS:
                    raw = np.array(ds.variables[v][:])  # (time, tile)
                    fill = ds.variables[v].getncattr("_FillValue")
                    valid = np.abs(raw - fill) > 1.0
                    contrib = np.where(valid, raw, 0.0)
                    sign = -1.0 if v == "CATDEF_INCR" else 1.0
                    cum_incr += sign * contrib.sum(axis=0)
                    if v == "WESNN1_INCR":
                        n_nonzero_windows += int((contrib != 0).any(axis=1).sum())
            n_files_total += 1

        print(f"Y{year}: processed {len(day_files)} days, running cum(increment mass) "
              f"area-weighted mean = {aw_mean(cum_incr):.4f} kg m-2", flush=True)

        shutil.rmtree(year_dir)

    print(f"\ntotal daily files processed: {n_files_total}")
    print(f"total 3-hourly windows with nonzero WESNN1_INCR: {n_nonzero_windows}")
    print(f"\ncum DA analysis-increment mass, area-weighted, all land: {aw_mean(cum_incr):.4f} kg m-2")
    print("gap (cum WCHANGE - true dStorage) from close_da_water_budget.py: -311.723 kg m-2")
    print(f"increment mass + gap = {aw_mean(cum_incr) + (-311.723):.4f} kg m-2  (should be near 0 if this explains the gap)")


if __name__ == "__main__":
    main()
