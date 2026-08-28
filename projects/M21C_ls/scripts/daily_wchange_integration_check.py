"""
Test whether integrating WCHANGELAND at DAILY resolution (instead of the
monthly-mean archive) closes the gap against the true restart-based OL
storage change over WY2001-WY2006.

Extracts only the needed tavg24_1d_lnd_Nt daily files, one tar (calendar
year) at a time, reads WCHANGELAND, accumulates cum(WCHANGE*86400), deletes
the extracted files, and moves to the next year to keep scratch usage low.

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import os
import shutil
import subprocess
import tarfile

import numpy as np
from netCDF4 import Dataset

N_TILE = 112573
SCRATCH = "/gpfsm/dnb34/amfox/scratch_daily_wchange"
TAR_ROOT = "/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/SMAP_EASEv2_M36_GLOBAL_cat_tars"

area_m2 = np.loadtxt(
    "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def",
    usecols=(6,), skiprows=1) * 1.0e6


def aw_mean(field):
    valid = np.isfinite(field)
    return float(np.average(field[valid], weights=area_m2[valid]))


# (year, months_needed) for the Oct2000-Sep2006 window
YEAR_MONTHS = {
    2000: [10, 11, 12],
    2001: list(range(1, 13)),
    2002: list(range(1, 13)),
    2003: list(range(1, 13)),
    2004: list(range(1, 13)),
    2005: list(range(1, 13)),
    2006: list(range(1, 10)),
}


def main():
    cum_wchange = np.zeros(N_TILE, dtype=np.float64)
    n_days_total = 0

    for year, months in YEAR_MONTHS.items():
        tarpath = os.path.join(TAR_ROOT, f"Y{year}.tar")
        year_dir = os.path.join(SCRATCH, f"Y{year}")
        os.makedirs(year_dir, exist_ok=True)

        with tarfile.open(tarpath, "r") as tf:
            members = [
                m for m in tf.getmembers()
                if "tavg24_1d_lnd_Nt" in m.name
                and any(f"/M{mo:02d}/" in m.name for mo in months)
            ]
            print(f"Y{year}: extracting {len(members)} daily tavg24 files...", flush=True)
            tf.extractall(path=SCRATCH, members=members)

        # process
        day_files = []
        for mo in months:
            month_dir = os.path.join(SCRATCH, f"Y{year}", f"M{mo:02d}")
            if os.path.isdir(month_dir):
                for fn in sorted(os.listdir(month_dir)):
                    if "tavg24_1d_lnd_Nt" in fn:
                        day_files.append(os.path.join(month_dir, fn))

        for path in day_files:
            with Dataset(path) as ds:
                wc = ds.variables["WCHANGELAND"][0, :]
                wc = np.ma.filled(wc.astype(np.float64), np.nan)
            cum_wchange += wc * 86400.0
            n_days_total += 1

        print(f"Y{year}: processed {len(day_files)} days, running cum(WCHANGE) area-weighted mean = "
              f"{aw_mean(cum_wchange):.4f} kg m-2", flush=True)

        # free scratch space before next year
        shutil.rmtree(year_dir)

    print(f"\ntotal days integrated: {n_days_total}")
    print(f"cum(WCHANGELAND), DAILY integration, area-weighted: {aw_mean(cum_wchange):.4f} kg m-2")
    print("cum(WCHANGELAND), MONTHLY integration (from earlier): -1.321 kg m-2")
    print("dStorage_true (restart endpoints), area-weighted    : -4.903 kg m-2")
    gap_daily = cum_wchange - (-4.903)
    print(f"\ngap (daily-integrated WCHANGE minus true dStorage), area-weighted: {aw_mean(gap_daily):.4f} kg m-2")
    print("gap using MONTHLY integration (from earlier): +3.582 kg m-2")


if __name__ == "__main__":
    main()
