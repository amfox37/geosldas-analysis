"""
Verify the PEATCLSM free-standing-water flux identity on the archived monthly
tavg24_1d_lnd_Nt products for LS_OLv8_M36 and LS_DAv8_M36.

On peat tiles (POROS >= PEATCLSM_POROS_THRESHOLD, catch_constants.F90), catchment.F90
computes:

    FSW_CHANGE = P - EVAP - RUNOFF - WCHANGE

as a store deliberately excluded from WTOT. This script tests, per tile per month,
whether

    PRECTOTCORRLAND - EVLAND - RUNSURFLAND - BASEFLOWLAND - WCHANGELAND
        ~= PEATCLSM_FSWCHANGE

on peat tiles, and that PEATCLSM_FSWCHANGE ~= 0 on non-peat tiles.

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import glob
import os
import numpy as np
from netCDF4 import Dataset

SOIL_PARAM = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/soil_param.dat"
POROS_COL = 6  # 0-indexed column 7 per the Discover brief

EXPERIMENTS = {
    "OL": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
    "DA": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg",
}

VARS = ["PRECTOTCORRLAND", "EVLAND", "RUNSURFLAND", "BASEFLOWLAND", "WCHANGELAND", "PEATCLSM_FSWCHANGE"]


def load_poros():
    data = np.loadtxt(SOIL_PARAM, usecols=(POROS_COL,))
    return data


def find_monthly_files(base):
    pattern = os.path.join(base, "Y*", "M*", "*tavg24_1d_lnd_Nt.monthly.*.nc4")
    files = sorted(glob.glob(pattern))
    return files


def load_vars(path):
    with Dataset(path) as ds:
        out = {}
        for v in VARS:
            raw = ds.variables[v][0, :]  # netCDF4 auto-masks using _FillValue
            out[v] = np.ma.filled(raw.astype(np.float64), np.nan)
    return out


def main():
    poros = load_poros()
    n_tile = poros.size
    peat_mask = poros >= 0.90
    print(f"tiles: {n_tile}, peat tiles (POROS>=0.90): {peat_mask.sum()} "
          f"({100*peat_mask.mean():.2f}%)")

    for label, base in EXPERIMENTS.items():
        files = find_monthly_files(base)
        print(f"\n=== {label}: {len(files)} monthly files ===")
        if not files:
            print("  no files found, skipping")
            continue

        n_months = len(files)
        peat_resid_sum = np.zeros(n_tile)
        peat_resid_sqsum = np.zeros(n_tile)
        n_valid_peat = np.zeros(n_tile, dtype=int)
        n_valid_fsw = np.zeros(n_tile, dtype=int)   # months where FSWCHANGE is not fill
        n_masked_fsw_nonpeat = np.zeros(n_tile, dtype=int)  # months FSWCHANGE fill on non-peat

        for f in files:
            v = load_vars(f)
            lhs = v["PRECTOTCORRLAND"] - v["EVLAND"] - v["RUNSURFLAND"] - v["BASEFLOWLAND"] - v["WCHANGELAND"]
            fsw = v["PEATCLSM_FSWCHANGE"]
            resid = lhs - fsw

            fsw_defined = ~np.isnan(fsw)
            n_valid_fsw += fsw_defined.astype(int)
            n_masked_fsw_nonpeat += ((~fsw_defined) & (~peat_mask)).astype(int)

            valid = peat_mask & fsw_defined & ~np.isnan(resid)
            n_valid_peat += valid.astype(int)
            peat_resid_sum += np.where(valid, resid, 0.0)
            peat_resid_sqsum += np.where(valid, resid**2, 0.0)

        with np.errstate(invalid="ignore", divide="ignore"):
            resid_mean = peat_resid_sum / n_valid_peat
            resid_rms = np.sqrt(peat_resid_sqsum / n_valid_peat)

        pm = peat_mask & (n_valid_peat > 0)

        print(f"  PEAT tiles (n={pm.sum()} of {peat_mask.sum()}): identity residual "
              f"(kg m-2 s-1), time-mean-of-tile:")
        print(f"    mean(resid)      = {np.nanmean(resid_mean[pm]):.3e}")
        print(f"    mean(|resid|)    = {np.nanmean(np.abs(resid_mean[pm])):.3e}")
        print(f"    max(|resid|)     = {np.nanmax(np.abs(resid_mean[pm])):.3e}")
        print(f"    mean(rms resid)  = {np.nanmean(resid_rms[pm]):.3e}")

        # FSWCHANGE is fill-valued off of peat in the archive (unlike the Fortran
        # initialization to 0.), so the non-peat check is: is it consistently masked?
        n_nonpeat = (~peat_mask).sum()
        always_masked = (n_masked_fsw_nonpeat == n_months).sum()
        print(f"  NON-PEAT tiles (n={n_nonpeat}): FSWCHANGE is fill-valued "
              f"(not literal 0) in the archive")
        print(f"    tiles masked in ALL {n_months} months: {always_masked} "
              f"({100*always_masked/n_nonpeat:.2f}%)")
        # any peat tiles unexpectedly masked, or non-peat tiles unexpectedly defined?
        peat_never_defined = peat_mask & (n_valid_fsw == 0)
        nonpeat_ever_defined = (~peat_mask) & (n_masked_fsw_nonpeat < n_months)
        print(f"    peat tiles never defined (n={peat_never_defined.sum()}), "
              f"non-peat tiles defined in >=1 month (n={nonpeat_ever_defined.sum()})")

        # relative to typical flux magnitude for context
        mid = load_vars(files[len(files)//2])
        typical_precip = np.nanmean(np.abs(mid["PRECTOTCORRLAND"][pm]))
        print(f"  (for scale: mean |PRECTOTCORRLAND| on peat tiles, mid-record month "
              f"= {typical_precip:.3e} kg m-2 s-1)")


if __name__ == "__main__":
    main()
