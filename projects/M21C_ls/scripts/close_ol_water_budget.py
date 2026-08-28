"""
Close the OL water budget over WY2001-WY2006 (Oct 2000 - Sep 2006), single
experiment, no DA differencing.

Tests, per tile, whether

    cum(P) - cum(EVLAND) - cum(RUNSURFLAND + BASEFLOWLAND)
        - dStorage_true - cum(PEATCLSM_FSWCHANGE) ~= 0

where dStorage_true is the restart-based TWLAND(2006-10-01) - TWLAND(2000-10-01)
(OLv8_twland_wy_endpoints_2000_2006.nc), and cum(PEATCLSM_FSWCHANGE) is the
6-year cumulative peat free-standing-water flux (OLv8_peat_fsw_2000_2024_compressed.nc).
Reports three progressively-corrected residuals, mirroring the brief's own
table: as published-style (no peat term, Sep-mean storage), + true storage
endpoint, + true storage endpoint + peat term.

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import glob
import os

import numpy as np
from netCDF4 import Dataset

N_TILE = 112573
SOIL_PARAM = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/soil_param.dat"
CATCHMENT_DEF = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def"

OL_CAT_ENS_AVG = "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg"
TWLAND_ENDPOINTS = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_twland_wy_endpoints_2000_2006.nc"
PEAT_FSW = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_peat_fsw_2000_2024_compressed.nc"
WATER_BUDGET_REF = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_water_budget_2000_2024_compressed.nc"

WY_START = "2000-10"  # Oct 2000
WY_END = "2006-09"    # Sep 2006, inclusive


def load_poros():
    return np.loadtxt(SOIL_PARAM, usecols=(6,))


def load_area_m2():
    return np.loadtxt(CATCHMENT_DEF, usecols=(6,), skiprows=1) * 1.0e6


def months_in_window():
    months = []
    y, m = 2000, 10
    while f"{y}-{m:02d}" != "2006-10":
        months.append((y, m))
        m += 1
        if m == 13:
            m = 1
            y += 1
    return months


def load_monthly_vars(path):
    with Dataset(path) as ds:
        out = {}
        for v in ["PRECTOTCORRLAND", "EVLAND", "RUNSURFLAND", "BASEFLOWLAND", "WCHANGELAND"]:
            raw = ds.variables[v][0, :]
            out[v] = np.ma.filled(raw.astype(np.float64), np.nan)
    return out


def days_in_month(year, month):
    if month == 12:
        nxt_y, nxt_m = year + 1, 1
    else:
        nxt_y, nxt_m = year, month + 1
    return (np.datetime64(f"{nxt_y}-{nxt_m:02d}-01") - np.datetime64(f"{year}-{month:02d}-01")).astype(int)


def main():
    poros = load_poros()
    peat_mask = poros >= 0.90
    area_m2 = load_area_m2()

    months = months_in_window()
    print(f"window: {months[0]} .. {months[-1]}  ({len(months)} months)")

    cum_P = np.zeros(N_TILE)
    cum_ET = np.zeros(N_TILE)
    cum_Runoff = np.zeros(N_TILE)
    cum_WCHANGE = np.zeros(N_TILE)

    for (y, m) in months:
        path = os.path.join(OL_CAT_ENS_AVG, f"Y{y}", f"M{m:02d}",
                             f"LS_OLv8_M36.tavg24_1d_lnd_Nt.monthly.{y}{m:02d}.nc4")
        v = load_monthly_vars(path)
        seconds = days_in_month(y, m) * 86400.0
        cum_P += v["PRECTOTCORRLAND"] * seconds
        cum_ET += v["EVLAND"] * seconds
        cum_Runoff += (v["RUNSURFLAND"] + v["BASEFLOWLAND"]) * seconds
        cum_WCHANGE += v["WCHANGELAND"] * seconds

    # true restart-based storage change
    with Dataset(TWLAND_ENDPOINTS) as ds:
        t = np.array(ds.variables["time"][:])  # days since 2000-06-01
        tw = np.array(ds.variables["TWLAND"][:, :])
    base = np.datetime64("2000-06-01")
    dates = base + t.astype("timedelta64[D]")
    i0 = int(np.where(dates == np.datetime64("2000-10-01"))[0][0])
    i1 = int(np.where(dates == np.datetime64("2006-10-01"))[0][0])
    dStorage_true = tw[i1, :] - tw[i0, :]

    # published-style approximation: Sep-mean TWLAND difference (2000-09 -> 2006-09)
    with Dataset(WATER_BUDGET_REF) as ds:
        t_ref = np.array(ds.variables["time"][:])
        tw_ref = np.array(ds.variables["TWLAND"][:, :])
    dates_ref = base + t_ref.astype("timedelta64[D]")
    j0 = int(np.where(dates_ref == np.datetime64("2000-09-01"))[0][0])
    j1 = int(np.where(dates_ref == np.datetime64("2006-09-01"))[0][0])
    dStorage_approx = np.ma.filled(tw_ref[j1, :], np.nan) - np.ma.filled(tw_ref[j0, :], np.nan)

    # cumulative peat FSW over the same window
    with Dataset(PEAT_FSW) as ds:
        t_fsw = np.array(ds.variables["time"][:])
        fsw = np.array(ds.variables["PEATCLSM_FSWCHANGE"][:, :])
    dates_fsw = base + t_fsw.astype("timedelta64[D]")
    keep = (dates_fsw >= np.datetime64("2000-10-01")) & (dates_fsw <= np.datetime64("2006-09-01"))
    fsw_window = fsw[keep, :]
    fsw_dates = dates_fsw[keep]
    seconds_per_month = np.array([
        days_in_month(int(str(d)[:4]), int(str(d)[5:7])) * 86400.0 for d in fsw_dates
    ])
    cum_FSW = (fsw_window * seconds_per_month[:, None]).sum(axis=0)

    # sanity: process-balance identity (should hold ~exactly for OL, no assim increments)
    process_diag = cum_P - cum_ET - cum_Runoff - cum_FSW - cum_WCHANGE

    def area_mean(field, mask=None):
        if mask is None:
            mask = np.ones(N_TILE, dtype=bool)
            valid = np.isfinite(field) & mask
        else:
            valid = np.isfinite(field) & mask
        return float(np.average(field[valid], weights=area_m2[valid]))

    residual_A = cum_P - cum_ET - cum_Runoff - dStorage_approx  # as-published style, no peat term
    residual_B = cum_P - cum_ET - cum_Runoff - dStorage_true    # + true storage endpoint, no peat term
    residual_C = cum_P - cum_ET - cum_Runoff - dStorage_true - cum_FSW  # + true storage + peat term

    print("\n=== OL water budget closure, WY2001-WY2006 (6-year cumulative, kg m-2), area-weighted, all 112,573 tiles ===")
    print(f"cum P                              : {area_mean(cum_P):10.3f}")
    print(f"cum ET                              : {area_mean(cum_ET):10.3f}")
    print(f"cum Runoff (surf+base)              : {area_mean(cum_Runoff):10.3f}")
    print(f"cum PEATCLSM_FSWCHANGE              : {area_mean(cum_FSW):10.3f}")
    print(f"dStorage, Sep-mean approx           : {area_mean(dStorage_approx):10.3f}")
    print(f"dStorage, true restart endpoints    : {area_mean(dStorage_true):10.3f}")
    print()
    print(f"residual A  (approx storage, no peat term)          : {area_mean(residual_A):10.3f} kg m-2  "
          f"({100*area_mean(residual_A)/area_mean(cum_P):.3f}% of P)")
    print(f"residual B  (true storage, no peat term)             : {area_mean(residual_B):10.3f} kg m-2  "
          f"({100*area_mean(residual_B)/area_mean(cum_P):.3f}% of P)")
    print(f"residual C  (true storage + peat term)               : {area_mean(residual_C):10.3f} kg m-2  "
          f"({100*area_mean(residual_C)/area_mean(cum_P):.3f}% of P)")
    print()
    print(f"process-balance identity P-ET-Runoff-FSW-WCHANGE     : "
          f"{area_mean(process_diag):.3e} kg m-2 (should be ~0, sanity check)")

    print("\n--- peat tiles only (n={}) ---".format(peat_mask.sum()))
    print(f"residual A : {area_mean(residual_A, peat_mask):10.3f} kg m-2")
    print(f"residual B : {area_mean(residual_B, peat_mask):10.3f} kg m-2")
    print(f"residual C : {area_mean(residual_C, peat_mask):10.3f} kg m-2")

    print("\n--- non-peat tiles only (n={}) ---".format((~peat_mask).sum()))
    print(f"residual A : {area_mean(residual_A, ~peat_mask):10.3f} kg m-2")
    print(f"residual B : {area_mean(residual_B, ~peat_mask):10.3f} kg m-2")
    print(f"residual C : {area_mean(residual_C, ~peat_mask):10.3f} kg m-2")

    print("\n=== isolating the leftover: integrated WCHANGELAND alone vs true restart storage change ===")
    # NOTE: WCHANGELAND already equals P-ET-Runoff-FSW by construction on peat tiles
    # (the identity we verified earlier), so FSW must NOT be added again here.
    gap = cum_WCHANGE - dStorage_true
    print(f"cum WCHANGELAND                      : {area_mean(cum_WCHANGE):10.3f} kg m-2")
    print(f"dStorage_true (restart endpoints)    : {area_mean(dStorage_true):10.3f} kg m-2")
    print(f"gap = integral - true                : {area_mean(gap):10.3f} kg m-2  "
          f"(all land, area-weighted)")
    print(f"gap, non-peat only                   : {area_mean(gap, ~peat_mask):10.3f} kg m-2")
    print(f"gap, peat only                        : {area_mean(gap, peat_mask):10.3f} kg m-2")
    print(f"(residual_C = process_diag + gap: {area_mean(process_diag):.3f} + {area_mean(gap):.3f} "
          f"= {area_mean(process_diag)+area_mean(gap):.3f}, cf. residual C above)")


if __name__ == "__main__":
    main()
