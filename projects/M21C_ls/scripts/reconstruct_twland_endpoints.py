"""
Reconstruct exact TWLAND at each water-year boundary (Oct 1, 2000-2006) for
OL and DA from catch_internal_rst restarts, per the Discover brief's Task 2
option (a):

    TWLAND = cdcr2/(1-wpwet) - catdef + rzexc + srfexc + capac + sum(wesnn)

Ensemble-mean over all 24 members. Compares the annual dStorage implied by
these true endpoints against the Sep/Oct two-month-mean approximation taken
from the existing *_water_budget_2000_2024_compressed.nc TWLAND field, for
each water year WY2001-WY2006 individually (not just the 6-year aggregate).

Run with the GEOSpyD full-path python3 (has netCDF4); no g5_modules needed.
"""

import numpy as np
from netCDF4 import Dataset

N_ENS = 24
N_TILE = 112573
YEARS = list(range(2000, 2007))  # Oct 1 of each year, 2000..2006

EXPERIMENTS = {
    "OL": {
        "restart_root": "/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs",
        "prefix": "LS_OLv8_M36",
        "water_budget_nc": "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/OLv8_water_budget_2000_2024_compressed.nc",
    },
    "DA": {
        "restart_root": "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs",
        "prefix": "LS_DAv8_M36",
        "water_budget_nc": "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/M21C_ls/output/monthly_flux_states/DAv8_water_budget_2000_2024_compressed.nc",
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

    twland = cdcr2 / (1.0 - wpwet) - catdef + rzexc + srfexc + capac + wesnn
    return twland


def twland_ens_mean(restart_root, prefix, year):
    acc = np.zeros(N_TILE, dtype=np.float64)
    for e in range(N_ENS):
        path = (f"{restart_root}/ens{e:04d}/Y{year}/M10/"
                f"{prefix}.catch_internal_rst.{year}1001_0000")
        acc += twland_one_member(path)
    return acc / N_ENS


def domain_mean(field, land_mask):
    return np.mean(field[land_mask])


def main():
    # land mask: finite lat/lon in the existing product marks valid land tiles
    with Dataset(EXPERIMENTS["OL"]["water_budget_nc"]) as ds:
        lat = np.array(ds.variables["lat"][:])
    land_mask = np.isfinite(lat)
    print(f"land tiles: {land_mask.sum()} of {N_TILE}")

    true_twland = {}  # [label][year] -> ensemble-mean TWLAND field
    for label, cfg in EXPERIMENTS.items():
        print(f"\n=== {label}: reconstructing TWLAND at Oct 1 of each year ===")
        true_twland[label] = {}
        for y in YEARS:
            tw = twland_ens_mean(cfg["restart_root"], cfg["prefix"], y)
            true_twland[label][y] = tw
            print(f"  Y{y}-10-01: domain-mean TWLAND = {domain_mean(tw, land_mask):.4f} kg m-2")

    # Sep/Oct-mean approximation from the monthly archive TWLAND
    approx_twland = {}  # [label][year] -> approx endpoint at "Oct 1 of year"
    for label, cfg in EXPERIMENTS.items():
        with Dataset(cfg["water_budget_nc"]) as ds:
            t = np.array(ds.variables["time"][:])  # days since 2000-06-01
            tw_all = ds.variables["TWLAND"][:, :]  # (time, tile), monthly means
        # time index for Sep of year y and Oct of year y
        base = np.datetime64("2000-06-01")
        dates = base + t.astype("timedelta64[D]")
        approx_twland[label] = {}
        for y in YEARS:
            sep_target = np.datetime64(f"{y}-09-01")
            oct_target = np.datetime64(f"{y}-10-01")
            i_sep = int(np.where(dates == sep_target)[0][0])
            i_oct = int(np.where(dates == oct_target)[0][0])
            approx = 0.5 * (np.ma.filled(tw_all[i_sep, :], np.nan)
                             + np.ma.filled(tw_all[i_oct, :], np.nan))
            approx_twland[label][y] = approx

    # annual water-year dStorage: TWLAND(Oct1,y+1) - TWLAND(Oct1,y), true vs approx
    print("\n=== annual (DA-OL) dStorage per water year: true restart-based vs Sep/Oct-mean approx ===")
    print(f"{'WY':6s} {'true dStorage':>16s} {'approx dStorage':>16s} {'diff':>10s}")
    for y in YEARS[:-1]:
        wy_label = f"WY{y+1}"

        true_dS_ol = domain_mean(true_twland["OL"][y+1] - true_twland["OL"][y], land_mask)
        true_dS_da = domain_mean(true_twland["DA"][y+1] - true_twland["DA"][y], land_mask)
        true_diff = true_dS_da - true_dS_ol

        approx_dS_ol = np.nanmean((approx_twland["OL"][y+1] - approx_twland["OL"][y])[land_mask])
        approx_dS_da = np.nanmean((approx_twland["DA"][y+1] - approx_twland["DA"][y])[land_mask])
        approx_diff = approx_dS_da - approx_dS_ol

        print(f"{wy_label:6s} {true_diff:16.4f} {approx_diff:16.4f} {true_diff-approx_diff:10.4f}")

    # six-year aggregate check (should reproduce the aggregate discussed in the brief)
    true_dS_ol_6yr = domain_mean(true_twland["OL"][2006] - true_twland["OL"][2000], land_mask)
    true_dS_da_6yr = domain_mean(true_twland["DA"][2006] - true_twland["DA"][2000], land_mask)
    approx_dS_ol_6yr = np.nanmean((approx_twland["OL"][2006] - approx_twland["OL"][2000])[land_mask])
    approx_dS_da_6yr = np.nanmean((approx_twland["DA"][2006] - approx_twland["DA"][2000])[land_mask])
    print("\n=== 6-year (WY2001-WY2006) aggregate ===")
    print(f"true   (DA-OL) dStorage : {true_dS_da_6yr - true_dS_ol_6yr:.4f} kg m-2")
    print(f"approx (DA-OL) dStorage : {approx_dS_da_6yr - approx_dS_ol_6yr:.4f} kg m-2")
    print(f"annualized true/approx  : {(true_dS_da_6yr - true_dS_ol_6yr)/6:.4f} / "
          f"{(approx_dS_da_6yr - approx_dS_ol_6yr)/6:.4f} kg m-2 yr-1")


if __name__ == "__main__":
    main()
