#!/usr/bin/env python3
"""
Verify H121 ASCAT z-score scaling by applying it to raw OFA obs and computing O-F.

The OLv7_M36_MULTI_type_13_H121 experiment is a monitor-mode OL run, so the 'obs'
field in each OFA file is the raw (unscaled) ASCAT SSM in degree-of-saturation %.
We apply our newly generated scaling parameters and compute:

    scaled_obs = m_mean + (m_std / o_std) * (raw_obs - o_mean)
    O-F = scaled_obs - fcst

If the scaling is working the mean O-F should tend toward zero.
"""
import argparse
from pathlib import Path
from datetime import date
import numpy as np
import netCDF4 as nc

# ---------- CLI ----------
parser = argparse.ArgumentParser(description="Verify z-score scaling O-F statistics.")
parser.add_argument("--stats-file", required=True, type=Path, help="Path to all_pentads stats nc4.")
parser.add_argument("--species", default="8,9,10", help="Comma-separated species IDs (default: 8,9,10).")
parser.add_argument("--species-names", default="", help="Comma-separated names matching --species order.")
parser.add_argument("--exp-path", default="/discover/nobackup/projects/land_da/hsaf_cdr_test", type=Path)
parser.add_argument("--exp-run",  default="OLv7_M36_MULTI_type_13_H121")
parser.add_argument("--domain",   default="SMAP_EASEv2_M36_GLOBAL")
parser.add_argument("--start-year", type=int, default=2017)
parser.add_argument("--end-year",   type=int, default=2020)
args = parser.parse_args()

STATS_FILE    = args.stats_file
ANA_BASE      = args.exp_path / args.exp_run / "output" / args.domain / "ana" / "ens_avg"
SPECIES_LIST  = [int(s) for s in args.species.split(",")]
if args.species_names:
    _names = args.species_names.split(",")
    SPECIES_NAMES = {sp: nm for sp, nm in zip(SPECIES_LIST, _names)}
else:
    SPECIES_NAMES = {sp: f"species_{sp}" for sp in SPECIES_LIST}
YEAR_RANGE = range(args.start_year, args.end_year + 1)

# ---------- load stats ----------
print("Loading scaling parameter file ...", flush=True)
with nc.Dataset(STATS_FILE) as ds:
    ll_lon = float(ds.variables['ll_lon'][:])
    ll_lat = float(ds.variables['ll_lat'][:])
    d_lon  = float(ds.variables['d_lon'][:])
    d_lat  = float(ds.variables['d_lat'][:])
    # shape: (pentad=73, lon=1440, lat=720)
    o_mean = ds.variables['o_mean'][:].filled(np.nan)
    o_std  = ds.variables['o_std'][:].filled(np.nan)
    m_mean = ds.variables['m_mean'][:].filled(np.nan)
    m_std  = ds.variables['m_std'][:].filled(np.nan)
    m_min  = ds.variables['m_min'][:].filled(np.nan)
    m_max  = ds.variables['m_max'][:].filled(np.nan)

N_LON, N_LAT = o_mean.shape[1], o_mean.shape[2]
print(f"Stats loaded: {o_mean.shape}  ll_lon={ll_lon} ll_lat={ll_lat} d={d_lon}", flush=True)

# ---------- helpers ----------
def doy(yr, mo, dy):
    return (date(yr, mo, dy) - date(yr, 1, 1)).days + 1

def pentad_idx(d):
    return min((d - 1) // 5, 72)   # 0-based, 0-72; DOY 366 clipped to pentad 72

# ---------- accumulators ----------
# [N, sum_omf, sum_omf2] per species per month
accum = {sp: {mo: np.zeros(3) for mo in range(1, 13)} for sp in SPECIES_LIST}
n_files = 0
n_obs_total = 0

# ---------- main loop ----------
for year in YEAR_RANGE:
    for month in range(1, 13):
        ofa_dir = ANA_BASE / f"Y{year}" / f"M{month:02d}"
        if not ofa_dir.exists():
            continue

        for ofa_file in sorted(ofa_dir.glob("*.nc4")):
            # parse date from filename token e.g. 20180601_0300z
            tok = ofa_file.stem.split('.')[-1]
            yr2, mo2, dy2 = int(tok[:4]), int(tok[4:6]), int(tok[6:8])
            pidx = pentad_idx(doy(yr2, mo2, dy2))

            with nc.Dataset(ofa_file) as ds:
                species  = ds.variables['species'][:].data.astype(int)
                obs_raw  = np.ma.filled(ds.variables['obs'][:],  np.nan).astype(np.float64)
                obs_fcst = np.ma.filled(ds.variables['fcst'][:], np.nan).astype(np.float64)
                obs_lon  = np.ma.filled(ds.variables['lon'][:],  np.nan).astype(np.float64)
                obs_lat  = np.ma.filled(ds.variables['lat'][:],  np.nan).astype(np.float64)

            for sp in SPECIES_LIST:
                mask = (species == sp) & np.isfinite(obs_raw) & np.isfinite(obs_fcst)
                if not np.any(mask):
                    continue

                raw  = obs_raw[mask]
                fcst = obs_fcst[mask]
                lon  = obs_lon[mask]
                lat  = obs_lat[mask]

                # ceiling convention (matches Fortran and Python stats generator)
                i_idx = np.ceil((lon - ll_lon) / d_lon).astype(int) - 1
                j_idx = np.ceil((lat - ll_lat) / d_lat).astype(int) - 1

                in_bounds = (i_idx >= 0) & (i_idx < N_LON) & (j_idx >= 0) & (j_idx < N_LAT)
                if not np.any(in_bounds):
                    continue

                ii = i_idx[in_bounds]
                jj = j_idx[in_bounds]

                om  = o_mean[pidx, ii, jj]
                os_ = o_std[ pidx, ii, jj]
                mm  = m_mean[pidx, ii, jj]
                ms  = m_std[ pidx, ii, jj]
                mlo = m_min[ii, jj]
                mhi = m_max[ii, jj]

                can_scale = (
                    (om > 0) & (mm > 0) & (os_ > 0)
                    & np.isfinite(om) & np.isfinite(os_) & np.isfinite(mm) & np.isfinite(ms)
                )
                if not np.any(can_scale):
                    continue

                ratio  = np.where(can_scale, ms / os_, np.nan)
                scaled = mm + ratio * (raw[in_bounds] - om)
                scaled = np.clip(scaled, mlo, mhi)
                scaled = np.where(can_scale, scaled, np.nan)

                omf   = scaled - fcst[in_bounds]
                valid = np.isfinite(omf)
                n     = int(np.sum(valid))
                if n > 0:
                    omf_v = omf[valid]
                    accum[sp][month][0] += n
                    accum[sp][month][1] += omf_v.sum()
                    accum[sp][month][2] += (omf_v ** 2).sum()
                    n_obs_total += n

            n_files += 1
            if n_files % 1000 == 0:
                print(f"  {n_files} files, {n_obs_total:,} obs ...", flush=True)

# ---------- report ----------
MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

print(f"\nProcessed {n_files} OFA files  |  {n_obs_total:,} scaled obs  |  years {YEAR_RANGE.start}-{YEAR_RANGE.stop-1}")

for sp in SPECIES_LIST:
    print(f"\n{SPECIES_NAMES[sp]}")
    print(f"  {'Month':<6} {'N':>8} {'Mean O-F':>12} {'Std O-F':>10}")
    print(f"  {'-'*40}")
    for mo in range(1, 13):
        N, S, S2 = accum[sp][mo]
        if N > 0:
            mean = S / N
            std  = np.sqrt(max(S2 / N - mean ** 2, 0.0))
            print(f"  {MONTHS[mo-1]:<6} {int(N):>8} {mean:>12.5f} {std:>10.5f}")

print("\n--- Overall ---")
print(f"{'Species':<25} {'N':>10} {'Mean O-F':>12} {'Std O-F':>10}")
print("-" * 60)
for sp in SPECIES_LIST:
    N  = sum(accum[sp][mo][0] for mo in range(1, 13))
    S  = sum(accum[sp][mo][1] for mo in range(1, 13))
    S2 = sum(accum[sp][mo][2] for mo in range(1, 13))
    if N > 0:
        mean = S / N
        std  = np.sqrt(max(S2 / N - mean ** 2, 0.0))
        print(f"{SPECIES_NAMES[sp]:<25} {int(N):>10} {mean:>12.5f} {std:>10.5f}")
