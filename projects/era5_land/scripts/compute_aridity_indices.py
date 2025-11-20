#!/usr/bin/env python3
"""
Compute aridity indices (AI, Budyko φ, CMI) from GEOS-LDAS daily tile files
over a YYYYY/Mmm directory tree, between a start and end date.

Assumes file names like: OLv8_M36_cd.tavg24_1d_lfs_Nt.20221222_1200z*.nc4
and daily variables: Tair, Qair, Wind, RefH, Psurf, SWdown, LWdown, HLWUP, RainfSnowf.
"""

import re
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

# ------------------- CONFIG -------------------
BASE = Path("/discover/nobackup/projects/land_da/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg")
START = pd.Timestamp("2018-08-01")
END   = pd.Timestamp("2024-06-30")

# file matching inside YYYYY/Mmm:
FILE_GLOB = "OLv8_M36_cd.tavg24_1d_*.*"   # keep broad; handles .nc/.nc4/.nc4.nc etc.

# chunking for performance (adjust for your memory)
CHUNKS = {"time": 64, "tile": 10000}

# Mask threshold to avoid divisions by ~0 in cold regions (annual totals)
EPS_MM = 1e-3

# Reference albedo for FAO-56 (grass)
ALPHA_REF = 0.23
# ---------------------------------------------

DATE_RE = re.compile(r"\.(\d{8})_")   # captures YYYYMMDD in e.g., ... .20221222_1200z ...

def collect_files(base: Path, start: pd.Timestamp, end: pd.Timestamp) -> list[str]:
    files = []
    for ydir in sorted(base.glob("Y*/")):
        for mdir in sorted(ydir.glob("M*/")):
            for f in mdir.glob(FILE_GLOB):
                m = DATE_RE.search(f.name)
                if not m:
                    continue
                dt = pd.to_datetime(m.group(1), format="%Y%m%d")
                if start <= dt <= end:
                    files.append(str(f))
    files.sort()
    if not files:
        raise SystemExit("No files found in range. Check BASE/START/END/FILE_GLOB.")
    return files

def main():
    files = collect_files(BASE, START, END)
    print(f"Found {len(files)} daily files between {START.date()} and {END.date()}.")

    # Open with xarray/dask
    ds = xr.open_mfdataset(
        files, combine="by_coords", parallel=True,
        decode_times=True, chunks=CHUNKS, engine="netcdf4"
    )

    FV = 1e15
    def clean(v): return ds[v].where(ds[v] < FV)

    # Inputs (daily means)
    Tair = clean("Tair") - 273.15               # °C
    qair = clean("Qair")                        # kg/kg
    wind_z = clean("Wind")                      # m/s at RefH
    zref = clean("RefH").fillna(2.0)            # m
    ps_kpa = (clean("Psurf") / 1000.0)          # Pa -> kPa
    sw_down = clean("SWdown")                   # W/m2
    lw_down = clean("LWdown")                   # W/m2
    lw_up   = clean("HLWUP")                    # W/m2
    P_rate  = clean("RainfSnowf")               # kg m-2 s-1 (== mm s-1)

    # Radiation → net Rn (MJ m-2 day-1)
    Rn_W = (1.0 - ALPHA_REF) * sw_down + lw_down - lw_up
    Rn_MJ = (Rn_W * 86400.0) / 1e6

    # Thermodynamics
    es = 0.6108 * np.exp(17.27 * Tair / (Tair + 237.3))      # kPa
    ea = (qair * ps_kpa) / (0.622 + 0.378 * qair)            # kPa
    vpd = (es - ea).clip(min=0.0)

    delta = 4098.0 * es / (Tair + 237.3) ** 2                # kPa / °C
    gamma = 0.000665 * ps_kpa                                 # kPa / °C

    # Wind at 2 m (FAO-56)
    z = xr.where(zref > 0.5, zref, 2.0)
    u2 = wind_z * (4.87 / np.log(67.8 * z - 5.42))

    # FAO-56 Penman–Monteith (mm day-1). Ground heat flux ~ 0 at daily.
    ET0_day = (
        0.408 * delta * (Rn_MJ) +
        gamma * (900.0 / (Tair + 273.0)) * u2 * vpd
    ) / (delta + gamma * (1.0 + 0.34 * u2))
    PET_day = ET0_day.clip(min=0.0)               # mm/day
    P_day   = (P_rate * 86400.0).clip(min=0.0)    # mm/day

    # --- Monthly & Annual aggregation ---
    P_mon   = P_day.resample(time="MS").sum()
    PET_mon = PET_day.resample(time="MS").sum()

    P_ann   = P_mon.resample(time="YS").sum()
    PET_ann = PET_mon.resample(time="YS").sum()

    # Indices (annual)
    AI_ann  = (P_ann / PET_ann).where(PET_ann > EPS_MM).rename("AI")
    PHI_ann = (PET_ann / P_ann).where(P_ann > EPS_MM).rename("Budyko_phi")
    CMI_ann = ((P_ann - PET_ann) / PET_ann.where(PET_ann > EPS_MM)).rename("CMI")

    # Climatologies (mean across available years)
    AI_clim  = AI_ann.mean("time").rename("AI_clim")
    PHI_clim = PHI_ann.mean("time").rename("Budyko_phi_clim")
    CMI_clim = CMI_ann.mean("time").rename("CMI_clim")

    # De Martonne (optional)
    T_ann = Tair.resample(time="YS").mean()
    IDM_ann = (P_ann / (T_ann + 10.0)).rename("DeMartonne")
    IDM_clim = IDM_ann.mean("time").rename("DeMartonne_clim")

    # Pack outputs
    out = xr.Dataset(
        data_vars=dict(
            AI=AI_ann, Budyko_phi=PHI_ann, CMI=CMI_ann,
            DeMartonne=IDM_ann,
            AI_clim=AI_clim, Budyko_phi_clim=PHI_clim, CMI_clim=CMI_clim,
            DeMartonne_clim=IDM_clim,
            P_annual=P_ann.rename("P_annual"),
            PET_annual=PET_ann.rename("PET_annual"),
        ),
        coords=dict(
            time=AI_ann.time,
            tile=ds["tile"] if "tile" in ds.coords else np.arange(ds.dims["tile"]),
            lat=ds["lat"], lon=ds["lon"]
        ),
        attrs=dict(
            description="Aridity indices from GEOS-LDAS daily tile output using FAO-56 Penman–Monteith",
            period_start=str(START.date()), period_end=str(END.date()),
            albedo_reference=ALPHA_REF, pet_method="FAO-56 Penman–Monteith",
            precip_var="RainfSnowf", note="Annual sums/means over resampled calendar years."
        )
    )

    # Add simple year coordinate
    out = out.assign_coords(year=("time", out["time.year"].values))

    # Write a single NetCDF
    outfile = BASE / f"aridity_indices_{START:%Y%m%d}_{END:%Y%m%d}.nc4"
    encoding = {v: {"zlib": True, "complevel": 3} for v in out.data_vars}
    out.to_netcdf(outfile, encoding=encoding)
    print(f"Wrote {outfile}")

if __name__ == "__main__":
    main()
