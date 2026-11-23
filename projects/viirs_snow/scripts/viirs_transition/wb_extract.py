#!/usr/bin/env python3
"""Analysis D, step 1: per-tile snow-year water-balance sums for one arm.

Reads the daily tavg24 files and the daily catch_progn_incr files of one arm (config/arms.json) and writes one
NetCDF with, per snow year (1 Aug - 31 Jul) and tile:

  annual sums (kg m-2): P, ET (EVLAND), RUNSURF, BASEFLOW, SNOWMELT (SMLAND), QINFIL, FSW (PEATCLSM_FSWCHANGE),
      WCHANGE (WCHANGELAND), I_snow (sum of WESNN1-3_INCR), I_snow_abs (sum over steps of |net layer incr|),
      n_incr (steps with a non-zero snow increment)
  annual state: SWE_max, SWE_end (31 Jul), n_days (days read)
  monthly (Aug..Jul): SWE, RZMC, SFMC, TWLAND, SCF (means); ET, RUNOFF, SNOWMELT, P, I_snow (totals)
  restart endpoints: TWLAND at 1 Aug of each year, from catch_internal_rst (single member):
      TWLAND = CDCR2/(1-WPWET) - CATDEF + RZEXC + SRFEXC + CAPAC + sum(WESNN)   [same as M21C_ls
      reconstruct_twland_endpoints.py]

Flux definitions follow projects/M21C_ls/scripts/water_year_snow_da_budget.py. Run with the GEOSpyD python3 -I.
"""
import argparse
import datetime as dt
import glob
import json
import os

import numpy as np
from netCDF4 import Dataset

REPO = "/discover/nobackup/projects/land_da/geosldas-analysis"
CFG = f"{REPO}/projects/viirs_snow/config/arms.json"
N_TILE = 112573
SEC = 86400.0
MONTHS = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]


def rd(ds, name):
    return np.ma.filled(ds.variables[name][:].astype(np.float64), np.nan).reshape(-1, N_TILE)


def twland_rst(path):
    with Dataset(path) as d:
        g = lambda v: np.ma.filled(d.variables[v][:].astype(np.float64), np.nan)
        return (g("CDCR2") / (1.0 - g("WPWET")) - g("CATDEF") + g("RZEXC") + g("SRFEXC") + g("CAPAC")
                + g("WESNN1") + g("WESNN2") + g("WESNN3"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sy", nargs="*", help="subset of snow-year ids (testing)")
    ap.add_argument("--rst-only", action="store_true", help="only (re)write TWLAND_rst into an existing --out file")
    ap.add_argument("--max-days", type=int, default=None, help="testing: only the first N days of each SY")
    a = ap.parse_args()
    cfg = json.load(open(CFG))
    arm = cfg["arms"][a.arm]
    exp = arm["exp_id"]
    root = f"{cfg['project_dir']}/{exp}/output/{cfg['domain']}"
    beg, end = dt.date.fromisoformat(arm["beg"]), dt.date.fromisoformat(arm["end"])
    sys_ = [s for s in cfg["snow_years"]
            if dt.date.fromisoformat(s["beg"]) >= beg and dt.date.fromisoformat(s["end"]) < end
            and (not a.sy or s["id"] in a.sy)]
    nsy = len(sys_)
    if a.rst_only:
        end_years = [dt.date.fromisoformat(x["beg"]).year for x in sys_] + [dt.date.fromisoformat(sys_[-1]["end"]).year]
        rst = np.stack([twland_rst(f"{root}/rs/ens0000/Y{y}/M08/{exp}.catch_internal_rst.{y}0801_0000") for y in end_years])
        with Dataset(a.out, "a") as o:
            o["TWLAND_rst"][:] = rst
            o["endpoint_year"][:] = end_years
            o.missing = json.dumps([m for m in json.loads(o.missing) if m[0] != "rst"])
        print("rst-only done", end_years)
        return
    z2 = lambda: np.zeros((nsy, N_TILE), np.float64)
    zm = lambda: np.zeros((nsy, 12, N_TILE), np.float64)
    ann = {k: z2() for k in ["P", "ET", "RUNSURF", "BASEFLOW", "SNOWMELT", "QINFIL", "FSW", "WCHANGE",
                              "I_snow", "I_snow_abs", "n_incr", "SWE_end"]}
    ann["SWE_max"] = np.full((nsy, N_TILE), -np.inf)
    nd = np.zeros(nsy, int)
    mmean = {k: zm() for k in ["SWE", "RZMC", "SFMC", "TWLAND", "SCF"]}
    mtot = {k: zm() for k in ["ET", "RUNOFF", "SNOWMELT", "P", "I_snow"]}
    mdays = np.zeros((nsy, 12), int)
    missing = []
    for i, s in enumerate(sys_):
        d0, d1 = dt.date.fromisoformat(s["beg"]), dt.date.fromisoformat(s["end"])
        day, k = d0, 0
        while day <= d1:
            if a.max_days and k >= a.max_days:
                break
            k += 1
            ymd = day.strftime("%Y%m%d")
            base = f"{root}/cat/ens0000/Y{day.year}/M{day.month:02d}/{exp}"
            fm, fi = f"{base}.tavg24_1d_lnd_Nt.{ymd}_1200z.nc4", f"{base}.catch_progn_incr.{ymd}.nc4"
            m = MONTHS.index(day.month)
            try:
                with Dataset(fm) as d:
                    P, ET = rd(d, "PRECTOTCORRLAND")[0] * SEC, rd(d, "EVLAND")[0] * SEC
                    RS, BF = rd(d, "RUNSURFLAND")[0] * SEC, rd(d, "BASEFLOWLAND")[0] * SEC
                    SM, QI = rd(d, "SMLAND")[0] * SEC, rd(d, "QINFILLAND")[0] * SEC
                    FSW, WC = rd(d, "PEATCLSM_FSWCHANGE")[0] * SEC, rd(d, "WCHANGELAND")[0] * SEC
                    SWE = rd(d, "SNOMASLAND")[0]
                    for n in ["RZMC", "SFMC", "TWLAND"]:
                        mmean[n][i, m] += rd(d, n)[0]
                    mmean["SCF"][i, m] += rd(d, "FRLANDSNO")[0]
            except (FileNotFoundError, OSError):
                missing.append(("tavg24", ymd, a.arm))
                day += dt.timedelta(days=1)
                continue
            Is = Ia = Ni = 0.0
            if a.arm != "OL":  # OL has no analysis increments
                try:
                    if os.path.exists(fi):
                        with Dataset(fi) as d:
                            inc = sum(rd(d, f"WESNN{j}_INCR") for j in (1, 2, 3))  # (8, tile)
                    else:  # first day of a run: only the 3-hourly files exist
                        parts = sorted(glob.glob(f"{base}.catch_progn_incr.{ymd}_*z.nc4"))
                        if not parts:
                            raise FileNotFoundError(fi)
                        inc = np.concatenate([sum(rd(Dataset(q), f"WESNN{j}_INCR") for j in (1, 2, 3)) for q in parts])
                        missing.append(("incr_3hourly_files_only", ymd, len(parts)))
                    inc = np.nan_to_num(inc)
                    Is, Ia, Ni = inc.sum(0), np.abs(inc).sum(0), (inc != 0).sum(0)
                except (FileNotFoundError, OSError):
                    missing.append(("incr", ymd, a.arm))
            for key, v in [("P", P), ("ET", ET), ("RUNSURF", RS), ("BASEFLOW", BF), ("SNOWMELT", SM),
                           ("QINFIL", QI), ("FSW", FSW), ("WCHANGE", WC), ("I_snow", Is), ("I_snow_abs", Ia),
                           ("n_incr", Ni)]:
                ann[key][i] += np.nan_to_num(v)
            ann["SWE_max"][i] = np.fmax(ann["SWE_max"][i], SWE)
            ann["SWE_end"][i] = SWE
            for key, v in [("ET", ET), ("RUNOFF", RS + BF), ("SNOWMELT", SM), ("P", P), ("I_snow", Is)]:
                mtot[key][i, m] += np.nan_to_num(v)
            mmean["SWE"][i, m] += SWE
            mdays[i, m] += 1
            nd[i] += 1
            day += dt.timedelta(days=1)
        print(f"{a.arm} {s['id']}: {nd[i]} days", flush=True)
    for n in mmean:
        mmean[n] /= np.maximum(mdays, 1)[:, :, None]
    ann["SWE_max"][~np.isfinite(ann["SWE_max"])] = np.nan
    # restart endpoints at 1 Aug
    end_years = [dt.date.fromisoformat(s["beg"]).year for s in sys_] + [dt.date.fromisoformat(sys_[-1]["end"]).year]
    rst = np.full((len(end_years), N_TILE), np.nan)
    for j, y in enumerate(end_years):
        p = f"{root}/rs/ens0000/Y{y}/M08/{exp}.catch_internal_rst.{y}0801_0000"
        if os.path.exists(p):
            rst[j] = twland_rst(p)
        else:
            missing.append(("rst", f"{y}0801", a.arm))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with Dataset(a.out, "w") as o:
        o.createDimension("sy", nsy); o.createDimension("month", 12); o.createDimension("tile", N_TILE)
        o.createDimension("endpoint", len(end_years))
        o.arm, o.exp_id, o.missing = a.arm, exp, json.dumps(missing)
        o.units_note = "annual/monthly totals in kg m-2; monthly means of state variables in native units"
        v = o.createVariable("sy_id", str, ("sy",)); v[:] = np.array([s["id"] for s in sys_], dtype=object)
        o.createVariable("n_days", "i4", ("sy",))[:] = nd
        o.createVariable("endpoint_year", "i4", ("endpoint",))[:] = end_years
        o.createVariable("month_order", "i4", ("month",))[:] = MONTHS
        o.createVariable("TWLAND_rst", "f4", ("endpoint", "tile"), zlib=True)[:] = rst
        for k, x in ann.items():
            o.createVariable(k, "f4", ("sy", "tile"), zlib=True)[:] = x
        for k, x in mmean.items():
            o.createVariable(f"m_{k}", "f4", ("sy", "month", "tile"), zlib=True)[:] = x
        for k, x in mtot.items():
            o.createVariable(f"m_{k}", "f4", ("sy", "month", "tile"), zlib=True)[:] = x
        for n in ["lat", "lon"]:
            with Dataset(fm) as d:
                o.createVariable(n, "f4", ("tile",))[:] = d.variables[n][:]
    print("missing:", missing)


if __name__ == "__main__":
    main()
