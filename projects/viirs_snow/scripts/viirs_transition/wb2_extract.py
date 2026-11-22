#!/usr/bin/env python3
"""Analysis D additions, step 1: fluxes on the input-gap / boundary days, summed per snow year, for one arm.

The annual sums from wb_extract.py include days on which some assimilated sensor had no input file (see plan,
"Input data gaps") and the first day of each run segment. This script re-reads only those days (daily tavg24 and
catch_progn_incr files) and stores their per-snow-year sums of P, ET, RUNOFF, SNOWMELT, FSW, I_snow so that
wb2_analyze.py can subtract them from every arm (common mask; storage cannot be adjusted this way, so only flux terms
are corrected). Same definitions as wb_extract.py. Run with the GEOSpyD python3 -I.
"""
import argparse, datetime as dt, glob, json, os
import numpy as np
from netCDF4 import Dataset

REPO = "/discover/nobackup/projects/land_da/geosldas-analysis"
CFG = f"{REPO}/projects/viirs_snow/config/arms.json"
N_TILE = 112573
SEC = 86400.0
# Aqua 2022-04-01..16 (missing source files; analysis A found 04-02..16), Terra 2022-10-11..22 and 2025-05-11..14,
# J2 2024-02-28, plus the segment boundary 2023-08-01.
GAPS = ([dt.date(2022, 4, 1) + dt.timedelta(d) for d in range(16)] + [dt.date(2022, 10, 11) + dt.timedelta(d) for d in range(12)]
        + [dt.date(2025, 5, 11) + dt.timedelta(d) for d in range(4)] + [dt.date(2024, 2, 28), dt.date(2023, 8, 1)])


def rd(ds, n):
    return np.ma.filled(ds.variables[n][:].astype(np.float64), np.nan).reshape(-1, N_TILE)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--arm", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = json.load(open(CFG)); arm = cfg["arms"][a.arm]; exp = arm["exp_id"]
    root = f"{cfg['project_dir']}/{exp}/output/{cfg['domain']}"
    beg, end = dt.date.fromisoformat(arm["beg"]), dt.date.fromisoformat(arm["end"])
    sys_ = [s for s in cfg["snow_years"] if dt.date.fromisoformat(s["beg"]) >= beg and dt.date.fromisoformat(s["end"]) < end]
    keys = ["P", "ET", "RUNOFF", "SNOWMELT", "FSW", "I_snow"]
    acc = {k: np.zeros((len(sys_), N_TILE)) for k in keys}
    nd = np.zeros(len(sys_), int); used = []
    for day in GAPS:
        for i, s in enumerate(sys_):
            if dt.date.fromisoformat(s["beg"]) <= day <= dt.date.fromisoformat(s["end"]):
                break
        else:
            continue
        ymd = day.strftime("%Y%m%d")
        base = f"{root}/cat/ens0000/Y{day.year}/M{day.month:02d}/{exp}"
        try:
            with Dataset(f"{base}.tavg24_1d_lnd_Nt.{ymd}_1200z.nc4") as d:
                v = {"P": rd(d, "PRECTOTCORRLAND")[0] * SEC, "ET": rd(d, "EVLAND")[0] * SEC,
                     "RUNOFF": (rd(d, "RUNSURFLAND")[0] + rd(d, "BASEFLOWLAND")[0]) * SEC,
                     "SNOWMELT": rd(d, "SMLAND")[0] * SEC, "FSW": rd(d, "PEATCLSM_FSWCHANGE")[0] * SEC}
        except OSError:
            continue
        Is = 0.0
        if a.arm != "OL":
            fi = f"{base}.catch_progn_incr.{ymd}.nc4"
            try:
                if os.path.exists(fi):
                    with Dataset(fi) as d:
                        Is = np.nan_to_num(sum(rd(d, f"WESNN{j}_INCR") for j in (1, 2, 3))).sum(0)
                else:
                    parts = sorted(glob.glob(f"{base}.catch_progn_incr.{ymd}_*z.nc4"))
                    Is = np.nan_to_num(np.concatenate([sum(rd(Dataset(q), f"WESNN{j}_INCR") for j in (1, 2, 3)) for q in parts])).sum(0)
            except (OSError, ValueError):
                pass
        v["I_snow"] = Is
        for k in keys:
            acc[k][i] += np.nan_to_num(v[k])
        nd[i] += 1; used.append(ymd)
    with Dataset(a.out, "w") as o:
        o.createDimension("sy", len(sys_)); o.createDimension("tile", N_TILE)
        o.arm, o.days = a.arm, json.dumps(used)
        o.createVariable("n_gap_days", "i4", ("sy",))[:] = nd
        for k in keys:
            o.createVariable(k, "f4", ("sy", "tile"), zlib=True)[:] = acc[k]
    print(a.arm, "gap days used:", len(used), nd)


if __name__ == "__main__":
    main()
