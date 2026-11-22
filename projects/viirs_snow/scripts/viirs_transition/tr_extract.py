#!/usr/bin/env python3
"""Analysis E step 1: daily snow state of the 5 arms around 1 Aug (the 2023 MODIS->VIIRS switch and the same calendar
window of 2024 and 2025), from the tavg24 daily files (SNOMASLAND = SWE kg m-2, FRLANDSNO = snow-covered fraction).

Per day and latitude band (S_hemisphere, N_0_30, N_30_45, N_45_60, N_60_90), for each arm:
    mean SWE, mean SCF, number of tiles with SWE > 1 kg m-2 (n_snow), mean of SWE over snow tiles is derived later.
Per day and band for each pair of arms (A - B) in PAIRS:  mean diff, rms diff (SWE, SCF), n tiles with |dSWE| > 1.
Per-tile SWE snapshots (2023 only, every 7 days from 07-31) for maps.
Output: output/transition/tr_<year>.nc (git-ignored).
Run: python3 -I tr_extract.py --year 2023 [--max-days N]
"""
import argparse, datetime as dt, json, os
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

REPO = "/discover/nobackup/projects/land_da/geosldas-analysis"
CFG = f"{REPO}/projects/viirs_snow/config/arms.json"
N_TILE = 112573
EDGES = [-90.0, 0.0, 30.0, 45.0, 60.0, 90.01]
BANDS = ["S_hemisphere", "N_0_30", "N_30_45", "N_45_60", "N_60_90"]
ARMS = ["OL", "T_to_J2", "TA_to_J1J2", "T_control", "TA_control"]
PAIRS = [("T_control", "T_to_J2"), ("TA_control", "TA_to_J1J2"),      # C1, C2 (same start state)
         ("T_to_J2", "OL"), ("TA_to_J1J2", "OL"), ("T_control", "OL"), ("TA_control", "OL"),
         ("TA_control", "T_control"), ("TA_to_J1J2", "T_to_J2")]
MM1 = 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--max-days", type=int)
    ap.add_argument("--out")
    a = ap.parse_args()
    cfg = json.load(open(CFG))
    d0, d1 = dt.date(a.year, 7, 15), dt.date(a.year, 12, 15)
    days = [d0 + dt.timedelta(i) for i in range((d1 - d0).days + 1)][: a.max_days]
    nd = len(days)
    root = {k: f"{cfg['project_dir']}/{v['exp_id']}/output/{cfg['domain']}" for k, v in cfg["arms"].items()}
    exp = {k: v["exp_id"] for k, v in cfg["arms"].items()}
    arm_beg = {k: dt.date.fromisoformat(v["beg"]) for k, v in cfg["arms"].items()}
    nb = len(BANDS)
    st = {v: np.full((nd, len(ARMS), nb), np.nan) for v in ["swe", "scf", "nsnow"]}
    pr = {v: np.full((nd, len(PAIRS), nb), np.nan) for v in ["dswe_mean", "dswe_rms", "dscf_mean", "dscf_rms", "n_dswe1"]}
    snaps, snap_days = [], []
    lat = None
    for i, day in enumerate(days):
        ymd = day.strftime("%Y%m%d")
        sw, sc = {}, {}
        for k in ARMS:
            if day < arm_beg[k]:
                continue
            f = f"{root[k]}/cat/ens0000/Y{day.year}/M{day.month:02d}/{exp[k]}.tavg24_1d_lnd_Nt.{ymd}_1200z.nc4"
            try:
                with Dataset(f) as d:
                    sw[k] = np.ma.filled(d["SNOMASLAND"][:].astype(np.float64), np.nan).reshape(-1)
                    sc[k] = np.ma.filled(d["FRLANDSNO"][:].astype(np.float64), np.nan).reshape(-1)
                    if lat is None:
                        lat = d["lat"][:].astype(np.float64)
                        bidx = np.digitize(lat, EDGES) - 1
            except (FileNotFoundError, OSError):
                pass
        if lat is None:
            continue
        for j, k in enumerate(ARMS):
            if k not in sw:
                continue
            for b in range(nb):
                m = bidx == b
                st["swe"][i, j, b] = np.nanmean(sw[k][m]); st["scf"][i, j, b] = np.nanmean(sc[k][m])
                st["nsnow"][i, j, b] = np.nansum(sw[k][m] > MM1)
        for j, (A, B) in enumerate(PAIRS):
            if A not in sw or B not in sw:
                continue
            ds, dc = sw[A] - sw[B], sc[A] - sc[B]
            for b in range(nb):
                m = bidx == b
                pr["dswe_mean"][i, j, b] = np.nanmean(ds[m]); pr["dswe_rms"][i, j, b] = np.sqrt(np.nanmean(ds[m] ** 2))
                pr["dscf_mean"][i, j, b] = np.nanmean(dc[m]); pr["dscf_rms"][i, j, b] = np.sqrt(np.nanmean(dc[m] ** 2))
                pr["n_dswe1"][i, j, b] = np.nansum(np.abs(ds[m]) > MM1)
        if a.year == 2023 and (day - dt.date(2023, 7, 31)).days % 7 == 0 and len(sw) == len(ARMS):
            snaps.append(np.stack([sw[k] for k in ARMS]).astype(np.float32)); snap_days.append(ymd)
        if i % 20 == 0:
            print(ymd, flush=True)
    out = a.out or f"{REPO}/projects/viirs_snow/output/transition/tr_{a.year}.nc"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with Dataset(out, "w") as o:
        o.createDimension("day", nd); o.createDimension("arm", len(ARMS)); o.createDimension("pair", len(PAIRS))
        o.createDimension("band", nb); o.createDimension("tile", N_TILE)
        o.arms, o.pairs, o.bands = json.dumps(ARMS), json.dumps(PAIRS), json.dumps(BANDS)
        o.createVariable("date", "i4", ("day",))[:] = [int(x.strftime("%Y%m%d")) for x in days]
        for k, v in st.items():
            o.createVariable(k, "f8", ("day", "arm", "band"))[:] = v
        for k, v in pr.items():
            o.createVariable(k, "f8", ("day", "pair", "band"))[:] = v
        if lat is not None:
            o.createVariable("lat", "f4", ("tile",))[:] = lat
        if snaps:
            o.createDimension("snap", len(snaps))
            o.createVariable("snap_date", "i4", ("snap",))[:] = [int(x) for x in snap_days]
            o.createVariable("snap_swe", "f4", ("snap", "arm", "tile"), zlib=True)[:] = np.stack(snaps)


if __name__ == "__main__":
    main()
