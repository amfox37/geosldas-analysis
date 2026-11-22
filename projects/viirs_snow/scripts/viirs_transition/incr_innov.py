#!/usr/bin/env python3
"""Analysis B (second-sensor innovations) of docs/viirs_transition_plan.md, from the 3-hourly ObsFcstAna files.

Per arm, per day x species (Aqua, Terra, J1, J2; mapped via obsparam_descr) x class x latitude band:
  class 0 = all tiles with an obs of that species that day
  class 1 = tiles with obs from BOTH sensors of the pair that day (Terra+Aqua for Aqua/Terra; J1+J2 for J1/J2)
Per cell: n, sum(obs-fcst), sum((obs-fcst)^2), sum|obs-fcst|, sum(obs-ana), sum(obs), n_assim (assim_flag==1).
One value per tile-day and species (the last record of the day; window-boundary duplicates are rare).
The forecast of a later step contains the update of an earlier one, so Aqua (13:30, step after Terra) in a TA arm
vs in the OL/T arms shows how much the Terra update already moved the state toward Aqua. J1/J2 are in the same step
(averaged into one obs per tile), so J1 vs J2 innovations only show how close the two sensors are.
Output: output/incr_stats/innov_<arm>.nc.   python3 -I incr_innov.py --arm TA_control [--beg D --end D --out F]
"""
import argparse, datetime as dt, json, time
from pathlib import Path
import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
BAND_EDGES = [-90.0, 0.0, 30.0, 45.0, 60.0, 90.01]
BAND_NAMES = ["S_hemisphere", "N_0_30", "N_30_45", "N_45_60", "N_60_90"]
NB = 5
VARS = ["n", "sum_omf", "sum_omf2", "sum_abs_omf", "sum_oma", "sum_obs", "n_assim"]


def dec(x):
    return x.tobytes().decode(errors="ignore").strip("\x00 ") if hasattr(x, "tobytes") else str(x).strip("\x00 ")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--arm", required=True); p.add_argument("--config", default=str(PROJECT / "config" / "arms.json"))
    p.add_argument("--beg"); p.add_argument("--end"); p.add_argument("--out")
    a = p.parse_args()
    cfg = json.load(open(a.config)); arm = cfg["arms"][a.arm]; exp = arm["exp_id"]; ntile = cfg["n_tiles"]
    ids = cfg["obsfcstana_species_ids"]
    names = [ids[k]["name"] for k in sorted(k for k in ids if k.isdigit())]
    prod = {sp["product"]: names.index(sp["name"]) for sp in cfg["species"].values() if sp["name"] in names}
    iA, iT, iJ1, iJ2 = (names.index(n) for n in ("Aqua", "Terra", "J1", "J2"))
    ns = len(names)
    root = Path(cfg["project_dir"]) / exp / "output" / cfg["domain"]
    pattern = cfg["output_layout"]["obsfcstana"].split("  (")[0]
    beg = dt.date.fromisoformat(a.beg or arm["beg"]); end = dt.date.fromisoformat(a.end or arm["end"])
    days = [beg + dt.timedelta(i) for i in range((end - beg).days + 1)]
    nd = len(days)
    out_a = np.zeros((nd, ns, 2, NB, len(VARS)), np.float64)
    tile_band = None
    t0 = time.time()
    for di, day in enumerate(days):
        ymd = f"{day:%Y%m%d}"
        have = np.zeros((ns, ntile), bool)
        obs = np.zeros((ns, ntile)); fc = np.zeros((ns, ntile)); an = np.zeros((ns, ntile)); asf = np.zeros((ns, ntile), bool)
        for hh in range(0, 24, 3):
            f = root / pattern.replace("<yyyy>", f"{day.year:04d}").replace("<mm>", f"{day.month:02d}") \
                .replace("<exp_id>", exp).replace("<yyyymmdd>", ymd).replace("<hh>", f"{hh:02d}")
            if not f.exists():
                continue
            d = netCDF4.Dataset(f); d.set_auto_mask(False)
            sp = d["species"][:]
            if len(sp) == 0:
                d.close(); continue
            tn = d["tilenum"][:] - 1
            o = d["obs"][:]; fo = d["fcst"][:]; ao = d["ana"][:]; af = d["assim_flag"][:]
            if tile_band is None:
                tile_band = np.full(ntile, -1, np.int8)
            lat = d["lat"][:]
            tile_band[tn] = (np.digitize(lat, BAND_EDGES) - 1).astype(np.int8)
            fids = d["obsparam_species_id"][:]; fdesc = [dec(x) for x in d["obsparam_descr"][:]]
            d.close()
            for lid, s in zip(fids, fdesc):
                if s not in prod:
                    continue
                m = sp == lid
                ci = prod[s]
                t = tn[m]
                have[ci, t] = True; obs[ci, t] = o[m]; fc[ci, t] = fo[m]; an[ci, t] = ao[m]; asf[ci, t] = af[m] == 1
        if tile_band is None:
            continue
        pair = {iT: have[iT] & have[iA], iA: have[iT] & have[iA], iJ1: have[iJ1] & have[iJ2], iJ2: have[iJ1] & have[iJ2]}
        for s in range(ns):
            for c, sel in ((0, have[s]), (1, pair[s])):
                m = sel & (tile_band >= 0) & np.isfinite(obs[s]) & np.isfinite(fc[s]) & (np.abs(fc[s]) < 1e10)
                if not m.any():
                    continue
                b = tile_band[m].astype(int); omf = obs[s][m] - fc[s][m]; oma = obs[s][m] - an[s][m]
                for k, w in enumerate((np.ones(m.sum()), omf, omf ** 2, np.abs(omf), oma, obs[s][m], asf[s][m].astype(float))):
                    out_a[di, s, c, :, k] = np.bincount(b, weights=w, minlength=NB)
        if di % 100 == 0:
            print(f"{exp} {ymd} {time.time() - t0:7.1f} s", flush=True)
    out = Path(a.out) if a.out else PROJECT / "output" / "incr_stats" / f"innov_{a.arm}.nc"
    out.parent.mkdir(parents=True, exist_ok=True)
    ds = netCDF4.Dataset(out, "w", format="NETCDF4")
    for n_, s_ in (("day", nd), ("species", ns), ("cls", 2), ("band", NB), ("var", len(VARS))):
        ds.createDimension(n_, s_)
    ds.setncattr("cls", "0 = all tiles with obs of the species; 1 = tiles with obs of both sensors of the pair that day")
    ds.setncattr("var_order", ",".join(VARS)); ds.setncattr("band_names", ",".join(BAND_NAMES)); ds.setncattr("species_names", ",".join(names))
    v = ds.createVariable("date", "i4", ("day",)); v[:] = [int(f"{x:%Y%m%d}") for x in days]
    v = ds.createVariable("stats", "f8", ("day", "species", "cls", "band", "var"), zlib=True, complevel=4); v[:] = out_a
    ds.close()
    print("wrote", out, f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
