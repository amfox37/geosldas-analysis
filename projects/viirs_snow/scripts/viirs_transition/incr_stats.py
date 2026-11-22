#!/usr/bin/env python3
"""Analysis B of docs/viirs_transition_plan.md: snow-increment structure per arm, from catch_progn_incr.

Reads the daily cat/ens0000 files (8 instantaneous 3-hourly steps, WESNN1..3_INCR in kg m-2) of one arm of
config/arms.json. Per tile and step the SWE increment is the sum over the 3 snow layers; a step "has an increment" if
|incr| >= THR (0.01 kg m-2; exact zeros are the no-obs/no-update case, tiny values are numerical residue). A step is at
the cap if |incr| >= CAP_TEST (0.999 x SCF_ANA_MAXINCRSWE = 5 kg m-2).
Writes output/incr_stats/incr_stats_<arm>.nc (git-ignored):
  per day x band x lon-sector (45 deg, 8 sectors):
     n_tiles_by_nsteps[day, 4, band, lon]  tiles with 0 / 1 / 2 / >=3 steps with an increment
     dbl_* : among tiles with >=2 increment steps, the first two: count, same-sign count, sum|first|, sum|second|,
             sum net (first+second), sum abs
     step_* [day, step, band, lon]: n_nonzero, sum net, sum abs, n_cap
     dbl_gap[day, 7]: step distance between first and second increment (1..7)
  per tile x snow year: net, abs, n_incr_steps, n_incr_days, n_double_days, n_cap_steps
  per day: n_steps_found
Run on Discover (GEOSpyD py3.13):  python3 -I incr_stats.py --arm TA_to_J1J2 [--beg D --end D --out F]
"""
import argparse, datetime as dt, json, time
from pathlib import Path
import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
THR = 0.01
CAP = 5.0
CAP_TEST = 0.999 * CAP
BAND_EDGES = [-90.0, 0.0, 30.0, 45.0, 60.0, 90.01]
BAND_NAMES = ["S_hemisphere", "N_0_30", "N_30_45", "N_45_60", "N_60_90"]
NB, NL = 5, 8


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--arm", required=True)
    p.add_argument("--config", default=str(PROJECT / "config" / "arms.json"))
    p.add_argument("--beg"); p.add_argument("--end"); p.add_argument("--out")
    a = p.parse_args()
    cfg = json.load(open(a.config))
    arm = cfg["arms"][a.arm]
    exp = arm["exp_id"]
    ntile = cfg["n_tiles"]
    syears = cfg["snow_years"]
    root = Path(cfg["project_dir"]) / exp / "output" / cfg["domain"] / "cat" / "ens0000"
    beg = dt.date.fromisoformat(a.beg or arm["beg"])
    end = dt.date.fromisoformat(a.end or arm["end"])
    days = [beg + dt.timedelta(i) for i in range((end - beg).days + 1)]
    nd = len(days)

    def sy_of(d):
        for i, s in enumerate(syears):
            if dt.date.fromisoformat(s["beg"]) <= d <= dt.date.fromisoformat(s["end"]):
                return i
        return -1

    ng = NB * NL
    ntn = np.zeros((nd, 4, ng), np.int32)
    dbl = np.zeros((nd, 6, ng), np.float64)   # count, same sign, sum|f|, sum|s|, sum net, sum abs
    sn = np.zeros((nd, 8, ng), np.int32)
    ssum = np.zeros((nd, 8, ng), np.float64)
    sabs = np.zeros((nd, 8, ng), np.float64)
    scap = np.zeros((nd, 8, ng), np.int32)
    gap = np.zeros((nd, 7), np.int32)
    nfound = np.zeros(nd, np.int8)
    nsy = len(syears)
    t_net = np.zeros((nsy, ntile), np.float64); t_abs = np.zeros_like(t_net)
    t_ns = np.zeros((nsy, ntile), np.int32); t_nd = np.zeros_like(t_ns); t_ndd = np.zeros_like(t_ns); t_cap = np.zeros_like(t_ns)
    grp = None
    lat = lon = None
    t0 = time.time()
    for di, day in enumerate(days):
        f = root / f"Y{day.year:04d}" / f"M{day.month:02d}" / f"{exp}.catch_progn_incr.{day:%Y%m%d}.nc4"
        if not f.exists():
            continue
        d = netCDF4.Dataset(f); d.set_auto_mask(False)
        if grp is None:
            lat = d["lat"][:].astype(np.float32); lon = d["lon"][:].astype(np.float32)
            band = np.digitize(lat, BAND_EDGES) - 1
            sec = np.clip(((lon + 180.0) // 45.0).astype(int), 0, NL - 1)
            grp = band * NL + sec
        w = sum(d[f"WESNN{i}_INCR"][:].astype(np.float64) for i in (1, 2, 3))
        d.close()
        w[~np.isfinite(w) | (np.abs(w) > 1e10)] = 0.0
        ns_found = w.shape[0]
        nfound[di] = ns_found
        nz = np.abs(w) >= THR
        cap = np.abs(w) >= CAP_TEST
        n = nz.sum(0)
        for k, m in enumerate((n == 0, n == 1, n == 2, n >= 3)):
            ntn[di, k] = np.bincount(grp[m], minlength=ng)
        for s in range(ns_found):
            sn[di, s] = np.bincount(grp[nz[s]], minlength=ng)
            ssum[di, s] = np.bincount(grp[nz[s]], weights=w[s][nz[s]], minlength=ng)
            sabs[di, s] = np.bincount(grp[nz[s]], weights=np.abs(w[s])[nz[s]], minlength=ng)
            scap[di, s] = np.bincount(grp[cap[s]], minlength=ng)
        m2 = n >= 2
        if m2.any():
            idx = np.flatnonzero(m2)
            order = np.argsort(~nz[:, idx], axis=0, kind="stable")      # nonzero steps first, in time order
            i1, i2 = order[0], order[1]
            fv = w[i1, idx]; sv = w[i2, idx]
            g = grp[idx]
            dbl[di, 0] = np.bincount(g, minlength=ng)
            dbl[di, 1] = np.bincount(g, weights=(fv * sv > 0).astype(float), minlength=ng)
            dbl[di, 2] = np.bincount(g, weights=np.abs(fv), minlength=ng)
            dbl[di, 3] = np.bincount(g, weights=np.abs(sv), minlength=ng)
            dbl[di, 4] = np.bincount(g, weights=fv + sv, minlength=ng)
            dbl[di, 5] = np.bincount(g, weights=np.abs(fv) + np.abs(sv), minlength=ng)
            gap[di] = np.bincount(np.clip(i2 - i1, 1, 7) - 1, minlength=7)
        sy = sy_of(day)
        if sy >= 0:
            wz = np.where(nz, w, 0.0)
            t_net[sy] += wz.sum(0); t_abs[sy] += np.abs(wz).sum(0)
            t_ns[sy] += n; t_nd[sy] += (n >= 1); t_ndd[sy] += (n >= 2); t_cap[sy] += cap.sum(0)
        if di % 100 == 0:
            print(f"{exp} {day} {time.time() - t0:7.1f} s", flush=True)

    out = Path(a.out) if a.out else PROJECT / "output" / "incr_stats" / f"incr_stats_{a.arm}.nc"
    out.parent.mkdir(parents=True, exist_ok=True)
    ds = netCDF4.Dataset(out, "w", format="NETCDF4")
    for n_, s_ in (("day", nd), ("nsteps_class", 4), ("dblvar", 6), ("step", 8), ("gap", 7), ("grp", ng), ("snowyear", nsy), ("tile", ntile)):
        ds.createDimension(n_, s_)
    z = dict(zlib=True, complevel=4)
    ds.setncattr("THR_kg_m-2", THR); ds.setncattr("CAP_TEST_kg_m-2", CAP_TEST)
    ds.setncattr("grp_definition", "grp = band*8 + lon_sector; bands " + ",".join(BAND_NAMES) + "; lon sector = floor((lon+180)/45) (0 = 180W..135W)")
    ds.setncattr("nsteps_classes", "0,1,2,>=3 steps with an increment")
    ds.setncattr("dblvar_order", "count,same_sign,sum_abs_first,sum_abs_second,sum_net_pair,sum_abs_pair")
    v = ds.createVariable("date", "i4", ("day",)); v[:] = [int(x.strftime("%Y%m%d")) for x in days]
    for name, arr, dims in (("n_tiles_by_nsteps", ntn, ("day", "nsteps_class", "grp")), ("dbl", dbl, ("day", "dblvar", "grp")),
                            ("step_n", sn, ("day", "step", "grp")), ("step_net", ssum, ("day", "step", "grp")),
                            ("step_abs", sabs, ("day", "step", "grp")), ("step_cap", scap, ("day", "step", "grp")),
                            ("dbl_gap", gap, ("day", "gap")), ("n_steps_found", nfound, ("day",)),
                            ("tile_net", t_net, ("snowyear", "tile")), ("tile_abs", t_abs, ("snowyear", "tile")),
                            ("tile_n_steps", t_ns, ("snowyear", "tile")), ("tile_n_days", t_nd, ("snowyear", "tile")),
                            ("tile_n_double_days", t_ndd, ("snowyear", "tile")), ("tile_n_cap", t_cap, ("snowyear", "tile"))):
        v = ds.createVariable(name, "f8" if arr.dtype == np.float64 and name.startswith(("dbl", "step_net", "step_abs")) else
                              ("f4" if name in ("tile_net", "tile_abs") else arr.dtype.str.replace("<", "")), dims, **z)
        v[:] = arr
    if lat is not None:
        for name, arr in (("lat", lat), ("lon", lon)):
            v = ds.createVariable(name, "f4", ("tile",)); v[:] = arr
    ds.close()
    print("wrote", out, f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
