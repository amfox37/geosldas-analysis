#!/usr/bin/env python3
"""Analysis A of docs/viirs_transition_plan.md: SCF observation availability per arm.

Reads the 3-hourly ObsFcstAna files (ana/ens_avg) of one arm of config/arms.json and writes one compact NetCDF
(output/obs_availability/obs_availability_<arm>.nc):

  per day x species x latitude band: n_obs, n_assim (assim_flag == 1), n_tiles (unique tiles with >= 1 obs)
  per day:                           n_steps (3-hourly files found), tiles with obs from both Terra and Aqua,
                                     both J1 and J2, and the unions (any MODIS, any VIIRS), tiles with Terra or Aqua only
  per month x species x tile:        days with at least one obs (int16)

Species are mapped through obsparam_descr of each file (never namelist indices): Aqua MYD10C1, Terra MOD10C1,
J1 VJ110C1, J2 VJ210C1 (names and products from config/arms.json). Obs are counted by the UTC date of the
file time stamp. A missing file or a species with no obs is not an error (e.g. the OL and VIIRS arms have no
monitor-only MODIS obs from 2026-01-01, because their MODIS path ends 2025-12-31); it just shows as zero.

Run on Discover (GEOSpyD py3.13, no g5_modules needed):
  python3 obs_availability.py --arm TA_to_J1J2                      # whole record of the arm
  python3 obs_availability.py --arm TA_to_J1J2 --beg 2024-01-01 --end 2024-01-31 --out /tmp/x.nc   # test
"""
import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]                                   # projects/viirs_snow
DEFAULT_CONFIG = PROJECT / "config" / "arms.json"
DEFAULT_OUT = PROJECT / "output" / "obs_availability"

# latitude bands (degrees north); band 0 is the whole southern hemisphere
BAND_EDGES = [(-90.0, 0.0), (0.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 90.01)]
BAND_NAMES = ["S_hemisphere", "N_0_30", "N_30_45", "N_45_60", "N_60_90"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arm", required=True, help="arm key in arms.json, e.g. TA_to_J1J2")
    p.add_argument("--config", default=str(DEFAULT_CONFIG))
    p.add_argument("--beg", help="first date YYYY-MM-DD (default: arm start)")
    p.add_argument("--end", help="last date YYYY-MM-DD, inclusive (default: arm end)")
    p.add_argument("--out", help="output NetCDF (default: output/obs_availability/obs_availability_<arm>.nc)")
    return p.parse_args()


def decode(x):
    return x.tobytes().decode(errors="ignore").strip("\x00 ") if hasattr(x, "tobytes") else str(x).strip("\x00 ")


def main():
    a = parse_args()
    cfg = json.load(open(a.config))
    arm = cfg["arms"][a.arm]
    exp_id = arm["exp_id"]

    # species order (canonical) and product -> index; from arms.json, not hard-coded
    ids = cfg["obsfcstana_species_ids"]
    names = [ids[k]["name"] for k in sorted(k for k in ids if k.isdigit())]          # Aqua Terra J1 J2
    prod_to_idx = {}
    for sp in cfg["species"].values():
        if sp["name"] in names:
            prod_to_idx[sp["product"]] = names.index(sp["name"])
    iA, iT, iJ1, iJ2 = (names.index(n) for n in ("Aqua", "Terra", "J1", "J2"))
    ns, nb, ntile = len(names), len(BAND_EDGES), cfg["n_tiles"]

    root = Path(cfg["project_dir"]) / exp_id / "output" / cfg["domain"]
    pattern = cfg["output_layout"]["obsfcstana"].split("  (")[0]      # drop the trailing description
    beg = dt.date.fromisoformat(a.beg or arm["beg"])
    end = dt.date.fromisoformat(a.end or arm["end"])
    days = [beg + dt.timedelta(d) for d in range((end - beg).days + 1)]
    months = sorted({(d.year, d.month) for d in days})
    mindex = {m: i for i, m in enumerate(months)}
    nd, nm = len(days), len(months)

    n_obs = np.zeros((nd, ns, nb), np.int32)
    n_assim = np.zeros((nd, ns, nb), np.int32)
    n_tiles = np.zeros((nd, ns, nb), np.int32)
    n_steps = np.zeros(nd, np.int8)
    both_TA = np.zeros(nd, np.int32)
    both_J = np.zeros(nd, np.int32)
    union_MODIS = np.zeros(nd, np.int32)
    union_VIIRS = np.zeros(nd, np.int32)
    only_T = np.zeros(nd, np.int32)       # Terra but not Aqua
    only_A = np.zeros(nd, np.int32)
    only_J1 = np.zeros(nd, np.int32)
    only_J2 = np.zeros(nd, np.int32)
    days_with_obs = np.zeros((nm, ns, ntile), np.int16)
    tile_lat = np.full(ntile, np.nan, np.float32)
    tile_lon = np.full(ntile, np.nan, np.float32)
    assim_seen = np.zeros((nd, ns), bool)                           # any assim_flag == 1 that day
    mapping_checked = {}
    t0 = time.time()
    n_files = 0

    for di, day in enumerate(days):
        ymd = day.strftime("%Y%m%d")
        seen = np.zeros((ns, ntile), bool)                          # tiles with obs that day, per species
        for hh in range(0, 24, 3):
            f = root / pattern.replace("<yyyy>", f"{day.year:04d}").replace("<mm>", f"{day.month:02d}") \
                .replace("<exp_id>", exp_id).replace("<yyyymmdd>", ymd).replace("<hh>", f"{hh:02d}")
            if not f.exists():
                continue
            n_steps[di] += 1
            n_files += 1
            d = netCDF4.Dataset(f)
            d.set_auto_mask(False)
            sp = d["species"][:]
            af = d["assim_flag"][:]
            tn = d["tilenum"][:] - 1
            lat = d["lat"][:]
            lon = d["lon"][:]
            fids = d["obsparam_species_id"][:]
            fdesc = [decode(x) for x in d["obsparam_descr"][:]]
            d.close()
            # local species id in this file -> canonical index, via the descriptor
            loc2can = {int(i): prod_to_idx[s] for i, s in zip(fids, fdesc) if s in prod_to_idx}
            key = tuple(sorted(loc2can.items()))
            if key not in mapping_checked:
                mapping_checked[key] = ymd
            if len(sp) == 0:
                continue
            can = np.full(len(sp), -1, np.int8)
            for lid, cidx in loc2can.items():
                can[sp == lid] = cidx
            band = np.full(len(sp), -1, np.int8)
            for bi, (lo_, hi_) in enumerate(BAND_EDGES):
                band[(lat >= lo_) & (lat < hi_)] = bi
            ok = (can >= 0) & (band >= 0)
            c, b, t, g = can[ok], band[ok], tn[ok], af[ok]
            np.add.at(n_obs[di], (c, b), 1)
            np.add.at(n_assim[di], (c, b), (g == 1).astype(np.int32))
            for s in range(ns):
                m = c == s
                if m.any():
                    seen[s, t[m]] = True
                    tile_lat[t[m]] = lat[ok][m]
                    tile_lon[t[m]] = lon[ok][m]
                    if (g[m] == 1).any():
                        assim_seen[di, s] = True
        # per-day unique tiles, per band (tiles sit in one band; use the tile latitude)
        if n_steps[di] == 0:
            continue
        for s in range(ns):
            idx = np.flatnonzero(seen[s])
            if idx.size:
                tl = tile_lat[idx]
                for bi, (lo_, hi_) in enumerate(BAND_EDGES):
                    n_tiles[di, s, bi] = int(((tl >= lo_) & (tl < hi_)).sum())
        days_with_obs[mindex[(day.year, day.month)]] += seen.astype(np.int16)
        both_TA[di] = int((seen[iT] & seen[iA]).sum())
        both_J[di] = int((seen[iJ1] & seen[iJ2]).sum())
        union_MODIS[di] = int((seen[iT] | seen[iA]).sum())
        union_VIIRS[di] = int((seen[iJ1] | seen[iJ2]).sum())
        only_T[di] = int((seen[iT] & ~seen[iA]).sum())
        only_A[di] = int((seen[iA] & ~seen[iT]).sum())
        only_J1[di] = int((seen[iJ1] & ~seen[iJ2]).sum())
        only_J2[di] = int((seen[iJ2] & ~seen[iJ1]).sum())
        if di % 100 == 0:
            print(f"{exp_id} {ymd}  {time.time() - t0:7.1f} s  files={n_files}", flush=True)

    if len(mapping_checked) > 1:
        print("NOTE: species descriptor mapping differs between files:", mapping_checked, file=sys.stderr)

    out = Path(a.out) if a.out else DEFAULT_OUT / f"obs_availability_{a.arm}.nc"
    out.parent.mkdir(parents=True, exist_ok=True)
    ds = netCDF4.Dataset(out, "w", format="NETCDF4")
    ds.createDimension("day", nd)
    ds.createDimension("species", ns)
    ds.createDimension("band", nb)
    ds.createDimension("month", nm)
    ds.createDimension("tile", ntile)
    z = dict(zlib=True, complevel=4)
    v = ds.createVariable("date", "i4", ("day",))
    v.units = "yyyymmdd (UTC date of the ObsFcstAna file time stamp)"
    v[:] = [int(d.strftime("%Y%m%d")) for d in days]
    v = ds.createVariable("month_yyyymm", "i4", ("month",))
    v[:] = [y * 100 + m for y, m in months]
    v = ds.createVariable("species_name", str, ("species",))
    v[:] = np.array(names, dtype=object)
    v = ds.createVariable("band_name", str, ("band",))
    v[:] = np.array(BAND_NAMES, dtype=object)
    for nm_, arr, dims, desc in (
        ("n_obs", n_obs, ("day", "species", "band"), "SCF obs in all 3-hourly ObsFcstAna files of the day (monitored or assimilated)"),
        ("n_assim", n_assim, ("day", "species", "band"), "obs with assim_flag == 1"),
        ("n_tiles", n_tiles, ("day", "species", "band"), "unique tiles with >= 1 obs that day"),
        ("n_steps", n_steps, ("day",), "3-hourly ObsFcstAna files found for the day (8 normally; first and last day of a run can have fewer)"),
        ("tiles_both_terra_aqua", both_TA, ("day",), "tiles with Terra and Aqua obs the same day"),
        ("tiles_both_j1_j2", both_J, ("day",), "tiles with J1 and J2 obs the same day"),
        ("tiles_any_modis", union_MODIS, ("day",), "tiles with Terra or Aqua obs"),
        ("tiles_any_viirs", union_VIIRS, ("day",), "tiles with J1 or J2 obs"),
        ("tiles_terra_only", only_T, ("day",), "tiles with Terra but not Aqua"),
        ("tiles_aqua_only", only_A, ("day",), "tiles with Aqua but not Terra"),
        ("tiles_j1_only", only_J1, ("day",), "tiles with J1 but not J2"),
        ("tiles_j2_only", only_J2, ("day",), "tiles with J2 but not J1"),
    ):
        v = ds.createVariable(nm_, arr.dtype, dims, **z)
        v[:] = arr
        v.long_name = desc
    v = ds.createVariable("assim_any", "i1", ("day", "species"), **z)
    v[:] = assim_seen.astype(np.int8)
    v.long_name = "1 if any obs of the species had assim_flag == 1 that day"
    v = ds.createVariable("days_with_obs", "i2", ("month", "species", "tile"), zlib=True, complevel=4,
                          chunksizes=(1, 1, ntile))
    v[:] = days_with_obs
    v.long_name = "days in the month with >= 1 obs for the tile and species (tile index = tilenum - 1)"
    v = ds.createVariable("tile_lat", "f4", ("tile",), **z)
    v[:] = tile_lat
    v.long_name = "latitude of tiles that ever had an obs (NaN otherwise)"
    v = ds.createVariable("tile_lon", "f4", ("tile",), **z)
    v[:] = tile_lon
    ds.arm = a.arm
    ds.exp_id = exp_id
    ds.source_root = str(root)
    ds.config = str(a.config)
    ds.note = ("Counts are of 3-hourly ObsFcstAna entries (obs mapped to tiles), by UTC date of the file. "
               "OL and VIIRS arms read monitor-only MODIS from gmao_ops, which ends 2025-12-31: zero MODIS obs "
               "from 2026-01-01 is expected for them (see arms.json known_issues).")
    ds.created = dt.datetime.now().isoformat(timespec="seconds")
    ds.close()
    print(f"wrote {out}  ({n_files} files in {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
