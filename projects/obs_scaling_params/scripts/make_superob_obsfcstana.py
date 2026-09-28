#!/usr/bin/env python
"""Write super-ob-style ldas_ObsFcstAna files from an existing (raw-mode) experiment.

Purpose: build a z-score scaling climatology that matches GEOSldas H SAF
super-obs (obs_param_nml%superob_grid_deg > 0) without rerunning the model.
The existing climatology is binned on raw-mode obs (one tile-mean obs per tile,
at the tile-mean location); super-obs are cell means placed at cell centres.

For each input ObsFcstAna file and each selected species, the obs records are
aggregated to a regular lat/lon grid (default 0.25 deg, same cell convention
as the GEOSldas super-ob reader).  Each occupied cell becomes one record with:

    lat/lon       = cell centre
    obs, fcst ... = unweighted mean over the tile records in the cell
    tilenum       = tilenum of the first record in the cell (not used by the
                    scaling statistics)

Approximations: the input obs are already tile means (the per-tile raw-obs
counts are not stored), and the model side is the mean of the tile
predictions in the cell, not the GEOSldas uniform-FOV prediction for a super-ob.

The output tree mirrors <exp>/output/<domain>/ana/ens_avg/Y*/M*/ with the same
file names, so scripts/run_scaling_params.py can be pointed at it unchanged.
Other species are dropped; the obsparam_* metadata variables are copied.
"""

from __future__ import annotations

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


AVG_VARS = ("obs", "obsvar", "fcst", "fcstvar", "ana", "anavar")


def aggregate_file(args: tuple[str, str, tuple[int, ...], float]) -> tuple[str, int, int]:
    src, dst, species, res = args
    with nc.Dataset(src) as d:
        sp = np.asarray(d["species"][:])
        m = np.isin(sp, species)
        cols = {v: np.ma.filled(d[v][:].astype(np.float64), np.nan)[m] for v in AVG_VARS}
        lat = np.asarray(d["lat"][:], dtype=np.float64)[m]
        lon = np.asarray(d["lon"][:], dtype=np.float64)[m]
        df = pd.DataFrame(cols)
        df["species"] = sp[m]
        df["tilenum"] = np.asarray(d["tilenum"][:])[m]
        df["assim_flag"] = np.asarray(d["assim_flag"][:])[m]
        n_lon, n_lat = int(round(360.0 / res)), int(round(180.0 / res))
        df["ci"] = np.minimum(np.floor((np.mod(lon + 180.0, 360.0)) / res).astype(np.int64), n_lon - 1)
        df["cj"] = np.minimum(np.floor((lat + 90.0) / res).astype(np.int64), n_lat - 1)
        agg = df.groupby(["species", "ci", "cj"], sort=True).agg(
            **{v: (v, "mean") for v in AVG_VARS},
            tilenum=("tilenum", "first"), assim_flag=("assim_flag", "max")).reset_index()
        agg["lon"] = -180.0 + (agg.ci + 0.5) * res
        agg["lat"] = -90.0 + (agg.cj + 0.5) * res

        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        tmp = dst + ".tmp"
        with nc.Dataset(tmp, "w", format="NETCDF4") as o:
            o.setncatts({a: d.getncattr(a) for a in d.ncattrs()})
            o.superob_note = (f"super-ob-style aggregation of {os.path.basename(src)} to {res:g}-deg "
                              "cells (make_superob_obsfcstana.py)")
            for name, dim in d.dimensions.items():
                o.createDimension(name, len(agg) if name == "n_obs" else
                                  (None if dim.isunlimited() else dim.size))
            for name, v in d.variables.items():
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                ov = o.createVariable(name, v.dtype, v.dimensions, fill_value=fill, zlib=True)
                ov.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != "_FillValue"})
                if v.dimensions == ("n_obs",):
                    vals = agg[name].to_numpy()
                    if np.issubdtype(v.dtype, np.floating) and fill is not None:
                        vals = np.where(np.isnan(vals), fill, vals)
                    ov[:] = vals.astype(v.dtype)
                else:
                    ov[:] = v[:]
    os.replace(tmp, dst)
    return src, int(m.sum()), len(agg)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--src-ana", type=Path, required=True, help=".../output/<domain>/ana/ens_avg")
    p.add_argument("--dst-ana", type=Path, required=True)
    p.add_argument("--species", default="12,13,14", help="ObsFcstAna species indices to aggregate")
    p.add_argument("--grid-res", type=float, default=0.25)
    p.add_argument("--start", default="2019-02", help="YYYY-MM")
    p.add_argument("--end", default="2020-08", help="YYYY-MM")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    a = p.parse_args()

    species = tuple(int(s) for s in a.species.split(","))
    months = pd.period_range(a.start, a.end, freq="M")
    jobs = []
    for mo in months:
        for f in sorted(glob.glob(str(a.src_ana / f"Y{mo.year}" / f"M{mo.month:02d}" / "*ldas_ObsFcstAna.*.nc4"))):
            dst = a.dst_ana / f"Y{mo.year}" / f"M{mo.month:02d}" / os.path.basename(f)
            if not os.path.exists(dst):
                jobs.append((f, str(dst), species, a.grid_res))
    print(f"{len(jobs)} files to aggregate with {a.workers} workers", flush=True)
    n_in = n_out = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, (src, ni, no) in enumerate(ex.map(aggregate_file, jobs, chunksize=4), 1):
            n_in += ni; n_out += no
            if k % 250 == 0 or k == len(jobs):
                print(f"  {k}/{len(jobs)} files, obs {n_in:,} -> super-obs {n_out:,} "
                      f"({n_out / max(n_in, 1):.2f})", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
