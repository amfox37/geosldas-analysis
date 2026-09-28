#!/usr/bin/env python
"""Fill gaps in a lat/lon z-score scaling-parameter file from the nearest valid cell.

GEOSldas (scale_obs_sfmc_zscore) rejects an obs when the scaling cell it falls
in has no valid statistics.  Climatologies built from ldas_ObsFcstAna files are
binned on the obs locations of the run they came from.  When those locations
do not cover every grid cell (e.g., raw-mode tile-mean H SAF locations on a
0.25 deg grid at high latitudes for CF0360), obs located elsewhere (super-obs
at cell centres, thinned obs at single grid points) are rejected.

For each pentad, every cell that fails the GEOSldas validity test (o_mean,
m_mean > 0; o_std, m_std >= 0; all finite) receives o_mean, o_std, m_mean, and
m_std from the nearest valid cell of the same pentad, if that cell is within
--max-dist-km (great-circle distance between cell centres).  m_min and m_max
(no pentad dimension) are filled the same way.  n_data is left unchanged (0 in
filled cells), and fill_distance_km records the fill distance (0 for original
cells, NaN where nothing was filled).

The input file is not modified.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import netCDF4 as nc
import numpy as np
from scipy.spatial import cKDTree


EARTH_RADIUS_KM = 6371.0
ZSCORE_VARS = ("o_mean", "o_std", "m_mean", "m_std")


def unit_vectors(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def chord_for_km(dist_km: float) -> float:
    return 2.0 * np.sin(dist_km / EARTH_RADIUS_KM / 2.0)


def km_from_chord(chord: np.ndarray) -> np.ndarray:
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))


def nearest_fill(valid: np.ndarray, xyz: np.ndarray, max_chord: float):
    """For invalid cells, index of the nearest valid cell within max_chord (flat indices)."""
    vflat = valid.ravel()
    if not vflat.any():
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([])
    tree = cKDTree(xyz[vflat])
    donors_all = np.flatnonzero(vflat)
    targets = np.flatnonzero(~vflat)
    dist, idx = tree.query(xyz[targets], k=1, distance_upper_bound=max_chord)
    ok = np.isfinite(dist)
    return targets[ok], donors_all[idx[ok]], km_from_chord(dist[ok])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", type=Path)
    parser.add_argument("outfile", type=Path)
    parser.add_argument("--max-dist-km", type=float, default=50.0)
    parser.add_argument("--dry-run", action="store_true",
                        help="report fill statistics without writing the output file")
    args = parser.parse_args()

    with nc.Dataset(args.infile) as src:
        ll_lon, ll_lat = float(src["ll_lon"][:]), float(src["ll_lat"][:])
        dlon, dlat = float(src["d_lon"][:]), float(src["d_lat"][:])
        n_lon, n_lat = src.dimensions["lon"].size, src.dimensions["lat"].size
        n_pentad = src.dimensions["pentad"].size
        data = {v: np.ma.filled(src[v][:].astype(np.float64), np.nan) for v in ZSCORE_VARS}
        m_min = np.ma.filled(src["m_min"][:].astype(np.float64), np.nan)
        m_max = np.ma.filled(src["m_max"][:].astype(np.float64), np.nan)

    lon_c = ll_lon + (np.arange(n_lon) + 0.5) * dlon
    lat_c = ll_lat + (np.arange(n_lat) + 0.5) * dlat
    lon2, lat2 = np.meshgrid(lon_c, lat_c, indexing="ij")          # (lon, lat), like the file
    xyz = unit_vectors(lat2.ravel(), lon2.ravel())
    max_chord = chord_for_km(args.max_dist_km)

    fill_dist = np.full((n_pentad, n_lon, n_lat), np.nan)
    total_filled = 0
    for p in range(n_pentad):
        om, os_, mm, ms = (data[v][p] for v in ZSCORE_VARS)
        with np.errstate(invalid="ignore"):
            valid = (np.isfinite(om) & np.isfinite(os_) & np.isfinite(mm) & np.isfinite(ms)
                     & (om > 0) & (mm > 0) & (os_ >= 0) & (ms >= 0))
        fd = fill_dist[p].ravel()
        fd[valid.ravel()] = 0.0
        targets, donors, dist_km = nearest_fill(valid, xyz, max_chord)
        for v in ZSCORE_VARS:
            arr = data[v][p].ravel()          # view into data[v][p]
            arr[targets] = arr[donors]
        fd[targets] = dist_km
        total_filled += targets.size
        if p % 12 == 0 or p == n_pentad - 1:
            print(f"pentad {p + 1:2d}: valid {int(valid.sum()):7d}, filled {targets.size:7d}"
                  + (f", fill distance median {np.median(dist_km):5.1f} km, max {dist_km.max():5.1f} km"
                     if targets.size else ""), flush=True)

    mm_valid = np.isfinite(m_min) & np.isfinite(m_max)
    t, d, dk = nearest_fill(mm_valid, xyz, max_chord)
    for arr in (m_min, m_max):
        flat = arr.ravel()
        flat[t] = flat[d]
    print(f"m_min/m_max: valid {int(mm_valid.sum())}, filled {t.size}")
    print(f"total pentad-cells filled: {total_filled:,}")

    if args.dry_run:
        return 0

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.infile, args.outfile)
    with nc.Dataset(args.outfile, "a") as out:
        for v in ZSCORE_VARS:
            out[v][:] = data[v]
        out["m_min"][:] = m_min
        out["m_max"][:] = m_max
        fdv = out.createVariable("fill_distance_km", "f4", ("pentad", "lon", "lat"),
                                 zlib=True, fill_value=np.float32(np.nan))
        fdv.long_name = ("great-circle distance to the cell the z-score stats were copied from "
                         "(0 = original stats, NaN = no stats)")
        fdv.units = "km"
        fdv[:] = fill_dist.astype(np.float32)
        out.gapfill_source = str(args.infile.resolve())
        out.gapfill_method = (f"nearest valid cell of the same pentad within {args.max_dist_km:g} km "
                              "(gapfill_scaling_params.py); n_data unchanged")
    print(f"Wrote {args.outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
