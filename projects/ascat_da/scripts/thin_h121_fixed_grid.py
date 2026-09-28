#!/usr/bin/env python
"""Write spatially thinned copies of H SAF ASCAT H121 files for GEOSldas.

Thinning is done before GEOSldas reads the observations, so the Fortran reader,
its QC, and the obs_param namelist (apart from %path) are unchanged.

H121 observations sit on a fixed 12.5 km Fibonacci grid identified by
location_id. A static keep-list of grid points is built once:

    one grid point per regular lat/lon cell (default 0.25 deg), the point
    closest to the cell centre (cos(lat)-scaled distance; ties go to the
    smaller location_id)

and applied to every file, so the same locations are kept at every overpass
and for every platform. Cells follow the GEOSldas H SAF super-ob convention
(i = floor((lon+180)/res), j = floor((lat+90)/res)), so thinned points line up
with the cells used by obs_param_nml%superob_grid_deg.

Selection happens before QC: if the kept point fails the reader's QC at an
overpass, that cell has no observation even if a neighbouring point passed.

The file lists (flistpath/Y*/M*/D*/<flistname>) contain bare file names, so
thinned copies keep the original names in <out-root>/<platform>/Y*/M*/ and
only obs_param_nml%path needs to change.

Outputs in <out-root>:
    <platform>/Y*/M*/<original file name>   thinned copies
    keep_list_grid<res>.csv                 kept location_id, lat, lon, cell
    thinning_summary_grid<res>.csv          per-file obs in/out
"""

from __future__ import annotations

import argparse
import os
import time
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


PLATFORMS = {
    "metop_a": "H121_METOPA.txt",
    "metop_b": "H121_H139_METOPB.txt",
    "metop_c": "H121_H139_METOPC.txt",
}

OBS_ROOT = Path("/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/H121")
FLIST_ROOT = Path("/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/flists")


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def month_subdir(date: datetime) -> Path:
    return Path(f"Y{date:%Y}") / f"M{date:%m}"


def list_files(flist_root: Path, obs_root: Path, platform: str,
               start: datetime, end: datetime) -> list[tuple[datetime, Path]]:
    """Files named in the GEOSldas daily file lists, as (list date, source path)."""
    files = []
    for date in iter_dates(start, end):
        flist = flist_root / f"Y{date:%Y}" / f"M{date:%m}" / f"D{date:%d}" / PLATFORMS[platform]
        if not flist.exists():
            print(f"  {platform} {date:%Y-%m-%d}: no file list {flist}", flush=True)
            continue
        for name in flist.read_text().split():
            src = obs_root / platform / month_subdir(date) / name
            if not src.exists():
                raise FileNotFoundError(f"listed in {flist} but missing: {src}")
            files.append((date, src))
    return files


def read_locations(path: Path) -> pd.DataFrame:
    with nc.Dataset(path) as ds:
        lid = ds["location_id"][:]
        lat = ds["latitude"][:]
        lon = ds["longitude"][:]
    ok = ~(np.ma.getmaskarray(lid) | np.ma.getmaskarray(lat) | np.ma.getmaskarray(lon))
    return pd.DataFrame({
        "location_id": np.asarray(lid[ok], dtype=np.int64),
        "lat": np.asarray(lat[ok], dtype=np.float64),
        "lon": np.asarray(lon[ok], dtype=np.float64),
    })


def build_location_table(sources: list[Path]) -> pd.DataFrame:
    """Unique location_id -> lat/lon over all files; checks the grid is static."""
    frames = []
    for kk, src in enumerate(sources, 1):
        frames.append(read_locations(src).drop_duplicates("location_id"))
        if kk % 50 == 0:
            print(f"  read locations from {kk}/{len(sources)} files", flush=True)
    loc = pd.concat(frames, ignore_index=True)
    spread = loc.groupby("location_id")[["lat", "lon"]].agg(lambda x: x.max() - x.min())
    max_spread = float(spread.to_numpy().max()) if len(spread) else 0.0
    if max_spread > 1.0e-5:
        raise ValueError(f"location_id maps to different lat/lon across files "
                         f"(max spread {max_spread:.2e} deg)")
    return loc.drop_duplicates("location_id").sort_values("location_id").reset_index(drop=True)


def select_keep_list(loc: pd.DataFrame, res: float) -> pd.DataFrame:
    """One location per res x res lat/lon cell: the one closest to the cell centre."""
    n_lon = int(round(360.0 / res))
    n_lat = int(round(180.0 / res))
    lon = np.mod(loc["lon"].to_numpy() + 180.0, 360.0) - 180.0
    lat = loc["lat"].to_numpy()
    i = np.minimum(np.floor((lon + 180.0) / res).astype(np.int64), n_lon - 1)
    j = np.minimum(np.floor((lat + 90.0) / res).astype(np.int64), n_lat - 1)
    centre_lon = -180.0 + (i + 0.5) * res
    centre_lat = -90.0 + (j + 0.5) * res
    dlon = (lon - centre_lon) * np.cos(np.deg2rad(lat))
    dlat = lat - centre_lat
    cand = loc.assign(cell_i=i, cell_j=j,
                      dist_deg=np.sqrt(dlon * dlon + dlat * dlat))
    keep = (cand.sort_values(["cell_j", "cell_i", "dist_deg", "location_id"])
                .drop_duplicates(["cell_j", "cell_i"], keep="first"))
    return keep.sort_values("location_id").reset_index(drop=True)


def write_subset(src: Path, dst: Path, keep_ids: np.ndarray, note: str) -> tuple[int, int]:
    """Copy src to dst keeping only obs whose location_id is in keep_ids.

    Raw (packed) values, fill values, attributes, and compression are preserved.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".tmp")
    with nc.Dataset(src) as s:
        s.set_auto_maskandscale(False)
        lid = s["location_id"][:]
        mask = np.isin(lid, keep_ids)
        n_in, n_out = int(lid.size), int(mask.sum())
        with nc.Dataset(tmp, "w", format=s.data_model) as d:
            atts = {a: s.getncattr(a) for a in s.ncattrs()}
            atts["history"] = (f"{atts.get('history', '')}\n"
                               f"{datetime.now():%Y-%m-%d %H:%M:%S} {note}").lstrip()
            d.setncatts(atts)
            for name, dim in s.dimensions.items():
                size = n_out if name == "obs" else (None if dim.isunlimited() else dim.size)
                d.createDimension(name, size)
            for name, v in s.variables.items():
                filt = v.filters() or {}
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                out = d.createVariable(
                    name, v.dtype, v.dimensions, fill_value=fill,
                    zlib=bool(filt.get("zlib", False)),
                    complevel=int(filt.get("complevel", 4) or 4),
                    shuffle=bool(filt.get("shuffle", False)),
                )
                out.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != "_FillValue"})
                out.set_auto_maskandscale(False)
                data = v[:]
                out[:] = data[mask] if v.dimensions == ("obs",) else data
    os.replace(tmp, dst)
    with nc.Dataset(dst) as d:
        d.set_auto_maskandscale(False)
        if not np.array_equal(d["location_id"][:], lid[mask]):
            raise RuntimeError(f"verification failed for {dst}")
    return n_in, n_out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2019-04-30"),
                        help="first file-list date (include the day before the run starts)")
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2019-05-03"))
    parser.add_argument("--grid-res", type=float, default=0.25, help="cell size [deg]")
    parser.add_argument("--platforms", nargs="+", default=list(PLATFORMS), choices=list(PLATFORMS))
    parser.add_argument("--obs-root", type=Path, default=OBS_ROOT)
    parser.add_argument("--flist-root", type=Path, default=FLIST_ROOT)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--keep-list-in", type=Path, default=None,
                        help="reuse an existing keep-list CSV (e.g. to extend the period "
                             "with identical locations) instead of building one")
    parser.add_argument("--select-only", action="store_true",
                        help="build and report the keep-list without writing files")
    parser.add_argument("--force", action="store_true", help="overwrite existing outputs")
    args = parser.parse_args()

    res = args.grid_res
    if abs(round(360.0 / res) * res - 360.0) > 1e-6 or abs(round(180.0 / res) * res - 180.0) > 1e-6:
        parser.error("--grid-res must divide 180 and 360")
    tag = f"grid{res:g}"

    print(f"Obs root:   {args.obs_root}")
    print(f"Flist root: {args.flist_root}")
    print(f"Out root:   {args.out_root}")
    print(f"Dates:      {args.start_date:%Y-%m-%d} to {args.end_date:%Y-%m-%d}")
    print(f"Grid res:   {res} deg")

    t0 = time.perf_counter()
    files = {p: list_files(args.flist_root, args.obs_root, p, args.start_date, args.end_date)
             for p in args.platforms}
    for p, fl in files.items():
        print(f"  {p}: {len(fl)} files", flush=True)

    if args.keep_list_in is not None:
        keep = pd.read_csv(args.keep_list_in)
        print(f"Keep-list:  {len(keep):,} locations from {args.keep_list_in}")
    else:
        loc = build_location_table([src for fl in files.values() for _, src in fl])
        keep = select_keep_list(loc, res)
        n_cells = len(keep)
        print(f"Locations:  {len(loc):,} unique grid points")
        print(f"Keep-list:  {n_cells:,} cells occupied -> {n_cells:,} points kept "
              f"({100.0 * n_cells / len(loc):.1f}% of grid points)")
        dist_km = keep["dist_deg"] * 111.2
        print(f"Distance of kept point to cell centre [km]: median {dist_km.median():.1f}, "
              f"95th pct {dist_km.quantile(0.95):.1f}, max {dist_km.max():.1f}")
        per_cell = loc.shape[0] / n_cells
        print(f"Grid points per occupied cell: {per_cell:.2f} (mean)")
        if not args.select_only:
            args.out_root.mkdir(parents=True, exist_ok=True)
            keep_path = args.out_root / f"keep_list_{tag}.csv"
            if keep_path.exists() and not args.force:
                raise FileExistsError(f"{keep_path} exists (use --force or --keep-list-in)")
            keep.to_csv(keep_path, index=False, float_format="%.6f")
            print(f"Wrote {keep_path}")

    if args.select_only:
        return 0

    keep_ids = np.sort(keep["location_id"].to_numpy(dtype=np.int64))
    note = (f"thinned by thin_h121_fixed_grid.py: one location_id per {res:g}-deg lat/lon "
            f"cell (closest to cell centre), {len(keep_ids)} locations in keep-list")

    rows = []
    for p, fl in files.items():
        for date, src in fl:
            dst = args.out_root / p / month_subdir(date) / src.name
            if dst.exists() and not args.force:
                print(f"  exists, skipping {dst.name}", flush=True)
                continue
            n_in, n_out = write_subset(src, dst, keep_ids, note)
            rows.append({"platform": p, "date": f"{date:%Y-%m-%d}", "file": src.name,
                         "obs_in": n_in, "obs_out": n_out})
        done = [r for r in rows if r["platform"] == p]
        if done:
            n_in = sum(r["obs_in"] for r in done)
            n_out = sum(r["obs_out"] for r in done)
            print(f"  {p}: wrote {len(done)} files, obs {n_in:,} -> {n_out:,} "
                  f"({100.0 * n_out / max(n_in, 1):.1f}% kept)", flush=True)

    if rows:
        summary_path = args.out_root / f"thinning_summary_{tag}.csv"
        summary = pd.DataFrame(rows)
        if summary_path.exists() and not args.force:
            summary = (pd.concat([pd.read_csv(summary_path), summary], ignore_index=True)
                         .drop_duplicates(["platform", "file"], keep="last"))
        summary.sort_values(["platform", "file"]).to_csv(summary_path, index=False)
        print(f"Wrote {summary_path}")

    print(f"Total elapsed: {(time.perf_counter() - t0) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
