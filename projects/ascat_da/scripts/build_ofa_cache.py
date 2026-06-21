#!/usr/bin/env python
"""Build a compact GEOSldas ObsFcstAna ASCAT cache.

The cache is one row per analysis date/cycle/species/tile. It is intended as
the first analysis object in the Legacy-vs-H121 notebook, before reproducing
super-obs from raw BUFR/H121 files.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


PRODUCTS = {
    "legacy": {"Metop-A": 9, "Metop-B": 10, "Metop-C": 11},
    "h121": {"Metop-A": 14, "Metop-B": 15, "Metop-C": 16},
}
SPECIES_META = {
    species: {"product": product, "platform": platform}
    for product, platforms in PRODUCTS.items()
    for platform, species in platforms.items()
}


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def filename_cycle(path: Path) -> int:
    token = [p for p in path.name.replace(".", "_").split("_") if p.endswith("z")][0]
    return int(token[:2]) // 3


def read_ofa_file(path: Path, date: datetime) -> pd.DataFrame:
    cycle = filename_cycle(path)
    species_keep = np.array(sorted(SPECIES_META), dtype=np.int16)
    with nc.Dataset(path) as ds:
        species = ds.variables["species"][:]
        mask = np.isin(species, species_keep)
        if not np.any(mask):
            return pd.DataFrame()

        obs = ds.variables["obs"][:][mask].astype("float32") * 100.0
        fcst = ds.variables["fcst"][:][mask].astype("float32") * 100.0
        ana = ds.variables["ana"][:][mask].astype("float32") * 100.0
        assim = ds.variables["assim_flag"][:][mask].astype("int16")
        sp = species[mask].astype("int16")

        frame = pd.DataFrame(
            {
                "date": date.strftime("%Y-%m-%d"),
                "cycle": np.int16(cycle),
                "window": np.int8(cycle),
                "species": sp,
                "tilenum": ds.variables["tilenum"][:][mask].astype("int64"),
                "lat": ds.variables["lat"][:][mask].astype("float32"),
                "lon": ds.variables["lon"][:][mask].astype("float32"),
                "obs_pct": obs,
                "fcst_pct": fcst,
                "ana_pct": ana,
                "innov_pct": obs - fcst,
                "incr_pct": ana - fcst,
                "assim_flag": assim,
            }
        )
    frame["product"] = frame["species"].map(lambda s: SPECIES_META[int(s)]["product"])
    frame["platform"] = frame["species"].map(lambda s: SPECIES_META[int(s)]["platform"])
    return frame


def collapse_ofa_rows(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return rows
    rows = rows.copy()
    rows["_assim_count"] = (rows["assim_flag"] != 0).astype("int16")
    return (
        rows.groupby(["date", "cycle", "window", "product", "platform", "species", "tilenum"], as_index=False)
        .agg(
            lat=("lat", "mean"),
            lon=("lon", "mean"),
            obs_pct=("obs_pct", "mean"),
            fcst_pct=("fcst_pct", "mean"),
            ana_pct=("ana_pct", "mean"),
            innov_pct=("innov_pct", "mean"),
            incr_pct=("incr_pct", "mean"),
            n_ofa_obs=("obs_pct", "size"),
            assim_count=("_assim_count", "sum"),
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2020-06-01"))
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2020-06-10"))
    parser.add_argument(
        "--ofa-dir",
        type=Path,
        default=Path(
            "/Users/amfox/Desktop/ASCAT_SSM_CDR/"
            "hsaf_cdr_test_DAv8_M36_202006_innov/output/SMAP_EASEv2_M36_GLOBAL/"
            "ana/ens_avg/Y2020/M06"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("projects/ascat_da/.cache/ofa"),
    )
    parser.add_argument("--version", default="v1")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for date in iter_dates(args.start_date, args.end_date):
        files = sorted(args.ofa_dir.glob(f"*{date:%Y%m%d}*.nc4"))
        print(f"{date:%Y-%m-%d}: {len(files)} OFA files", flush=True)
        day_frames = [read_ofa_file(path, date) for path in files]
        day_frames = [frame for frame in day_frames if not frame.empty]
        if day_frames:
            frames.append(pd.concat(day_frames, ignore_index=True))

    if not frames:
        raise RuntimeError("No OFA ASCAT rows found")

    cache = collapse_ofa_rows(pd.concat(frames, ignore_index=True))
    cache["assim_frac"] = cache["assim_count"] / cache["n_ofa_obs"]

    tag = f"{args.start_date:%Y%m%d}_{args.end_date:%Y%m%d}_{args.version}"
    pkl = args.out_dir / f"ofa_ascat_tile_cycle_{tag}.pkl"
    csv = args.out_dir / f"ofa_ascat_tile_cycle_summary_{tag}.csv"
    cache.to_pickle(pkl)

    summary = (
        cache.groupby(["product", "platform", "species"], as_index=False)
        .agg(
            tile_cycles=("tilenum", "size"),
            total_ofa_obs=("n_ofa_obs", "sum"),
            mean_obs_pct=("obs_pct", "mean"),
            mean_fcst_pct=("fcst_pct", "mean"),
            mean_innov_pct=("innov_pct", "mean"),
            mean_incr_pct=("incr_pct", "mean"),
            mean_assim_frac=("assim_frac", "mean"),
        )
    )
    summary.to_csv(csv, index=False)
    print(f"Wrote {len(cache):,} tile/cycle rows: {pkl}")
    print(f"Wrote summary: {csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
