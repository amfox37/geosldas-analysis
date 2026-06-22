#!/usr/bin/env python
"""Build global ASCAT Legacy/H121 0.25-degree grid super-ob caches.

Same raw-read/QC pipeline as build_global_superob_cache.py, but observations
are binned onto a regular 0.25-deg lat/lon grid (cell-mean, like
compare_legacy_bufr_vs_H121.ipynb) instead of looked up against a GEOS
tilecoord. Used to check whether spatial patterns found in tile/cycle
super-obs (e.g. the Legacy/H121 high-latitude bias) are an artifact of the
GEOS tile super-obbing step or a property of the raw observations.

Each output row is one: date + product + platform + grid cell + GEOS cycle.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
ASCAT_ROOT = REPO_ROOT / "projects" / "ascat_da"

if str(ASCAT_ROOT) not in sys.path:
    sys.path.insert(0, str(ASCAT_ROOT))

from lib.qc import QC_DEFAULT_BUFR, QC_DEFAULT_H121  # noqa: E402
from lib.readers import read_bufr, read_h121  # noqa: E402


GRID_RES = 0.25

PRODUCTS = {
    "legacy": {
        "base": "legacy_bufr",
        "reader": "bufr",
        "platforms": {
            "Metop-A": {"subdir": "metop_a", "bufr_prefix": "M02-ASCA-ASCSMO02-NA-5.0-"},
            "Metop-B": {"subdir": "metop_b", "bufr_prefix": "M01-ASCA-ASCSMO02-NA-5.0-"},
            "Metop-C": {"subdir": "metop_c", "bufr_prefix": "M03-ASCA-ASCSMO02-NA-5.0-"},
        },
    },
    "h121": {
        "base": "H121",
        "reader": "h121",
        "platforms": {
            "Metop-A": {"subdir": "metop_a"},
            "Metop-B": {"subdir": "metop_b"},
            "Metop-C": {"subdir": "metop_c"},
        },
    },
}


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def month_dir(base: Path, subdir: str, date: datetime) -> Path:
    return base / subdir / f"Y{date:%Y}" / f"M{date:%m}"


def cache_path(cache_dir: Path, date: datetime, version: str) -> Path:
    return cache_dir / f"ascat_grid025_superobs_{date:%Y%m%d}_{version}.pkl"


def read_raw(product: str, cfg: dict[str, str], obs_root: Path, date: datetime) -> dict:
    meta = PRODUCTS[product]
    base = obs_root / str(meta["base"])
    data_dir = month_dir(base, cfg["subdir"], date)
    if meta["reader"] == "bufr":
        return read_bufr(str(data_dir), date, cfg["bufr_prefix"], domain=None, qc=QC_DEFAULT_BUFR)
    if meta["reader"] == "h121":
        return read_h121(str(data_dir), date, domain=None, qc=QC_DEFAULT_H121)
    raise ValueError(f"Unknown reader for product {product!r}")


def grid_super_obs(raw: dict, res: float = GRID_RES) -> pd.DataFrame:
    if len(raw["ssm"]) == 0:
        return pd.DataFrame(columns=["grid_lat", "grid_lon", "cycle", "ssm", "n_obs"])
    grid_lat = (np.floor(raw["lat"] / res) * res + res / 2.0).astype("float32")
    grid_lon = (np.floor(raw["lon"] / res) * res + res / 2.0).astype("float32")
    df = pd.DataFrame({
        "grid_lat": grid_lat,
        "grid_lon": grid_lon,
        "cycle": raw["cycle"].astype("int16"),
        "ssm": raw["ssm"].astype("float32"),
    })
    return (
        df.groupby(["grid_lat", "grid_lon", "cycle"], as_index=False)
        .agg(ssm=("ssm", "mean"), n_obs=("ssm", "size"))
    )


def build_day(date: datetime, obs_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    summary_rows = []
    for product, meta in PRODUCTS.items():
        for platform, cfg in meta["platforms"].items():
            t0 = time.perf_counter()
            raw = read_raw(product, cfg, obs_root, date)
            gridded = grid_super_obs(raw)
            elapsed = time.perf_counter() - t0
            gridded = gridded.copy()
            gridded.insert(0, "date", date.strftime("%Y-%m-%d"))
            gridded.insert(1, "product", product)
            gridded.insert(2, "platform", platform)
            gridded["raw_obs_read"] = len(raw["ssm"])
            frames.append(gridded)
            summary_rows.append({
                "date": date.strftime("%Y-%m-%d"),
                "product": product,
                "platform": platform,
                "raw_obs_read": len(raw["ssm"]),
                "grid_cells": len(gridded),
                "elapsed_s": round(elapsed, 2),
            })
            print(
                f"  {date:%Y-%m-%d} {product:6s} {platform}: "
                f"raw={len(raw['ssm']):,} grid_cells={len(gridded):,} elapsed={elapsed:.1f}s",
                flush=True,
            )

    day_frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    return day_frame, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2020-06-01"))
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2020-06-10"))
    parser.add_argument(
        "--obs-root", type=Path,
        default=Path("/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample"),
    )
    parser.add_argument(
        "--cache-dir", type=Path,
        default=ASCAT_ROOT / ".cache" / "grid025_superobs",
    )
    parser.add_argument("--version", default="grid025_v1_default_qc")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    obs_root = args.obs_root.expanduser()
    cache_dir = args.cache_dir.expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"Obs root: {obs_root}")
    print(f"Cache dir: {cache_dir}")
    print(f"Date range: {args.start_date:%Y-%m-%d} to {args.end_date:%Y-%m-%d}")
    print(f"Version: {args.version}")

    all_summary = []
    overall_t0 = time.perf_counter()
    for date in iter_dates(args.start_date, args.end_date):
        out = cache_path(cache_dir, date, args.version)
        if out.exists() and not args.force:
            print(f"{date:%Y-%m-%d}: cache exists, skipping {out.name}", flush=True)
            continue
        print(f"{date:%Y-%m-%d}: building 0.25-deg grid super-ob cache", flush=True)
        t0 = time.perf_counter()
        day_frame, summary = build_day(date, obs_root)
        day_frame.to_pickle(out)
        summary["cache_file"] = out.name
        summary["day_elapsed_s"] = round(time.perf_counter() - t0, 2)
        all_summary.append(summary)
        print(f"{date:%Y-%m-%d}: wrote {len(day_frame):,} rows to {out} "
              f"in {time.perf_counter() - t0:.1f}s", flush=True)

    if all_summary:
        summary_all = pd.concat(all_summary, ignore_index=True)
        summary_path = cache_dir / f"ascat_grid025_superobs_summary_{args.version}.csv"
        if summary_path.exists() and not args.force:
            existing = pd.read_csv(summary_path)
            summary_all = pd.concat([existing, summary_all], ignore_index=True)
            summary_all = summary_all.drop_duplicates(["date", "product", "platform"], keep="last")
            summary_all = summary_all.sort_values(["date", "product", "platform"])
        summary_all.to_csv(summary_path, index=False)
        print(f"Wrote summary: {summary_path}")

    print(f"Total elapsed: {(time.perf_counter() - overall_t0) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
