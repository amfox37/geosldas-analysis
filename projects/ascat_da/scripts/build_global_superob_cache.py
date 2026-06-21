#!/usr/bin/env python
"""Build global ASCAT Legacy/H121 GEOS tile/cycle super-ob caches.

This reads the raw observation files from the local Discover-style sample,
applies the same notebook QC and GEOS-centered cycle assignment, maps each
accepted observation to the experiment M36 tilecoord, then writes compact
derived super-ob tables. Raw individual retrievals are not cached.

Each output row is one:

    date + product + platform + tilenum + GEOS cycle

with the mean SSM super-ob value and simple diagnostics for the observations
that went into that tile/cycle.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
ASCAT_ROOT = REPO_ROOT / "projects" / "ascat_da"
COMMON_IO = REPO_ROOT / "common" / "python" / "io"

for path in (ASCAT_ROOT, COMMON_IO):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from read_GEOSldas import read_tilecoord, read_tilegrids  # noqa: E402
from lib.qc import QC_DEFAULT_BUFR, QC_DEFAULT_H121  # noqa: E402
from lib.readers import read_bufr, read_h121  # noqa: E402
from lib.superob import build_tile_lookup, form_super_obs  # noqa: E402


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
    return cache_dir / f"ascat_global_superobs_{date:%Y%m%d}_{version}.pkl"


def read_raw(product: str, cfg: dict[str, str], obs_root: Path, date: datetime) -> dict:
    meta = PRODUCTS[product]
    base = obs_root / str(meta["base"])
    data_dir = month_dir(base, cfg["subdir"], date)
    if meta["reader"] == "bufr":
        return read_bufr(
            str(data_dir),
            date,
            cfg["bufr_prefix"],
            domain=None,
            qc=QC_DEFAULT_BUFR,
        )
    if meta["reader"] == "h121":
        return read_h121(str(data_dir), date, domain=None, qc=QC_DEFAULT_H121)
    raise ValueError(f"Unknown reader for product {product!r}")


def superobs_to_frame(date: datetime, product: str, platform: str, raw: dict, superobs: dict) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "date": date.strftime("%Y-%m-%d"),
            "product": product,
            "platform": platform,
            "tilenum": superobs["tilenum"].astype("int64"),
            "cycle": superobs["cycle"].astype("int16"),
            "window": superobs["window"].astype("int8"),
            "i_indg": superobs["ij"][:, 0].astype("int16"),
            "j_indg": superobs["ij"][:, 1].astype("int16"),
            "lat": superobs["lat"].astype("float32"),
            "lon": superobs["lon"].astype("float32"),
            "ssm_pct": superobs["ssm"].astype("float32"),
            "n_obs": superobs["count"].astype("int16"),
            "ssm_std_pct": superobs["ssm_std"].astype("float32"),
            "ssm_min_pct": superobs["ssm_min"].astype("float32"),
            "ssm_max_pct": superobs["ssm_max"].astype("float32"),
            "raw_obs_read": len(raw["ssm"]),
        }
    )


def build_day(date: datetime, obs_root: Path, tile_coord: dict, tile_grid: dict, lookup: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    summary_rows = []
    for product, meta in PRODUCTS.items():
        for platform, cfg in meta["platforms"].items():
            t0 = time.perf_counter()
            raw = read_raw(product, cfg, obs_root, date)
            superobs = form_super_obs(
                raw["lat"],
                raw["lon"],
                raw["ssm"],
                window=raw["window"],
                cycle=raw["cycle"],
                tile_coord=tile_coord,
                tile_grid=tile_grid,
                lookup=lookup,
            )
            elapsed = time.perf_counter() - t0
            frame = superobs_to_frame(date, product, platform, raw, superobs)
            frames.append(frame)
            summary_rows.append(
                {
                    "date": date.strftime("%Y-%m-%d"),
                    "product": product,
                    "platform": platform,
                    "raw_obs_read": len(raw["ssm"]),
                    "superobs": len(frame),
                    "total_superob_inputs": int(frame["n_obs"].sum()) if len(frame) else 0,
                    "elapsed_s": round(elapsed, 2),
                }
            )
            print(
                f"  {date:%Y-%m-%d} {product:6s} {platform}: "
                f"raw={len(raw['ssm']):,} superobs={len(frame):,} "
                f"elapsed={elapsed:.1f}s",
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
        "--obs-root",
        type=Path,
        default=Path("/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample"),
        help="Root containing legacy_bufr/, H121/, and tilecoord/.",
    )
    parser.add_argument(
        "--tile-base",
        type=Path,
        default=None,
        help="Path stem for .ldas_tilecoord.bin and .ldas_tilegrids.bin.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=ASCAT_ROOT / ".cache" / "global_superobs",
        help="Directory for per-day super-ob pickle files and summary CSV.",
    )
    parser.add_argument("--version", default="geos_cycle_global_v1")
    parser.add_argument("--force", action="store_true", help="Rebuild existing daily cache files.")
    args = parser.parse_args()

    obs_root = args.obs_root.expanduser()
    cache_dir = args.cache_dir.expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)

    tile_base = args.tile_base
    if tile_base is None:
        tile_base = obs_root / "tilecoord" / "hsaf_cdr_test_DAv8_M36_202006_innov"
    tile_base = tile_base.expanduser()

    print(f"Obs root: {obs_root}")
    print(f"Tile base: {tile_base}")
    print(f"Cache dir: {cache_dir}")
    print(f"Date range: {args.start_date:%Y-%m-%d} to {args.end_date:%Y-%m-%d}")
    print(f"Version: {args.version}")

    tile_coord = read_tilecoord(str(tile_base) + ".ldas_tilecoord.bin")
    _, tile_grid = read_tilegrids(str(tile_base) + ".ldas_tilegrids.bin")
    lookup = build_tile_lookup(tile_coord)

    all_summary = []
    overall_t0 = time.perf_counter()
    for date in iter_dates(args.start_date, args.end_date):
        out = cache_path(cache_dir, date, args.version)
        if out.exists() and not args.force:
            print(f"{date:%Y-%m-%d}: cache exists, skipping {out.name}", flush=True)
            continue
        print(f"{date:%Y-%m-%d}: building global super-ob cache", flush=True)
        t0 = time.perf_counter()
        day_frame, summary = build_day(date, obs_root, tile_coord, tile_grid, lookup)
        day_frame.to_pickle(out)
        summary["cache_file"] = out.name
        summary["day_elapsed_s"] = round(time.perf_counter() - t0, 2)
        all_summary.append(summary)
        print(
            f"{date:%Y-%m-%d}: wrote {len(day_frame):,} rows to {out} "
            f"in {time.perf_counter() - t0:.1f}s",
            flush=True,
        )

    if all_summary:
        summary_all = pd.concat(all_summary, ignore_index=True)
        summary_path = cache_dir / f"ascat_global_superobs_summary_{args.version}.csv"
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
