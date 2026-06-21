#!/usr/bin/env python
"""Compare raw ASCAT superobs with GEOSldas ObsFcstAna obs.

The notebook compares Legacy BUFR and H121 raw observations after assigning
them to GEOS M36 tiles. This diagnostic asks a narrower question:

    Do raw Legacy BUFR and H121 superobs reproduce the ASCAT species stored
    in ObsFcstAna files for the same GEOS tilenum/cycle?

Raw observations are assigned to centered GEOSldas 3-hour analysis cycles.
OFA ASCAT obs are stored as fraction saturation, so they are multiplied by
100 before comparison with BUFR percent saturation.
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path

import netCDF4 as nc
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
ASCAT_ROOT = REPO_ROOT / "projects" / "ascat_da"
COMMON_IO = REPO_ROOT / "common" / "python" / "io"

for path in (ASCAT_ROOT, COMMON_IO):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from read_GEOSldas import read_tilecoord, read_tilegrids  # noqa: E402
from lib.qc import QC_DEFAULT_BUFR, QC_DEFAULT_H121  # noqa: E402
from lib.readers import read_bufr, read_h121  # noqa: E402
from lib.superob import form_super_obs  # noqa: E402


PRODUCTS = {
    "legacy": {
        "base": "legacy_bufr",
        "platforms": {
            "Metop-A": {"bufr_prefix": "M02-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_a", "species": 9},
            "Metop-B": {"bufr_prefix": "M01-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_b", "species": 10},
            "Metop-C": {"bufr_prefix": "M03-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_c", "species": 11},
        },
    },
    "h121": {
        "base": "H121",
        "platforms": {
            "Metop-A": {"subdir": "metop_a", "species": 14},
            "Metop-B": {"subdir": "metop_b", "species": 15},
            "Metop-C": {"subdir": "metop_c", "species": 16},
        },
    },
}


def parse_domain(text: str) -> tuple[float, float, float, float]:
    vals = [float(x.strip()) for x in text.split(",")]
    if len(vals) != 4:
        raise argparse.ArgumentTypeError("domain must be lat0,lon0,lat1,lon1")
    return tuple(vals)  # type: ignore[return-value]


def month_dir(base: Path, subdir: str, date: datetime) -> Path:
    return base / subdir / f"Y{date:%Y}" / f"M{date:%m}"


def inbox(lat: np.ndarray, lon: np.ndarray, domain: tuple[float, float, float, float]) -> np.ndarray:
    lat0, lon0, lat1, lon1 = domain
    return (lat >= lat0) & (lat <= lat1) & (lon >= lon0) & (lon <= lon1)


def filename_window(path: Path, shift: int) -> int:
    """Return GEOSldas cycle from OFA filename timestamp plus integer shift."""
    token = [p for p in path.name.replace(".", "_").split("_") if p.endswith("z")][0]
    hour = int(token[:2])
    return (hour // 3 + shift) % 8


def read_ofa_species(
    ofa_dir: Path,
    date: datetime,
    species: int,
    domain: tuple[float, float, float, float],
    window_shift: int,
) -> dict[tuple[int, int], list[dict[str, float | int | str]]]:
    """Read OFA obs grouped by (tilenum, shifted GEOSldas cycle)."""
    out: dict[tuple[int, int], list[dict[str, float | int | str]]] = {}
    for fpath in sorted(ofa_dir.glob(f"*{date:%Y%m%d}*.nc4")):
        win = filename_window(fpath, window_shift)
        with nc.Dataset(fpath) as ds:
            sp = ds.variables["species"][:]
            lat = ds.variables["lat"][:]
            lon = ds.variables["lon"][:]
            tile = ds.variables["tilenum"][:]
            obs = ds.variables["obs"][:]
            assim = ds.variables["assim_flag"][:]
            mask = (sp == species) & inbox(lat, lon, domain)
            for t, la, lo, ob, af in zip(tile[mask], lat[mask], lon[mask], obs[mask], assim[mask]):
                key = (int(t), int(win))
                out.setdefault(key, []).append(
                    {
                        "tilenum": int(t),
                        "window": int(win),
                        "obs_pct": float(ob) * 100.0,
                        "lat": float(la),
                        "lon": float(lo),
                        "assim_flag": int(af),
                        "file": fpath.name,
                    }
                )
    return out


def collapse_ofa(groups: dict[tuple[int, int], list[dict[str, float | int | str]]]) -> dict[tuple[int, int], dict[str, float | int | str]]:
    """Average duplicate OFA entries for the same tile/window, if any."""
    collapsed = {}
    for key, rows in groups.items():
        obs = np.array([float(r["obs_pct"]) for r in rows])
        collapsed[key] = {
            "tilenum": key[0],
            "window": key[1],
            "obs_pct": float(obs.mean()),
            "count": len(rows),
            "files": ";".join(sorted({str(r["file"]) for r in rows})),
        }
    return collapsed


def summarize(
    raw: dict[tuple[int, int], dict[str, float | int]],
    ofa: dict[tuple[int, int], dict[str, float | int | str]],
) -> dict[str, float | int]:
    common = sorted(set(raw) & set(ofa))
    raw_only = set(raw) - set(ofa)
    ofa_only = set(ofa) - set(raw)
    if not common:
        return {
            "common": 0,
            "raw_total": len(raw),
            "ofa_total": len(ofa),
            "raw_only": len(raw_only),
            "ofa_only": len(ofa_only),
            "raw_min": np.nan,
            "raw_max": np.nan,
            "raw_mean": np.nan,
            "ofa_min": np.nan,
            "ofa_max": np.nan,
            "ofa_mean": np.nan,
            "bias": np.nan,
            "rmse": np.nan,
            "mae": np.nan,
            "median_abs": np.nan,
            "maxabs": np.nan,
        }
    raw_vals = np.array([float(raw[k]["ssm_pct"]) for k in common])
    ofa_vals = np.array([float(ofa[k]["obs_pct"]) for k in common])
    diff = ofa_vals - raw_vals
    return {
        "common": len(common),
        "raw_total": len(raw),
        "ofa_total": len(ofa),
        "raw_only": len(raw_only),
        "ofa_only": len(ofa_only),
        "raw_min": float(raw_vals.min()),
        "raw_max": float(raw_vals.max()),
        "raw_mean": float(raw_vals.mean()),
        "ofa_min": float(ofa_vals.min()),
        "ofa_max": float(ofa_vals.max()),
        "ofa_mean": float(ofa_vals.mean()),
        "bias": float(diff.mean()),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "mae": float(np.mean(np.abs(diff))),
        "median_abs": float(np.median(np.abs(diff))),
        "maxabs": float(np.max(np.abs(diff))),
    }


def format_summary(summary: dict[str, float | int]) -> str:
    if int(summary["common"]) == 0:
        return (
            f"matched=0/{summary['ofa_total']} OFA; "
            f"raw-only={summary['raw_only']} ofa-only={summary['ofa_only']}"
        )
    return (
        f"matched={summary['common']}/{summary['ofa_total']} OFA "
        f"raw-only={summary['raw_only']} ofa-only={summary['ofa_only']} "
        f"range raw={float(summary['raw_min']):.1f}-{float(summary['raw_max']):.1f}% "
        f"ofa={float(summary['ofa_min']):.1f}-{float(summary['ofa_max']):.1f}% "
        f"bias={float(summary['bias']):+.3f}% "
        f"rmse={float(summary['rmse']):.3f}% "
        f"mae={float(summary['mae']):.3f}% "
        f"medabs={float(summary['median_abs']):.3f}% "
        f"maxabs={float(summary['maxabs']):.3f}%"
    )


def read_raw_product(
    product: str,
    cfg: dict[str, int | str],
    obs_root: Path,
    date: datetime,
    domain: tuple[float, float, float, float],
) -> dict[str, np.ndarray]:
    """Read one platform/product raw obs with notebook-compatible QC."""
    base = obs_root / str(PRODUCTS[product]["base"])
    subdir = str(cfg["subdir"])
    if product == "legacy":
        return read_bufr(
            str(month_dir(base, subdir, date)),
            date,
            str(cfg["bufr_prefix"]),
            domain=domain,
            qc=QC_DEFAULT_BUFR,
        )
    if product == "h121":
        return read_h121(
            str(month_dir(base, subdir, date)),
            date,
            domain=domain,
            qc=QC_DEFAULT_H121,
        )
    raise ValueError(f"Unknown product: {product}")


def superobs_to_keyed_rows(super_obs: dict[str, np.ndarray]) -> dict[tuple[int, int], dict[str, float | int]]:
    """Convert form_super_obs output to rows keyed by (tilenum, cycle)."""
    return {
        (int(t), int(cyc)): {
            "tilenum": int(t),
            "window": int(w),
            "cycle": int(cyc),
            "ssm_pct": float(v),
            "raw_count": int(c),
            "lat": float(la),
            "lon": float(lo),
        }
        for t, w, cyc, v, c, la, lo in zip(
            super_obs["tilenum"],
            super_obs["window"],
            super_obs["cycle"],
            super_obs["ssm"],
            super_obs["count"],
            super_obs["lat"],
            super_obs["lon"],
        )
    }


def compare_one(
    product: str,
    platform: str,
    cfg: dict[str, int | str],
    obs_root: Path,
    ofa_dir: Path,
    date: datetime,
    domain: tuple[float, float, float, float],
    shifts: list[int],
    tile_coord: dict,
    tile_grid_domain: dict,
) -> tuple[list[dict[str, float | int | str]], dict[str, float | int | str]]:
    raw_obs = read_raw_product(product, cfg, obs_root, date, domain)
    raw_so = form_super_obs(
        raw_obs["lat"],
        raw_obs["lon"],
        raw_obs["ssm"],
        window=raw_obs["window"],
        cycle=raw_obs["cycle"],
        tile_coord=tile_coord,
        tile_grid=tile_grid_domain,
    )
    raw = superobs_to_keyed_rows(raw_so)
    species = int(cfg["species"])

    print(f"\n{product} {platform} species {species}: raw tile-cycles={len(raw)}")
    summaries = {}
    ofa_by_shift = {}
    for shift in shifts:
        ofa_groups = read_ofa_species(ofa_dir, date, species, domain, shift)
        ofa = collapse_ofa(ofa_groups)
        ofa_by_shift[shift] = ofa
        summaries[shift] = summarize(raw, ofa)
        print(f"  OFA cycle shift {shift:+d}: {format_summary(summaries[shift])}")

    best_shift = max(
        shifts,
        key=lambda s: (
            int(summaries[s]["common"]),
            -float(summaries[s]["rmse"]) if int(summaries[s]["common"]) else -np.inf,
        ),
    )
    best_ofa = ofa_by_shift[best_shift]
    best_summary = summaries[best_shift]
    common = sorted(set(raw) & set(best_ofa))
    print(f"  best shift: {best_shift:+d}")

    csv_rows: list[dict[str, float | int | str]] = []
    if common:
        rows = []
        for key in common:
            raw_val = float(raw[key]["ssm_pct"])
            ofa_val = float(best_ofa[key]["obs_pct"])
            rows.append((abs(ofa_val - raw_val), key, raw_val, ofa_val))
        rows.sort(reverse=True)
        print("  worst matched residuals:")
        for _absdiff, key, raw_val, ofa_val in rows[:5]:
            print(
                f"    tilenum={key[0]} cycle={key[1]} window={raw[key]['window']} "
                f"raw={raw_val:.3f}% ofa={ofa_val:.3f}% "
                f"diff={ofa_val - raw_val:+.3f}%"
            )
        for _, key, raw_val, ofa_val in rows:
            csv_rows.append(
                {
                    "product": product,
                    "platform": platform,
                    "species": species,
                    "window_shift": best_shift,
                    "tilenum": key[0],
                    "cycle": key[1],
                    "window": raw[key]["window"],
                    "raw_superob_pct": raw_val,
                    "ofa_obs_pct": ofa_val,
                    "diff_pct": ofa_val - raw_val,
                    "raw_count": raw[key]["raw_count"],
                    "ofa_count": best_ofa[key]["count"],
                    "ofa_files": best_ofa[key]["files"],
                }
            )
    summary_row: dict[str, float | int | str] = {
        "date": date.strftime("%Y-%m-%d"),
        "product": product,
        "platform": platform,
        "species": species,
        "best_shift": best_shift,
        **best_summary,
    }
    return csv_rows, summary_row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="2020-06-02", help="Date to compare, YYYY-MM-DD")
    parser.add_argument(
        "--domain",
        type=parse_domain,
        default=(35.0, -99.0, 38.0, -96.0),
        help="Bounding box as lat0,lon0,lat1,lon1",
    )
    parser.add_argument(
        "--obs-root",
        type=Path,
        default=Path("/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample"),
        help="Root containing legacy_bufr/, H121/, and tilecoord/",
    )
    parser.add_argument(
        "--ofa-dir",
        type=Path,
        default=Path(
            "/Users/amfox/Desktop/ASCAT_SSM_CDR/"
            "hsaf_cdr_test_DAv8_M36_202006_innov/output/SMAP_EASEv2_M36_GLOBAL/"
            "ana/ens_avg/Y2020/M06"
        ),
        help="Directory containing ObsFcstAna nc4 files.",
    )
    parser.add_argument(
        "--tile-base",
        type=Path,
        default=None,
        help="Path stem for .ldas_tilecoord.bin and .ldas_tilegrids.bin",
    )
    parser.add_argument(
        "--window-shifts",
        default="0",
        help="Comma-separated diagnostic shifts applied to OFA filename cycles.",
    )
    parser.add_argument(
        "--products",
        default="legacy,h121",
        help="Comma-separated products to compare: legacy,h121",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV path for matched rows using the best shift per platform.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=None,
        help="Optional CSV path for one summary row per product/platform.",
    )
    args = parser.parse_args()

    date = datetime.strptime(args.date, "%Y-%m-%d")
    obs_root = args.obs_root.expanduser()
    ofa_dir = args.ofa_dir.expanduser()
    shifts = [int(x.strip()) for x in args.window_shifts.split(",") if x.strip()]
    products = [x.strip() for x in args.products.split(",") if x.strip()]

    tile_base = args.tile_base
    if tile_base is None:
        tile_base = obs_root / "tilecoord" / "hsaf_cdr_test_DAv8_M36_202006_innov"
    tile_base = tile_base.expanduser()

    print(f"Date: {date:%Y-%m-%d}")
    print(f"Domain: {args.domain}")
    print(f"Obs root: {obs_root}")
    print(f"OFA dir: {ofa_dir}")
    print(f"Products: {products}")
    print(f"Window shifts tested: {shifts}")

    tile_coord = read_tilecoord(str(tile_base) + ".ldas_tilecoord.bin")
    _, tile_grid_domain = read_tilegrids(str(tile_base) + ".ldas_tilegrids.bin")

    csv_rows = []
    summary_rows = []
    for product in products:
        if product not in PRODUCTS:
            raise ValueError(f"Unknown product {product!r}; choose from {sorted(PRODUCTS)}")
        for platform, cfg in PRODUCTS[product]["platforms"].items():
            rows, summary = compare_one(
                product,
                platform,
                cfg,
                obs_root,
                ofa_dir,
                date,
                args.domain,
                shifts,
                tile_coord,
                tile_grid_domain,
            )
            csv_rows.extend(rows)
            summary_rows.append(summary)

    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()) if csv_rows else [])
            if csv_rows:
                writer.writeheader()
                writer.writerows(csv_rows)
        print(f"\nWrote matched diagnostics: {args.csv}")

    if args.summary_csv is not None:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()) if summary_rows else [])
            if summary_rows:
                writer.writeheader()
                writer.writerows(summary_rows)
        print(f"Wrote summary diagnostics: {args.summary_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
