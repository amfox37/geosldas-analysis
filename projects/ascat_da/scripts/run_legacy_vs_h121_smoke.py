#!/usr/bin/env python
"""Smoke-test Legacy BUFR vs H121 obs processing.

This is a lightweight companion to
``projects/ascat_da/notebooks/legacy_vs_h121_obs.ipynb``. It avoids plotting
and tests the parts most likely to break:

* Legacy BUFR and H121 file discovery/reading.
* QC filtering.
* GEOSldas tilecoord/tilegrids reading.
* Raw-ob assignment to GEOS M36 tile numbers.
* Windowed super-ob formation and matched Legacy/H121 tile-window statistics.

Run with the regrid environment, for example:

    /Users/amfox/mamba/envs/regrid/bin/python \
        projects/ascat_da/scripts/run_legacy_vs_h121_smoke.py
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

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
from lib.superob import form_super_obs, match_tiles  # noqa: E402


PLATFORMS = {
    "Metop-A": {"bufr_prefix": "M02-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_a"},
    "Metop-B": {"bufr_prefix": "M01-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_b"},
    "Metop-C": {"bufr_prefix": "M03-ASCA-ASCSMO02-NA-5.0-", "subdir": "metop_c"},
}

def parse_domain(text: str) -> tuple[float, float, float, float]:
    vals = [float(x.strip()) for x in text.split(",")]
    if len(vals) != 4:
        raise argparse.ArgumentTypeError("domain must be lat0,lon0,lat1,lon1")
    return tuple(vals)  # type: ignore[return-value]


def month_dir(base: Path, subdir: str, date: datetime) -> Path:
    return base / subdir / f"Y{date:%Y}" / f"M{date:%m}"


def stats_line(legacy: np.ndarray, h121: np.ndarray) -> str:
    if len(legacy) == 0:
        return "matched=0"
    bias = float(np.mean(h121 - legacy))
    rmse = float(np.sqrt(np.mean((h121 - legacy) ** 2)))
    corr = float(np.corrcoef(legacy, h121)[0, 1]) if len(legacy) > 1 else np.nan
    return f"matched={len(legacy)} bias={bias:+.2f}% rmse={rmse:.2f}% r={corr:.3f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="2020-06-02", help="Date to read, YYYY-MM-DD")
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
        "--tile-base",
        type=Path,
        default=None,
        help="Path stem for .ldas_tilecoord.bin and .ldas_tilegrids.bin",
    )
    args = parser.parse_args()

    date = datetime.strptime(args.date, "%Y-%m-%d")
    obs_root = args.obs_root.expanduser()
    tile_base = args.tile_base
    if tile_base is None:
        tile_base = obs_root / "tilecoord" / "hsaf_cdr_test_DAv8_M36_202006_innov"
    tile_base = tile_base.expanduser()

    tilecoord_path = Path(f"{tile_base}.ldas_tilecoord.bin")
    tilegrids_path = Path(f"{tile_base}.ldas_tilegrids.bin")

    print(f"Date: {date:%Y-%m-%d}")
    print(f"Domain: {args.domain}")
    print(f"Obs root: {obs_root}")
    print(f"Tilecoord: {tilecoord_path}")
    print(f"Tilegrids: {tilegrids_path}")

    tile_coord = read_tilecoord(str(tilecoord_path))
    _, tile_grid_domain = read_tilegrids(str(tilegrids_path))
    print(
        "Tile metadata: "
        f"{tile_coord['N_tile']} tiles, "
        f"{tile_grid_domain['gridtype'].strip()} "
        f"{tile_grid_domain['N_lon']}x{tile_grid_domain['N_lat']}"
    )

    bufr_base = obs_root / "legacy_bufr"
    h121_base = obs_root / "H121"
    failures = []

    for platform, cfg in PLATFORMS.items():
        legacy = read_bufr(
            str(month_dir(bufr_base, cfg["subdir"], date)),
            date,
            cfg["bufr_prefix"],
            domain=args.domain,
            qc=QC_DEFAULT_BUFR,
        )
        h121 = read_h121(
            str(month_dir(h121_base, cfg["subdir"], date)),
            date,
            domain=args.domain,
            qc=QC_DEFAULT_H121,
        )

        legacy_so = form_super_obs(
            legacy["lat"],
            legacy["lon"],
            legacy["ssm"],
            window=legacy["window"],
            cycle=legacy["cycle"],
            tile_coord=tile_coord,
            tile_grid=tile_grid_domain,
        )
        h121_so = form_super_obs(
            h121["lat"],
            h121["lon"],
            h121["ssm"],
            window=h121["window"],
            cycle=h121["cycle"],
            tile_coord=tile_coord,
            tile_grid=tile_grid_domain,
        )
        matched_legacy, matched_h121 = match_tiles(legacy_so, h121_so)

        print(
            f"{platform}: "
            f"raw Legacy={len(legacy['ssm'])} H121={len(h121['ssm'])}; "
            f"tile-windows Legacy={len(legacy_so['ssm'])} H121={len(h121_so['ssm'])}; "
            f"{stats_line(matched_legacy, matched_h121)}"
        )

        if len(legacy["ssm"]) == 0 or len(h121["ssm"]) == 0:
            failures.append(f"{platform}: missing raw obs")
        if len(matched_legacy) == 0:
            failures.append(f"{platform}: no matched tile-windows")

    if failures:
        print("FAIL:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print("ASCAT smoke test OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
