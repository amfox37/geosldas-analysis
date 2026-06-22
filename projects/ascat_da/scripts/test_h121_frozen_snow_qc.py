#!/usr/bin/env python
"""Test: does screening H121's own snow/frozen-soil probability fields change
the Legacy-vs-H121 high-latitude bias?

GEOSldas does not use H121's `snow_cover_probability` / `frozen_soil_probability`
fields -- it screens frozen/snow conditions with its own land-model state
instead. Those two H121 fields are themselves derived from ERA5 climatology
(see their `comment` attributes), not real-time retrievals, so this is a test
of an *additional* climatological screen layered on top of the current default
QC, not a replacement for model-based screening.

H121 only -- Legacy has no analogous fields and is QC-unaffected, so the
existing Legacy v6 tile/cycle cache is reused for matching in the notebook.
Output cache version: geos_cycle_global_v7_frozen_snow_test (H121 rows only).
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
ASCAT_ROOT = REPO_ROOT / "projects" / "ascat_da"
COMMON_IO = REPO_ROOT / "common" / "python" / "io"

for path in (ASCAT_ROOT, COMMON_IO):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from read_GEOSldas import read_tilecoord, read_tilegrids  # noqa: E402
from lib.qc import QC_DEFAULT_H121, apply_h121_qc  # noqa: E402
from lib.readers import _hour_to_geos_cycle  # noqa: E402
from lib.superob import build_tile_lookup, form_super_obs  # noqa: E402


PLATFORMS = {"Metop-A": "metop_a", "Metop-B": "metop_b", "Metop-C": "metop_c"}
SNOW_MAX = 50.0    # reject if snow_cover_probability >= this
FROZEN_MAX = 50.0  # reject if frozen_soil_probability >= this


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def _load(ds, name, nan_fill=False):
    v = ds.variables[name][:]
    if np.ma.is_masked(v) or isinstance(v, np.ma.MaskedArray):
        out = v.filled(v.fill_value).astype(float)
        if nan_fill:
            out[np.ma.getmaskarray(v)] = np.nan
        return out
    return np.asarray(v, float)


def read_h121_with_frozen_snow(h121_dir: Path, date: datetime) -> dict:
    out = {k: [] for k in ("lat", "lon", "ssm", "cycle", "window")}
    pattern = os.path.join(str(h121_dir), f"*_{date:%Y%m%d}*.nc")
    for fpath in sorted(glob.glob(pattern)):
        ds = nc.Dataset(fpath)
        try:
            lat = ds.variables["latitude"][:].filled(np.nan)
            lon = ds.variables["longitude"][:].filled(np.nan)
            ssm = ds.variables["surface_soil_moisture"][:].filled(np.nan)

            sf = _load(ds, "surface_flag")
            pf = _load(ds, "processing_flag")
            cf = _load(ds, "correction_flag")
            tc = _load(ds, "topographic_complexity")
            wf = _load(ds, "wetland_fraction")
            subsfc = _load(ds, "subsurface_scattering_probability", nan_fill=True)
            sens = _load(ds, "surface_soil_moisture_sensitivity", nan_fill=True)
            bsflag = _load(ds, "backscatter40_flag")
            snow_prob = _load(ds, "snow_cover_probability", nan_fill=True)
            frozen_prob = _load(ds, "frozen_soil_probability", nan_fill=True)

            subsfc_for_qc = np.where(np.isnan(subsfc), 0.0, subsfc)
            sens_for_qc = np.where(np.isnan(sens), -np.inf, sens)
            arrays = dict(ssm=ssm, sf=sf, pf=pf, cf=cf, tc=tc, wf=wf,
                          subsfc=subsfc_for_qc, sens=sens_for_qc, bsflag=bsflag)
            mask = apply_h121_qc(arrays, QC_DEFAULT_H121)

            # additional test screen: reject high snow/frozen-soil probability
            # (missing -> NaN -> comparison False -> not rejected by this screen)
            mask &= ~(snow_prob >= SNOW_MAX)
            mask &= ~(frozen_prob >= FROZEN_MAX)

            try:
                t_var = ds.variables["time"]
                times = nc.num2date(t_var[:], getattr(t_var, "units", ""))
                hour_frac = np.array([t.hour + t.minute / 60. for t in times], dtype=float)
            except Exception:
                fname = os.path.basename(fpath)
                tok = [p for p in fname.split("_") if "T" in p and len(p) >= 13][0]
                hour_frac = np.full(len(lat), int(tok[9:11]) + int(tok[11:13]) / 60., dtype=float)

            out["lat"].extend(lat[mask].tolist())
            out["lon"].extend(lon[mask].tolist())
            out["ssm"].extend(ssm[mask].tolist())
            cyc, win = _hour_to_geos_cycle(hour_frac[mask], interval="right_closed")
            out["cycle"].extend(cyc.tolist())
            out["window"].extend(win.tolist())
        finally:
            ds.close()
    return {k: np.array(v) for k, v in out.items()}


def cache_path(cache_dir: Path, date: datetime, version: str) -> Path:
    return cache_dir / f"ascat_global_superobs_{date:%Y%m%d}_{version}.pkl"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2020-06-01"))
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2020-06-10"))
    parser.add_argument("--obs-root", type=Path, default=Path("/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample"))
    parser.add_argument("--tile-base", type=Path, default=None)
    parser.add_argument("--cache-dir", type=Path, default=ASCAT_ROOT / ".cache" / "global_superobs")
    parser.add_argument("--version", default="geos_cycle_global_v7_frozen_snow_test")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    obs_root = args.obs_root.expanduser()
    cache_dir = args.cache_dir.expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)

    tile_base = args.tile_base or (obs_root / "tilecoord" / "hsaf_cdr_test_DAv8_M36_202006_innov")
    tile_base = tile_base.expanduser()

    print(f"Snow/frozen-soil probability test thresholds: snow>={SNOW_MAX}%, frozen>={FROZEN_MAX}% rejected")
    tile_coord = read_tilecoord(str(tile_base) + ".ldas_tilecoord.bin")
    _, tile_grid = read_tilegrids(str(tile_base) + ".ldas_tilegrids.bin")
    lookup = build_tile_lookup(tile_coord)

    overall_t0 = time.perf_counter()
    for date in iter_dates(args.start_date, args.end_date):
        out = cache_path(cache_dir, date, args.version)
        if out.exists() and not args.force:
            print(f"{date:%Y-%m-%d}: cache exists, skipping {out.name}", flush=True)
            continue
        frames = []
        for platform, subdir in PLATFORMS.items():
            t0 = time.perf_counter()
            h121_dir = obs_root / "H121" / subdir / f"Y{date:%Y}" / f"M{date:%m}"
            raw = read_h121_with_frozen_snow(h121_dir, date)
            superobs = form_super_obs(
                raw["lat"], raw["lon"], raw["ssm"],
                window=raw["window"], cycle=raw["cycle"],
                tile_coord=tile_coord, tile_grid=tile_grid, lookup=lookup,
            )
            frame = pd.DataFrame({
                "date": date.strftime("%Y-%m-%d"),
                "product": "h121",
                "platform": platform,
                "tilenum": superobs["tilenum"].astype("int64"),
                "cycle": superobs["cycle"].astype("int16"),
                "window": superobs["window"].astype("int8"),
                "lat": superobs["lat"].astype("float32"),
                "lon": superobs["lon"].astype("float32"),
                "ssm_pct": superobs["ssm"].astype("float32"),
                "n_obs": superobs["count"].astype("int16"),
                "ssm_std_pct": superobs["ssm_std"].astype("float32"),
                "ssm_min_pct": superobs["ssm_min"].astype("float32"),
                "ssm_max_pct": superobs["ssm_max"].astype("float32"),
                "raw_obs_read": len(raw["ssm"]),
            })
            frames.append(frame)
            print(f"  {date:%Y-%m-%d} {platform}: raw={len(raw['ssm']):,} "
                  f"superobs={len(frame):,} elapsed={time.perf_counter() - t0:.1f}s", flush=True)

        day_frame = pd.concat(frames, ignore_index=True)
        day_frame.to_pickle(out)
        print(f"{date:%Y-%m-%d}: wrote {len(day_frame):,} rows to {out}", flush=True)

    print(f"Total elapsed: {(time.perf_counter() - overall_t0) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
