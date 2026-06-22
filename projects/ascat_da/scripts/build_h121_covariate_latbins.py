#!/usr/bin/env python
"""Latitude-binned H121 QC covariates (sensitivity, subsurface scattering).

Checks whether the Legacy-vs-H121 high-latitude SSM bias correlates with
degraded H121 retrieval quality at high latitude (low backscatter
sensitivity / high subsurface-scattering probability over boreal
forest/tundra), independent of which QC thresholds are applied. Obs are
screened with a *minimal* QC (ssm range, water-flag, processing-flag only —
no sens_min/subsfc_max/bsflag screening) so the covariates reflect what
the obs look like before any of those filters act on them.

Output is small (one row per date/platform/lat-bin) so per-day raw reads
are aggregated in place rather than cached at the per-obs level.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import netCDF4 as nc


REPO_ROOT = Path(__file__).resolve().parents[3]
ASCAT_ROOT = REPO_ROOT / "projects" / "ascat_da"

if str(ASCAT_ROOT) not in sys.path:
    sys.path.insert(0, str(ASCAT_ROOT))

from lib.qc import QC_DEFAULT_H121, apply_h121_qc  # noqa: E402

QC_MINIMAL_H121 = {**QC_DEFAULT_H121, "sens_min": None, "subsfc_max": None, "bsflag_bad_bits": None}

LAT_BIN_EDGES = np.arange(-60, 91, 10)
PLATFORMS = {"Metop-A": "metop_a", "Metop-B": "metop_b", "Metop-C": "metop_c"}


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def read_day_platform(h121_dir: Path, date: datetime) -> dict:
    out = {k: [] for k in ("lat", "ssm", "sens", "subsfc")}
    pattern = f"*_{date.strftime('%Y%m%d')}*.nc"
    for fpath in sorted(h121_dir.glob(pattern)):
        ds = nc.Dataset(fpath)
        try:
            lat = ds.variables["latitude"][:].filled(np.nan)
            ssm = ds.variables["surface_soil_moisture"][:].filled(np.nan)

            def _load(name, nan_fill=False):
                # lib.readers.read_h121's `_load` uses `.data`, which returns the
                # *unscaled* raw fill value for masked entries (e.g. -2147483648
                # instead of the scaled ~-214.7) -- harmless there because
                # sens_min/subsfc filters reject those rows anyway, but it corrupts
                # any mean computed over them once those filters are disabled.
                v = ds.variables[name][:]
                if np.ma.is_masked(v) or isinstance(v, np.ma.MaskedArray):
                    out = v.filled(v.fill_value).astype(float)
                    if nan_fill:
                        out[np.ma.getmaskarray(v)] = np.nan
                    return out
                return np.asarray(v, float)

            sf = _load("surface_flag")
            pf = _load("processing_flag")
            cf = _load("correction_flag")
            tc = _load("topographic_complexity")
            wf = _load("wetland_fraction")
            subsfc = _load("subsurface_scattering_probability", nan_fill=True)
            sens = _load("surface_soil_moisture_sensitivity", nan_fill=True)
            bsflag = _load("backscatter40_flag")
            # QC screening treats "not computed" subsfc as benign (0% contamination
            # risk); keep that convention for the mask, but track true NaNs
            # separately for the covariate means below.
            subsfc_for_qc = np.where(np.isnan(subsfc), 0.0, subsfc)

            arrays = dict(ssm=ssm, sf=sf, pf=pf, cf=cf, tc=tc, wf=wf, subsfc=subsfc_for_qc, sens=np.where(np.isnan(sens), -np.inf, sens), bsflag=bsflag)
            mask = apply_h121_qc(arrays, QC_MINIMAL_H121)

            out["lat"].extend(lat[mask].tolist())
            out["ssm"].extend(ssm[mask].tolist())
            out["sens"].extend(sens[mask].tolist())
            out["subsfc"].extend(subsfc[mask].tolist())
        finally:
            ds.close()
    return {k: np.array(v) for k, v in out.items()}


def bin_day(date: datetime, platform: str, obs_root: Path) -> pd.DataFrame:
    h121_dir = obs_root / "H121" / PLATFORMS[platform] / f"Y{date:%Y}" / f"M{date:%m}"
    data = read_day_platform(h121_dir, date)
    if len(data["lat"]) == 0:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    df["lat_bin"] = pd.cut(df["lat"], bins=LAT_BIN_EDGES)
    agg = df.groupby("lat_bin", observed=True).agg(
        n=("ssm", "size"),
        ssm_mean=("ssm", "mean"),
        sens_mean=("sens", lambda s: np.nanmean(s) if s.notna().any() else np.nan),
        sens_n=("sens", lambda s: s.notna().sum()),
        subsfc_mean=("subsfc", lambda s: np.nanmean(s) if s.notna().any() else np.nan),
        subsfc_n=("subsfc", lambda s: s.notna().sum()),
    ).reset_index()
    agg.insert(0, "date", date.strftime("%Y-%m-%d"))
    agg.insert(1, "platform", platform)
    return agg


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2020-06-01"))
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2020-06-10"))
    parser.add_argument("--obs-root", type=Path, default=Path("/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample"))
    parser.add_argument("--out", type=Path, default=ASCAT_ROOT / ".cache" / "h121_covariate_latbins.csv")
    args = parser.parse_args()

    obs_root = args.obs_root.expanduser()
    rows = []
    t0 = time.perf_counter()
    for date in iter_dates(args.start_date, args.end_date):
        for platform in PLATFORMS:
            dt0 = time.perf_counter()
            agg = bin_day(date, platform, obs_root)
            rows.append(agg)
            print(f"{date:%Y-%m-%d} {platform}: {len(agg)} lat bins, {time.perf_counter() - dt0:.1f}s", flush=True)

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {len(out):,} rows to {args.out}")
    print(f"Total elapsed: {(time.perf_counter() - t0) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
