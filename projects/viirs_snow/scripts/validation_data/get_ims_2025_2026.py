#!/usr/bin/env python3
"""Fetch IMS 24 km daily snow cover for 2025-01-01..2026-07-31 and write one NetCDF per year.

2020-2024 already exist as yearly IMS-on-M36 files in projects/IMS/output/ (ims_snowcover_24km_<YYYY>_on_m36_nearest.nc4,
built from raw files that are not kept, except 2024 in projects/IMS/output/gz). This script adds 2025 and 2026 using
the functions of projects/IMS/scripts/download_ims_ascii_to_nc.py (imported, not modified):
  raw .asc.gz  ->  <out>/gz (temporary)  ->  <out>/ims_snowcover_24km_<YYYY>.nc4  (variable ims_snowcover, time,y,x)
The M36 nearest-neighbour regridding is then done with projects/IMS/scripts/regrid_ims_to_m36_nearest.py (see README).

Also spot-checks the downloader against the existing 2024 raw files (byte-identical download and identical parsed grid).

Run (login node, GEOSpyD py3.13):
  python get_ims_2025_2026.py --out /discover/nobackup/projects/land_da/VIIRS_testing/validation_data/ims
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import requests

IMS_SCRIPTS = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/scripts")
EXISTING_GZ = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/gz")
BASE = "https://noaadata.apps.nsidc.org/NOAA/G02156"


def load_ims_module():
    spec = importlib.util.spec_from_file_location("ims_dl", IMS_SCRIPTS / "download_ims_ascii_to_nc.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--years", type=int, nargs="+", default=[2025, 2026])
    ap.add_argument("--end-2026-doy", type=int, default=212, help="212 = 2026-07-31")
    ap.add_argument("--sleep", type=float, default=0.2)
    ap.add_argument("--keep-gz", action="store_true")
    a = ap.parse_args()
    ims = load_ims_module()
    out = Path(a.out); gz = out / "gz"; gz.mkdir(parents=True, exist_ok=True)
    nx = ny = 1024
    lat, lon = ims.read_latlon(IMS_SCRIPTS / "ims_24km_latlon.nc4", nx, ny)
    sess = requests.Session()
    report = {"years": {}, "spot_check_2024": []}

    # --- spot check against existing 2024 raw files (download to a separate cache dir)
    chk = out / "gz_spotcheck"; chk.mkdir(exist_ok=True)
    for doy in (15, 120, 200, 300):
        fn = ims.build_filename(2024, doy, "24km", "v1.3")
        ims.download_file_with_retries(sess, f"{BASE}/24km/2024/{fn}", chk / fn, 60, 4, 2.0)
        same_bytes = md5(chk / fn) == md5(EXISTING_GZ / fn)
        g_new = ims.parse_ims_ascii_gz(chk / fn, nx, ny); g_old = ims.parse_ims_ascii_gz(EXISTING_GZ / fn, nx, ny)
        report["spot_check_2024"].append({"file": fn, "identical_bytes": bool(same_bytes),
                                          "identical_grid": bool(np.array_equal(g_new, g_old))})
        (chk / fn).unlink()
    chk.rmdir()

    for year in a.years:
        end = a.end_2026_doy if year == 2026 else ims.days_in_year(year)
        n = end
        of = out / f"ims_snowcover_24km_{year}.nc4"
        ds, var = ims.init_output_netcdf(of, n, ny, nx, year, lat, lon)
        ok, missing, failed, miss_days = 0, 0, 0, []
        try:
            for ti, doy in enumerate(range(1, end + 1)):
                d = ims.day_from_year_doy(year, doy)
                fn = ims.build_filename(year, doy, "24km", "v1.3")
                ds["time"][ti] = float((d - ims.date(year, 1, 1)).days)
                try:
                    if not (gz / fn).exists():
                        ims.download_file_with_retries(sess, f"{BASE}/24km/{year}/{fn}", gz / fn, 60, 4, 2.0)
                        time.sleep(a.sleep)
                    var[ti, :, :] = ims.parse_ims_ascii_gz(gz / fn, nx, ny)
                    ok += 1
                except FileNotFoundError:
                    var[ti, :, :] = np.nan; missing += 1; miss_days.append(d.isoformat())
                except Exception as e:  # parse failure
                    var[ti, :, :] = np.nan; failed += 1; miss_days.append(f"{d.isoformat()} ({type(e).__name__})")
                if (doy % 60) == 0:
                    print(f"{year} doy {doy}/{end}", flush=True)
        finally:
            ds.close()
        report["years"][str(year)] = {"file": str(of), "days_expected": n, "days_ok": ok,
                                      "days_missing_404": missing, "days_failed_parse": failed, "missing_days": miss_days}
        print(report["years"][str(year)], flush=True)

    (out / "ims_fetch_report.json").write_text(json.dumps(report, indent=2))
    if not a.keep_gz:
        pass  # raw .gz are removed only after verification by the caller (verify_validation_data.py)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
