#!/usr/bin/env python3
"""
Download the IMS 4 km daily snow/ice ASCII files of all study days (periods.py) and the 4 km lat/lon grids,
for Part A3 (docs/viirs_obs_study_plan.md). Public NSIDC archive G02156 (no login).

Output (git-ignored): <data>/viirs_obs_study/ims/<YYYY>/ims<YYYY><DDD>_00UTC_4km_v1.3.asc.gz,
                      <data>/viirs_obs_study/ims/imslat_4km.bin.gz, imslon_4km.bin.gz
Existing files are skipped, so the script can be rerun.

usage: get_ims.py [P1 P2 ...]      (default: all periods)
"""

import sys
import time
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from periods import DATA_DIR, PERIODS, days  # noqa: E402

ROOT = 'https://noaadata.apps.nsidc.org/NOAA/G02156'
OUT = DATA_DIR.parent / 'ims'


def ims_file(day):
    return OUT / f'{day.year}' / f'ims{day.year}{day.timetuple().tm_yday:03d}_00UTC_4km_v1.3.asc.gz'


def fetch(session, url, path):
    if path.exists() and path.stat().st_size > 0:
        return 'have'
    path.parent.mkdir(parents=True, exist_ok=True)
    for k in range(3):
        try:
            r = session.get(url, timeout=120)
            if r.status_code == 200:
                path.write_bytes(r.content)
                return 'ok'
            status = f'HTTP {r.status_code}'
        except requests.RequestException as e:
            status = str(e)
        time.sleep(2 ** k)
    return f'FAILED ({status})'


def main():
    keys = sys.argv[1:] or list(PERIODS)
    s = requests.Session()
    for name in ('imslat_4km.bin.gz', 'imslon_4km.bin.gz'):
        print(name, fetch(s, f'{ROOT}/metadata/{name}', OUT / name))
    for k in keys:
        for d in days(k):
            p = ims_file(d)
            print(k, d, fetch(s, f'{ROOT}/4km/{d.year}/{p.name}', p), f'{p.stat().st_size / 1e6:.1f} MB' if p.exists() else '')


if __name__ == '__main__':
    main()
