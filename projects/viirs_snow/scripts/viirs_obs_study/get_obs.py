#!/usr/bin/env python3
"""
Download the VIIRS and MODIS daily CMG files of all study periods (periods.py).

Uses ../viirs_scf_reader_test/get_snow_cmg_obs.py (NSIDC via CMR; needs Earthdata login in ~/.netrc).
Files that already exist are skipped, so the script can be rerun after an interruption.
At the end it writes <DATA_DIR>/manifest.txt: one line per expected file with size and original
file name (LocalGranuleID for VIIRS), or MISSING.

usage: get_obs.py [P1 P2 ...]      (default: all periods)
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from periods import DATA_DIR, PERIODS, cmg_file, days, products  # noqa: E402

GET = HERE.parent / 'viirs_scf_reader_test' / 'get_snow_cmg_obs.py'


def granule_id(path):
    """Original file name stored in the file (VIIRS LocalGranuleID; MODIS LOCALGRANULEID in CoreMetadata.0)."""
    import netCDF4
    try:
        with netCDF4.Dataset(path) as ds:
            if path.suffix == '.h5':
                return getattr(ds, 'LocalGranuleID', '')
            meta = getattr(ds, 'CoreMetadata.0', '')
            i = meta.find('LOCALGRANULEID')
            j = meta.find('END_OBJECT', i)
            seg = meta[i:j] if i >= 0 and j > i else ''
            q = seg.split('"')
            return q[1] if len(q) >= 3 else ''
    except OSError:
        return 'UNREADABLE'


def main():
    periods = sys.argv[1:] or list(PERIODS)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for p in periods:
        first, last, why = PERIODS[p]
        print(f'== {p} {first}..{last} ({why}): {" ".join(products(p))}', flush=True)
        subprocess.run([sys.executable, str(GET), str(DATA_DIR), first, last, *products(p)], check=True)

    lines, n_missing = [], 0
    for p in PERIODS:
        for prod in products(p):
            for d in days(p):
                f = cmg_file(prod, d)
                if f.exists():
                    lines.append(f'{p} {prod} {d} {f.stat().st_size:10d} {granule_id(f)}')
                else:
                    lines.append(f'{p} {prod} {d} MISSING')
                    n_missing += 1
    (DATA_DIR / 'manifest.txt').write_text('\n'.join(lines) + '\n')
    print(f'manifest: {len(lines)} files expected, {n_missing} missing -> {DATA_DIR / "manifest.txt"}')


if __name__ == '__main__':
    main()
