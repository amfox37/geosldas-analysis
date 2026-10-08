"""Fetch a small subset of the H SAF-only cycling runs (super-obs vs thinning) for local plotting.

Run on Discover with GEOSpyD:
    /usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3 fetch_cycling_superob_vs_thin.py

Writes $NOBACKUP/cf0360_cycling_sob_vs_thin/ and a .tar.gz of it. Copy the tarball to
~/Desktop/M21C_land_local/ on the Mac and unpack it there.

Per experiment, day and window:
  - the ens_avg ObsFcstAna file (copied as is)
  - SFMC_FCST and SFMC_ANA (ensemble mean, tile space) from inst3_1d_lndfcstana_Nt, as float32 .npy
Plus the raw H121 files (Metop-A/B/C) for the days not already on the Mac.
"""
import datetime as dt
import glob
import os
import shutil
import tarfile
from pathlib import Path

import numpy as np
import netCDF4 as nc

ROOT = Path('/discover/nobackup/projects/land_da/M21C_testing')
EXPS = {
    'sob025fp': 'M21C_test_CF0360_sob025fp_HSAF_xc03125',
    'thin25km': 'M21C_test_CF0360_thin25km_HSAF_xc03125',
    'OL_monitor': 'M21C_test_CF0360_OL_monitor_HSAF',
}
DAYS = [dt.date(2019, 5, 1), dt.date(2019, 5, 10)]
H121_ROOT = Path('/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/H121')
H121_DAYS = [dt.date(2019, 5, 10)]          # 2019-05-01 is already local
SATS = ['metop_a', 'metop_b', 'metop_c']

OUT = Path(os.environ.get('NOBACKUP', '.')) / 'cf0360_cycling_sob_vs_thin'
OUT.mkdir(parents=True, exist_ok=True)

missing = []
for key, exp in EXPS.items():
    d_out = OUT / key
    d_out.mkdir(exist_ok=True)
    base = ROOT / exp / 'output' / 'CF0360x6C_GLOBAL'
    for day in DAYS:
        for h in range(0, 24, 3):
            t = dt.datetime(day.year, day.month, day.day, h)
            tag = f'{t:%Y%m%d_%H}00z'
            ymd = f'Y{t:%Y}/M{t:%m}'
            ofa = glob.glob(str(base / 'ana' / 'ens_avg' / ymd / f'*ObsFcstAna.{tag}.nc4'))
            if ofa:
                shutil.copy2(ofa[0], d_out / Path(ofa[0]).name)
            else:
                missing.append(f'{key} ObsFcstAna {tag}')
            fa = glob.glob(str(base / 'cat' / 'ens_avg' / ymd / f'*inst3_1d_lndfcstana_Nt.{tag}.nc4'))
            if not fa:
                missing.append(f'{key} lndfcstana {tag}')
                continue
            with nc.Dataset(fa[0]) as d:
                for v in ['SFMC_FCST', 'SFMC_ANA']:
                    if v in d.variables:
                        np.save(d_out / f'{v.lower()}_{tag}.npy', np.asarray(d[v][0, :], np.float32))
                    else:
                        missing.append(f'{key} {v} {tag} (variables: {list(d.variables)[:20]})')
    print(f'{key}: {len(list(d_out.iterdir()))} files')

# raw H121: files whose sensing start falls in [day - 3 h, day + 24 h), which covers the 00z-21z windows
for day in H121_DAYS:
    t0 = dt.datetime(day.year, day.month, day.day) - dt.timedelta(hours=3)
    for sat in SATS:
        d_out = OUT / 'H121' / sat
        d_out.mkdir(parents=True, exist_ok=True)
        n = 0
        for d in {t0.date(), day}:
            for f in glob.glob(str(H121_ROOT / sat / f'Y{d:%Y}' / f'M{d:%m}' / f'*_{d:%Y%m%d}??0000_*.nc')):
                start = dt.datetime.strptime(Path(f).name.split('_')[-6], '%Y%m%d%H%M%S')
                if t0 <= start < t0 + dt.timedelta(hours=27):
                    shutil.copy2(f, d_out / Path(f).name)
                    n += 1
        print(f'H121 {sat} {day}: {n} files')

print('\nmissing:' if missing else '\nnothing missing')
for m in missing:
    print('  ', m)

tar = OUT.with_suffix('.tar.gz')
with tarfile.open(tar, 'w:gz') as tf:
    tf.add(OUT, arcname=OUT.name)
print(f'\nwrote {tar} ({tar.stat().st_size / 1e6:.0f} MB)')
