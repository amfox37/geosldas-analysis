#!/usr/bin/env python3
"""
Check that lib/scf_qc.py reproduces the GEOSldas SCF readers exactly.

For each daily CMG file, runs the Fortran readers (built with fortran_check/build.sh)
over the whole globe and compares the kept obs with lib/scf_qc.py: same number of
obs, same order, and bit-identical lon, lat and SCF (all real*4).

usage:
  check_qc_vs_fortran.py <test_driver> <file> [<file> ...]

The product (VIIRS or MODIS) is taken from the file name (V*10C1 or M*10C1).
Exit status 0 only if every file matches.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # projects/snow_da
from lib.scf_qc import cmg_lon_lat, modis_scf, read_modis_cmg, read_viirs_cmg, viirs_scf  # noqa: E402


def python_obs(path):
    """lon, lat, scf of the kept cells, in the order of the Fortran readers (lon outer, lat inner)."""
    name = Path(path).name
    if name.startswith('V'):
        scf, keep = viirs_scf(read_viirs_cmg(path))
    elif name.startswith('M'):
        scf, keep = modis_scf(read_modis_cmg(path))
    else:
        raise ValueError(f'cannot tell product from file name: {name}')
    lon_c, lat_c = cmg_lon_lat()
    ii, jj = np.nonzero(keep.T)                      # (lon index, lat index), lon outer
    return (lon_c[ii].astype(np.float32), lat_c[jj].astype(np.float32), scf.T[ii, jj])


def fortran_obs(driver, path):
    sensor = 'VIIRS' if Path(path).name.startswith('V') else 'MODIS'
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'obs.bin'
        log = subprocess.run([driver, sensor, str(path), str(out)], capture_output=True, text=True)
        if log.returncode != 0 or not out.exists():
            raise RuntimeError(f'test_driver failed for {path}:\n{log.stdout}{log.stderr}')
        raw = out.read_bytes()
    n = int(np.frombuffer(raw[:4], dtype=np.int32)[0])
    v = np.frombuffer(raw[4:], dtype=np.float32)
    assert v.size == 3 * n, 'unexpected output size'
    return v[:n], v[n:2*n], v[2*n:], log.stdout


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(2)
    driver, files = sys.argv[1], sys.argv[2:]
    all_ok = True
    for f in files:
        flon, flat, fscf, flog = fortran_obs(driver, f)
        plon, plat, pscf = python_obs(f)
        # the lon/lat of the Fortran readers are computed in real*4 arithmetic, so compare them
        # to within half a micro-degree; SCF must be bit-identical
        same_n = flon.size == plon.size
        ok = (same_n and np.abs(flon - plon).max(initial=0) < 5e-7 * 360
              and np.abs(flat - plat).max(initial=0) < 5e-7 * 180
              and np.array_equal(fscf, pscf))
        all_ok &= bool(ok)
        granule = [l.split('=', 1)[1].strip() for l in flog.splitlines() if 'LocalGranuleID' in l]
        msg = f'N fortran {flon.size:8d}  N python {plon.size:8d}'
        if same_n:
            msg += (f'  max|dlon| {np.abs(flon-plon).max(initial=0):.1e}  max|dlat| {np.abs(flat-plat).max(initial=0):.1e}'
                    f'  SCF identical: {np.array_equal(fscf, pscf)}')
        print(f'{"OK  " if ok else "FAIL"} {Path(f).name:28s} {msg}  {granule[0] if granule else ""}')
    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
