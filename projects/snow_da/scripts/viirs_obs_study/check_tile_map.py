#!/usr/bin/env python3
"""
Build (or load) the CMG -> tile maps for M36 and CF0360 and check them against GEOSldas.

1. M36: every tile center falls in its own EASE cell (i_indg, j_indg).
2. End to end, per tile space: the SCF obs in GEOSldas ObsFcstAna (OFA) files are
   reproduced from the raw CMG files with lib/scf_qc.py (reader QC) + lib/obs2tile.py:
   same tile-mean SCF (to 1e-5) and same tile-mean obs lon/lat (to 1e-3 deg). The
   reader averages the lon/lat of the kept cells, so the lon/lat test checks the exact
   cell set of each tile. Each OFA obs is compared with the daily files of the day
   before, of, and after the window (the reader picks the file by local overpass time)
   and counts as matched if one of them agrees.
   - M36: LS_DAv8_M36_as_test3, MYD10C1, 2020-05-06 (local run)
   - CF0360: M21C_test_CF0360_jpss_5d4d636_OL (monitors MYD10C1, MOD10C1, VJ110C1,
     VNP10C1), 2020-01-09 and -10, copied from Discover to data/viirs_obs_study/validation_cf0360/
3. CF0360: tiles that can receive CMG cells at all, and cells per tile.

usage: check_tile_map.py
"""

import glob
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from lib import obs2tile, scf_qc  # noqa: E402
from periods import DATA_DIR  # noqa: E402

VAL = DATA_DIR.parent
CACHE = VAL / 'tilemap'

M36 = dict(
    kind='M36',
    tilecoord=Path('/Users/amfox/Desktop/LS_DAv8_M36_as_test3/output/SMAP_EASEv2_M36_GLOBAL/rc_out/'
                   'LS_DAv8_M36_as_test3.ldas_tilecoord.bin'),
    ofa=('/Users/amfox/Desktop/LS_DAv8_M36_as_test3/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y2020/M05/'
         '*ObsFcstAna.20200506_*.nc4'),
    cmg=VAL / 'validation_m36',
    products=('MYD10C1',),
)

CF_RUN = VAL / 'validation_cf0360' / 'M21C_test_CF0360_jpss_5d4d636_OL' / 'output' / 'CF0360x6C_GLOBAL'
# tile space of the 5d4d636 runs (copied from run A2 on Discover; same experiment setup as the OL run)
CF_TILECOORD = (VAL / 'validation_cf0360' / 'tiles' / 'M21C_test_CF0360_jpss_5d4d636_A2' / 'output' /
                'CF0360x6C_GLOBAL' / 'rc_out' / 'M21C_test_CF0360_jpss_5d4d636_A2.ldas_tilecoord.bin')
CF0360 = dict(
    kind='CF0360',
    tilecoord=CF_TILECOORD,
    ofa=str(CF_RUN / 'ana/ens_avg/Y2020/M01/*ObsFcstAna.20200109_*.nc4'),
    ofa2=str(CF_RUN / 'ana/ens_avg/Y2020/M01/*ObsFcstAna.20200110_*.nc4'),
    cmg=VAL / 'validation_cf0360' / 'cmg',
    products=('MYD10C1', 'MOD10C1', 'VJ110C1', 'VNP10C1'),
)


def check_m36_centers(tc):
    r, s = obs2tile.ease2_rs(tc['com_lon'], tc['com_lat'])
    i, j = obs2tile._fortran_nint(r), obs2tile._fortran_nint(s)
    bad = (i != tc['i_indg']) | (j != tc['j_indg'])
    print(f'M36 tile centers in own EASE cell: {tc["N_tile"] - bad.sum()} / {tc["N_tile"]}')
    return not bad.any()


def read_ofa(files, products):
    """OFA SCF obs per (file, product): tile (0-based), obs, lon, lat. Files are closed after reading."""
    import netCDF4
    out = {}
    for fn in files:
        with netCDF4.Dataset(fn) as d:
            ids = {str(s): int(i) for i, s in zip(d['obsparam_species_id'][:], d['obsparam_descr'][:])}
            species = np.asarray(d['species'][:])
            cols = {k: np.asarray(d[k][:]) for k in ('tilenum', 'obs', 'lon', 'lat')}
        for prod in products:
            if prod not in ids:
                continue
            sel = species == ids[prod]
            out[(fn, prod)] = (cols['tilenum'][sel] - 1, cols['obs'][sel], cols['lon'][sel], cols['lat'][sel])
    return out


def raw_tile_means(tile_map, n_tile, prod, cmg_dir):
    """Per daily file: (n cells, mean SCF, mean lon, mean lat) per tile, with the reader QC."""
    lon2d, lat2d = np.meshgrid(*obs2tile._cmg_lonlat32())
    lon2d, lat2d = lon2d.astype(np.float32), lat2d.astype(np.float32)
    out = {}
    for f in sorted(Path(cmg_dir, prod).glob('*/*.h*')):
        if prod.startswith('V'):
            scf, keep = scf_qc.viirs_scf(scf_qc.read_viirs_cmg(f))
        else:
            scf, keep = scf_qc.modis_scf(scf_qc.read_modis_cmg(f))
        n, mean = obs2tile.tile_mean(tile_map, scf, keep, n_tile)
        _, mlon = obs2tile.tile_mean(tile_map, lon2d, keep, n_tile)
        _, mlat = obs2tile.tile_mean(tile_map, lat2d, keep, n_tile)
        out[f.name.split('.')[1][1:]] = (n, mean, mlon, mlat)      # key YYYYDDD
    return out


def check_ofa(cfg, tile_map, tc):
    files = sorted(glob.glob(cfg['ofa'])) + (sorted(glob.glob(cfg['ofa2'])) if 'ofa2' in cfg else [])
    if not files:
        print(f'{cfg["kind"]}: no ObsFcstAna files found ({cfg["ofa"]}); end-to-end check skipped')
        return None
    ofa = read_ofa(files, cfg['products'])
    all_ok = True
    for prod in cfg['products']:
        raw = raw_tile_means(tile_map, tc['N_tile'], prod, cfg['cmg'])
        n_obs = n_ok = n_extra = 0
        for fn in files:
            if (fn, prod) not in ofa:
                continue
            t, obs, olon, olat = ofa[(fn, prod)]
            stamp = datetime.strptime(Path(fn).name.split('.')[-2], '%Y%m%d_%H%Mz')
            keys = [(stamp + timedelta(days=k)).strftime('%Y%j') for k in (-1, 0, 1)]
            ok = np.zeros(t.size, bool)
            for key in keys:
                if key not in raw:
                    continue
                n, mean, mlon, mlat = raw[key]
                ok |= ((n[t] > 0) & (np.abs(mean[t] - obs) < 1e-5) &
                       (np.abs(mlon[t] - olon) < 1e-3) & (np.abs(mlat[t] - olat) < 1e-3))
            n_obs += t.size
            n_ok += ok.sum()
        print(f'{cfg["kind"]} {prod}: OFA obs reproduced {n_ok} / {n_obs} ({100 * n_ok / max(n_obs, 1):.2f}%)')
        all_ok &= n_ok >= 0.9999 * n_obs       # 1 obs per product per 2 days is the known residual
    return all_ok


def check_cf_counts(cf, tc):
    n = np.bincount(cf[cf >= 0], minlength=tc['N_tile'])
    print(f'CF0360 tiles: {tc["N_tile"]}; with >= 1 CMG cell: {(n > 0).sum()} '
          f'({100 * (n > 0).mean():.2f}%); cells per tile (with cells) median {np.median(n[n > 0]):.0f}, '
          f'max {n.max()}')
    print(f'CF0360 land CMG cells mapped to a tile: {(cf >= 0).sum()}')


def main():
    m36, tc36 = obs2tile.cmg_tile_map(M36['tilecoord'], 'M36', CACHE)
    ok1 = check_m36_centers(tc36)
    ok2 = check_ofa(M36, m36, tc36)
    print(f'CF0360 tilecoord: {CF0360["tilecoord"]}')
    cf, tccf = obs2tile.cmg_tile_map(CF0360['tilecoord'], 'CF0360', CACHE)
    check_cf_counts(cf, tccf)
    ok3 = check_ofa(CF0360, cf, tccf)
    print('ALL OK (>= 99.99% of OFA obs reproduced)' if ok1 and ok2 and ok3 else 'CHECK FAILED (see counts above)')


if __name__ == '__main__':
    main()
