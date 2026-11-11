"""
Map 0.05 deg CMG cells to GEOSldas tiles, as the SCF obs readers do.

Mirrors ``get_tile_num_for_obs`` -> ``get_tile_num_from_latlon``
(GEOSldas_GridComp, ``LDAS_TileCoordRoutines.F90``; checked at commit 6d0ed32), with
FOV = 0 as for the SCF species:

- the search grid is the perturbation grid (``get_pert_grid``, ``LDAS_PertRoutines.F90``);
- **EASEv2 (M36):** the pert grid is the EASE grid; an obs goes to the single tile of
  the EASE cell that contains it (nint of the 0-based EASE row/column), none if the
  cell has no land tile;
- **cube-sphere (CF0360):** the pert grid is the global lat/lon grid
  ``DE001440xPE001080`` (4N x 3N cells for CF N; 0.25 x 1/6 deg; edges at -180/-90);
  each tile is listed in the pert cell that contains its center; an obs goes to the
  tile of its pert cell with the smallest |dlon| + |dlat| to the tile center (first
  tile wins a tie);
- in both cases the obs is dropped if it lies outside the bounding box
  (``min/max_lon/lat``) of the chosen tile.

The readers then average all kept obs of a tile (one obs per tile and window).

The map depends only on the tile space, so it is built once per tilecoord file and
cached as an int32 array (3600, 7200) of 0-based tile indices into the tilecoord
arrays (-1 = no tile).

CMG cell centers, distances and pert-cell indices are computed in float32 as in the
Fortran (cell centers with a fused multiply-add, see _cmg_lonlat32). The EASE row/column
(MAPL_ease_convert) is computed in float64 with pyproj.
"""

import struct
from pathlib import Path

import numpy as np

from .scf_qc import CMG_N_LAT, CMG_N_LON, cmg_lon_lat

NO_TILE = -1

# EASE-Grid 2.0 global (EPSG:6933), as in MAPL_ease_convert
EASE2_M36 = dict(n_lon=964, n_lat=406, cell=36032.220840584, x_min=-17367530.44516138,
                 y_max=7314540.830638505)


# ── tile files (little-endian Fortran sequential, as written by GEOSldas) ─────

def read_tilecoord(fname):
    """tile_coord fields of a GEOSldas ``*.ldas_tilecoord.bin`` as numpy arrays."""
    fields = ['tile_id', 'typ', 'pfaf', 'com_lon', 'com_lat', 'min_lon', 'max_lon',
              'min_lat', 'max_lat', 'i_indg', 'j_indg', 'frac_cell', 'frac_pfaf', 'area', 'elev']
    ints = {'tile_id', 'typ', 'pfaf', 'i_indg', 'j_indg'}
    tc = {}
    with open(fname, 'rb') as f:
        f.read(4)
        nt = struct.unpack('<i', f.read(4))[0]
        f.read(4)
        for name in fields:
            f.read(4)
            tc[name] = np.frombuffer(f.read(nt * 4), dtype='<i4' if name in ints else '<f4').copy()
            f.read(4)
    tc['N_tile'] = nt
    return tc


# ── CMG cell centers in float32 ───────────────────────────────────────────────

def _cmg_lonlat32():
    """
    CMG cell centers as the readers compute them in real*4 (read_JPSS_SCF_hdf /
    read_MODIS_SCF_hdf): lon_c = (CMG_ll_lon + 0.5*CMG_dlon) + k*CMG_dlon,
    lat_c = (CMG_ur_lat - 0.5*CMG_dlat) - k*CMG_dlat, with the multiply-add fused
    (one rounding to real*4), as the Intel compiler on Discover does. This rounding
    matters: tile bounding-box edges often lie exactly on CMG cell centers (the BCS
    raster is 1/120 deg). Checked: CF0360 ObsFcstAna 455427/455428 VJ110C1 obs and M36
    17827/17827 MYD10C1 obs reproduced; plain real*4 (two roundings) gives 95.6% on CF0360.
    """
    f = np.float32
    d = np.float64(f(0.05))                 # real*4 0.05
    c_lon = np.float64(f(-180.) + f(0.5) * f(0.05))
    c_lat = np.float64(f(90.) - f(0.5) * f(0.05))
    # float64 holds c + k*d exactly (24-bit operands), so one cast = one rounding (FMA)
    lon = (c_lon + np.arange(CMG_N_LON) * d).astype(f)
    lat = (c_lat - np.arange(CMG_N_LAT) * d).astype(f)
    return lon, lat


# ── EASEv2 ────────────────────────────────────────────────────────────────────

def ease2_rs(lon, lat, grid=EASE2_M36):
    """0-based EASEv2 column r and row s (float) of lon/lat, as MAPL_ease_convert."""
    from pyproj import Transformer
    tr = Transformer.from_crs('EPSG:4326', 'EPSG:6933', always_xy=True)
    x, y = tr.transform(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
    r = (x - grid['x_min']) / grid['cell'] - 0.5
    s = (grid['y_max'] - y) / grid['cell'] - 0.5
    return r, s


def _fortran_nint(x):
    """Fortran nint: round half away from zero."""
    return (np.sign(x) * np.floor(np.abs(x) + 0.5)).astype(np.int64)


def build_map_ease(tc, grid=EASE2_M36):
    """CMG cell -> 0-based tile index for an EASEv2 tile space (one tile per EASE cell)."""
    lookup = np.full((grid['n_lat'], grid['n_lon']), NO_TILE, dtype=np.int32)
    if np.unique(tc['j_indg'].astype(np.int64) * grid['n_lon'] + tc['i_indg']).size != tc['N_tile']:
        raise ValueError('more than one tile per EASE cell')
    lookup[tc['j_indg'], tc['i_indg']] = np.arange(tc['N_tile'], dtype=np.int32)

    lon32, lat32 = _cmg_lonlat32()
    out = np.full((CMG_N_LAT, CMG_N_LON), NO_TILE, dtype=np.int32)
    lon2 = np.broadcast_to(lon32[None, :], (1, CMG_N_LON))
    for j in range(CMG_N_LAT):                      # row by row keeps memory small
        r, s = ease2_rs(lon2[0], np.full(CMG_N_LON, lat32[j]), grid)
        i, k = _fortran_nint(r), _fortran_nint(s)
        i = np.clip(i, 0, grid['n_lon'] - 1)
        k = np.clip(k, 0, grid['n_lat'] - 1)
        out[j] = lookup[k, i]
    _bbox_check(out, tc, lon32, lat32)
    return out


# ── cube-sphere: lat/lon pert grid ────────────────────────────────────────────

def pert_grid_cs(n_cs):
    """Global lat/lon pert grid of a CF<n_cs> tile space (get_pert_grid + LDAS_create_grid_g)."""
    n_lon, n_lat = 4 * n_cs, 3 * n_cs
    return dict(name=f'DE{n_lon:06d}xPE{n_lat:06d}', n_lon=n_lon, n_lat=n_lat,
                ll_lon=np.float32(-180.), ll_lat=np.float32(-90.),
                dlon=np.float32(360.) / np.float32(n_lon), dlat=np.float32(180.) / np.float32(n_lat))


def _pert_ij(pg, lon32, lat32):
    """1-based pert-grid indices, as get_ij_ind_from_latlon for LatLon (i_dir = j_dir = +1)."""
    i = np.ceil((lon32 - pg['ll_lon']) / pg['dlon']).astype(np.int64)
    j = np.ceil((lat32 - pg['ll_lat']) / pg['dlat']).astype(np.int64)
    return i, j


def build_map_cs(tc, n_cs=360):
    """CMG cell -> 0-based tile index for a cube-sphere tile space (via the lat/lon pert grid)."""
    pg = pert_grid_cs(n_cs)
    ti, tj = _pert_ij(pg, tc['com_lon'], tc['com_lat'])
    tkey = (tj - 1) * pg['n_lon'] + (ti - 1)
    order = np.argsort(tkey, kind='stable')          # tiles of a cell in tile order
    n_cell = pg['n_lon'] * pg['n_lat']
    count = np.bincount(tkey, minlength=n_cell)
    start = np.concatenate([[0], np.cumsum(count)[:-1]])

    lon32, lat32 = _cmg_lonlat32()
    out = np.full((CMG_N_LAT, CMG_N_LON), NO_TILE, dtype=np.int32)
    ci, _ = _pert_ij(pg, lon32, np.zeros(1, np.float32))
    ci = np.clip(ci, 1, pg['n_lon'])
    for j in range(CMG_N_LAT):
        lat = lat32[j]
        cj = int(np.clip(np.ceil((lat - pg['ll_lat']) / pg['dlat']), 1, pg['n_lat']))
        key = (cj - 1) * pg['n_lon'] + (ci - 1)
        cnt, st = count[key], start[key]
        best = np.full(CMG_N_LON, NO_TILE, dtype=np.int64)
        dmin = np.full(CMG_N_LON, np.float32(1e10), dtype=np.float32)
        for k in range(int(cnt.max()) if cnt.size else 0):
            has = cnt > k
            t = order[st[has] + k]
            d = (np.abs(lon32[has] - tc['com_lon'][t]) + np.abs(lat - tc['com_lat'][t])).astype(np.float32)
            better = d < dmin[has]
            idx = np.flatnonzero(has)[better]
            dmin[idx] = d[better]
            best[idx] = t[better]
        out[j] = best
    _bbox_check(out, tc, lon32, lat32)
    return out


def _bbox_check(out, tc, lon32, lat32):
    """Drop obs outside the bounding box of the chosen tile (FOV = 0), in place."""
    for j in range(CMG_N_LAT):
        t = out[j]
        ok = t >= 0
        tt = t[ok]
        lo = lon32[ok]
        inside = ((lo >= tc['min_lon'][tt]) & (lo <= tc['max_lon'][tt]) &
                  (lat32[j] >= tc['min_lat'][tt]) & (lat32[j] <= tc['max_lat'][tt]))
        row = out[j]
        row[np.flatnonzero(ok)[~inside]] = NO_TILE


# ── cache ─────────────────────────────────────────────────────────────────────

def land_tiles(tc):
    """
    Land (catchment, typ 100) tiles only, in file order.

    The tilecoord files of the M21C CF0360 runs also list lake (typ 19) and land-ice
    (typ 20) tiles after the land tiles; the land model and its obs readers work on the
    land tiles only (ObsFcstAna tilenum <= number of land tiles), so only those can
    receive obs.
    """
    land = tc['typ'] == 100
    out = {k: v[land] for k, v in tc.items() if isinstance(v, np.ndarray)}
    out['N_tile'] = int(land.sum())
    return out


def cmg_tile_map(tilecoord_file, kind, cache_dir):
    """
    Cached CMG -> tile map over the land tiles. kind = 'M36' (EASEv2 M36) or 'CF0360'.

    Returns (map int32 (3600, 7200), land tilecoord dict); map values index the land tiles.
    """
    tc = land_tiles(read_tilecoord(tilecoord_file))
    cache = Path(cache_dir) / f'cmg2tile_v4_{kind}_{Path(tilecoord_file).stem}.npy'   # v4: FMA real*4 centers, land tiles
    if cache.exists():
        return np.load(cache), tc
    m = build_map_ease(tc) if kind == 'M36' else build_map_cs(tc, 360)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache, m)
    return m, tc


def tile_mean(tile_map, values, keep, n_tile):
    """Per-tile number of kept cells and mean value (NaN where none), as the readers average."""
    t = tile_map[keep]
    v = values[keep]
    ok = t >= 0
    n = np.bincount(t[ok], minlength=n_tile)
    s = np.bincount(t[ok], weights=v[ok].astype(np.float64), minlength=n_tile)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return n, mean
