"""
Map 0.05 deg CMG cells to GEOSldas tiles, as the SCF obs readers do.

Mirrors ``get_tile_num_for_obs`` -> ``get_tile_num_from_latlon``
(GEOSldas_GridComp, ``LDAS_TileCoordRoutines.F90``; M36 checked at 6d0ed32, cube-sphere
follows PR #252 at e36655c), with FOV = 0 as for the SCF species:

- **EASEv2 (M36):** an obs goes to the single tile of the EASE cell that contains it (nint of
  the 0-based EASE row/column), none if the cell has no land tile;
- **cube-sphere (CF0360):** as GEOSldas_GridComp PR #252 (``get_tile_num_from_latlon_c3``):
  the cube-sphere cell of the obs (``MAPL_GetGlobalHorzIJIndex``, reimplemented in ``cs_ij``)
  and, if no tile bounding box in it contains the obs, the cells of 8 points at +/- half a
  cell; among the tiles whose bounding box contains the obs, the smallest normalized
  distance (dx/hx)**2 + (dy/hy)**2 wins. FOV = 0 only. (The earlier lat/lon pert-grid
  search, which lost ~45% of land CMG cells, is gone.)

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


# ── cube-sphere: cube-sphere grid cells (GEOSldas_GridComp PR #252) ───────────

_SQR2 = 1.41421356237310
_ALPHA = 0.615479708670387          # asin(1/sqrt(3))
_SHIFT = 0.174532925199433          # pi/18, the "Fuji" shift of the GEOS cube-sphere grid
_EPS = 2.220446049250313e-16        # epsilon(1.0d0)


def cs_ij(lon_deg, lat_deg, n_cs):
    """
    1-based global cube-sphere cell (i, j) of lon/lat [deg]: i = 1..N, j = 1..6N (six faces
    stacked in j), as MAPL_GetGlobalHorzIJIndex (unstretched gnomonic equiangular grid).
    lon/lat are converted to radians in float64, as get_ij_ind_from_latlon_c3 does.
    """
    lon = np.asarray(lon_deg, dtype=np.float64) * (np.pi / 180.) + _SHIFT
    lat = np.asarray(lat_deg, dtype=np.float64) * (np.pi / 180.)
    xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    xyz = xyz / np.abs(xyz).max(axis=0)
    x, y, z = xyz
    dalpha = 2. * _ALPHA / n_cs

    def idx(u, v):
        return (np.ceil((np.arctan(u / _SQR2) + _ALPHA) / dalpha).astype(np.int64),
                np.ceil((np.arctan(v / _SQR2) + _ALPHA) / dalpha).astype(np.int64))

    # faces tested in MAPL's order, so edge points go to the lowest face
    tests = [(np.abs(x - 1) <= _EPS, y, z), (np.abs(y - 1) <= _EPS, -x, z),
             (np.abs(z - 1) <= _EPS, -x, -y), (np.abs(x + 1) <= _EPS, -z, -y),
             (np.abs(y + 1) <= _EPS, -z, x), (np.abs(z + 1) <= _EPS, y, x)]
    ii = np.full(lon.shape, -1, dtype=np.int64)
    jj = np.full(lon.shape, -1, dtype=np.int64)
    done = np.zeros(lon.shape, dtype=bool)
    for face, (m, u, v) in enumerate(tests):
        m = m & ~done
        if m.any():
            i, j = idx(u[m], v[m])
            j = np.where(j == 0, 1, j)
            j = np.where(j == n_cs + 1, n_cs, j)
            ii[m], jj[m] = i, j + n_cs * face
            done |= m
    ii = np.where(ii == 0, 1, ii)
    ii = np.where(ii == n_cs + 1, n_cs, ii)
    return ii, jj


def _cs_tile_lists(tc, n_cs):
    """CSR lists of the tiles of each cube-sphere cell (tile order), from tilecoord i/j_indg."""
    key = (tc['j_indg'].astype(np.int64) - 1) * n_cs + (tc['i_indg'].astype(np.int64) - 1)
    order = np.argsort(key, kind='stable')
    count = np.bincount(key, minlength=n_cs * 6 * n_cs)
    start = np.concatenate([[0], np.cumsum(count)[:-1]])
    return order, count, start


def _best_in_cells(lon, lat, ci, cj, tc, lists, n_cs, best, dmin):
    """
    get_tile_num_in_bbox_c3 for one cell per point: among the tiles of the cell whose bounding
    box contains the point, keep the one with the smallest (dx/hx)^2 + (dy/hy)^2 (strictly
    smaller than dmin). Updates best, dmin in place.
    """
    order, count, start = lists
    key = (cj - 1) * n_cs + (ci - 1)
    cnt, st = count[key], start[key]
    f = np.float32
    for k in range(int(cnt.max()) if cnt.size else 0):
        has = np.flatnonzero(cnt > k)
        t = order[st[has] + k]
        lo, la = lon[has], lat[has]
        inside = ((tc['min_lon'][t] <= lo) & (tc['max_lon'][t] >= lo) &
                  (tc['min_lat'][t] <= la) & (tc['max_lat'][t] >= la))
        hx = np.maximum(f(0.5) * (tc['max_lon'][t] - tc['min_lon'][t]), f(1e-6))
        hy = np.maximum(f(0.5) * (tc['max_lat'][t] - tc['min_lat'][t]), f(1e-6))
        d = ((np.abs(lo - tc['com_lon'][t]) / hx) ** 2 + (np.abs(la - tc['com_lat'][t]) / hy) ** 2).astype(f)
        better = inside & (d < dmin[has])
        idx = has[better]
        dmin[idx] = d[better]
        best[idx] = t[better]


def build_map_cs(tc, n_cs=360):
    """
    CMG cell -> 0-based tile index for a cube-sphere tile space, FOV = 0
    (get_tile_num_from_latlon_c3, GEOSldas_GridComp PR #252):

    1. tiles in the cube-sphere cell that contains the obs: among those whose bounding box
       contains it, the smallest (dx/hx)^2 + (dy/hy)^2;
    2. if none: the same over the cells of the 8 points at +/- half a cell (0.5*dlon, dlon =
       360/(4N); lon offset divided by max(cos(lat), 0.1)); the max_dist lattice and the
       nearest-center fallback (step 3) need max_dist > 0 and do nothing for FOV = 0.
    No tile if both fail. The bounding-box/max_dist drop is then satisfied by construction.
    """
    f = np.float32
    lists = _cs_tile_lists(tc, n_cs)
    lon32, lat32 = _cmg_lonlat32()
    dlon = f(360.) / f(4 * n_cs)
    out = np.full((CMG_N_LAT, CMG_N_LON), NO_TILE, dtype=np.int32)
    for j in range(CMG_N_LAT):
        lat = np.full(CMG_N_LON, lat32[j], dtype=f)
        best = np.full(CMG_N_LON, NO_TILE, dtype=np.int64)
        dmin = np.full(CMG_N_LON, f(1e10), dtype=f)
        ci, cj = cs_ij(lon32, lat, n_cs)
        _best_in_cells(lon32, lat, ci, cj, tc, lists, n_cs, best, dmin)
        todo = np.flatnonzero(best < 0)
        if todo.size:
            lo, la = lon32[todo], lat[todo]
            dlat8 = f(0.5) * dlon
            dlon8 = dlat8 / np.maximum(np.cos(f(np.pi) / f(180.) * la), f(0.1))
            bt = np.full(todo.size, NO_TILE, dtype=np.int64)
            dm = np.full(todo.size, f(1e10), dtype=f)
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    nlo = (lo + f(di) * dlon8).astype(f)
                    nla = (la + f(dj) * dlat8).astype(f)
                    ni, nj = cs_ij(nlo, nla, n_cs)
                    # search cells of the neighbor points, test the bounding box at the obs
                    _best_in_cells(lo, la, ni, nj, tc, lists, n_cs, bt, dm)
            best[todo] = bt
        out[j] = best
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
    cache = Path(cache_dir) / f'cmg2tile_v5_{kind}_{Path(tilecoord_file).stem}.npy'   # v5: CF via cube-sphere cells (PR #252); v4 was the old pert-grid search
    if cache.exists():
        return np.load(cache), tc
    m = build_map_ease(tc) if kind == 'M36' else build_map_cs(tc, int(kind[2:]))
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
