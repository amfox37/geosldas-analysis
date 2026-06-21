"""
M36 EASE-Grid tile assignment and super-obs formation.

Uses GEOSldas tilecoord/tilegrids information to reproduce the relevant
pieces of:
  LDAS_TileCoordRoutines.F90  — get_tile_num_from_latlon  (cell index lookup)
  clsm_ensupd_read_obs.F90    — super-obs averaging loop  (arithmetic mean per tile)
"""

import numpy as np
from pyproj import Proj

# ── EASE-Grid 2.0 M36 constants ───────────────────────────────────────────────
EASE2     = Proj('+proj=cea +lat_ts=30 +datum=WGS84 +lon_0=0 +units=m')
M36_CELL  = 36032.220840584          # metres
X_ORIGIN, Y_ORIGIN = EASE2(-180, 90)


def latlon_to_ij(lat, lon):
    """Lat/lon → 0-based M36 EASE cell indices (col i, row j)."""
    x, y = EASE2(np.asarray(lon, float), np.asarray(lat, float))
    i = np.floor((x - X_ORIGIN) / M36_CELL).astype(int)
    j = np.floor((Y_ORIGIN - y) / M36_CELL).astype(int)
    return i, j


def ij_to_corners(i, j):
    """Return (lons, lats) of 4 CCW corners of M36 cell (i, j)."""
    x_lo = X_ORIGIN + i * M36_CELL
    x_hi = X_ORIGIN + (i + 1) * M36_CELL
    y_hi = Y_ORIGIN - j * M36_CELL
    y_lo = Y_ORIGIN - (j + 1) * M36_CELL
    lons, lats = EASE2([x_lo, x_hi, x_hi, x_lo],
                       [y_lo, y_lo, y_hi, y_hi], inverse=True)
    return np.array(lons), np.array(lats)


def _as_clean_gridtype(tile_grid):
    return str(tile_grid.get('gridtype', '')).strip()


def build_tile_lookup(tile_coord):
    """Build a GEOS M36 cell -> tile lookup from tilecoord arrays.

    For EASEv2 M36 land DA domains there should be at most one land tile per
    grid cell. If duplicates exist, keep the largest frac_cell, then the
    smallest tile_id, matching the representative-tile convention used
    elsewhere in this repo.
    """
    i_indg = np.asarray(tile_coord['i_indg'], dtype=np.int64)
    j_indg = np.asarray(tile_coord['j_indg'], dtype=np.int64)
    tile_id = np.asarray(tile_coord['tile_id'], dtype=np.int64)
    frac_cell = np.asarray(tile_coord.get('frac_cell', np.ones_like(tile_id)), dtype=float)

    order = np.lexsort((tile_id, -frac_cell, j_indg, i_indg))
    lookup = {}
    for idx in order:
        key = (int(i_indg[idx]), int(j_indg[idx]))
        if key not in lookup:
            lookup[key] = idx
    return lookup


def obs_to_tiles(lat, lon, tile_coord, tile_grid=None, lookup=None):
    """Assign obs locations to GEOS tile numbers using tilecoord grid cells.

    Returns
    -------
    dict with arrays:
        tilenum     GEOS tile number / tile_id for each obs, -1 if no land tile
        tile_index  zero-based index into tile_coord arrays, -1 if no land tile
        ij          M36 EASE cell indices from the obs lat/lon
    """
    if tile_grid is not None:
        gridtype = _as_clean_gridtype(tile_grid)
        if 'EASEv2_M36' not in gridtype:
            raise ValueError(f'Only EASEv2_M36 tile grids are supported here, got {gridtype!r}')

    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    i_arr, j_arr = latlon_to_ij(lat, lon)

    if lookup is None:
        lookup = build_tile_lookup(tile_coord)

    tile_id = np.asarray(tile_coord['tile_id'], dtype=np.int64)
    out_tile = np.full(len(lat), -1, dtype=np.int64)
    out_idx = np.full(len(lat), -1, dtype=np.int64)

    for n, key in enumerate(zip(i_arr, j_arr)):
        idx = lookup.get((int(key[0]), int(key[1])))
        if idx is not None:
            out_idx[n] = idx
            out_tile[n] = tile_id[idx]

    return {'tilenum': out_tile, 'tile_index': out_idx, 'ij': np.column_stack([i_arr, j_arr])}


def form_super_obs(lat, lon, ssm, window=None, cycle=None, tile_coord=None, tile_grid=None, lookup=None):
    """Reproduce GEOSldas super-obs: arithmetic mean of all QC'd obs per tile.

    If window is provided, obs are grouped by (tile_i, tile_j, window) so that
    each assimilation window gets its own super-obs, matching GEOSldas behaviour.

    Parameters
    ----------
    lat, lon, ssm : array-like
        Obs locations and SSM values (% saturation).
    window : array-like of int (0–7) or None
        GEOSldas assimilation window index per obs. If None, all obs are pooled.
    cycle : array-like of int or None
        GEOSldas cycle index relative to the requested data date. Use this
        when available to avoid mixing same-day 0000z with next-day 0000z.
    tile_coord, tile_grid : dict or None
        GEOSldas tilecoord/tilegrids structures. If provided, obs are assigned
        to GEOS tile numbers via tilecoord before averaging. If omitted, the
        fallback grouping is by M36 i/j cell from the EASE constants above.
    lookup : dict or None
        Optional precomputed output from build_tile_lookup(tile_coord). Supplying
        this avoids rebuilding the same tile lookup in repeated global calls.

    Returns
    -------
    dict with arrays:
        tilenum — GEOS tile number per super-ob (-1 for fallback mode)
        ij      — (N, 2) M36 cell indices [col, row]
        cycle   — GEOSldas cycle per super-ob (or -1 if not split)
        window  — window index per super-ob (or -1 if not split)
        ssm     — mean SSM per super-ob
        lat     — mean lat per super-ob
        lon     — mean lon per super-ob
        count   — N obs contributing to each super-ob
        ssm_std — within-super-ob standard deviation
        ssm_min — minimum raw SSM contributing to the super-ob
        ssm_max — maximum raw SSM contributing to the super-ob
    """
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    ssm = np.asarray(ssm, float)

    empty = dict(
        tilenum=np.array([], dtype=int),
        ij=np.empty((0, 2), dtype=int),
        cycle=np.array([], dtype=int),
        window=np.array([], dtype=int),
        ssm=np.array([], dtype=float),
        lat=np.array([], dtype=float),
        lon=np.array([], dtype=float),
        count=np.array([], dtype=int),
        ssm_std=np.array([], dtype=float),
        ssm_min=np.array([], dtype=float),
        ssm_max=np.array([], dtype=float),
    )
    if len(ssm) == 0:
        return empty

    if tile_coord is not None:
        assigned = obs_to_tiles(lat, lon, tile_coord, tile_grid=tile_grid, lookup=lookup)
        valid = assigned['tilenum'] > 0
        lat = lat[valid]
        lon = lon[valid]
        ssm = ssm[valid]
        tile_arr = assigned['tilenum'][valid]
        ij_arr = assigned['ij'][valid]
        if len(ssm) == 0:
            return empty
        if window is not None:
            window = np.asarray(window, int)[valid]
        if cycle is not None:
            cycle = np.asarray(cycle, int)[valid]
    else:
        i_arr, j_arr = latlon_to_ij(lat, lon)
        tile_arr = np.full(len(lat), -1, dtype=np.int64)
        ij_arr = np.column_stack([i_arr, j_arr])

    if cycle is not None:
        cyc_arr = np.asarray(cycle, int)
        keys = np.column_stack([tile_arr, cyc_arr])
    elif window is not None:
        win_arr = np.asarray(window, int)
        keys = np.column_stack([tile_arr, win_arr])
    else:
        keys = tile_arr[:, None]

    unique_keys, inv = np.unique(keys, axis=0, return_inverse=True)
    n = len(unique_keys)

    out_count = np.bincount(inv, minlength=n)
    out_ssm = np.bincount(inv, weights=ssm, minlength=n) / out_count
    out_lat = np.bincount(inv, weights=lat, minlength=n) / out_count
    out_lon = np.bincount(inv, weights=lon, minlength=n) / out_count

    out_sumsq = np.bincount(inv, weights=ssm * ssm, minlength=n)
    out_var = np.maximum(out_sumsq / out_count - out_ssm * out_ssm, 0.0)
    out_std = np.sqrt(out_var)

    order = np.argsort(inv, kind='stable')
    starts = np.r_[0, np.cumsum(out_count)[:-1]]
    sorted_ssm = ssm[order]
    out_min = np.minimum.reduceat(sorted_ssm, starts)
    out_max = np.maximum.reduceat(sorted_ssm, starts)
    out_ij = ij_arr[order[starts]].astype(int)

    out_tile = unique_keys[:, 0].astype(int)
    out_cycle = unique_keys[:, 1].astype(int) if (cycle is not None or window is not None) else np.full(n, -1)
    if cycle is not None:
        out_win = out_cycle % 8
    elif window is not None:
        out_win = out_cycle
    else:
        out_win = np.full(n, -1)

    return dict(tilenum=out_tile, ij=out_ij, cycle=out_cycle, window=out_win, ssm=out_ssm,
                lat=out_lat, lon=out_lon, count=out_count, ssm_std=out_std,
                ssm_min=out_min, ssm_max=out_max)


def match_tiles(so_a, so_b):
    """Find super-obs in so_a and so_b that share the same tile (and window).

    Both inputs are dicts as returned by form_super_obs.

    Returns
    -------
    (vals_a, vals_b) — paired SSM arrays for matched tiles/windows
    """
    time_key = 'cycle' if 'cycle' in so_a and 'cycle' in so_b else 'window'
    key_a = {tuple(row): so_a['ssm'][k]
             for k, row in enumerate(np.column_stack([so_a['tilenum'], so_a[time_key]]))}
    key_b = {tuple(row): so_b['ssm'][k]
             for k, row in enumerate(np.column_stack([so_b['tilenum'], so_b[time_key]]))}

    common = set(key_a.keys()) & set(key_b.keys())
    if not common:
        return np.array([]), np.array([])

    vals_a = np.array([key_a[k] for k in common])
    vals_b = np.array([key_b[k] for k in common])
    return vals_a, vals_b
