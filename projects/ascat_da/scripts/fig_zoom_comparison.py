"""
Zoomed 1-degree box: Legacy BUFR vs H121 obs + reproduced GEOSldas super-obs.

Steps reproduced from Fortran:
  1. QC filter obs (same flags as read_obs_sm_ASCAT_EUMET / read_obs_sm_ASCAT_HSAF)
  2. Assign each obs to M36 EASE tile via grid-cell index lookup  (get_tile_num_from_latlon)
  3. Average obs within each tile  (super-obs formation loop)
  4. Validate against OFA monitor-run `obs` field
"""

import numpy as np
import glob, os
from datetime import datetime, timedelta

import eccodes as ecc
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Proj

# ── EASE-Grid 2.0 M36 constants ──────────────────────────────────────────────
EASE2     = Proj('+proj=cea +lat_ts=30 +datum=WGS84 +lon_0=0 +units=m')
M36_CELL  = 36032.220840584                      # metres
X_ORIGIN, Y_ORIGIN = EASE2(-180, 90)            # upper-left corner of global grid


def latlon_to_ease_ij(lat, lon):
    """Obs lat/lon → 0-based M36 EASE cell index (i=col, j=row).
    Reproduces get_ij_ind_from_latlon in LDAS_TileCoordRoutines.F90."""
    x, y = EASE2(np.asarray(lon, float), np.asarray(lat, float))
    i = np.floor((x - X_ORIGIN) / M36_CELL).astype(int)
    j = np.floor((Y_ORIGIN - y) / M36_CELL).astype(int)
    return i, j


def ease_ij_to_corners(i, j):
    """Return lon/lat of 4 corners (CCW) of M36 cell (i, j)."""
    x_lo = X_ORIGIN + i * M36_CELL
    x_hi = X_ORIGIN + (i + 1) * M36_CELL
    y_hi = Y_ORIGIN - j * M36_CELL
    y_lo = Y_ORIGIN - (j + 1) * M36_CELL
    lons, lats = EASE2([x_lo, x_hi, x_hi, x_lo],
                        [y_lo, y_lo, y_hi, y_hi], inverse=True)
    return np.array(lons), np.array(lats)


def compute_super_obs(lat, lon, ssm):
    """Reproduce GEOSldas super-obs formation: arithmetic mean per M36 tile.

    Returns:
        tile_ij  — (N_tiles, 2) array of (i, j) indices
        so_ssm   — mean SSM per tile (same units as input ssm)
        so_lat   — mean lat per tile
        so_lon   — mean lon per tile
        so_cnt   — N obs per tile
        obs_tile — per-obs tile index into tile_ij rows (-1 = outside domain)
    """
    i_arr, j_arr = latlon_to_ease_ij(lat, lon)
    ij = np.stack([i_arr, j_arr], axis=1)
    unique_ij, inv = np.unique(ij, axis=0, return_inverse=True)
    n = len(unique_ij)
    so_ssm = np.array([ssm[inv == k].mean() for k in range(n)])
    so_lat = np.array([lat[inv == k].mean() for k in range(n)])
    so_lon = np.array([lon[inv == k].mean() for k in range(n)])
    so_cnt = np.array([(inv == k).sum()      for k in range(n)])
    return unique_ij, so_ssm, so_lat, so_lon, so_cnt, inv


# ── Paths ─────────────────────────────────────────────────────────────────────
BUFR_BASE = '/Users/amfox/Desktop/ASCAT_SSM_CDR/legacy_bufr'
H121_BASE = '/Users/amfox/Desktop/ASCAT_SSM_CDR/H121'
OUT_DIR   = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'report'))

DATE = datetime(2020, 1, 1)

PLATFORMS = {
    'Metop-A': {'bufr_prefix': 'M02-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_a'},
    'Metop-B': {'bufr_prefix': 'M01-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_b'},
    'Metop-C': {'bufr_prefix': 'M03-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_c'},
}

LAT0, LON0 = 35.0, -99.0
LAT1, LON1 = LAT0 + 3.0, LON0 + 3.0


# ── Readers ───────────────────────────────────────────────────────────────────
def read_bufr_qcd(bufr_dir, date, prefix):
    KEYS = {'ssm': 'surfaceSoilMoisture', 'smpf': 'soilMoistureProcessingFlag',
            'smcf': 'soilMoistureCorrectionFlag', 'tpcx': 'topographicComplexity',
            'iwfr': 'inundationAndWetlandFraction', 'lat': 'latitude', 'lon': 'longitude'}
    lats, lons, ssms = [], [], []
    for f in sorted(glob.glob(os.path.join(bufr_dir, f'{prefix}{date.strftime("%Y%m%d")}*.bfr'))):
        with open(f, 'rb') as fh:
            bufr = ecc.codes_bufr_new_from_file(fh)
            while bufr is not None:
                ecc.codes_set(bufr, 'unpack', 1)
                try:
                    arrays = {k: ecc.codes_get_array(bufr, v) for k, v in KEYS.items()}
                    n = len(arrays['lat'])
                    for k in arrays:
                        if len(arrays[k]) == 1:
                            arrays[k] = np.repeat(arrays[k], n)
                    q = ((arrays['ssm'] >= 0) & (arrays['ssm'] <= 100) &
                         (arrays['smpf'] == 0) &
                         ((arrays['smcf'] == 0) | (arrays['smcf'] == 4)) &
                         (arrays['tpcx'] >= 0) & (arrays['tpcx'] <= 10) &
                         (arrays['iwfr'] >= 0) & (arrays['iwfr'] <= 10))
                    lats.extend(arrays['lat'][q].tolist())
                    lons.extend(arrays['lon'][q].tolist())
                    ssms.extend(arrays['ssm'][q].tolist())
                except Exception:
                    pass
                ecc.codes_release(bufr)
                bufr = ecc.codes_bufr_new_from_file(fh)
    return np.array(lats), np.array(lons), np.array(ssms)


def read_h121_qcd(h121_dir, date):
    lats, lons, ssms = [], [], []
    for f in sorted(glob.glob(os.path.join(h121_dir, f'*_{date.strftime("%Y%m%d")}*.nc'))):
        ds = nc.Dataset(f)
        lat = ds.variables['latitude'][:].filled(np.nan)
        lon = ds.variables['longitude'][:].filled(np.nan)
        ssm = ds.variables['surface_soil_moisture'][:].filled(np.nan)
        pf  = ds.variables['processing_flag'][:].filled(255)
        cf  = ds.variables['correction_flag'][:].filled(255)
        wf  = ds.variables['wetland_fraction'][:].filled(-128)
        tc  = ds.variables['topographic_complexity'][:].filled(-128)
        ds.close()
        q = ((ssm >= 0) & (ssm <= 100) & (pf == 0) &
             ((cf == 0) | (cf == 4)) & (wf >= 0) & (wf <= 10) & (tc >= 0) & (tc <= 10))
        lats.extend(lat[q].tolist())
        lons.extend(lon[q].tolist())
        ssms.extend(ssm[q].tolist())
    return np.array(lats), np.array(lons), np.array(ssms)


def inbox(lat, lon):
    return (lat >= LAT0) & (lat <= LAT1) & (lon >= LON0) & (lon <= LON1)


def deduplicate(lat, lon, ssm):
    """Same DGG node visited by multiple passes → mean SSM, scaled marker size."""
    keys = np.round(np.stack([lat, lon], axis=1), 4)
    unique_keys, inv = np.unique(keys, axis=0, return_inverse=True)
    n = len(unique_keys)
    mean_ssm = np.array([ssm[inv == k].mean() for k in range(n)])
    counts   = np.array([(inv == k).sum()      for k in range(n)])
    return unique_keys[:, 0], unique_keys[:, 1], mean_ssm, counts



# ── Load per-platform raw obs + compute super-obs ────────────────────────────
plat_data = {}
for plat, cfg in PLATFORMS.items():
    print(f"Loading {plat}...")
    bla, blo, bss = read_bufr_qcd(f'{BUFR_BASE}/{cfg["subdir"]}/2020/01',
                                   DATE, cfg['bufr_prefix'])
    hla, hlo, hss = read_h121_qcd(f'{H121_BASE}/{cfg["subdir"]}/2020/01', DATE)
    bq = inbox(bla, blo)
    hq = inbox(hla, hlo)

    # clip to box
    bla, blo, bss = bla[bq], blo[bq], bss[bq]
    hla, hlo, hss = hla[hq], hlo[hq], hss[hq]

    # super-obs (Python reproduction of GEOSldas tile assignment + averaging)
    b_ij, b_so, b_so_lat, b_so_lon, b_cnt, b_inv = compute_super_obs(bla, blo, bss)
    h_ij, h_so, h_so_lat, h_so_lon, h_cnt, h_inv = compute_super_obs(hla, hlo, hss)

    plat_data[plat] = dict(
        bufr_lat=bla, bufr_lon=blo, bufr_ssm=bss,
        h121_lat=hla, h121_lon=hlo, h121_ssm=hss,
        b_ij=b_ij, b_so=b_so, b_so_lat=b_so_lat, b_so_lon=b_so_lon, b_cnt=b_cnt,
        h_ij=h_ij, h_so=h_so, h_so_lat=h_so_lat, h_so_lon=h_so_lon, h_cnt=h_cnt,
    )
    print(f'  Legacy: {len(bla)} obs → {len(b_ij)} tiles')
    print(f'  H121:   {len(hla)} obs → {len(h_ij)} tiles')



# ── M36 grid line drawing ─────────────────────────────────────────────────────
def draw_m36_grid(ax, lat0, lon0, lat1, lon1):
    margin = 0.6
    xs = [EASE2(lon0-margin, lat0-margin)[0], EASE2(lon1+margin, lat0-margin)[0],
          EASE2(lon0-margin, lat1+margin)[0], EASE2(lon1+margin, lat1+margin)[0]]
    ys = [EASE2(lon0-margin, lat0-margin)[1], EASE2(lon1+margin, lat0-margin)[1],
          EASE2(lon0-margin, lat1+margin)[1], EASE2(lon1+margin, lat1+margin)[1]]
    x_lo, x_hi = min(xs)-M36_CELL, max(xs)+M36_CELL
    y_lo, y_hi = min(ys)-M36_CELL, max(ys)+M36_CELL
    i_lo = int(np.floor((x_lo - X_ORIGIN) / M36_CELL))
    i_hi = int(np.ceil( (x_hi - X_ORIGIN) / M36_CELL))
    j_lo = int(np.floor((Y_ORIGIN - y_hi) / M36_CELL))
    j_hi = int(np.ceil( (Y_ORIGIN - y_lo) / M36_CELL))
    x_bounds = X_ORIGIN + np.arange(i_lo, i_hi+1) * M36_CELL
    y_bounds = Y_ORIGIN - np.arange(j_lo, j_hi+1) * M36_CELL
    x_samp = np.linspace(x_lo, x_hi, 50)
    y_samp = np.linspace(y_lo, y_hi, 50)
    kw = dict(color='#888', lw=0.6, transform=ccrs.PlateCarree(), zorder=3)
    for xb in x_bounds:
        lo, la = EASE2(np.full(50, xb), y_samp, inverse=True)
        ax.plot(lo, la, **kw)
    for yb in y_bounds:
        lo, la = EASE2(x_samp, np.full(50, yb), inverse=True)
        ax.plot(lo, la, **kw)


# ── Plot ──────────────────────────────────────────────────────────────────────
PLATNAMES = list(PLATFORMS.keys())
NCOLS = 5
pad   = 0.05
extent = [LON0-pad, LON1+pad, LAT0-pad, LAT1+pad]
VMIN, VMAX = 0, 80
CMAP = plt.cm.YlOrBr
NORM = mcolors.Normalize(VMIN, VMAX)

fig, axes = plt.subplots(
    len(PLATNAMES), NCOLS, figsize=(21, 4.2 * len(PLATNAMES)),
    subplot_kw={'projection': ccrs.PlateCarree()}
)


def setup_ax(ax, row, col):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#efefef')
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.5,
                   edgecolor='#bbb', facecolor='none')
    draw_m36_grid(ax, LAT0, LON0, LAT1, LON1)
    is_left   = (col == 0)
    is_bottom = (row == len(PLATNAMES) - 1)
    gl = ax.gridlines(draw_labels=True, linewidth=0, x_inline=False, y_inline=False)
    gl.top_labels    = False
    gl.bottom_labels = is_bottom
    gl.left_labels   = is_left
    gl.right_labels  = False


def draw_super_obs_tiles(ax, ij_arr, so_vals, cnt_arr, raw_lon, raw_lat,
                          label_fontsize=7):
    for k in range(len(ij_arr)):
        ci, cj = ij_arr[k]
        clons, clats = ease_ij_to_corners(ci, cj)
        color = CMAP(NORM(so_vals[k]))
        poly = mpatches.Polygon(
            np.column_stack([clons, clats]),
            facecolor=color, edgecolor='#444', linewidth=0.7, alpha=0.85,
            transform=ccrs.PlateCarree(), zorder=4
        )
        ax.add_patch(poly)
        ax.text(np.mean(clons), np.mean(clats), str(cnt_arr[k]),
                ha='center', va='center', fontsize=label_fontsize,
                fontweight='bold', color='k',
                transform=ccrs.PlateCarree(), zorder=6)
    ax.scatter(raw_lon, raw_lat, color='k', s=6, alpha=0.4,
               transform=ccrs.PlateCarree(), zorder=5)


col_titles = ['Both overlaid',
              'Legacy BUFR raw', 'Legacy super-obs',
              'H121 CDR raw',    'H121 super-obs']
for j, t in enumerate(col_titles):
    axes[0, j].set_title(t, fontsize=10, fontweight='bold', pad=4)

for i, plat in enumerate(PLATNAMES):
    d = plat_data[plat]
    bla, blo, bss = d['bufr_lat'], d['bufr_lon'], d['bufr_ssm']
    hla, hlo, hss = d['h121_lat'], d['h121_lon'], d['h121_ssm']

    bla_u, blo_u, bss_u, b_dup = deduplicate(bla, blo, bss)
    hla_u, hlo_u, hss_u, h_dup = deduplicate(hla, hlo, hss)
    B_BASE, H_BASE = 45, 15
    b_sz = B_BASE * b_dup
    h_sz = H_BASE * h_dup

    for j in range(NCOLS):
        setup_ax(axes[i, j], i, j)
    axes[i, 0].text(-0.18, 0.5, plat, transform=axes[i, 0].transAxes,
                    fontsize=10, fontweight='bold', va='center', rotation=90)

    # ── col 0: overlay ────────────────────────────────────────────────────
    axes[i, 0].scatter(hlo_u, hla_u, color='#e07b00', s=h_sz, marker='o',
                       alpha=0.75, transform=ccrs.PlateCarree(), zorder=4,
                       label=f'H121 ({len(hla)} obs, {len(hla_u)} nodes)')
    axes[i, 0].scatter(blo_u, bla_u, color='white',    s=b_sz*1.6, marker='s',
                       transform=ccrs.PlateCarree(), zorder=5)
    axes[i, 0].scatter(blo_u, bla_u, color='steelblue', s=b_sz, marker='s',
                       transform=ccrs.PlateCarree(), zorder=6,
                       label=f'Legacy ({len(bla)} obs, {len(bla_u)} nodes)')
    if i == 0:
        axes[i, 0].legend(loc='upper left', fontsize=7, framealpha=0.9)

    # ── col 1: Legacy raw obs ─────────────────────────────────────────────
    axes[i, 1].scatter(blo_u, bla_u, c=bss_u, cmap=CMAP,
                       vmin=VMIN, vmax=VMAX, s=b_sz, marker='s',
                       transform=ccrs.PlateCarree(), zorder=5)
    axes[i, 1].set_title(f'{len(bla)} obs → {len(d["b_ij"])} tiles',
                          fontsize=9, pad=2)

    # ── col 2: Legacy super-obs tiles ─────────────────────────────────────
    draw_super_obs_tiles(axes[i, 2],
                         d['b_ij'], d['b_so'], d['b_cnt'], blo, bla)
    axes[i, 2].set_title(f'{len(d["b_ij"])} Legacy super-obs', fontsize=9, pad=2)

    # ── col 3: H121 raw obs ───────────────────────────────────────────────
    axes[i, 3].scatter(hlo_u, hla_u, c=hss_u, cmap=CMAP,
                       vmin=VMIN, vmax=VMAX, s=h_sz, marker='o',
                       transform=ccrs.PlateCarree(), zorder=5)
    axes[i, 3].set_title(f'{len(hla)} obs → {len(d["h_ij"])} tiles',
                          fontsize=9, pad=2)

    # ── col 4: H121 super-obs tiles ───────────────────────────────────────
    draw_super_obs_tiles(axes[i, 4],
                         d['h_ij'], d['h_so'], d['h_cnt'], hlo, hla)
    axes[i, 4].set_title(f'{len(d["h_ij"])} H121 super-obs', fontsize=9, pad=2)


# ── Colourbar ─────────────────────────────────────────────────────────────────
fig.subplots_adjust(right=0.90, hspace=0.08, wspace=0.07)
cax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
sm  = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
cb  = fig.colorbar(sm, cax=cax)
cb.set_label('SSM (% saturation)', fontsize=10)

fig.suptitle(
    f'ASCAT obs & GEOSldas super-obs  |  2020-01-01  |  '
    f'{LAT0:.0f}–{LAT1:.0f}°N, {LON0:.0f}–{LON1:.0f}°E',
    fontsize=12, y=1.01
)

out_path = os.path.join(OUT_DIR, 'fig_zoom_comparison.png')
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out_path}")
