"""
Scatter plot: Legacy BUFR super-obs vs H121 CDR super-obs on matched M36 tiles.
One panel per Metop platform + one combined panel.
"""

import numpy as np
import glob, os
from datetime import datetime

import eccodes as ecc
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pyproj import Proj

EASE2     = Proj('+proj=cea +lat_ts=30 +datum=WGS84 +lon_0=0 +units=m')
M36_CELL  = 36032.220840584
X_ORIGIN, Y_ORIGIN = EASE2(-180, 90)

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

COLORS = {'Metop-A': '#1f77b4', 'Metop-B': '#ff7f0e', 'Metop-C': '#2ca02c'}


def latlon_to_ease_ij(lat, lon):
    x, y = EASE2(np.asarray(lon, float), np.asarray(lat, float))
    i = np.floor((x - X_ORIGIN) / M36_CELL).astype(int)
    j = np.floor((Y_ORIGIN - y) / M36_CELL).astype(int)
    return i, j


def inbox(lat, lon):
    return (lat >= LAT0) & (lat <= LAT1) & (lon >= LON0) & (lon <= LON1)


def compute_super_obs(lat, lon, ssm):
    i_arr, j_arr = latlon_to_ease_ij(lat, lon)
    ij = np.stack([i_arr, j_arr], axis=1)
    unique_ij, inv = np.unique(ij, axis=0, return_inverse=True)
    n = len(unique_ij)
    so_ssm = np.array([ssm[inv == k].mean() for k in range(n)])
    so_cnt = np.array([(inv == k).sum()      for k in range(n)])
    return unique_ij, so_ssm, so_cnt


def read_bufr_qcd(bufr_dir, date, prefix):
    lats, lons, ssms = [], [], []
    for f in sorted(glob.glob(os.path.join(bufr_dir,
                              f'{prefix}{date.strftime("%Y%m%d")}*.bfr'))):
        with open(f, 'rb') as fh:
            bufr = ecc.codes_bufr_new_from_file(fh)
            while bufr is not None:
                ecc.codes_set(bufr, 'unpack', 1)
                try:
                    KEYS = {'ssm': 'surfaceSoilMoisture',
                            'smpf': 'soilMoistureProcessingFlag',
                            'smcf': 'soilMoistureCorrectionFlag',
                            'tpcx': 'topographicComplexity',
                            'iwfr': 'inundationAndWetlandFraction',
                            'lat': 'latitude', 'lon': 'longitude'}
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
    for f in sorted(glob.glob(os.path.join(h121_dir,
                              f'*_{date.strftime("%Y%m%d")}*.nc'))):
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
             ((cf == 0) | (cf == 4)) &
             (wf >= 0) & (wf <= 10) & (tc >= 0) & (tc <= 10))
        lats.extend(lat[q].tolist())
        lons.extend(lon[q].tolist())
        ssms.extend(ssm[q].tolist())
    return np.array(lats), np.array(lons), np.array(ssms)


# ── Load data and find matched tiles per platform ─────────────────────────────
matched = {}   # plat → {'legacy': array, 'h121': array, 'n_leg': array, 'n_h121': array}

for plat, cfg in PLATFORMS.items():
    print(f"Loading {plat}...")
    bla, blo, bss = read_bufr_qcd(f'{BUFR_BASE}/{cfg["subdir"]}/2020/01',
                                   DATE, cfg['bufr_prefix'])
    hla, hlo, hss = read_h121_qcd(f'{H121_BASE}/{cfg["subdir"]}/2020/01', DATE)

    bq = inbox(bla, blo);  bla, blo, bss = bla[bq], blo[bq], bss[bq]
    hq = inbox(hla, hlo);  hla, hlo, hss = hla[hq], hlo[hq], hss[hq]

    b_ij, b_so, b_cnt = compute_super_obs(bla, blo, bss)
    h_ij, h_so, h_cnt = compute_super_obs(hla, hlo, hss)

    # index h121 tiles by (i,j) for fast lookup
    h_lookup = {(h_ij[k, 0], h_ij[k, 1]): (h_so[k], h_cnt[k])
                for k in range(len(h_ij))}

    leg_vals, h_vals, leg_n, h_n = [], [], [], []
    for k in range(len(b_ij)):
        key = (b_ij[k, 0], b_ij[k, 1])
        if key in h_lookup:
            leg_vals.append(b_so[k])
            h_vals.append(h_lookup[key][0])
            leg_n.append(b_cnt[k])
            h_n.append(h_lookup[key][1])

    matched[plat] = dict(
        legacy=np.array(leg_vals),
        h121=np.array(h_vals),
        n_leg=np.array(leg_n),
        n_h121=np.array(h_n),
    )
    print(f"  matched tiles: {len(leg_vals)}")


# ── Plot ───────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True, sharex=True)

VMIN, VMAX = 30, 100
ref = np.array([VMIN, VMAX])


def annotate(ax, x, y, color, label):
    if len(x) < 2:
        return
    r, _, _, _, _ = stats.linregress(x, y)
    bias = np.mean(y - x)
    rmse = np.sqrt(np.mean((y - x) ** 2))
    r_val = np.corrcoef(x, y)[0, 1]
    ax.scatter(x, y, color=color, s=30, alpha=0.75, label=label, zorder=3)
    txt = f'n={len(x)}\nbias={bias:+.1f}%\nRMSE={rmse:.1f}%\nr={r_val:.2f}'
    ax.text(0.04, 0.97, txt, transform=ax.transAxes,
            fontsize=8, va='top', color=color,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=2))


all_leg, all_h121 = [], []

for j, (plat, d) in enumerate(matched.items()):
    ax = axes[j]
    ax.plot(ref, ref, 'k--', lw=0.8, zorder=1)
    annotate(ax, d['legacy'], d['h121'], COLORS[plat], plat)
    ax.set_xlim(VMIN, VMAX)
    ax.set_ylim(VMIN, VMAX)
    ax.set_title(plat, fontsize=11)
    ax.set_xlabel('Legacy BUFR super-obs (% sat)', fontsize=9)
    if j == 0:
        ax.set_ylabel('H121 CDR super-obs (% sat)', fontsize=9)
    ax.set_aspect('equal')
    all_leg.extend(d['legacy'].tolist())
    all_h121.extend(d['h121'].tolist())

# combined panel
ax = axes[3]
ax.plot(ref, ref, 'k--', lw=0.8, zorder=1)
for plat, d in matched.items():
    ax.scatter(d['legacy'], d['h121'], color=COLORS[plat],
               s=30, alpha=0.7, label=plat, zorder=3)
all_leg = np.array(all_leg)
all_h121 = np.array(all_h121)
if len(all_leg) >= 2:
    bias = np.mean(all_h121 - all_leg)
    rmse = np.sqrt(np.mean((all_h121 - all_leg) ** 2))
    r_val = np.corrcoef(all_leg, all_h121)[0, 1]
    txt = f'n={len(all_leg)}\nbias={bias:+.1f}%\nRMSE={rmse:.1f}%\nr={r_val:.2f}'
    ax.text(0.04, 0.97, txt, transform=ax.transAxes,
            fontsize=8, va='top', color='k',
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=2))
ax.set_title('All platforms', fontsize=11)
ax.set_xlabel('Legacy BUFR super-obs (% sat)', fontsize=9)
ax.legend(fontsize=8, loc='lower right')
ax.set_aspect('equal')

fig.suptitle(
    f'M36 tile super-obs: Legacy BUFR vs H121 CDR  |  2020-01-01  |  '
    f'{LAT0:.0f}–{LAT1:.0f}°N, {LON0:.0f}–{LON1:.0f}°E',
    fontsize=11, y=1.01
)
fig.tight_layout()

out_path = os.path.join(OUT_DIR, 'fig_tile_scatter.png')
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out_path}")
