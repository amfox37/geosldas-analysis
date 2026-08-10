#!/usr/bin/env python3

"""
Stage 2 (portable -- no Discover paths, no netCDF4): reconstruct the
footprint peat fraction for every candidate obs near the 6 boundary tiles,
using ONLY the CSV data files in this same directory (produced by
extract_footprint_data.py), and plot:

  (a) per-tile local map: nearby model (catchment) tiles coloured by
      porosity, PEATCLSM peat/mineral threshold marked, obs footprint
      centres overlaid coloured by kept/lost.
  (b) reconstructed footprint peat fraction vs kept/lost, with the 0.10
      QC threshold marked, across all 6 tiles combined.
  (c) one paired footprint map per target tile, comparing the nearest actual
      retained and rejected observations over the local porosity grid.
  (d) tile 3823 actual observation time series: all baseline candidates in
      grey, retained observations in green, and rejected observations in red.

This script hard-asserts that every candidate footprint is reconstructable,
that every reconstructed decision agrees with the OFA-derived outcome, and
that every observation lies on the expected side of the 0.10 threshold.

Run from anywhere:
  python3 plot_footprint_case_study.py
Requires: numpy, pandas, matplotlib (no netCDF4, no cluster paths).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data')
FIGURE_DIR = os.path.join(HERE, 'figures')
OUTPUT_DIR = os.path.join(HERE, 'output')
os.makedirs(FIGURE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

df_tiles = pd.read_csv(os.path.join(DATADIR, 'local_tile_geometry_porosity.csv'))
df_obs = pd.read_csv(os.path.join(DATADIR, 'obs_candidates.csv'))
const = pd.read_csv(os.path.join(DATADIR, 'constants.csv')).set_index('name')['value']

FOV_KM = const['FOV_km']
FAC_SEARCH = const['fac_search_FOV_km']
POROS_THRESHOLD = const['PEATCLSM_POROS_THRESHOLD']
MAX_PEAT_FRAC = const['ASCAT_max_peat_frac']
EARTH_RADIUS_KM = const['MAPL_RADIUS_km']

TARGET_TILES = sorted(df_tiles['target_tile'].unique())

# ----------------------------------------------------------------------------------
# footprint peat-fraction reconstruction (see extract_footprint_data.py docstring
# for the exact Fortran lines this replicates)
# ----------------------------------------------------------------------------------


def dist_km2deg(dist_km, lat_deg):
    """dist_km2deg() in clsm_ensupd_upd_routines.F90 line 6206."""
    dist_y_deg = dist_km * (180.0 / np.pi) / EARTH_RADIUS_KM
    dist_x_deg = dist_y_deg / np.cos(np.deg2rad(lat_deg))
    return dist_x_deg, dist_y_deg


def footprint_peat_fraction(obs_lon, obs_lat, tiles):
    """tiles: DataFrame with com_lon, com_lat, area, poros (one target tile's
    neighborhood). Returns (peatfrac, n_tiles_in_ellipse, weights_df)."""

    # search radius = fac_search_FOV_km * FOV, then rescaled ndst2 reduces to
    # normalizing distance by the FOV itself (not the 2x-expanded search radius)
    fov_x_deg, fov_y_deg = dist_km2deg(FOV_KM, obs_lat)
    r_x_deg, r_y_deg = dist_km2deg(FAC_SEARCH * FOV_KM, obs_lat)

    dlon = tiles['com_lon'].values - obs_lon
    dlat = tiles['com_lat'].values - obs_lat

    ndst2_search = (dlon / r_x_deg) ** 2 + (dlat / r_y_deg) ** 2
    in_ellipse = ndst2_search <= 1.0

    if not np.any(in_ellipse):
        return np.nan, 0, None

    ndst2_fov = (dlon[in_ellipse] / fov_x_deg) ** 2 + (dlat[in_ellipse] / fov_y_deg) ** 2
    gauss_w = np.exp(-0.5 * ndst2_fov)
    area_w = tiles['area'].values[in_ellipse]
    weight = gauss_w * area_w

    poros = tiles['poros'].values[in_ellipse]
    peat_ind = (poros >= POROS_THRESHOLD).astype(float)

    wsum = weight.sum()
    if wsum <= 0:
        return np.nan, int(in_ellipse.sum()), None

    peatfrac = float(np.sum(weight * peat_ind) / wsum)

    wdf = pd.DataFrame({
        'com_lon': tiles['com_lon'].values[in_ellipse],
        'com_lat': tiles['com_lat'].values[in_ellipse],
        'poros': poros,
        'peat_indicator': peat_ind,
        'weight': weight,
    })
    return peatfrac, int(in_ellipse.sum()), wdf


def centers_to_edges(centers):
    """Infer plotting-cell edges from sorted grid-cell centers."""
    centers = np.asarray(centers, dtype=float)
    if centers.size < 2:
        raise ValueError('At least two centers are required to infer grid edges')
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def nearest_outcome_pair(obs_t):
    """Return the geographically closest actual kept/lost observation pair."""
    kept = obs_t[obs_t['kept']].reset_index(drop=True)
    lost = obs_t[~obs_t['kept']].reset_index(drop=True)
    if kept.empty or lost.empty:
        raise ValueError('A paired map requires at least one kept and one lost observation')

    mean_lat = float(obs_t['lat'].mean())
    dx_km = (
        np.deg2rad(kept['lon'].to_numpy()[:, None] - lost['lon'].to_numpy()[None, :])
        * EARTH_RADIUS_KM * np.cos(np.deg2rad(mean_lat))
    )
    dy_km = (
        np.deg2rad(kept['lat'].to_numpy()[:, None] - lost['lat'].to_numpy()[None, :])
        * EARTH_RADIUS_KM
    )
    distance_km = np.hypot(dx_km, dy_km)
    kept_index, lost_index = np.unravel_index(np.argmin(distance_km), distance_km.shape)
    return kept.iloc[kept_index], lost.iloc[lost_index], float(distance_km[kept_index, lost_index])


def plot_footprint_panel(ax, obs_row, tiles_t, outcome):
    """Plot one actual observation footprint over the local porosity grid."""
    grid = tiles_t.drop_duplicates('tile_id').copy()
    lon_centers = np.sort(grid['com_lon'].unique())
    lat_centers = np.sort(grid['com_lat'].unique())
    lon_edges = centers_to_edges(lon_centers)
    lat_edges = centers_to_edges(lat_centers)
    poros_grid = (
        grid.pivot(index='com_lat', columns='com_lon', values='poros')
        .reindex(index=lat_centers, columns=lon_centers)
        .to_numpy()
    )

    mesh = ax.pcolormesh(
        lon_edges, lat_edges, poros_grid, cmap='YlGnBu', vmin=0.3, vmax=1.0,
        shading='flat', edgecolors='0.70', linewidth=0.55, zorder=0,
    )

    # Ring every PEATCLSM grid cell and the OFA owner tile explicitly.
    for row_index, lat_center in enumerate(lat_centers):
        for column_index, lon_center in enumerate(lon_centers):
            if np.isfinite(poros_grid[row_index, column_index]) and poros_grid[row_index, column_index] >= POROS_THRESHOLD:
                ax.add_patch(Rectangle(
                    (lon_edges[column_index], lat_edges[row_index]),
                    lon_edges[column_index + 1] - lon_edges[column_index],
                    lat_edges[row_index + 1] - lat_edges[row_index],
                    fill=False, edgecolor='#d62728', linewidth=1.8, zorder=2,
                ))

    owner = grid[grid['tile_id'] == int(obs_row['tilenum'])]
    if len(owner) != 1:
        raise ValueError(f"Expected one owner tile row for tilenum {int(obs_row['tilenum'])}, found {len(owner)}")
    owner_lon = float(owner['com_lon'].iloc[0])
    owner_lat = float(owner['com_lat'].iloc[0])
    owner_col = int(np.argmin(np.abs(lon_centers - owner_lon)))
    owner_row = int(np.argmin(np.abs(lat_centers - owner_lat)))
    ax.add_patch(Rectangle(
        (lon_edges[owner_col], lat_edges[owner_row]),
        lon_edges[owner_col + 1] - lon_edges[owner_col],
        lat_edges[owner_row + 1] - lat_edges[owner_row],
        fill=False, edgecolor='black', linewidth=2.6, zorder=3,
    ))
    ax.scatter(owner_lon, owner_lat, marker='s', s=48, facecolor='white', edgecolor='black', zorder=6)

    obs_lon = float(obs_row['lon'])
    obs_lat = float(obs_row['lat'])
    peatfrac, _, weights = footprint_peat_fraction(obs_lon, obs_lat, grid)
    if weights is None:
        raise ValueError(f'No contributing footprint tiles for {obs_row.timestamp}')

    weight_fraction = weights['weight'] / weights['weight'].sum()
    ax.scatter(
        weights['com_lon'], weights['com_lat'],
        s=80 + 260 * weight_fraction, facecolors='none', edgecolors='#cc00cc',
        linewidths=2.0, zorder=5,
    )
    for (_, cell), fraction in zip(weights.iterrows(), weight_fraction):
        ax.annotate(
            f'{100.0 * fraction:.0f}%', (cell['com_lon'], cell['com_lat']),
            xytext=(0, -18), textcoords='offset points', ha='center', va='top',
            fontsize=8, color='#7a007a', zorder=7,
        )

    fov_x_deg, fov_y_deg = dist_km2deg(FOV_KM, obs_lat)
    search_x_deg, search_y_deg = dist_km2deg(FAC_SEARCH * FOV_KM, obs_lat)
    outcome_color = '#17823b' if outcome == 'kept' else '#c51b1d'
    ax.add_patch(Ellipse(
        (obs_lon, obs_lat), 2.0 * fov_x_deg, 2.0 * fov_y_deg,
        fill=False, edgecolor=outcome_color, linewidth=2.2, linestyle='-', zorder=4,
    ))
    ax.add_patch(Ellipse(
        (obs_lon, obs_lat), 2.0 * search_x_deg, 2.0 * search_y_deg,
        fill=False, edgecolor=outcome_color, linewidth=1.8, linestyle='--', zorder=4,
    ))
    ax.scatter(
        obs_lon, obs_lat, marker='*', s=180, facecolor=outcome_color,
        edgecolor='black', linewidth=0.7, zorder=8,
    )

    platform = {8: 'Metop-A', 9: 'Metop-B'}.get(int(obs_row['species']), f"species {int(obs_row['species'])}")
    relation = '<' if peatfrac < MAX_PEAT_FRAC else '>='
    ax.set_title(
        f"{outcome.upper()}: {obs_row['timestamp']} {platform}\n"
        f"peat fraction={peatfrac:.3f} {relation} {MAX_PEAT_FRAC:.2f}"
    )
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(lon_edges[0], lon_edges[-1])
    ax.set_ylim(lat_edges[0], lat_edges[-1])
    ax.set_aspect(1.0 / np.cos(np.deg2rad(obs_lat)), adjustable='box')
    return mesh


# ----------------------------------------------------------------------------------
# compute peat fraction for every candidate obs
# ----------------------------------------------------------------------------------

peatfracs = np.full(len(df_obs), np.nan)
n_in_ellipse = np.zeros(len(df_obs), dtype=int)

for t in TARGET_TILES:
    tiles_t = df_tiles[df_tiles['target_tile'] == t]
    mask = df_obs['tilenum'] == t
    for i in np.nonzero(mask.values)[0]:
        pf, n, _ = footprint_peat_fraction(df_obs['lon'].iloc[i], df_obs['lat'].iloc[i], tiles_t)
        peatfracs[i] = pf
        n_in_ellipse[i] = n

df_obs['peatfrac'] = peatfracs
df_obs['n_tiles_in_ellipse'] = n_in_ellipse
df_obs['predicted_reject'] = df_obs['peatfrac'] >= MAX_PEAT_FRAC

out_csv = os.path.join(OUTPUT_DIR, 'footprint_case_study_results.csv')
df_obs.to_csv(out_csv, index=False)
print(f'wrote {out_csv}')

# ----------------------------------------------------------------------------------
# confusion matrix: predicted (peatfrac>=0.10 => reject) vs actual (kept/lost)
# ----------------------------------------------------------------------------------

valid = df_obs['peatfrac'].notna()
actual_reject = ~df_obs['kept']
pred_reject = df_obs['predicted_reject']

tp = int(np.sum(valid & actual_reject & pred_reject))     # rejected, predicted reject
tn = int(np.sum(valid & ~actual_reject & ~pred_reject))   # kept, predicted keep
fp = int(np.sum(valid & ~actual_reject & pred_reject))    # kept, predicted reject (disagreement)
fn = int(np.sum(valid & actual_reject & ~pred_reject))    # rejected, predicted keep (disagreement)

print(f'\nconfusion matrix (n={valid.sum()} candidates with reconstructable footprint):')
print(f'  predicted reject & actually lost (agree) : {tp}')
print(f'  predicted keep   & actually kept (agree) : {tn}')
print(f'  predicted reject & actually kept (DISAGREE): {fp}')
print(f'  predicted keep   & actually lost (DISAGREE): {fn}')
print(f'  agreement rate: {100*(tp+tn)/valid.sum():.2f}%')
if (~valid).sum() > 0:
    print(f'  ({(~valid).sum()} candidates had no tile in the search ellipse -- excluded)')

# ----------------------------------------------------------------------------------
# hard validation assertions -- see module docstring
# ----------------------------------------------------------------------------------

n_unreconstructable = int((~valid).sum())
if n_unreconstructable > 0:
    bad = df_obs.loc[~valid, ['timestamp', 'species', 'tilenum', 'lon', 'lat']]
    print(f'\nUNRECONSTRUCTABLE candidates (no tile in search ellipse):\n{bad.to_string()}')
assert n_unreconstructable == 0, (
    f'{n_unreconstructable}/{len(df_obs)} candidates have no reconstructable footprint '
    '(no tile fell inside the search ellipse) -- see printed rows above')
print(f'\n[ASSERT PASSED] all {len(df_obs)} candidates have a reconstructable footprint')

disagreements = df_obs.loc[valid & (pred_reject != actual_reject)]
if len(disagreements) > 0:
    print(
        '\nDISAGREEMENTS between reconstructed QC and actual outcome:\n'
        f'{disagreements[["timestamp", "species", "tilenum", "kept", "peatfrac"]].to_string()}'
    )
assert len(disagreements) == 0, (
    f'{len(disagreements)} candidates disagree between reconstructed QC and actual '
    'outcome -- see printed rows above')
print(
    '[ASSERT PASSED] reconstructed QC agrees with actual outcome for all '
    f'{int(valid.sum())} reconstructable candidates (100% agreement)'
)

wrong_side_kept = df_obs.loc[valid & df_obs['kept'] & (df_obs['peatfrac'] >= MAX_PEAT_FRAC)]
wrong_side_lost = df_obs.loc[valid & ~df_obs['kept'] & (df_obs['peatfrac'] < MAX_PEAT_FRAC)]
if len(wrong_side_kept) > 0 or len(wrong_side_lost) > 0:
    print(
        f'\nWRONG-SIDE-OF-THRESHOLD rows: kept-but->=threshold={len(wrong_side_kept)}, '
        f'lost-but-<threshold={len(wrong_side_lost)}'
    )
assert len(wrong_side_kept) == 0 and len(wrong_side_lost) == 0, (
    'some observations fall on the wrong side of the 0.10 threshold relative to their '
    'kept/lost outcome -- see printed rows above')
print(
    f'[ASSERT PASSED] no observation lies on the wrong side of the {MAX_PEAT_FRAC} '
    'threshold (all kept < threshold, all lost >= threshold)'
)

# ----------------------------------------------------------------------------------
# (a) per-tile local maps: porosity + obs footprints coloured by kept/lost
# ----------------------------------------------------------------------------------

fig, axes = plt.subplots(2, 3, figsize=(19, 12))
axes = axes.ravel()

for i, t in enumerate(TARGET_TILES):
    ax = axes[i]
    tiles_t = df_tiles[df_tiles['target_tile'] == t]
    obs_t = df_obs[df_obs['tilenum'] == t]

    sc = ax.scatter(tiles_t['com_lon'], tiles_t['com_lat'], c=tiles_t['poros'],
                     s=np.clip(tiles_t['area'] / tiles_t['area'].max() * 300, 20, 300),
                     cmap='YlGnBu', vmin=0.3, vmax=1.0, edgecolors='0.4', linewidths=0.4,
                     zorder=1)
    peat_mask = tiles_t['poros'] >= POROS_THRESHOLD
    ax.scatter(tiles_t.loc[peat_mask, 'com_lon'], tiles_t.loc[peat_mask, 'com_lat'],
               facecolors='none', edgecolors='red', s=60, linewidths=1.3,
               label=f'poros>={POROS_THRESHOLD} (peat)', zorder=2)

    lost = obs_t[~obs_t['kept']]
    kept = obs_t[obs_t['kept']]
    ax.scatter(lost['lon'], lost['lat'], marker='x', c='black', s=35,
               label=f'lost obs (n={len(lost)})', zorder=3)
    ax.scatter(kept['lon'], kept['lat'], marker='*', c='lime', s=140,
               edgecolors='black', linewidths=0.6, label=f'kept obs (n={len(kept)})', zorder=4)

    ax.set_title(f'tile {t}: local porosity + obs footprints')
    ax.set_xlabel('lon')
    ax.set_ylabel('lat')
    ax.legend(fontsize=7, loc='best')
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('POROS [m3/m3]')

plt.tight_layout()
map_out = os.path.join(FIGURE_DIR, 'footprint_case_study_maps.png')
fig.savefig(map_out, dpi=130)
plt.close(fig)
print(f'wrote {map_out}')

# ----------------------------------------------------------------------------------
# (b) footprint peat fraction vs kept/lost, all 6 tiles combined
# ----------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 6))
d = df_obs[valid]
kept_pf = d.loc[d['kept'], 'peatfrac']
lost_pf = d.loc[~d['kept'], 'peatfrac']

rng = np.random.default_rng(0)
ax.scatter(rng.normal(0, 0.02, len(lost_pf)), lost_pf, c='red', s=14, alpha=0.4,
           label=f'lost (n={len(lost_pf)})')
ax.scatter(rng.normal(1, 0.02, len(kept_pf)) + 1, kept_pf, c='green', s=30, alpha=0.8,
           label=f'kept (n={len(kept_pf)})')
ax.axhline(MAX_PEAT_FRAC, color='k', ls='--', lw=1,
           label=f'QC threshold ({MAX_PEAT_FRAC})')
ax.set_xticks([0, 2])
ax.set_xticklabels(['lost', 'kept'])
ax.set_xlim(-0.5, 2.5)
ax.set_ylabel('reconstructed footprint peat fraction')
ax.set_title('Footprint peat fraction vs actual retain/reject outcome\n'
              f'(6 boundary tiles combined, agreement={100*(tp+tn)/valid.sum():.1f}%)')
ax.legend()
ax.grid(True, ls='--', alpha=0.3)
plt.tight_layout()
scatter_out = os.path.join(FIGURE_DIR, 'footprint_case_study_threshold.png')
fig.savefig(scatter_out, dpi=130)
plt.close(fig)
print(f'wrote {scatter_out}')

# ----------------------------------------------------------------------------------
# (c) one paired footprint map per target tile: nearest actual kept/lost locations
# ----------------------------------------------------------------------------------

ellipse_dir = os.path.join(FIGURE_DIR, 'footprint_ellipse_maps')
os.makedirs(ellipse_dir, exist_ok=True)

legend_handles = [
    Line2D([0], [0], color='0.2', lw=2.2, ls='-', label=f'{FOV_KM:g} km Gaussian scale'),
    Line2D([0], [0], color='0.2', lw=1.8, ls='--', label=f'{FAC_SEARCH * FOV_KM:g} km contribution cutoff'),
    Line2D([0], [0], marker='*', color='none', markerfacecolor='0.35', markeredgecolor='black',
           markersize=12, label='OFA super-observation centre'),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='white', markeredgecolor='black',
           markersize=8, label='M36 owner tile centre'),
    Line2D([0], [0], marker='o', color='none', markerfacecolor='none', markeredgecolor='#cc00cc',
           markersize=10, label='contributing tile (label=weight)'),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='none', markeredgecolor='#d62728',
           markersize=10, label=f'PEATCLSM cell (POROS>={POROS_THRESHOLD:g})'),
]

for target_tile in TARGET_TILES:
    tiles_t = df_tiles[df_tiles['target_tile'] == target_tile]
    obs_t = df_obs[df_obs['tilenum'] == target_tile]
    kept_row, lost_row, separation_km = nearest_outcome_pair(obs_t)

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.0))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.84, bottom=0.27, wspace=0.22)
    kept_mesh = plot_footprint_panel(axes[0], kept_row, tiles_t, 'kept')
    plot_footprint_panel(axes[1], lost_row, tiles_t, 'lost')

    cbar_ax = fig.add_axes([0.18, 0.15, 0.64, 0.035])
    cbar = fig.colorbar(kept_mesh, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('POROS (m3 m-3)')
    fig.legend(
        handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, 0.02),
        ncol=3, fontsize=8, frameon=True,
    )
    fig.suptitle(
        f'Tile {target_tile}: nearest retained/rejected OFA footprint pair '
        f'(centres {separation_km:.1f} km apart)',
        fontsize=13, y=0.965,
    )

    ellipse_out = os.path.join(ellipse_dir, f'tile_{target_tile}_footprint_ellipses.png')
    fig.savefig(ellipse_out, dpi=150)
    plt.close(fig)
    print(f'wrote {ellipse_out}')

# ----------------------------------------------------------------------------------
# (d) tile 3823 actual observation time series
# ----------------------------------------------------------------------------------

TS_TILE = 3823
d_ts = df_obs[df_obs['tilenum'] == TS_TILE].copy()
d_ts['date'] = pd.to_datetime(d_ts['timestamp'], format='%Y%m%d_%H%M')
d_ts = d_ts.sort_values('date')
lost_ts = d_ts[~d_ts['kept']]
kept_ts = d_ts[d_ts['kept']]

fig, ax = plt.subplots(figsize=(13, 5.5))
ax.scatter(
    d_ts['date'], d_ts['obs'], c='0.75', s=22, zorder=1,
    label=f'all FOV12.5 candidates (n={len(d_ts)})',
)
ax.scatter(
    lost_ts['date'], lost_ts['obs'], c='red', marker='x', s=40, zorder=2,
    label=f'rejected by peat-QC (n={len(lost_ts)})',
)
ax.scatter(
    kept_ts['date'], kept_ts['obs'], c='green', marker='o', s=55,
    edgecolors='black', linewidths=0.6, zorder=3,
    label=f'retained by peat-QC (n={len(kept_ts)})',
)
ax.set_ylabel('ASCAT H SAF soil moisture observation (m3 m-3)')
ax.set_xlabel('Date (2015)')
ax.set_title(
    f'Tile {TS_TILE}: actual OFA observation time series, species 8 and 9 combined\n'
    f'{len(d_ts)} FOV12.5 candidates; {len(kept_ts)} retained under peat-QC '
    f'({100 * len(kept_ts) / len(d_ts):.1f}% survival)'
)
ax.legend(loc='upper right')
ax.grid(True, ls='--', alpha=0.3)
fig.autofmt_xdate()
plt.tight_layout()
ts_out = os.path.join(FIGURE_DIR, f'footprint_case_study_timeseries_{TS_TILE}.png')
fig.savefig(ts_out, dpi=130)
plt.close(fig)
print(f'wrote {ts_out}')

print('\nDONE', flush=True)
