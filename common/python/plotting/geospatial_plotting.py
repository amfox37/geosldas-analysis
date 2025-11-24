# geospatial_plotting.py

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import get_cmap
import numpy as np
import xarray as xr
import re
import warnings
from pathlib import Path
from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

# === Predefined region bounds ===
REGION_BOUNDS = {
    'global':       [-180, 180, -60, 85],
    'north_america': [-125, -66.5, 24.5, 49.5],
    'australia':    [112, 154, -39, -10],
    'tibet':        [73, 105, 26, 40],
    'south_america': [0, 45, -36, -10],
    'conus':        [-125, -66.0, 24.0, 50.0],
    'cygnss':     [-140, 180, -40, 40]
}

# === Utilities ===
def format_number(num):
    if abs(num) < 0.01:
        return f"{num:.4f}"
    elif abs(num) < 1.0:
        return f"{num:.2f}"
    return f"{num:.3g}"

def summarize_array_stats(array, units='na', title=''):
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])
    if 'Relative $\\Delta$ StdDev' in title:
        return f"Mean = {mean:.1f}±{std:.1f} {units}"
    return f"Mean = {format_number(mean)}±{format_number(std)} {units}"

def load_ease_grid(ease_path=None):
    """
    Load EASE grid from bundled binaries. If ease_path is None, use the
    packaged files under common/python/plotting/ease_grids.
    """
    if ease_path is None:
        here = Path(__file__).resolve().parent
        ease_path = here / "ease_grids"
    ease_path = Path(ease_path)
    lats = np.fromfile(ease_path / 'EASE2_M36km.lats.964x406x1.double', dtype=np.float64).reshape((406,964))
    lons = np.fromfile(ease_path / 'EASE2_M36km.lons.964x406x1.double', dtype=np.float64).reshape((406,964))
    return lats, lons

def build_ease_grid_mapping(array, lats_row, lons_col):
    grid = np.full((lats_row.size, lons_col.size), np.nan, dtype=np.float64)
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if not np.isnan(array[i, 0]):
            grid[row, col] = array[i, 0]
    return grid

def colorbar_info(array):
    abmm = np.nanmax(np.abs(array[:, 0]))
    if np.nanmin(array[:, 0]) < 0:
        return -abmm, abmm, 'RdBu'
    return np.nanmin(array[:, 0]), np.nanmax(array[:, 0]), 'viridis'

def setup_map_features(ax):  # Land = grey, Ocean = white
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=False,
                      linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax.add_feature(cfeature.LAND, facecolor='lightgray')  # Draw land above
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=2)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=2)

def setup_colorbar(sc, ax, cmin, cmax, units):
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04)
    cbar.set_ticks(np.arange(cmin, cmax+1e-9, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(f'({units})', fontsize=12)

def save_plot(title, output_dir='./plots', fmt='png', dpi=600):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    clean_title = re.sub('[^0-9a-zA-Z]+', '_', title)
    savename = output_path / f"{clean_title}.{fmt}"
    print(f"Saving figure as {savename}")
    plt.savefig(savename, dpi=dpi, bbox_inches='tight', format=fmt)

# --- Discrete colormap + norm from bin edges ---

def make_discrete_cmap(edges, base_cmap="Spectral"):
    """
    edges: 1D array of bin edges (length K+1)
    returns (ListedColormap, BoundaryNorm)
    """
    edges = np.asarray(edges)
    n_bins = len(edges) - 1
    cmap_base = get_cmap(base_cmap)
    # sample one color per bin (evenly across the base colormap)
    colors = [cmap_base(i / max(n_bins-1, 1)) for i in range(n_bins)]
    cmap = ListedColormap(colors, name=f"{base_cmap}_discrete_{n_bins}")
    norm = BoundaryNorm(edges, ncolors=cmap.N, clip=False)
    return cmap, norm

def setup_colorbar_discrete(ax, norm, cmap, edges, label, orientation="horizontal"):
    """
    Discrete colorbar with \"≤ edge\" labels like your legend.
    """
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(mappable, ax=ax, orientation=orientation, pad=0.05, fraction=0.04)
    # put ticks at the upper edge of each bin (skip the first lower bound)
    tick_locs = edges[1:]
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels([f"≤ {v:g}" for v in tick_locs])
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)  # Set xticklabels to vertical
    cbar.set_label(label, fontsize=12)
    return cbar
    
# === Plotting Entry Point ===
def plot_region(array, region_bounds, ease_path='../test_data', 
                saveflag=False, meanflag=False, plot_title='regional_plot', 
                units='na', cmin=None, cmax=None, cmap=None, norm=None,
                output_dir='./plots', save_fmt='png', save_dpi=600, 
                star_lon=None, star_lat=None, overlay_points=None,
                discrete_edges=None, base_cmap="Spectral"):  # returns fig, ax in Jupyter
    lon_min, lon_max, lat_min, lat_max = region_bounds
    lats, lons = load_ease_grid(ease_path)
    lats_row, lons_col = lats[:,1], lons[1,:]
    grid = build_ease_grid_mapping(array, lats_row, lons_col) 

    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0

    textstr = summarize_array_stats(array, units, plot_title)

    # === Color handling ===
    if discrete_edges is not None:
        # Discrete Spectral with your bin edges
        ai_cmap, ai_norm = make_discrete_cmap(discrete_edges, base_cmap)
        cmap, norm = ai_cmap, ai_norm
        cmin = cmax = None  # ignored when using norm
    else:
        # Continuous fallback (original behavior)
        if cmin is None or cmax is None:
            cmin, cmax, cmap = colorbar_info(array)
        cmap = plt.get_cmap(cmap if cmap else ('RdBu_r' if cmin < 0 else 'viridis'), 20).copy()
        norm = None

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    setup_map_features(ax)

    sc = ax.pcolormesh(
        lons, lats, grid, transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm, vmin=None if norm is not None else cmin,
        vmax=None if norm is not None else cmax, shading="auto"
    )

    if overlay_points is not None:
        overlay_vals = overlay_points[:, 0]
        overlay_lons = overlay_points[:, 1]
        overlay_lats = overlay_points[:, 2]
        ax.scatter(
            overlay_lons, overlay_lats, c=overlay_vals,
            cmap=cmap, norm=norm,
            vmin=None if norm is not None else cmin,
            vmax=None if norm is not None else cmax,
            s=12, linewidths=0.2, edgecolor='black',
            transform=ccrs.PlateCarree(), zorder=3
        )

    if star_lon and star_lat:
        ax.plot(star_lon, star_lat, 'r*', markersize=10, transform=ccrs.PlateCarree())

    # Colorbar
    if discrete_edges is not None:
        setup_colorbar_discrete(sc.axes, norm, cmap, discrete_edges, f'({units})')
    else:
        setup_colorbar(sc, ax, cmin, cmax, units)

    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        save_plot(plot_title, output_dir=output_dir, fmt=save_fmt, dpi=save_dpi)

    return fig, ax

# === Example Usage ===
if __name__ == '__main__':
    # Example: load or generate dummy data for testing
    # array format: (value, lon, lat)
    array = np.array([
        [0.5, 100.0, 35.0],
        [0.6, 101.0, 36.0],
        [0.4, 102.0, 34.0]
    ])

    region = 'tibet'
    plot_region(array,
                region_bounds=REGION_BOUNDS[region],
                plot_title='Tibetan Plateau Test Plot',
                units='Test Units',
                meanflag=True,
                saveflag=False)
