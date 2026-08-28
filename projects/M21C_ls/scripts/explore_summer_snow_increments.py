#!/usr/bin/env python3
"""Where does snow DA add snow water during summer?

Figure 15's monthly climatology shows non-trivial positive snow-DA input in
Jun-Sep over a mask that is supposed to exclude perennial snow. This maps the
summer input and breaks it down by latitude, elevation, and longitude sector.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np, xarray as xr, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import TwoSlopeNorm
from pyproj import CRS, Transformer

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "common/python/io"))
from read_GEOSldas import read_tilecoord

DATA = Path("/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2")
# EASEv2 M36 global grid. Verified against the tilecoord cell centres to 1e-4 deg.
SCALE, NCOL, NROW = 36032.220840584, 964, 406

d = xr.open_dataset(ROOT / "output/monthly_synthesis_diagnostics/water_year_snow_da_budget"
                    / "water_year_tile_budgets.nc")
tid = d.tile.values.astype(int)
lat, lon, area = d.lat.values, d.lon.values, d.area.values
tc = read_tilecoord(str(DATA / "LS_OLv8_M36.ldas_tilecoord.bin"))
elev = np.asarray(tc["elev"])[tid]
icol = np.asarray(tc["i_indg"])[tid].astype(int)
irow = np.asarray(tc["j_indg"])[tid].astype(int)
sn = d["snow_net_monthly"].values                      # (wy, month, tile), Oct..Sep


def grid_edges() -> tuple[np.ndarray, np.ndarray]:
    """Cell-edge longitudes and latitudes of the EASEv2 M36 global grid.

    The grid is a cylindrical equal-area projection, so columns carry constant
    longitude and rows constant latitude: the mesh is regular in lon/lat (though
    unevenly spaced in lat) and can be handed straight to pcolormesh.

    Equal area means cell shape distorts away from the 30N standard parallel:
    cells run about 41 x 31 km at the equator but 4.5 x 285 km at 83.8N. The
    radial slivers this produces at high latitude in a polar projection are the
    true footprints, not a plotting artefact.
    """
    transformer = Transformer.from_crs(
        CRS.from_proj4("+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 "
                       "+ellps=WGS84 +datum=WGS84 +units=m"),
        CRS.from_epsg(4326), always_xy=True)
    x = (np.arange(NCOL + 1) - 0.5 - (NCOL - 1) / 2.0) * SCALE
    y = ((NROW - 1) / 2.0 - (np.arange(NROW + 1) - 0.5)) * SCALE
    lon_edges, _ = transformer.transform(x, np.zeros_like(x))
    _, lat_edges = transformer.transform(np.zeros_like(y), y)
    return lon_edges, lat_edges


def to_grid(values: np.ndarray) -> np.ma.MaskedArray:
    """Scatter per-tile values onto the native grid; off-domain cells stay masked."""
    field = np.full((NROW, NCOL), np.nan, dtype="float64")
    field[irow, icol] = values
    return np.ma.masked_invalid(field)


LON_EDGES, LAT_EDGES = grid_edges()
TOP = int(np.flatnonzero(LAT_EDGES > 18.0).max())      # keep the mapped hemisphere only

panels = [("Jun", [8]), ("Jul", [9]), ("Aug", [10]), ("Sep", [11])]
fields = {name: np.nansum(sn[:, idx, :], axis=1).mean(axis=0) for name, idx in panels}

fig = plt.figure(figsize=(13.5, 12.2))
proj = ccrs.NorthPolarStereo(central_longitude=0)
vmax = 12.0
theta = np.linspace(0, 2 * np.pi, 200)
boundary = mpath.Path(np.column_stack([0.5 + 0.5 * np.cos(theta), 0.5 + 0.5 * np.sin(theta)]))
for k, (name, _) in enumerate(panels):
    ax = fig.add_subplot(2, 2, k + 1, projection=proj)
    ax.set_extent([-180, 180, 22, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.93", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
    v = fields[name]
    mesh = ax.pcolormesh(LON_EDGES, LAT_EDGES[:TOP + 1], to_grid(v)[:TOP],
                         cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax),
                         transform=ccrs.PlateCarree(), shading="flat", rasterized=True)
    ax.coastlines(linewidth=0.35, color="0.35")
    dom = np.nansum(v * area) / np.nansum(area)
    ax.set_title(f"({'abcd'[k]}) {name}   domain mean {dom:.2f} kg m$^{{-2}}$",
                 fontsize=11.5, fontweight="bold")
    ax.set_boundary(boundary, transform=ax.transAxes)
cb = fig.colorbar(mesh, ax=fig.axes, orientation="horizontal", fraction=0.035,
                  pad=0.03, extend="both", shrink=0.55)
cb.set_label("Mean monthly snow-DA input, WY2001–WY2006 (kg m$^{-2}$)", fontsize=10.5)
fig.suptitle("Summer and early-autumn snow-DA input over the seasonal-snow domain",
             fontsize=13.5, y=0.95)
out = ROOT / "output/monthly_synthesis_diagnostics/figures/summer_snow_da_input_maps.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out.relative_to(ROOT))

# ---- breakdowns ----
jja = np.nansum(sn[:, [8, 9, 10], :], axis=1).mean(axis=0)
tot_pos = np.nansum(np.where(jja > 0, jja, 0) * area)
print("\nElevation breakdown of JJA positive input")
print(f"  {'elev band (m)':>16s} {'share':>8s} {'n tiles':>9s} {'mean lat':>9s}")
for lo, hi in ((0, 250), (250, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 9000)):
    m = (elev >= lo) & (elev < hi)
    if not m.any():
        continue
    share = 100 * np.nansum(np.where(jja[m] > 0, jja[m], 0) * area[m]) / tot_pos
    print(f"  {lo:6d}-{hi:6d} {share:7.1f}% {m.sum():9d} {np.average(lat[m], weights=area[m]):8.1f}")

print("\nLongitude sector, 50-75N only")
sectors = [("N America", -170, -60), ("Greenland/Atl", -60, 0),
           ("Europe/W Russia", 0, 60), ("C Siberia", 60, 120), ("E Siberia", 120, 180)]
band = (lat >= 50) & (lat < 75)
tot_band = np.nansum(np.where(jja[band] > 0, jja[band], 0) * area[band])
for name, lo, hi in sectors:
    m = band & (lon >= lo) & (lon < hi)
    share = 100 * np.nansum(np.where(jja[m] > 0, jja[m], 0) * area[m]) / tot_band
    print(f"  {name:16s} {share:5.1f}%  ({m.sum():5d} tiles)")


# ---- 90E discontinuity diagnostic ----
# A sharp meridional step in the analysis increments sits on the Tibetan Plateau
# at exactly 90.00E, a grid column edge. MODIS SCF is read per assimilation window
# in the longitude band where local solar time matches the overpass, so with
# 3-hourly windows those band edges fall on multiples of 45 degrees.
band = (lat >= 28) & (lat <= 38)
jja = np.nansum(sn[:, [8, 9, 10], :], axis=1).mean(axis=0)
ol_swe = d["ol_snow_mass_monthly"].values[:, [8, 9, 10], :].mean(axis=(0, 1))
cols = np.arange(700, 750)
col_lon = -180.0 + cols * (360.0 / NCOL) + 0.5 * (360.0 / NCOL)
prof = {n: np.array([np.nanmean(v[band & (icol == c)]) if (band & (icol == c)).sum() else np.nan
                     for c in cols]) for n, v in (("incr", jja), ("ol", ol_swe))}

fig, ax = plt.subplots(1, 2, figsize=(14, 5.2),
                       gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.28))
axm = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree()); ax[0].remove()
axm.set_extent([72, 108, 25, 42], ccrs.PlateCarree())
axm.add_feature(cfeature.LAND, facecolor="0.94")
axm.coastlines(lw=0.4); axm.add_feature(cfeature.BORDERS, lw=0.3, edgecolor="0.6")
mesh = axm.pcolormesh(LON_EDGES, LAT_EDGES[:TOP + 1], to_grid(jja)[:TOP], cmap="RdBu_r",
                      norm=TwoSlopeNorm(vcenter=0, vmin=-60, vmax=60),
                      transform=ccrs.PlateCarree())
axm.plot([90, 90], [25, 42], color="k", lw=1.2, ls="--", transform=ccrs.PlateCarree())
axm.text(90.5, 41.0, "90°E", transform=ccrs.PlateCarree(), fontsize=10, fontweight="bold")
gl = axm.gridlines(draw_labels=True, lw=0.3, color="0.75"); gl.top_labels = gl.right_labels = False
axm.set_title("(a) JJA snow-DA input, Tibetan Plateau", loc="left", fontweight="bold")
fig.colorbar(mesh, ax=axm, orientation="horizontal", pad=0.08, fraction=0.05,
             extend="both", label="kg m$^{-2}$")

ax[1].plot(col_lon, prof["incr"], "-o", ms=3.5, color="#c0392b", label="JJA snow-DA input")
ax[1].axvline(90.0, color="k", ls="--", lw=1.2)
ax[1].set_xlabel("longitude (°E)"); ax[1].set_ylabel("JJA snow-DA input (kg m$^{-2}$)", color="#c0392b")
ax[1].tick_params(axis="y", labelcolor="#c0392b")
ax[1].grid(color="0.9", lw=0.7); ax[1].set_axisbelow(True)
twin = ax[1].twinx()
twin.plot(col_lon, prof["ol"], "-s", ms=3.5, color="#2b6cb0", label="OL snow mass")
twin.set_ylabel("OL JJA snow mass (kg m$^{-2}$)", color="#2b6cb0")
twin.tick_params(axis="y", labelcolor="#2b6cb0")
ax[1].set_title("(b) Zonal profile, 28–38°N", loc="left", fontweight="bold")
ax[1].text(0.03, 0.30, "increment steps at the meridian;\nmodel background does not",
           transform=ax[1].transAxes, va="top", fontsize=9.5,
           bbox=dict(fc="white", ec="0.8", boxstyle="round,pad=0.4"))
out = ROOT / "output/monthly_synthesis_diagnostics/figures/snow_da_90E_discontinuity.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out.relative_to(ROOT))
