#!/usr/bin/env python3
"""Make the draft figures for the CYGNSS L1 paper outline (paper/paper_outline.md).

Writes PNGs to paper/figures/. Those PNGs are tracked in git so that the outline renders on GitHub.

Inputs (all local, gitignored under output/):
  output/cygl1_fixedop_2020_2022_stats/         O-F stats bundle (from cygl1_fixedop_2020_2022_stats.tar.gz)
  output/ismn_fixedop_20200101_20221231/        ISMN station skill + per-tile L1 counts
                                                (from cygl1_ismn_fixedop_2020_2022_inputs.tar.gz)
  output/operator_story_20201116_0000z/         schematic, from notebooks/cygl1_operator_story_figure.ipynb
Figs. 4 and 5 are drawn from tables in docs/cygl1_obs_quality_vs_qc_report.md and
docs/cygl1_obs_error_correlation_report.md, because their inputs are on Discover.

The primary L1 arm is DA_L1_full_xc015_coh040216_fixedop: coherency filter 0.40-2.16, errstd 2.75 dB, xcorr 0.15.

Run with: /Users/amfox/mamba/envs/regrid/bin/python paper/make_paper_figures.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyBboxPatch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from PIL import Image

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "paper" / "figures"
STATS = PROJ / "output" / "cygl1_fixedop_2020_2022_stats"
ISMN = PROJ / "output" / "ismn_fixedop_20200101_20221231"
STORY = PROJ / "output" / "operator_story_20201116_0000z" / "cygl1_operator_story_20201116_0000z_samp163383_ch3_err275.png"

DPI = 150
DOMAIN = (-118.0, -106.0, 29.0, 40.0)
L1_AREA_MIN_OBS = 100

# Arms: tag in the stats bundle, label, colour (Okabe-Ito).
ARMS = {
    "full_xc015_coh040216_fixedop": ("L1", "#D55E00"),
    "l3_fixedop": ("L3", "#0072B2"),
    "smap_fixedop": ("SMAP", "#009E73"),
}
ISMN_RUNS = {"L1coh": ("L1", "#D55E00"), "L3": ("L3", "#0072B2"), "SMAP": ("SMAP", "#009E73")}
# Spring degradation windows (O-F report §3), inclusive month ranges.
SPRING = [("2020-05", "2020-06"), ("2021-04", "2021-06"), ("2022-02", "2022-05")]
GROUPS = {"SMOS Tb": [0, 1, 2, 3], "SMAP Tb": [4, 5, 6, 7], "ASCAT SM": [8, 9, 10], "CYGNSS L3 SM": [11], "CYGNSS L1": [12]}

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "figure.titlesize": 11, "axes.spines.top": False, "axes.spines.right": False,
})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path.relative_to(PROJ)}")


def map_axes(fig, pos):
    ax = fig.add_subplot(*pos, projection=ccrs.PlateCarree())
    ax.set_extent(DOMAIN, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor="0.35")
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.6, edgecolor="0.25")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor="0.25")
    return ax


def _edges(centers):
    c = np.unique(np.round(np.asarray(centers, dtype=float), 4))
    mid = 0.5 * (c[:-1] + c[1:])
    return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def tile_squares(ax, lon, lat, values, color=None, **kw):
    """Draw each M36 tile as its EASE cell rectangle (edges at midpoints between tile centres)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    # Edges from the full 909-tile grid, so any subset of tiles gets the same cells.
    all_tiles = pd.read_csv(STATS / "tile_coords.csv")
    xc, xe = _edges(all_tiles.com_lon)
    yc, ye = _edges(all_tiles.com_lat)
    ix = np.searchsorted(xc, np.round(lon, 4))
    iy = np.searchsorted(yc, np.round(lat, 4))
    verts = [[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])] for i, j in zip(ix, iy)]
    kw.pop("s", None)
    vmin, vmax = kw.pop("vmin", None), kw.pop("vmax", None)
    pc = PolyCollection(verts, transform=ccrs.PlateCarree(), linewidths=0, **kw)
    if color is not None:
        pc.set_facecolor(color)
    else:
        pc.set_array(np.ma.masked_invalid(np.asarray(values, float)))
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def shade_spring(ax):
    for a, b in SPRING:
        ax.axvspan(pd.Timestamp(a + "-01"), pd.Timestamp(b + "-01") + pd.offsets.MonthEnd(1), color="0.9", zorder=0, lw=0)


def load_tiles():
    tiles = pd.read_csv(STATS / "tile_coords.csv")
    counts = pd.read_csv(ISMN / "tile_l1_obs_count_2020_2022.csv")
    tiles = tiles.merge(counts, left_on="tile_axis_index", right_on="tile_index")
    return tiles


# ---------------------------------------------------------------------------------------------
# Fig. 1: domain, L1 coverage, in-situ stations
# ---------------------------------------------------------------------------------------------
def fig01_domain():
    tiles = load_tiles()
    st = pd.read_csv(ISMN / "ismn_skill_stations.csv")
    st = st[(st.run == "OL") & (st.domain == "surface")]
    st = st.merge(tiles[["tile_axis_index", "n_l1_obs"]], left_on="tile_index", right_on="tile_axis_index")
    scored = st[st.n_l1_obs >= L1_AREA_MIN_OBS]
    other = st[st.n_l1_obs < L1_AREA_MIN_OBS]

    fig = plt.figure(figsize=(6.4, 5.6))
    ax = map_axes(fig, (1, 1, 1))
    has = tiles.n_l1_obs > 0
    tile_squares(ax, tiles.com_lon[~has], tiles.com_lat[~has], None, color="0.88", s=16)
    sc = tile_squares(ax, tiles.com_lon[has], tiles.com_lat[has], np.log10(tiles.n_l1_obs[has]), cmap="viridis",
                      vmin=1, vmax=np.log10(tiles.n_l1_obs.max()), s=16)
    ax.scatter(other.lon, other.lat, marker="o", s=10, facecolor="none", edgecolor="0.3", linewidths=0.6,
               transform=ccrs.PlateCarree(), zorder=5)
    markers = {"SCAN": "^", "USCRN": "D", "SNOTEL": "o", "SOILSCAPE": "*"}
    for net, m in markers.items():
        s = scored[scored.network == net]
        ax.scatter(s.lon, s.lat, marker=m, s=38 if m == "*" else 22, facecolor="white", edgecolor="k", linewidths=0.8,
                   transform=ccrs.PlateCarree(), zorder=6, label=f"{net} ({len(s)})")
    ax.scatter([], [], marker="o", s=10, facecolor="none", edgecolor="0.3", label=f"not scored ({len(other)})")
    ax.legend(loc="lower left", frameon=True, framealpha=0.9, title="ISMN surface stations", title_fontsize=8)
    cb = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label("CYGNSS L1 obs per tile, 2020–2022 (log$_{10}$)")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.6", alpha=0.6, xlocs=range(-118, -105, 3),
                      ylocs=range(29, 41, 2))
    gl.top_labels = gl.right_labels = False
    n_any, n_area = int(has.sum()), int((tiles.n_l1_obs >= L1_AREA_MIN_OBS).sum())
    ax.set_title(f"Domain (909 M36 tiles): {n_any} tiles with L1 obs, {n_area} with ≥{L1_AREA_MIN_OBS} (L1 area)")
    save(fig, "fig01_domain_coverage.png")


# ---------------------------------------------------------------------------------------------
# Fig. 2: operator schematic (made by the notebook; downscaled here for the document)
# ---------------------------------------------------------------------------------------------
def fig02_schematic():
    im = Image.open(STORY).convert("RGB")
    w = 1600
    im = im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)
    OUT.mkdir(parents=True, exist_ok=True)
    im.save(OUT / "fig02_operator_schematic.png", optimize=True)
    print("wrote paper/figures/fig02_operator_schematic.png")


# ---------------------------------------------------------------------------------------------
# Fig. 3: processing chain
# ---------------------------------------------------------------------------------------------
def fig03_pipeline():
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.6)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec="0.3", lw=0.8))
        ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="top", fontsize=9, fontweight="bold")
        ax.text(x + w / 2, y + h - 0.42, body, ha="center", va="top", fontsize=7.5, linespacing=1.35)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="0.25", lw=1.0))

    off, on = "#fdf0e6", "#e8f1f8"
    ax.text(0.05, 3.5, "Offline preprocessing (once per observation)", fontsize=9, color="#9a3b00", va="top")
    ax.text(0.05, 1.72, "GEOSldas EnKF (every 3-h cycle, every ensemble member)", fontsize=9, color="#004a80", va="top")
    y1, h = 1.95, 1.35
    box(0.05, y1, 2.2, h, "CYGNSS L1 v3.2", "DDMs over Arizona + 200 km\nland QC flags\ncoherency ratio 0.40–2.16", off)
    box(2.55, y1, 2.4, h, "IGOT forward model", "Copernicus DEM (~30 m)\nWAF kernel K$_p$, pixel factor σ$_p$\n3×5 crop around the peak", off)
    box(5.25, y1, 2.2, h, "Tile coefficients", "C$_t$ = s Σ$_{p∈t}$ K$_p$σ$_p$A$_t$\nvegetation A$_t$ from τ$_t$\nsparse, per M36 tile", off)
    box(7.75, y1, 2.2, h, "Daily product", "best obs per tile\nper 3-h window\nNetCDF (C$_t$, IG/JG, y)", off)
    y2 = 0.05
    box(7.75, y2, 2.2, h, "Reader", "obs in the 3-h window\nowner tile + time stamp\nz-score scaling", on)
    box(5.25, y2, 2.2, h, "Operator", "R$_t$ from SFMC, clay,\nporosity, incidence\nH = 10 log$_{10}$ Σ C$_t$R$_t$", on)
    box(2.55, y2, 2.4, h, "EnKF update", "errstd 2.75 dB\nGaussian R, xcorr 0.15°\nlocalization 1.25°", on)
    box(0.05, y2, 2.2, h, "Catchment state", "surface / root-zone\nexcess increments\n(SRFEXC, RZEXC)", on)
    for x0, x1 in [(2.25, 2.55), (4.95, 5.25), (7.45, 7.75)]:
        arrow(x0, y1 + h / 2, x1, y1 + h / 2)
    arrow(8.85, y1, 8.85, y2 + h)
    for x0, x1 in [(7.75, 7.45), (5.25, 4.95), (2.55, 2.25)]:
        arrow(x0, y2 + h / 2, x1, y2 + h / 2)
    save(fig, "fig03_processing_chain.png")


# ---------------------------------------------------------------------------------------------
# Fig. 4: L1 quality vs coherency ratio (QC report decile table, Jan-Jun 2020, fixed build)
# ---------------------------------------------------------------------------------------------
def fig04_coherency():
    edges = [0.07, 0.40, 0.53, 0.66, 0.81, 0.99, 1.22, 1.48, 1.78, 2.16, 10.17]
    inn_mean = [-0.40, -0.83, -0.78, -0.73, -0.64, -0.48, -0.32, -0.08, 0.36, 1.35]
    inn_std = [3.96, 2.81, 2.45, 2.51, 2.27, 2.16, 2.20, 2.24, 2.43, 3.26]
    r_desroz = [15.43, 8.24, 6.29, 6.50, 5.22, 4.61, 4.64, 4.75, 5.70, 11.74]
    r_assumed = [9.34, 8.78, 8.30, 8.26, 7.62, 7.33, 7.00, 6.59, 6.07, 5.88]
    r_smap = [-0.158, -0.248, -0.277, -0.255, -0.252, -0.267, -0.226, -0.209, -0.172, -0.088]
    r_l3 = [0.201, 0.227, 0.237, 0.222, 0.251, 0.284, 0.273, 0.280, 0.253, 0.175]
    x = np.arange(10)
    labels = [f"{a:.2f}–{b:.2f}" for a, b in zip(edges[:-1], edges[1:])]

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
    for ax in axes:
        ax.axvspan(-0.5, 0.5, color="0.9", lw=0)
        ax.axvspan(8.5, 9.5, color="0.9", lw=0)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=7)
        ax.set_xlabel("coherency_ratio decile")
    ax = axes[0]
    ax.plot(x, inn_std, "o-", color="k", label="std")
    ax.plot(x, inn_mean, "s--", color="0.5", label="mean")
    ax.axhline(0, color="0.7", lw=0.6)
    ax.set_ylabel("O − F$_{OL}$ (dB)")
    ax.set_title("(a) Innovations vs open loop")
    ax.legend()
    ax = axes[1]
    ax.plot(x, r_desroz, "o-", color="#D55E00", label="Desroziers estimate")
    ax.plot(x, r_assumed, "s--", color="0.4", label="assumed (scaled)")
    ax.set_ylabel("obs-error variance (dB$^2$)")
    ax.set_title("(b) Observation-error variance")
    ax.legend()
    ax = axes[2]
    ax.plot(x, np.abs(r_smap), "o-", color="#009E73", label="|r| with SMAP Tb innov.")
    ax.plot(x, r_l3, "s-", color="#0072B2", label="r with CYGNSS L3 innov.")
    ax.set_ylabel("correlation")
    ax.set_title("(c) Shared soil-moisture signal")
    ax.legend(loc="lower center")
    fig.suptitle("CYGNSS L1 quality vs coherency ratio (Jan–Jun 2020, ~86,000 obs); shaded deciles are removed by the screen")
    save(fig, "fig04_coherency_quality.png")


# ---------------------------------------------------------------------------------------------
# Fig. 5: L1 innovation correlation vs distance (error-correlation report, OL 2020-2022)
# ---------------------------------------------------------------------------------------------
def fig05_error_corr():
    d = np.array([0.025, 0.075, 0.125, 0.175, 0.25, 0.35, 0.45, 0.55])
    allp = np.array([0.566, 0.396, 0.317, 0.236, 0.202, 0.157, 0.162, 0.126])
    lo = np.array([0.542, 0.374, 0.288, 0.217, 0.186, 0.144, 0.149, 0.118])
    hi = np.array([0.589, 0.414, 0.347, 0.257, 0.215, 0.170, 0.176, 0.137])
    same = np.array([0.576, 0.435, 0.385, 0.290, 0.293, 0.240, 0.264, 0.240])
    diff = np.array([0.397, 0.209, 0.218, 0.144, 0.133, 0.112, 0.117, 0.091])
    # Model-implied innovation correlation for the R used (Gaussian, xcorr 0.15 deg), report §2.7.
    dm = np.array([0.032, 0.067, 0.127, 0.174, 0.257, 0.357, 0.450, 0.551])
    g015 = np.array([0.975, 0.908, 0.727, 0.564, 0.323, 0.167, 0.116, 0.096])

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.fill_between(d, lo, hi, color="0.8", lw=0)
    ax.plot(d, allp, "o-", color="k", label="all pairs (90% CI)")
    ax.plot(d, same, "^-", color="#D55E00", label="same track")
    ax.plot(d, diff, "v-", color="#0072B2", label="different spacecraft")
    ax.plot(dm, g015, "--", color="0.45", label="implied by R used (Gaussian, 0.15°)")
    ax.set_xlabel("separation (°)")
    ax.set_ylabel("correlation of L1 innovations (O − F$_{OL}$)")
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 0.6)
    ax.legend(frameon=False)
    ax.set_title("L1 innovation correlation vs distance (open loop, 2020–2022)")
    save(fig, "fig05_error_correlation.png")


# ---------------------------------------------------------------------------------------------
# Fig. 6: monthly O-F stdv change vs OL, 2020-2022
# ---------------------------------------------------------------------------------------------
def monthly_table():
    a = pd.read_csv(STATS / "summary" / "month_arm_omf_stdv_pct_vs_OL_2020.csv", dtype={"period": str})
    b = pd.read_csv(STATS / "summary" / "month_arm_omf_stdv_pct_vs_OL_2021_2022.csv", dtype={"period": str})
    m = pd.concat([a, b])
    m = m[m.period.str.len() == 6].copy()
    m["time"] = pd.to_datetime(m.period, format="%Y%m") + pd.Timedelta(days=14)
    return m


def fig06_monthly():
    m = monthly_table()
    panels = [("SMAP", "SMAP Tb"), ("SMOS", "SMOS Tb"), ("ASCAT", "ASCAT soil moisture"), ("L3", "CYGNSS L3 soil moisture")]
    fig, axes = plt.subplots(4, 1, figsize=(8.5, 8.4), sharex=True, constrained_layout=True)
    for ax, (col, title) in zip(axes, panels):
        shade_spring(ax)
        ax.axhline(0, color="0.5", lw=0.7)
        for tag in ["full_xc015_coh040216_fixedop", "l3_fixedop"]:
            lab, c = ARMS[tag]
            s = m[m.arm == tag].sort_values("time")
            own = " (assimilated)" if (tag == "l3_fixedop" and col == "L3") else ""
            ax.plot(s.time, s[col], "o-", ms=3, lw=1.3, color=c, label=f"{lab}{own}")
        ax.set_ylabel("Δ O−F stdv (%)")
        ax.set_title(title, loc="left")
        ax.legend(loc="lower left", ncol=2, frameon=False)
    fig.suptitle("Monthly change in O−F standard deviation vs the open loop (negative = better); shaded = spring windows")
    save(fig, "fig06_monthly_omf_change.png")


# ---------------------------------------------------------------------------------------------
# Fig. 7: per-tile maps of O-F stdv change, pooled 2020-2022
# ---------------------------------------------------------------------------------------------
def tile_pct(tag, group, nmin=50):
    da = xr.open_dataset(STATS / "stats" / f"temporal_stats_DA_paired_{tag}_20200101_20221231.nc4")
    ol = xr.open_dataset(STATS / "stats" / f"temporal_stats_OL_paired_monitor_xmask_{tag}_20200101_20221231.nc4")
    idx = GROUPS[group]
    d = da.OmF_stdv.values[:, idx]
    o = ol.OmF_stdv.values[:, idx]
    n = da.N_data.values[:, idx]
    pct = np.where((n >= nmin) & (o > 0), 100.0 * (d - o) / o, np.nan)
    with np.errstate(all="ignore"):
        return np.nanmean(pct, axis=1)


def fig07_maps():
    tiles = load_tiles()
    cols = ["SMAP Tb", "ASCAT SM", "CYGNSS L3 SM"]
    rows = ["full_xc015_coh040216_fixedop", "l3_fixedop"]
    fig = plt.figure(figsize=(11, 6.4))
    vmax = 10
    for i, tag in enumerate(rows):
        for j, grp in enumerate(cols):
            ax = map_axes(fig, (2, 3, i * 3 + j + 1))
            v = tile_pct(tag, grp)
            sc = tile_squares(ax, tiles.com_lon, tiles.com_lat, v, cmap="RdBu_r", vmin=-vmax, vmax=vmax, s=9)
            lab = ARMS[tag][0]
            own = " (assimilated)" if (tag == "l3_fixedop" and grp == "CYGNSS L3 SM") else ""
            ax.set_title(f"{lab} DA: {grp}{own}", fontsize=9)
            ax.text(0.02, 0.03, f"tile median {np.nanmedian(v):+.1f}%", transform=ax.transAxes, fontsize=7,
                    bbox=dict(fc="white", ec="none", alpha=0.8))
    cax = fig.add_axes([0.25, 0.05, 0.5, 0.02])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("Δ O−F stdv vs open loop, 2020–2022 (%; negative = better)")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.12, wspace=0.05, hspace=0.12)
    save(fig, "fig07_omf_change_maps.png")


# ---------------------------------------------------------------------------------------------
# Fig. 8: increment information: corr(dF, O-F_OL) and alpha_opt for SMAP Tb
# ---------------------------------------------------------------------------------------------
def fig08_noise_gain():
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 4.8), sharex=True, constrained_layout=True)
    for tag in ["full_xc015_coh040216_fixedop", "l3_fixedop"]:
        lab, c = ARMS[tag]
        ng = pd.read_csv(STATS / "summary" / f"noise_gain_{tag}_202001_202212.csv")
        ng = ng[(ng.grp == "SMAP") & (ng.month.str.len() == 7)].copy()
        ng["time"] = pd.to_datetime(ng.month) + pd.Timedelta(days=14)
        axes[0].plot(ng.time, ng["corr"], "o-", ms=3, lw=1.3, color=c, label=lab)
        axes[1].plot(ng.time, ng.alpha_opt, "o-", ms=3, lw=1.3, color=c, label=lab)
    for ax in axes:
        shade_spring(ax)
    axes[0].axhline(0, color="0.5", lw=0.7)
    axes[0].set_ylabel("corr(ΔF, O − F$_{OL}$)")
    axes[0].set_title("(a) Do the forecast changes point toward the SMAP Tb observations?", loc="left")
    axes[0].legend(frameon=False, ncol=2)
    axes[1].axhline(1, color="0.5", lw=0.7, ls="--")
    axes[1].axhline(0, color="0.5", lw=0.7)
    axes[1].set_ylabel("α$_{opt}$")
    axes[1].set_title("(b) Optimal increment scaling (1 = right size, <1 too large, ≈0 no information)", loc="left")
    save(fig, "fig08_increment_information.png")


# ---------------------------------------------------------------------------------------------
# Fig. 9 + Table 4: in-situ skill over the L1 area, tile-cluster bootstrap
# ---------------------------------------------------------------------------------------------
def ismn_deltas(nboot=4000, seed=1):
    st = pd.read_csv(ISMN / "ismn_skill_stations.csv")
    counts = pd.read_csv(ISMN / "tile_l1_obs_count_2020_2022.csv").set_index("tile_index").n_l1_obs
    st = st[st.tile_index.map(counts) >= L1_AREA_MIN_OBS]
    rng = np.random.default_rng(seed)
    rows = []
    for dom in ["surface", "rz"]:
        s = st[st.domain == dom]
        for metric, sign in [("R", 1), ("anomR", 1), ("ubRMSE", 1)]:
            wide = s.pivot_table(index=["station_key", "tile_index"], columns="run", values=metric).dropna()
            for comp, base in [("L1coh", "OL"), ("L3", "OL"), ("SMAP", "OL"), ("L1coh", "L3")]:
                dlt = (wide[comp] - wide[base]).reset_index()
                dlt.columns = ["station_key", "tile_index", "d"]
                by_tile = dlt.groupby("tile_index").d.agg(["sum", "count"])
                t = by_tile.index.values
                boot = np.empty(nboot)
                for k in range(nboot):
                    pick = rng.choice(t, size=t.size, replace=True)
                    sub = by_tile.loc[pick]
                    boot[k] = sub["sum"].sum() / sub["count"].sum()
                better = (dlt.d < 0).mean() if metric == "ubRMSE" else (dlt.d > 0).mean()
                rows.append(dict(domain=dom, metric=metric, comp=comp, base=base, n_st=len(dlt), n_tile=t.size,
                                 mean=dlt.d.mean(), lo=np.percentile(boot, 2.5), hi=np.percentile(boot, 97.5),
                                 better=better))
    return pd.DataFrame(rows)


def fig09_insitu(res):
    metrics = [("R", "ΔR"), ("anomR", "Δ anomaly R"), ("ubRMSE", "Δ ubRMSE (10$^{-3}$ m$^3$ m$^{-3}$)")]
    fig, axes = plt.subplots(2, 3, figsize=(10, 4.6), constrained_layout=True)
    comps = ["L1coh", "L3", "SMAP"]
    for i, dom in enumerate(["surface", "rz"]):
        for j, (met, lab) in enumerate(metrics):
            ax = axes[i, j]
            ax.axvline(0, color="0.5", lw=0.7)
            for k, comp in enumerate(comps):
                r = res[(res.domain == dom) & (res.metric == met) & (res.comp == comp) & (res.base == "OL")].iloc[0]
                f = 1e3 if met == "ubRMSE" else 1.0
                name, c = ISMN_RUNS[comp]
                sig = (r.lo > 0) or (r.hi < 0)
                ax.errorbar(r["mean"] * f, k, xerr=[[(r["mean"] - r.lo) * f], [(r.hi - r["mean"]) * f]], fmt="o",
                            color=c, mfc=c if sig else "white", ms=6, capsize=3, lw=1.3)
            ax.set_yticks(range(3))
            ax.set_yticklabels([ISMN_RUNS[c][0] for c in comps] if j == 0 else [])
            ax.invert_yaxis()
            if i == 1:
                ax.set_xlabel(lab)
            if j == 0:
                n = res[(res.domain == dom) & (res.metric == "R")].iloc[0]
                ax.set_ylabel(f"{'surface' if dom == 'surface' else 'root zone'}\n{n.n_st} stations, {n.n_tile} tiles")
    fig.suptitle("In-situ skill change vs open loop over the L1 area, 2020–2022 (95% tile-cluster bootstrap CI; filled = CI excludes 0)")
    save(fig, "fig09_insitu_skill.png")


def print_table4(res):
    print("\nTable 4 (paste into the outline):")
    print("| domain | metric | L1 − OL | L3 − OL | SMAP − OL | L1 − L3 |")
    print("|---|---|---|---|---|---|")
    for dom in ["surface", "rz"]:
        for met in ["R", "anomR", "ubRMSE"]:
            cells = []
            for comp, base in [("L1coh", "OL"), ("L3", "OL"), ("SMAP", "OL"), ("L1coh", "L3")]:
                r = res[(res.domain == dom) & (res.metric == met) & (res.comp == comp) & (res.base == base)].iloc[0]
                f, fmt = (1e3, "{:+.2f}") if met == "ubRMSE" else (1.0, "{:+.3f}")
                sig = "\\*" if (r.lo > 0 or r.hi < 0) else ""
                cells.append(f"{fmt.format(r['mean'] * f)} [{fmt.format(r.lo * f)}, {fmt.format(r.hi * f)}]{sig}")
            unit = " (10⁻³)" if met == "ubRMSE" else ""
            print(f"| {'surface' if dom == 'surface' else 'root zone'} | {met}{unit} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    fig01_domain()
    fig02_schematic()
    fig03_pipeline()
    fig04_coherency()
    fig05_error_corr()
    fig06_monthly()
    fig07_maps()
    fig08_noise_gain()
    res = ismn_deltas()
    fig09_insitu(res)
    print_table4(res)
