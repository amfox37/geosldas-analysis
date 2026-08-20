#!/usr/bin/env python3
"""CYGNSS L1 observation-error sweep: OmF standard deviation by experiment.

Compares the scaled open-loop baseline against three CYGNSS L1 assimilation
runs with full, half and quarter observation-error variance, over calendar
2020 on the CONUS-southwest ("AZ") box. CYGNSS L1 (species 12) is the only
assimilated species; everything else is monitored, so its OmF is an
independent check on whether the assimilation helps.

Follows the conventions in plot_az_omf_summary.py: species are grouped into
observing systems, statistics are N_data-weighted, and tiles with fewer than
--nmin observations are excluded from those weighted means.
"""
from __future__ import annotations
from pathlib import Path
import argparse, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import TwoSlopeNorm, Normalize
from netCDF4 import Dataset
import pickle
from pyproj import CRS, Transformer

PROJECT = Path(__file__).resolve().parent.parent
REPO = PROJECT.parents[1]
sys.path.insert(0, str(REPO / "common/python/io"))
from read_GEOSldas import read_tilecoord

STATS = PROJECT / "output" / "stats_output"
TILECOORD = PROJECT / "OLv8_M36_all_sensors_AZ_describe" / "OLv8_M36_all_sensors_AZ.ldas_tilecoord.bin"
OUT = PROJECT / "output" / "R_sweep_figures"

EXPERIMENTS = [
    ("OL baseline", "OL_scaled_baseline", "#333333", "-"),
    ("full R", "DA_cygl1assim_fullR", "#4c78a8", "-"),
    ("half R", "DA_cygl1assim_halfR", "#f58518", "-"),
    ("quarter R", "DA_cygl1assim_quarterR", "#c0392b", "-"),
]
GROUPS = [
    ("SMOS", (0, 1, 2, 3), "K"),
    ("SMAP", (4, 5, 6, 7), "K"),
    ("ASCAT", (8, 9, 10), "%"),
    ("CYGNSS L3", (11,), "m$^3$ m$^{-3}$"),
    ("CYGNSS L1*", (12,), "dB"),
]
EXTENT = (-118.6, -105.2, 28.6, 40.4)
SCALE, NCOL, NROW = 36032.220840584, 964, 406


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nmin", type=int, default=20,
                   help="minimum per-tile observation count for weighted means")
    return p.parse_args()


def load_temporal(tag):
    with Dataset(STATS / f"temporal_stats_{tag}_allspecies_20200101_20201231.nc4") as nc:
        d = {k: np.asarray(nc[k][:], dtype=float) for k in nc.variables}
    for k in d:
        d[k][~np.isfinite(d[k])] = np.nan
    return d


def load_monthly(tag):
    with open(STATS / f"spatial_stats_{tag}_allspecies_202001_202012.pkl", "rb") as f:
        return pickle.load(f)


def group_mean(values, counts, indices, nmin):
    """N_data-weighted mean across a group's species, per tile."""
    v = values[:, list(indices)]
    w = counts[:, list(indices)]
    ok = np.isfinite(v) & np.isfinite(w) & (w >= nmin)
    num = np.where(ok, v * w, 0.0).sum(axis=1)
    den = np.where(ok, w, 0.0).sum(axis=1)
    return np.divide(num, den, out=np.full(den.shape, np.nan), where=den > 0)


def grid_edges():
    tr = Transformer.from_crs(
        CRS.from_proj4("+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"),
        CRS.from_epsg(4326), always_xy=True)
    x = (np.arange(NCOL + 1) - 0.5 - (NCOL - 1) / 2.0) * SCALE
    y = ((NROW - 1) / 2.0 - (np.arange(NROW + 1) - 0.5)) * SCALE
    lon, _ = tr.transform(x, np.zeros_like(x))
    _, lat = tr.transform(np.zeros_like(y), y)
    return lon, lat


def base_map(ax):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="0.94", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), lw=0.55, edgecolor="0.25", zorder=3)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), lw=0.45, edgecolor="0.35", zorder=3)
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_1_states_provinces_lakes",
                                                "50m", facecolor="none"),
                   lw=0.45, edgecolor="0.35", zorder=3)


def main():
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tc = read_tilecoord(str(TILECOORD))
    icol = np.asarray(tc["i_indg"]).astype(int)
    irow = np.asarray(tc["j_indg"]).astype(int)
    LON, LAT = grid_edges()
    c0, c1 = icol.min(), icol.max() + 1
    r0, r1 = irow.min(), irow.max() + 1

    temporal = {tag: load_temporal(tag) for _, tag, _, _ in EXPERIMENTS}
    monthly = {tag: load_monthly(tag) for _, tag, _, _ in EXPERIMENTS}

    def to_grid(v):
        f = np.full((NROW, NCOL), np.nan)
        f[irow, icol] = v
        return np.ma.masked_invalid(f[r0:r1, c0:c1])

    # ---------------- Figure 1: maps ----------------
    nrow, ncol = len(GROUPS), len(EXPERIMENTS)
    fig = plt.figure(figsize=(4.1 * ncol, 3.25 * nrow))
    for gi, (gname, idx, units) in enumerate(GROUPS):
        base = group_mean(temporal["OL_scaled_baseline"]["OmF_stdv"],
                          temporal["OL_scaled_baseline"]["N_data"], idx, args.nmin)
        diffs = {tag: group_mean(temporal[tag]["OmF_stdv"], temporal[tag]["N_data"], idx, args.nmin) - base
                 for _, tag, _, _ in EXPERIMENTS[1:]}
        dmax = np.nanpercentile(np.abs(np.concatenate(list(diffs.values()))), 98)
        dmax = max(dmax, 1e-9)
        for ei, (ename, tag, _, _) in enumerate(EXPERIMENTS):
            ax = fig.add_subplot(nrow, ncol, gi * ncol + ei + 1, projection=ccrs.PlateCarree())
            base_map(ax)
            if ei == 0:
                v, cmap, norm = base, "viridis", Normalize(*np.nanpercentile(base, [2, 98]))
            else:
                v, cmap, norm = diffs[tag], "RdBu_r", TwoSlopeNorm(vcenter=0, vmin=-dmax, vmax=dmax)
            m = ax.pcolormesh(LON[c0:c1 + 1], LAT[r0:r1 + 1], to_grid(v),
                              cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), shading="flat", zorder=1)
            if gi == 0:
                ax.set_title("O$-$F stdv, open loop" if ei == 0 else f"{ename} $-$ OL",
                             fontsize=10.5, fontweight="bold")
            if ei == 0:
                ax.text(-0.06, 0.5, f"{gname}\n({units})", transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=10, fontweight="bold")
            else:
                mv = np.nansum(v * np.isfinite(v)) / max(np.isfinite(v).sum(), 1)
                ax.text(0.02, 0.03, f"mean {mv:+.3g}", transform=ax.transAxes, fontsize=8,
                        bbox=dict(fc="white", ec="0.7", alpha=0.85, pad=1.5))
            cb = fig.colorbar(m, ax=ax, fraction=0.038, pad=0.02)
            cb.ax.tick_params(labelsize=7)
    fig.suptitle("CYGNSS L1 observation-error sweep: O$-$F standard deviation, 2020\n"
                 "red = degraded relative to open loop.  CYGNSS L1 is the only assimilated species.",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0.01, 0, 1, 0.975])
    p1 = OUT / "R_sweep_fig01_omf_stdv_maps.png"
    fig.savefig(p1, dpi=170, bbox_inches="tight"); plt.close(fig)
    print("wrote", p1.relative_to(PROJECT))

    # ---------------- Figure 2: monthly time series ----------------
    months = [str(m) for m in monthly["OL_scaled_baseline"]["date_vec"]]
    x = np.arange(len(months))
    fig, axes = plt.subplots(len(GROUPS), 2, figsize=(13.5, 2.5 * len(GROUPS)),
                             gridspec_kw=dict(width_ratios=[1.35, 1.0], wspace=0.22, hspace=0.35))
    for gi, (gname, idx, units) in enumerate(GROUPS):
        ax, axr = axes[gi]
        series = {}
        for ename, tag, colour, ls in EXPERIMENTS:
            n = np.asarray(monthly[tag]["N_data"])[:, list(idx)]
            v = np.asarray(monthly[tag]["OmF_stdv"])[:, list(idx)]
            ok = np.isfinite(v) & np.isfinite(n) & (n > 0)
            s = np.divide(np.where(ok, v * n, 0).sum(axis=1), np.where(ok, n, 0).sum(axis=1),
                          out=np.full(len(months), np.nan), where=np.where(ok, n, 0).sum(axis=1) > 0)
            series[tag] = s
            ax.plot(x, s, ls, color=colour, lw=1.8, marker="o", ms=3.5, label=ename)
        for ename, tag, colour, ls in EXPERIMENTS[1:]:
            axr.plot(x, 100 * (series[tag] / series["OL_scaled_baseline"] - 1), ls, color=colour,
                     lw=1.8, marker="o", ms=3.5, label=ename)
        axr.axhline(0, color="0.3", lw=0.9)
        for a in (ax, axr):
            a.set_xticks(x); a.set_xticklabels([m[4:] for m in months], fontsize=8)
            a.grid(color="0.92"); a.set_axisbelow(True)
        ax.set_ylabel(f"O$-$F stdv ({units})", fontsize=9)
        axr.set_ylabel("change vs OL (%)", fontsize=9)
        ax.set_title(f"({chr(97+gi)}) {gname}", loc="left", fontweight="bold", fontsize=11)
        if gi == 0:
            ax.legend(frameon=False, fontsize=8.5, ncols=4)
        if gi == len(GROUPS) - 1:
            ax.set_xlabel("month of 2020", fontsize=9); axr.set_xlabel("month of 2020", fontsize=9)
    fig.suptitle("CYGNSS L1 observation-error sweep: monthly O$-$F standard deviation, 2020",
                 fontsize=13, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    p2 = OUT / "R_sweep_fig02_omf_stdv_monthly.png"
    fig.savefig(p2, dpi=170, bbox_inches="tight"); plt.close(fig)
    print("wrote", p2.relative_to(PROJECT))


    # ---------------- Figure 3: observation counts ----------------
    DAYS = 366  # calendar 2020
    fig = plt.figure(figsize=(4.1 * len(GROUPS), 6.6))
    gs = fig.add_gridspec(2, len(GROUPS), height_ratios=[1.25, 1.0], hspace=0.30, wspace=0.18)
    for gi, (gname, idx, _units) in enumerate(GROUPS):
        ax = fig.add_subplot(gs[0, gi], projection=ccrs.PlateCarree())
        base_map(ax)
        counts = np.nansum(temporal["OL_scaled_baseline"]["N_data"][:, list(idx)], axis=1) / DAYS
        counts[counts <= 0] = np.nan
        m = ax.pcolormesh(LON[c0:c1 + 1], LAT[r0:r1 + 1], to_grid(counts), cmap="viridis",
                          norm=Normalize(0, np.nanpercentile(counts, 98)),
                          transform=ccrs.PlateCarree(), shading="flat", zorder=1)
        ax.set_title(f"({chr(97+gi)}) {gname}", fontsize=10.5, fontweight="bold", loc="left")
        ax.text(0.02, 0.03, f"mean {np.nanmean(counts):.2f} d$^{{-1}}$\n{np.isfinite(counts).sum()} tiles",
                transform=ax.transAxes, fontsize=8,
                bbox=dict(fc="white", ec="0.7", alpha=0.85, pad=1.5))
        cb = fig.colorbar(m, ax=ax, fraction=0.038, pad=0.02); cb.ax.tick_params(labelsize=7)
    axb = fig.add_subplot(gs[1, :])
    months = [str(mm) for mm in monthly["OL_scaled_baseline"]["date_vec"]]
    x = np.arange(len(months))
    for gname, idx, _u in GROUPS:
        tot = np.nansum(np.asarray(monthly["OL_scaled_baseline"]["N_data"])[:, list(idx)], axis=1)
        axb.plot(x, tot, "-o", ms=4, lw=1.8, label=gname)
    spread = max(abs(np.nansum(np.asarray(monthly[t]["N_data"])) /
                     np.nansum(np.asarray(monthly["OL_scaled_baseline"]["N_data"])) - 1)
                 for _, t, _, _ in EXPERIMENTS)
    axb.set_yscale("log"); axb.set_xticks(x); axb.set_xticklabels([mm[4:] for mm in months])
    axb.set_xlabel("month of 2020"); axb.set_ylabel("observations per month")
    axb.grid(color="0.92"); axb.set_axisbelow(True)
    axb.legend(frameon=False, ncols=5, fontsize=9)
    axb.set_title(f"(f) Monthly assimilated/monitored observation counts  "
                  f"(experiments agree to within {100*spread:.2f}%)",
                  loc="left", fontweight="bold", fontsize=11)
    fig.suptitle("CYGNSS L1 R sweep: observation counts, 2020 (maps show mean observations per day)",
                 fontsize=13, y=0.98)
    p3 = OUT / "R_sweep_fig03_observation_counts.png"
    fig.savefig(p3, dpi=170, bbox_inches="tight"); plt.close(fig)
    print("wrote", p3.relative_to(PROJECT))

    # ---------------- Figure 4: mean observation and forecast ----------------
    fig, axes = plt.subplots(len(GROUPS), 2, figsize=(13.5, 2.6 * len(GROUPS)),
                             gridspec_kw=dict(width_ratios=[1.0, 1.45], wspace=0.22, hspace=0.38),
                             subplot_kw=None)
    for gi, (gname, idx, units) in enumerate(GROUPS):
        axes[gi, 0].remove()
        ax = fig.add_subplot(len(GROUPS), 2, gi * 2 + 1, projection=ccrs.PlateCarree())
        base_map(ax)
        omean = group_mean(temporal["OL_scaled_baseline"]["O_mean"],
                           temporal["OL_scaled_baseline"]["N_data"], idx, args.nmin)
        m = ax.pcolormesh(LON[c0:c1 + 1], LAT[r0:r1 + 1], to_grid(omean), cmap="viridis",
                          norm=Normalize(*np.nanpercentile(omean, [2, 98])),
                          transform=ccrs.PlateCarree(), shading="flat", zorder=1)
        ax.set_title(f"({chr(97+gi)}) {gname}", fontsize=10.5, fontweight="bold", loc="left")
        cb = fig.colorbar(m, ax=ax, fraction=0.038, pad=0.02)
        cb.ax.tick_params(labelsize=7); cb.set_label(f"mean obs ({units})", fontsize=8)

        axt = axes[gi, 1]
        def wser(tag, key):
            n = np.asarray(monthly[tag]["N_data"])[:, list(idx)]
            v = np.asarray(monthly[tag][key])[:, list(idx)]
            ok = np.isfinite(v) & np.isfinite(n) & (n > 0)
            den = np.where(ok, n, 0).sum(axis=1)
            return np.divide(np.where(ok, v * n, 0).sum(axis=1), den,
                             out=np.full(len(months), np.nan), where=den > 0)
        axt.plot(x, wser("OL_scaled_baseline", "O_mean"), "-", color="black", lw=2.4,
                 marker="s", ms=4, label="observation", zorder=5)
        for ename, tag, colour, _ls in EXPERIMENTS:
            axt.plot(x, wser(tag, "F_mean"), "--", color=colour, lw=1.7, marker="o", ms=3.2,
                     label=f"forecast, {ename}")
        axt.set_xticks(x); axt.set_xticklabels([mm[4:] for mm in months], fontsize=8)
        axt.grid(color="0.92"); axt.set_axisbelow(True)
        axt.set_ylabel(f"mean ({units})", fontsize=9)
        if gi == 0:
            h, lb = axt.get_legend_handles_labels()
            fig.legend(h, lb, loc="upper center", ncols=5, frameon=False, fontsize=9,
                       bbox_to_anchor=(0.5, 0.955))
        if gi == len(GROUPS) - 1:
            axt.set_xlabel("month of 2020", fontsize=9)
    fig.suptitle("CYGNSS L1 R sweep: mean observation and forecast, 2020\n"
                 "observations are identical across experiments; only the forecast responds to assimilation",
                 fontsize=12.5, y=1.005)
    p4 = OUT / "R_sweep_fig04_mean_obs_and_forecast.png"
    fig.savefig(p4, dpi=170, bbox_inches="tight"); plt.close(fig)
    print("wrote", p4.relative_to(PROJECT))


    # ---------------- Figure 5: O-A, is the analysis doing anything? ----------------
    fig = plt.figure(figsize=(4.1 * len(GROUPS), 8.4))
    gs = fig.add_gridspec(2, len(GROUPS), height_ratios=[1.0, 1.15], hspace=0.34, wspace=0.30)
    for gi, (gname, idx, units) in enumerate(GROUPS):
        ax = fig.add_subplot(gs[0, gi])
        xb = np.arange(len(EXPERIMENTS)); w = 0.38
        f = [np.nanmean(group_mean(temporal[t]["OmF_stdv"], temporal[t]["N_data"], idx, args.nmin))
             for _, t, _, _ in EXPERIMENTS]
        a = [np.nanmean(group_mean(temporal[t]["OmA_stdv"], temporal[t]["N_data"], idx, args.nmin))
             for _, t, _, _ in EXPERIMENTS]
        ax.bar(xb - w/2, f, w, color="0.72", edgecolor="0.3", lw=0.5, label="O$-$F (background)")
        ax.bar(xb + w/2, a, w, color=[c for _, _, c, _ in EXPERIMENTS], edgecolor="0.3", lw=0.5,
               label="O$-$A (analysis)")
        for k in range(len(EXPERIMENTS)):
            d = 100 * (a[k] / f[k] - 1)
            ax.annotate(f"{d:+.1f}%", (xb[k] + w/2, a[k]), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=7.6,
                        color="#1a7f37" if d < 0 else "#c0392b", fontweight="bold")
        ax.set_xticks(xb); ax.set_xticklabels([e for e, _, _, _ in EXPERIMENTS], fontsize=8, rotation=20)
        ax.set_ylabel(f"stdv ({units})", fontsize=9)
        lo, hi = min(min(f), min(a)), max(max(f), max(a))
        ax.set_ylim(lo - 0.18 * (hi - lo), hi + 0.22 * (hi - lo))
        ax.grid(axis="y", color="0.92"); ax.set_axisbelow(True)
        ax.set_title(f"({chr(97+gi)}) {gname}", loc="left", fontweight="bold", fontsize=10.5)
        if gi == 0:
            ax.legend(frameon=False, fontsize=7.6, loc="upper left")

        axm = fig.add_subplot(gs[1, gi], projection=ccrs.PlateCarree())
        base_map(axm)
        tag = EXPERIMENTS[-1][1]
        eff = (group_mean(temporal[tag]["OmA_stdv"], temporal[tag]["N_data"], idx, args.nmin)
               - group_mean(temporal[tag]["OmF_stdv"], temporal[tag]["N_data"], idx, args.nmin))
        lim = max(np.nanpercentile(np.abs(eff), 98), 1e-12)
        mm = axm.pcolormesh(LON[c0:c1 + 1], LAT[r0:r1 + 1], to_grid(eff), cmap="RdBu_r",
                            norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim),
                            transform=ccrs.PlateCarree(), shading="flat", zorder=1)
        axm.set_title(f"({chr(97+len(GROUPS)+gi)}) {gname}", loc="left", fontweight="bold", fontsize=10)
        axm.text(0.02, 0.03, f"mean {np.nanmean(eff):+.3g}", transform=axm.transAxes, fontsize=8,
                 bbox=dict(fc="white", ec="0.7", alpha=0.85, pad=1.5))
        cb = fig.colorbar(mm, ax=axm, fraction=0.038, pad=0.02)
        cb.ax.tick_params(labelsize=7); cb.set_label(f"O$-$A $-$ O$-$F ({units})", fontsize=7.5)
    fig.suptitle("Is the analysis doing anything?  O$-$A versus O$-$F standard deviation, 2020\n"
                 "top: domain means, percentage is O$-$A relative to O$-$F (green = analysis fits better).  "
                 "bottom: O$-$A minus O$-$F at quarter R, blue = analysis pulled toward that sensor",
                 fontsize=12.5, y=1.005)
    p5 = OUT / "R_sweep_fig05_oma_vs_omf.png"
    fig.savefig(p5, dpi=170, bbox_inches="tight"); plt.close(fig)
    print("wrote", p5.relative_to(PROJECT))


if __name__ == "__main__":
    main()
