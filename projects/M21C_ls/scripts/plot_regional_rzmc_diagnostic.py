#!/usr/bin/env python3
"""Diagnostic figures for the regional RZMC DA-OL transition analysis."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs, cartopy.feature as cfeature

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from trend_breakpoint_series import MonthlySeriesLoader  # noqa: E402
from m21c_periods import load_period_frames  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "output" / "regional_rzmc_transitions"
FIGDIR = ROOT / "docs" / "regional_rzmc_diagnostic_figures"
REG = json.loads((ROOT/"config"/"regional_rzmc_regions.json").read_text())["regions"]
PCOL = ["#e8e2d0","#dce6f0","#dcecdc","#f6e2d2","#e6dcee","#d8ece8","#f6dce6","#dce4f0","#ece8d8"]


def fig_masks():
    d = xr.open_dataset(OUT.parent/"trends_breakpoints"/"RZMC_delta_valid_land_trend_statistics.nc")
    slope = d["slope"].values*1000.0
    sig = d["significant_fdr"].values.astype(bool)
    with MonthlySeriesLoader() as L:
        lat = np.asarray(L.lat.values); lon = np.asarray(L.lon.values)
    fig = plt.figure(figsize=(13.0, 6.2))
    ax = fig.add_subplot(1,1,1, projection=ccrs.Robinson())
    ax.set_global(); ax.set_extent([-180,180,-60,84], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="0.45")
    ok = np.isfinite(slope)
    s = ax.scatter(lon[ok], lat[ok], c=slope[ok], s=0.6, marker="s",
                   cmap="RdBu_r", vmin=-0.76, vmax=0.76,
                   transform=ccrs.PlateCarree(), rasterized=True)
    m = ok & sig
    ax.scatter(lon[m], lat[m], c="k", s=0.06, marker=".", alpha=0.35,
               transform=ccrs.PlateCarree(), rasterized=True)
    for r in REG:
        b = r["bounds"]
        if b is None: continue
        ax.add_patch(Rectangle((b["lon_min"], b["lat_min"]),
                               b["lon_max"]-b["lon_min"], b["lat_max"]-b["lat_min"],
                               transform=ccrs.PlateCarree(), facecolor="none",
                               edgecolor="#111111", linewidth=1.8, zorder=6))
        ax.text(b["lon_min"]+1.0, b["lat_max"]-2.5, r["label"], transform=ccrs.PlateCarree(),
                fontsize=8.5, fontweight="bold", va="top", zorder=7,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.3", lw=0.6, alpha=0.9))
    cb = fig.colorbar(s, ax=ax, orientation="horizontal", pad=0.04, shrink=0.6, extend="both")
    cb.set_label(r"RZMC DA $-$ OL trend ($\times10^{-3}$ m$^3$ m$^{-3}$ yr$^{-1}$)")
    ax.set_title("Fixed diagnostic regions on the Fig. 16 RZMC DA $-$ OL trend field\n"
                 "(black stipple = significant after FDR; boxes are pre-specified, not trend-selected)",
                 fontsize=10.5)
    fig.savefig(FIGDIR/"regional_masks_on_trend_field.png", dpi=200, bbox_inches="tight")
    plt.close(fig); print("wrote regional_masks_on_trend_field.png")


def fig_series():
    ds = xr.open_dataset(OUT/"regional_rzmc_monthly.nc")
    t = pd.DatetimeIndex(ds.time.values)
    tab = pd.read_csv(OUT/"regional_breakpoint_table.csv")
    _, fine, _, _ = load_period_frames()
    label_of = {r["region_id"]: r["label"] for r in REG}
    order = [r["region_id"] for r in REG]

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.4), sharex=True)
    for k, rid in enumerate(order):
        ax = axes.flat[k]
        v = ds["delta"].sel(region=rid).values
        for j, row in enumerate(fine.itertuples()):
            ax.axvspan(row.start, row.end, color=PCOL[j % len(PCOL)], zorder=0, lw=0)
            ax.text(row.start + (row.end-row.start)/2, 0.015, row.period_id,
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=6.2, color="0.30", zorder=6)
        ax.axhline(0, color="0.35", lw=0.7, zorder=1)
        ax.plot(t, v, color="#b0483f", lw=0.55, alpha=0.55, zorder=2)
        rm = pd.Series(v, index=t).rolling(12, center=True, min_periods=12).mean()
        ax.plot(rm.index, rm.values, color="#8c1d16", lw=1.9, zorder=3)
        sub = tab[(tab.region_id == rid) & (tab.series == "delta") & tab.detected.notna()]
        for _, r in sub.iterrows():
            ax.axvline(pd.Timestamp(r.detected), color="#1f4e79", lw=1.7, ls="--", zorder=4)
            ax.text(pd.Timestamp(r.detected), 0.965, pd.Timestamp(r.detected).strftime("%Y-%m"),
                    transform=ax.get_xaxis_transform(), rotation=90, fontsize=6.6,
                    color="#1f4e79", ha="right", va="top", zorder=5)
        olsub = tab[(tab.region_id == rid) & (tab.series == "ol") & tab.detected.notna()]
        for _, r in olsub.iterrows():
            ax.axvline(pd.Timestamp(r.detected), color="#6b6b6b", lw=1.4, ls=":", zorder=4)
        n = int(pd.read_csv(OUT/"regional_support_audit.csv").set_index("region_id").loc[rid,"n_tiles"])
        ax.set_title(f"({chr(97+k)}) {label_of[rid]}  (n = {n:,} tiles)", fontsize=10)
        ax.grid(axis="y", color="0.88", lw=0.6); ax.set_axisbelow(True)
        if k % 3 == 0: ax.set_ylabel(r"RZMC DA $-$ OL (m$^3$ m$^{-3}$)")
        ax.xaxis.set_major_locator(mdates.YearLocator(4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    h = [plt.Line2D([],[],color="#b0483f",lw=0.8,alpha=0.6),
         plt.Line2D([],[],color="#8c1d16",lw=2.0),
         plt.Line2D([],[],color="#1f4e79",lw=1.7,ls="--"),
         plt.Line2D([],[],color="#6b6b6b",lw=1.4,ls=":")]
    fig.legend(h, ["monthly DA $-$ OL","12-month running mean",
                   "accepted changepoint (DA $-$ OL)","accepted changepoint (OL control)"],
               loc="lower center", ncols=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5,-0.005))
    fig.suptitle("Regional area-weighted RZMC DA $-$ OL with blind changepoints, June 2000 – May 2024",
                 fontsize=12.5, y=0.98)
    fig.tight_layout(rect=[0,0.045,1,0.955])
    fig.savefig(FIGDIR/"regional_rzmc_series_changepoints.png", dpi=200, bbox_inches="tight")
    plt.close(fig); print("wrote regional_rzmc_series_changepoints.png")




def fig_membership():
    """Map the tiles actually selected by each region (definition figure)."""
    from matplotlib.lines import Line2D
    with MonthlySeriesLoader() as L:
        lat = np.asarray(L.lat.values); lon = np.asarray(L.lon.values)
        vl = np.asarray(L.mask("valid_land").values).astype(bool)
    cols = {"aus": "#d62728", "safr": "#2ca02c", "nafr_me": "#ff7f0e",
            "wna": "#1f77b4", "neur": "#9467bd"}
    fig = plt.figure(figsize=(13, 6.0))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global(); ax.set_extent([-180, 180, -60, 84], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="0.5")
    ax.scatter(lon[vl], lat[vl], c="0.86", s=0.35, marker="s",
               transform=ccrs.PlateCarree(), rasterized=True)
    handles, labels = [], []
    for r in REG:
        b = r["bounds"]
        if b is None:
            continue
        m = (vl & (lat >= b["lat_min"]) & (lat <= b["lat_max"])
             & (lon >= b["lon_min"]) & (lon <= b["lon_max"]))
        ax.scatter(lon[m], lat[m], c=cols[r["region_id"]], s=0.5, marker="s",
                   transform=ccrs.PlateCarree(), rasterized=True)
        handles.append(Line2D([], [], marker="s", ls="", color=cols[r["region_id"]], ms=6))
        labels.append(f"{r['label']}  (n={int(m.sum()):,})")
    ax.legend(handles, labels, loc="lower left", frameon=False, fontsize=8.5,
              ncols=2, bbox_to_anchor=(0.0, -0.13))
    ax.set_title("Tiles actually selected by each region "
                 "(grey = all other valid land)", fontsize=11)
    fig.savefig(FIGDIR/"regional_tile_membership.png", dpi=200, bbox_inches="tight")
    plt.close(fig); print("wrote regional_tile_membership.png")

if __name__ == "__main__":
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig_masks(); fig_membership(); fig_series()
