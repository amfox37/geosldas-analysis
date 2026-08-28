#!/usr/bin/env python3
"""Map OmF standard-deviation changes for a CYGNSS L1 density experiment.

The primary metric follows the land-sweeper manuscript convention:

    100 * (sigma_OL - sigma_DA) / sigma_OL

Positive values indicate lower OmF variability in DA. Observation systems with
multiple species are combined with observation-count weighting after requiring
the requested minimum count in both experiments.
"""
from __future__ import annotations

import csv
import argparse
import calendar
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from netCDF4 import Dataset
import numpy as np


PROJECT = Path(__file__).resolve().parent.parent
REPO = PROJECT.parents[1]
sys.path.insert(0, str(REPO / "common/python/io"))
from read_GEOSldas import read_tilecoord


STATS = PROJECT / "output/thinning_expts"
OL_FILE = STATS / "temporal_stats_OL_paired_monitor_20200101_20201031.nc4"
TILECOORD = (
    PROJECT
    / "OLv8_M36_all_sensors_AZ_describe"
    / "OLv8_M36_all_sensors_AZ.ldas_tilecoord.bin"
)
OUT = STATS / "figures"
NMIN = 10
PERCENT_LIMIT = 20.0
EXTENT = (-118.6, -105.2, 28.6, 40.4)


@dataclass(frozen=True)
class Group:
    name: str
    indices: tuple[int, ...]
    units: str


GROUPS = (
    Group("SMOS", (0, 1, 2, 3), "K"),
    Group("SMAP", (4, 5, 6, 7), "K"),
    Group("ASCAT", (8, 9, 10), "m$^3$ m$^{-3}$"),
    Group("CYGNSS L3", (11,), "m$^3$ m$^{-3}$"),
    Group("CYGNSS L1", (12,), "dB"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experiment",
        choices=("intermediate", "dense"),
        default="intermediate",
        help="CYGNSS L1 density experiment to compare with OL",
    )
    parser.add_argument(
        "--ol-support",
        choices=("standard", "xmask"),
        default="standard",
        help="use the full OL monitor or the intermediate retained-event mask",
    )
    parser.add_argument(
        "--end-month",
        default="202010",
        help="final statistics month in YYYYMM format (default: 202010)",
    )
    return parser.parse_args()


def period_settings(end_month: str) -> tuple[str, str, str]:
    try:
        end = datetime.strptime(end_month, "%Y%m")
    except ValueError as exc:
        raise ValueError("--end-month must use YYYYMM format") from exc
    if end < datetime(2020, 1, 1):
        raise ValueError("--end-month cannot precede 202001")
    last_day = calendar.monthrange(end.year, end.month)[1]
    file_period = f"20200101_{end:%Y%m}{last_day:02d}"
    month_period = f"202001_{end:%Y%m}"
    label = f"January 2020-{end:%B %Y}"
    return file_period, month_period, label


def load_stats(path: Path) -> dict[str, np.ndarray]:
    with Dataset(path) as nc:
        data = {
            name: np.ma.filled(nc[name][:], np.nan).astype(float)
            for name in ("OmF_stdv", "N_data")
        }
    return data


def common_group_values(
    ol: dict[str, np.ndarray], da: dict[str, np.ndarray], group: Group
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    idx = list(group.indices)
    ol_std = ol["OmF_stdv"][:, idx]
    da_std = da["OmF_stdv"][:, idx]
    ol_n = ol["N_data"][:, idx]
    da_n = da["N_data"][:, idx]
    common = (
        np.isfinite(ol_std)
        & np.isfinite(da_std)
        & np.isfinite(ol_n)
        & np.isfinite(da_n)
        & (ol_n >= NMIN)
        & (da_n >= NMIN)
    )

    def weighted(values: np.ndarray, counts: np.ndarray) -> np.ndarray:
        denominator = np.where(common, counts, 0.0).sum(axis=1)
        return np.divide(
            np.where(common, values * counts, 0.0).sum(axis=1),
            denominator,
            out=np.full(denominator.shape, np.nan),
            where=denominator > 0,
        )

    ol_group = weighted(ol_std, ol_n)
    da_group = weighted(da_std, da_n)
    percent = 100.0 * (ol_group - da_group) / ol_group
    common_count = common.sum(axis=1)
    return ol_group, da_group, percent, common_count


def base_map(ax) -> None:
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="0.94", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), lw=0.5, edgecolor="0.25")
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), lw=0.4, edgecolor="0.35")
    ax.add_feature(
        cfeature.NaturalEarthFeature(
            "cultural", "admin_1_states_provinces_lakes", "50m", facecolor="none"
        ),
        lw=0.4,
        edgecolor="0.4",
    )


def scatter_map(ax, lon, lat, values, cmap, norm):
    valid = np.isfinite(values)
    return ax.scatter(
        lon[valid],
        lat[valid],
        c=values[valid],
        s=21,
        marker="s",
        linewidths=0,
        cmap=cmap,
        norm=norm,
        transform=ccrs.PlateCarree(),
        zorder=1,
    )


def aggregate_percent(
    ol: dict[str, np.ndarray], da: dict[str, np.ndarray], group: Group
) -> float:
    idx = list(group.indices)
    common = (
        np.isfinite(ol["OmF_stdv"][:, idx])
        & np.isfinite(da["OmF_stdv"][:, idx])
        & (ol["N_data"][:, idx] >= NMIN)
        & (da["N_data"][:, idx] >= NMIN)
    )

    def pooled(data: dict[str, np.ndarray]) -> float:
        values = data["OmF_stdv"][:, idx]
        counts = data["N_data"][:, idx]
        return np.sum(np.where(common, values * counts, 0.0)) / np.sum(
            np.where(common, counts, 0.0)
        )

    ol_value = pooled(ol)
    da_value = pooled(da)
    return 100.0 * (ol_value - da_value) / ol_value


def main() -> None:
    args = parse_args()
    experiment = args.experiment
    label = experiment.capitalize()
    file_period, month_period, period_label = period_settings(args.end_month)
    if args.ol_support == "xmask" and experiment != "intermediate":
        raise ValueError("xmask OL support is available only for the intermediate experiment")
    ol_file = (
        STATS / f"temporal_stats_OL_paired_monitor_xmask_intermediate_{file_period}.nc4"
        if args.ol_support == "xmask"
        else STATS / f"temporal_stats_OL_paired_monitor_{file_period}.nc4"
    )
    ol_label = "xmask-paired open loop" if args.ol_support == "xmask" else "open loop"
    output_tag = f"{experiment}_vs_ol" + ("_xmask" if args.ol_support == "xmask" else "")
    if args.end_month != "202010":
        output_tag += f"_{month_period}"
    da_file = STATS / f"temporal_stats_DA_paired_{experiment}_{file_period}.nc4"
    OUT.mkdir(parents=True, exist_ok=True)
    ol = load_stats(ol_file)
    da = load_stats(da_file)
    tilecoord = read_tilecoord(str(TILECOORD))
    lon = np.asarray(tilecoord["com_lon"], dtype=float)
    lat = np.asarray(tilecoord["com_lat"], dtype=float)

    results = {}
    summary_rows = []
    for group in GROUPS:
        ol_group, da_group, percent, common_count = common_group_values(ol, da, group)
        results[group.name] = (ol_group, da_group, percent)
        valid = np.isfinite(percent)
        summary_rows.append(
            {
                "system": group.name,
                "common_tiles": int(valid.sum()),
                "aggregate_percent_reduction": aggregate_percent(ol, da, group),
                "tile_mean_percent_reduction": float(np.nanmean(percent)),
                "tile_median_percent_reduction": float(np.nanmedian(percent)),
                "fraction_tiles_improved": float(np.mean(percent[valid] > 0)),
                "common_species_tile_cells": int(common_count.sum()),
            }
        )

    summary_path = OUT / f"{output_tag}_omf_stdv_summary.csv"
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    # Standard relative-reduction figure.
    fig, axes = plt.subplots(
        3,
        2,
        figsize=(11.5, 12.2),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    axes = axes.ravel()
    norm = TwoSlopeNorm(vmin=-PERCENT_LIMIT, vcenter=0.0, vmax=PERCENT_LIMIT)
    for panel, (group, row, ax) in enumerate(zip(GROUPS, summary_rows, axes)):
        base_map(ax)
        percent = results[group.name][2]
        image = scatter_map(ax, lon, lat, percent, "RdBu", norm)
        ax.set_title(f"({chr(97 + panel)}) {group.name}", loc="left", fontweight="bold")
        ax.text(
            0.02,
            0.03,
            f"aggregate {row['aggregate_percent_reduction']:+.2f}%\n"
            f"{row['common_tiles']} common tiles",
            transform=ax.transAxes,
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9, "pad": 2},
        )
    axes[-1].axis("off")
    colorbar_axis = axes[-1].inset_axes((0.03, 0.42, 0.94, 0.10))
    colorbar = fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label(
        "OmF standard-deviation reduction (%)\n"
        f"100 x (OL - {experiment}) / OL\n"
        "positive = lower variability",
        fontsize=9,
    )
    fig.suptitle(
        f"{label} CYGNSS L1 experiment versus {ol_label}, {period_label}\n"
        f"common species-level support, N_data >= {NMIN}",
        fontsize=14,
        y=0.98,
    )
    fig.subplots_adjust(top=0.91, bottom=0.04, hspace=0.24, wspace=0.10)
    percent_path = OUT / f"{output_tag}_omf_stdv_percent_maps.png"
    fig.savefig(percent_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Side-by-side values make large percentages and support differences auditable.
    fig = plt.figure(figsize=(14.5, 17.0))
    for row_index, group in enumerate(GROUPS):
        ol_group, da_group, percent = results[group.name]
        finite = np.concatenate((ol_group[np.isfinite(ol_group)], da_group[np.isfinite(da_group)]))
        lo, hi = np.nanpercentile(finite, (2, 98))
        absolute_norm = Normalize(vmin=lo, vmax=hi)
        for column, (values, title, cmap, map_norm) in enumerate(
            (
                (ol_group, "OL OmF stddev", "viridis", absolute_norm),
                (da_group, f"{label} OmF stddev", "viridis", absolute_norm),
                (percent, "Relative reduction", "RdBu", norm),
            )
        ):
            ax = fig.add_subplot(5, 3, row_index * 3 + column + 1, projection=ccrs.PlateCarree())
            base_map(ax)
            image = scatter_map(ax, lon, lat, values, cmap, map_norm)
            if row_index == 0:
                ax.set_title(title, fontweight="bold", fontsize=11)
            if column == 0:
                ax.text(
                    -0.08,
                    0.5,
                    group.name,
                    transform=ax.transAxes,
                    rotation=90,
                    va="center",
                    ha="center",
                    fontweight="bold",
                )
            colorbar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.045, pad=0.04)
            colorbar.set_label("%" if column == 2 else group.units, fontsize=8)
            colorbar.ax.tick_params(labelsize=7)
    fig.suptitle(
        f"OmF standard deviation: {ol_label} and {experiment} CYGNSS L1 experiment\n"
        f"January-October 2020; common species-level support, N_data >= {NMIN}",
        fontsize=14,
        y=0.995,
    )
    fig.subplots_adjust(top=0.95, bottom=0.02, hspace=0.32, wspace=0.12)
    comparison_path = OUT / f"{output_tag}_omf_stdv_comparison_maps.png"
    fig.savefig(comparison_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {percent_path.relative_to(PROJECT)}")
    print(f"wrote {comparison_path.relative_to(PROJECT)}")
    print(f"wrote {summary_path.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
