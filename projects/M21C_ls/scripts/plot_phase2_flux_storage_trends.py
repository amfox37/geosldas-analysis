#!/usr/bin/env python3
"""Plot focused Phase 2 ET, runoff, and storage trend/transition results."""

from __future__ import annotations

import shutil
import warnings
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from plot_trends_breakpoints_manuscript_figures import (
    MAP_LAT_MIN,
    add_map_background,
    build_tile_grid,
    load_periods,
    nice_limit,
    panel_label,
    segmented_diverging_scale,
    shade_periods,
    tile_values_to_grid,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = (
    PROJECT_ROOT / "output" / "trends_breakpoints" / "phase2_flux_storage"
)
TREND_DIR = RESULT_DIR / "tile_trends"
OUTPUT_DIR = RESULT_DIR / "figures"
DOC_FIG_DIR = PROJECT_ROOT / "docs" / "trends_breakpoints_report_figures"

VARIABLES = {
    "EVLAND": {
        "label": "Evapotranspiration",
        "units": r"kg m$^{-2}$ month$^{-1}$ yr$^{-1}$",
    },
    "TOTAL_RUNOFF": {
        "label": "Total runoff",
        "units": r"kg m$^{-2}$ month$^{-1}$ yr$^{-1}$",
    },
    "TWLAND": {
        "label": "Total land-water storage",
        "units": r"kg m$^{-2}$ yr$^{-1}$",
    },
}
SERIES = ("ol", "da", "delta")
SERIES_LABELS = {"ol": "OL", "da": "DA", "delta": "DA - OL"}


def _save_review_figure(fig: plt.Figure, name: str) -> tuple[Path, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_FIG_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / f"{name}.png"
    review = DOC_FIG_DIR / output.name
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in create_collection",
            category=RuntimeWarning,
        )
        fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    shutil.copy2(output, review)
    plt.close(fig)
    return output, review


def _load_trend(variable: str, series: str) -> xr.Dataset:
    path = TREND_DIR / f"{variable}_{series}_valid_land_trend_statistics.nc"
    if not path.exists():
        raise FileNotFoundError(path)
    with xr.open_dataset(path) as dataset:
        return dataset[["slope", "significant_fdr"]].load()


def plot_trend_maps() -> tuple[Path, Path]:
    fig = plt.figure(figsize=(13.6, 10.8))
    grid_spec = fig.add_gridspec(
        6,
        3,
        height_ratios=[1.0, 0.09, 1.0, 0.09, 1.0, 0.09],
        left=0.055,
        right=0.992,
        bottom=0.055,
        top=0.97,
        hspace=0.12,
        wspace=0.02,
    )
    plate = ccrs.PlateCarree()
    panel_index = 0
    for row_index, (variable, cfg) in enumerate(VARIABLES.items()):
        fields = {series: _load_trend(variable, series) for series in SERIES}
        lat = np.asarray(fields["ol"]["lat"].values)
        lon = np.asarray(fields["ol"]["lon"].values)
        tile_grid = build_tile_grid(lat, lon)
        visible = lat >= MAP_LAT_MIN
        state_values = np.concatenate(
            [np.asarray(fields[series]["slope"].values)[visible] for series in ("ol", "da")]
        )
        delta_values = np.asarray(fields["delta"]["slope"].values)[visible]
        state_vmax = nice_limit(state_values)
        delta_vmax = nice_limit(delta_values)
        state_cmap, state_norm, _ = segmented_diverging_scale(state_vmax)
        delta_cmap, delta_norm, _ = segmented_diverging_scale(delta_vmax)
        axes: list[plt.Axes] = []
        meshes = []

        for column, series in enumerate(SERIES):
            ax = fig.add_subplot(grid_spec[2 * row_index, column], projection=ccrs.Robinson())
            axes.append(ax)
            add_map_background(ax)
            cmap, norm = (
                (state_cmap, state_norm)
                if series != "delta"
                else (delta_cmap, delta_norm)
            )
            slope = np.asarray(fields[series]["slope"].values)
            mesh = ax.pcolormesh(
                tile_grid["lon_edges"],
                tile_grid["lat_edges"],
                tile_values_to_grid(slope, tile_grid),
                transform=plate,
                cmap=cmap,
                norm=norm,
                shading="auto",
                rasterized=True,
                zorder=1,
            )
            meshes.append(mesh)
            significant = (
                np.asarray(fields[series]["significant_fdr"].values, dtype=bool)
                & visible
            )
            ax.scatter(
                lon[significant],
                lat[significant],
                s=0.42,
                c="black",
                alpha=0.38,
                linewidths=0,
                transform=plate,
                rasterized=True,
                zorder=2,
            )
            if row_index == 0:
                ax.set_title(SERIES_LABELS[series], fontweight="bold", pad=5)
            panel_label(ax, f"({chr(97 + panel_index)})")
            panel_index += 1

        axes[0].text(
            -0.075,
            0.5,
            cfg["label"],
            transform=axes[0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )
        state_axis = fig.add_subplot(grid_spec[2 * row_index + 1, :2])
        delta_axis = fig.add_subplot(grid_spec[2 * row_index + 1, 2])
        state_bar = fig.colorbar(
            meshes[0],
            cax=state_axis,
            orientation="horizontal",
            extend="both",
            ticks=np.linspace(-state_vmax, state_vmax, 5),
        )
        delta_bar = fig.colorbar(
            meshes[2],
            cax=delta_axis,
            orientation="horizontal",
            extend="both",
            ticks=np.linspace(-delta_vmax, delta_vmax, 5),
        )
        state_bar.set_label(f"OL and DA trend ({cfg['units']})", labelpad=2)
        delta_bar.set_label(f"DA - OL trend ({cfg['units']})", labelpad=2)
        state_bar.ax.tick_params(length=2.5, pad=1)
        delta_bar.ax.tick_params(length=2.5, pad=1)

    return _save_review_figure(fig, "phase2_flux_storage_trend_maps")


def plot_transition_series() -> tuple[Path, Path]:
    periods = load_periods()
    with xr.open_dataset(RESULT_DIR / "interrupted_series_monthly.nc") as source:
        monthly = source.load()
    detections = pd.read_csv(RESULT_DIR / "changepoint_detections.csv", parse_dates=["break_date"])
    accepted = detections[detections["accepted_detection"].astype(bool)]
    coefficients = pd.read_csv(RESULT_DIR / "interrupted_series_coefficients.csv")

    fig, axes = plt.subplots(3, 1, figsize=(12.7, 8.2), sharex=True)
    colors = {"EVLAND": "#238b7b", "TOTAL_RUNOFF": "#15616d", "TWLAND": "#b44b4b"}
    time = pd.DatetimeIndex(monthly.time.values)
    for index, (ax, (variable, cfg)) in enumerate(zip(axes, VARIABLES.items())):
        series_id = f"{variable}_delta_valid_land__delta"
        shade_periods(ax, periods)
        observed = pd.Series(
            monthly["observed"].sel(series=series_id).values.astype(float), index=time
        )
        fitted = pd.Series(
            monthly["fitted"].sel(series=series_id).values.astype(float), index=time
        )
        ax.plot(
            time,
            observed.rolling(12, center=True, min_periods=6).mean(),
            color=colors[variable],
            linewidth=1.25,
            label="12-month mean",
            zorder=3,
        )
        ax.plot(
            time,
            fitted.rolling(12, center=True, min_periods=6).mean(),
            color="black",
            linewidth=1.0,
            label="Interrupted-series fit",
            zorder=4,
        )
        current_breaks = accepted[
            (accepted["variable"] == variable)
            & (accepted["source_series"] == "delta")
        ]
        for break_date in current_breaks["break_date"]:
            ax.axvline(break_date, color="#d9902f", linewidth=1.2, zorder=5)
        p6 = coefficients[
            (coefficients["variable"] == variable)
            & (coefficients["source_series"] == "delta")
            & (coefficients["coefficient"] == "level_change_P6")
        ].iloc[0]
        ax.set_title(
            f"{cfg['label']} DA - OL   |   P6 level change {p6['estimate']:+.3g}; "
            f"q={p6['q_value']:.3f}",
            loc="left",
            fontweight="bold",
            pad=5,
        )
        ax.set_ylabel(cfg["units"].replace(" yr$^{-1}$", ""))
        ax.axhline(0, color="0.35", linewidth=0.7, zorder=2)
        ax.grid(axis="y", color="0.88", linewidth=0.7)
        panel_label(ax, f"({chr(97 + index)})")
    axes[0].legend(loc="lower left", ncol=2, frameon=False)
    axes[-1].set_xlabel("Year")
    return _save_review_figure(fig, "phase2_flux_storage_transition_series")


def main() -> int:
    for output, review in (plot_trend_maps(), plot_transition_series()):
        print(f"Wrote {output}")
        print(f"Wrote {review}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
