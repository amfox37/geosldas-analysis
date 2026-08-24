#!/usr/bin/env python3
"""Render final M21C trend figures (Fig. 16 and the supplementary trend maps).

This is a plotting-only workflow. It reads the accepted Phase 1 trend,
products and does not refit any analysis. Figure 17 is generated separately by
`plot_regional_rzmc_periods.py`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import xarray as xr


PROJECT = Path(__file__).resolve().parents[1]
REPO = PROJECT.parents[1]
RESULT_DIR = PROJECT / "output" / "trends_breakpoints"
OUTPUT_DIR = PROJECT / "output" / "paper_figures"
DOC_FIG_DIR = PROJECT / "docs" / "paper_figures"
REPORT_PATH = PROJECT / "docs" / "trends_breakpoints_manuscript_figure_report.md"
REGISTRY_PATH = PROJECT / "config" / "observing_system_registry.json"

MAP_LAT_MIN = -60.0
LAND_FACE = "0.90"
ROBUST_PERCENTILE = 98.0
ZERO_HALF_WIDTH_FRACTION_OF_FULL_RANGE = 0.02
SERIES_ORDER = ("ol", "da", "delta")
SERIES_LABELS = {"ol": "OL", "da": "DA", "delta": "DA - OL"}
MAIN_TREND_DELTA_SCALES = {"RZMC": "separate", "FRLANDSNO": "state"}
SNOW_TREND_DELTA_SCALES = {"SNOMASLAND": "state", "SNODPLAND": "state"}

PERIOD_COLORS = {
    "P1": "#f2efe6",
    "P2": "#e7f0f7",
    "P3": "#eef4e8",
    "P4": "#f7eadf",
    "P5": "#ebe7f3",
    "P6": "#e8f3f1",
    "P7": "#f5edf3",
    "P8": "#edf0f5",
    "P9": "#f2f0e8",
}

VARIABLES = {
    "PRECTOTCORRLAND": {
        "mask": "valid_land",
        "label": "Precipitation",
        "units": r"kg m$^{-2}$ month$^{-1}$ yr$^{-1}$",
        "multiplier": 1.0,
    },
    "SFMC": {
        "mask": "valid_land",
        "label": "Surface soil moisture",
        "units": r"m$^3$ m$^{-3}$ yr$^{-1}$",
        "multiplier": 1000.0,
    },
    "RZMC": {
        "mask": "valid_land",
        "label": "Root-zone soil moisture",
        "units": r"m$^3$ m$^{-3}$ yr$^{-1}$",
        "multiplier": 1000.0,
    },
    "SNOMASLAND": {
        "mask": "seasonal_snow",
        "label": "Snow mass",
        "units": r"kg m$^{-2}$ yr$^{-1}$",
        "multiplier": 1.0,
    },
    "SNODPLAND": {
        "mask": "seasonal_snow",
        "label": "Snow depth",
        "units": r"m yr$^{-1}$",
        "multiplier": 1000.0,
    },
    "FRLANDSNO": {
        "mask": "seasonal_snow",
        "label": "Snow-cover fraction",
        "units": r"SCF yr$^{-1}$",
        "multiplier": 1000.0,
    },
}

P6_SERIES = {
    "RZMC_delta_valid_land__delta": {
        "label": "Valid land",
        "group": "rzmc",
        "expected": (0.001021, 0.000488, 0.001558),
    },
    "RZMC_delta_warm_snowfree_monthly__delta": {
        "label": "Warm snow-free",
        "group": "rzmc",
        "expected": (0.001267, 0.000436, 0.002037),
    },
    "SFMC_INC_RMS_value_valid_land__value": {
        "label": "SFMC increment RMS",
        "group": "increment",
        "expected": (0.001107, 0.000677, 0.001505),
    },
    "RZMC_INC_RMS_value_valid_land__value": {
        "label": "RZMC increment RMS",
        "group": "increment",
        "expected": (0.000120, 0.000077, 0.000162),
    },
    "soil_water_abs_activity_value_valid_land__value": {
        "label": "Absolute activity",
        "group": "activity",
        "expected": (10.612, 7.056, 14.715),
    },
    "soil_water_net_approx_value_valid_land__value": {
        "label": "Signed net correction",
        "group": "activity",
        "expected": (1.530, 0.589, 2.467),
    },
}

TIMESERIES = {
    "RZMC_delta_valid_land__delta": ("RZMC DA - OL", "#c44e52"),
    "soil_water_abs_activity_value_valid_land__value": ("Soil-water absolute activity", "#7b6ea8"),
    "SFMC_INC_RMS_value_valid_land__value": ("SFMC increment RMS", "#238b7b"),
    "RZMC_INC_RMS_value_valid_land__value": ("RZMC increment RMS", "#d9a32f"),
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 9.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
)


@dataclass(frozen=True)
class FigureAsset:
    name: str
    png: Path
    pdf: Path
    caption: str
    sources: tuple[str, ...]
    provenance: str


def load_periods() -> pd.DataFrame:
    registry = json.loads(REGISTRY_PATH.read_text())
    periods = pd.DataFrame(registry["fine_periods"])
    periods["start"] = pd.to_datetime(periods["start"])
    periods["end"] = pd.to_datetime(periods["end"])
    if periods["period_id"].tolist() != [f"P{i}" for i in range(1, 10)]:
        raise AssertionError("Observing-system registry does not contain ordered P1-P9 periods")
    if periods.loc[periods["period_id"] == "P6", "start"].item() != pd.Timestamp("2015-04-01"):
        raise AssertionError("P6 does not begin in April 2015")
    for left, right in zip(periods.iloc[:-1].itertuples(), periods.iloc[1:].itertuples()):
        if left.end + pd.Timedelta(days=1) != right.start:
            raise AssertionError(f"Registry periods {left.period_id} and {right.period_id} are not contiguous")
    return periods


def centers_to_edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    middle = 0.5 * (centers[:-1] + centers[1:])
    return np.concatenate(
        ([centers[0] - 0.5 * np.diff(centers)[0]], middle, [centers[-1] + 0.5 * np.diff(centers)[-1]])
    )


def build_tile_grid(lat: np.ndarray, lon: np.ndarray, decimals: int = 6) -> dict[str, np.ndarray | tuple[int, int]]:
    lat_round = np.round(np.asarray(lat), decimals)
    lon_round = np.round(np.asarray(lon), decimals)
    lat_unique = np.unique(lat_round)
    lon_unique = np.unique(lon_round)
    i_lat = np.searchsorted(lat_unique, lat_round)
    i_lon = np.searchsorted(lon_unique, lon_round)
    pair_code = i_lat.astype(np.int64) * lon_unique.size + i_lon
    if np.unique(pair_code).size != pair_code.size:
        raise ValueError("Duplicate rounded tile centers prevent grid reconstruction")
    lon_edges, lat_edges = np.meshgrid(centers_to_edges(lon_unique), centers_to_edges(lat_unique))
    return {
        "i_lat": i_lat,
        "i_lon": i_lon,
        "shape": (lat_unique.size, lon_unique.size),
        "lon_edges": lon_edges,
        "lat_edges": lat_edges,
    }


def tile_values_to_grid(values: np.ndarray, grid: dict[str, object]) -> np.ndarray:
    output = np.full(grid["shape"], np.nan, dtype=np.float32)
    output[grid["i_lat"], grid["i_lon"]] = np.asarray(values, dtype=np.float32)
    return output


def nice_limit(values: np.ndarray, percentile: float = ROBUST_PERCENTILE) -> float:
    finite = np.abs(np.asarray(values, dtype=float))
    finite = finite[np.isfinite(finite)]
    raw = float(np.percentile(finite, percentile))
    if not np.isfinite(raw) or raw <= 0:
        return 1.0
    power = 10.0 ** np.floor(np.log10(raw))
    return float(np.ceil(raw / power * 5.0) / 5.0 * power)


def segmented_diverging_scale(vmax: float) -> tuple[ListedColormap, BoundaryNorm, np.ndarray]:
    zero_half_width = 2.0 * vmax * ZERO_HALF_WIDTH_FRACTION_OF_FULL_RANGE
    bounds = np.concatenate(
        (
            np.linspace(-vmax, -zero_half_width, 6),
            [zero_half_width],
            np.linspace(zero_half_width, vmax, 6)[1:],
        )
    )
    base = plt.get_cmap("RdBu_r")
    colors = (
        [base(value) for value in np.linspace(0.05, 0.40, 5)]
        + [(1, 1, 1, 1)]
        + [base(value) for value in np.linspace(0.60, 0.95, 5)]
    )
    cmap = ListedColormap(colors)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap, BoundaryNorm(bounds, cmap.N, clip=False), bounds


def load_trend_fields(variable: str) -> tuple[np.ndarray, np.ndarray, dict[str, dict[str, np.ndarray]]]:
    cfg = VARIABLES[variable]
    fields: dict[str, dict[str, np.ndarray]] = {}
    reference_lat = reference_lon = None
    for series in SERIES_ORDER:
        path = RESULT_DIR / f"{variable}_{series}_{cfg['mask']}_trend_statistics.nc"
        with xr.open_dataset(path) as dataset:
            lat = np.asarray(dataset["lat"].values)
            lon = np.asarray(dataset["lon"].values)
            if reference_lat is None:
                reference_lat, reference_lon = lat, lon
            elif not (np.array_equal(reference_lat, lat) and np.array_equal(reference_lon, lon)):
                raise AssertionError(f"Tile coordinates differ for {variable}/{series}")
            fields[series] = {
                "slope": np.asarray(dataset["slope"].values, dtype=float),
                "significant": np.asarray(dataset["significant_fdr"].values, dtype=bool),
            }
            if dataset["significant_fdr"].attrs.get("long_name", "").find("FDR") < 0:
                raise AssertionError(f"{path.name} does not identify significant_fdr as FDR inference")
    assert reference_lat is not None and reference_lon is not None
    return reference_lat, reference_lon, fields


def add_map_background(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor=LAND_FACE, edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="white", edgecolor="none", zorder=0)
    ax.coastlines(resolution="110m", linewidth=0.4, color="0.30", zorder=3)
    ax.set_extent([-180, 180, MAP_LAT_MIN, 90], crs=ccrs.PlateCarree())


def panel_label(ax: plt.Axes, text: str, x: float = 0.015, y: float = 0.985) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontweight="bold",
        fontsize=9.5,
        zorder=8,
        bbox={"boxstyle": "square,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.80},
    )


def display_units(variable: str) -> str:
    cfg = VARIABLES[variable]
    if cfg["multiplier"] == 1000.0:
        return rf"Trend ($\times 10^{{-3}}$ {cfg['units']})"
    return f"Trend ({cfg['units']})"


def save_figure(fig: plt.Figure, basename: str) -> tuple[Path, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_FIG_DIR.mkdir(parents=True, exist_ok=True)
    png = OUTPUT_DIR / f"{basename}.png"
    pdf = OUTPUT_DIR / f"{basename}.pdf"
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in create_collection",
            category=RuntimeWarning,
        )
        fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    shutil.copy2(png, DOC_FIG_DIR / png.name)
    plt.close(fig)
    return png, pdf


def plot_trend_rows(
    row_variables: list[str],
    basename: str,
    delta_scales: dict[str, str] | None = None,
) -> tuple[Path, Path]:
    delta_scales = {} if delta_scales is None else dict(delta_scales)
    invalid = {variable: mode for variable, mode in delta_scales.items() if mode not in {"separate", "state"}}
    if invalid:
        raise ValueError(f"delta_scales values must be 'separate' or 'state': {invalid}")
    n_rows = len(row_variables)
    fig = plt.figure(figsize=(13.6, 3.45 * n_rows + 0.55))
    height_ratios: list[float] = []
    for _ in row_variables:
        height_ratios.extend([1.0, 0.095])
    grid_spec = fig.add_gridspec(
        2 * n_rows,
        3,
        height_ratios=height_ratios,
        left=0.055,
        right=0.992,
        bottom=0.065,
        top=0.955,
        hspace=0.12,
        wspace=0.02,
    )
    plate = ccrs.PlateCarree()
    panel_index = 0

    for row_index, variable in enumerate(row_variables):
        cfg = VARIABLES[variable]
        lat, lon, fields = load_trend_fields(variable)
        tile_grid = build_tile_grid(lat, lon)
        visible = lat >= MAP_LAT_MIN
        multiplier = float(cfg["multiplier"])
        state_values = np.concatenate(
            [fields[name]["slope"][visible] * multiplier for name in ("ol", "da")]
        )
        delta_values = fields["delta"]["slope"][visible] * multiplier
        state_vmax = nice_limit(state_values)
        row_delta_scale = delta_scales.get(variable, "separate")
        delta_vmax = state_vmax if row_delta_scale == "state" else nice_limit(delta_values)
        state_cmap, state_norm, _ = segmented_diverging_scale(state_vmax)
        delta_cmap, delta_norm, _ = segmented_diverging_scale(delta_vmax)
        axes: list[plt.Axes] = []
        meshes = []

        for column, series in enumerate(SERIES_ORDER):
            ax = fig.add_subplot(grid_spec[2 * row_index, column], projection=ccrs.Robinson())
            axes.append(ax)
            add_map_background(ax)
            cmap, norm = (state_cmap, state_norm) if series != "delta" else (delta_cmap, delta_norm)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="invalid value encountered in create_collection",
                    category=RuntimeWarning,
                )
                mesh = ax.pcolormesh(
                    tile_grid["lon_edges"],
                    tile_grid["lat_edges"],
                    tile_values_to_grid(fields[series]["slope"] * multiplier, tile_grid),
                    transform=plate,
                    cmap=cmap,
                    norm=norm,
                    shading="auto",
                    rasterized=True,
                    zorder=1,
                )
            meshes.append(mesh)
            significant = fields[series]["significant"] & visible
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
        cax_state = fig.add_subplot(grid_spec[2 * row_index + 1, :2])
        cax_delta = fig.add_subplot(grid_spec[2 * row_index + 1, 2])
        state_bar = fig.colorbar(
            meshes[0],
            cax=cax_state,
            orientation="horizontal",
            extend="both",
            ticks=np.linspace(-state_vmax, state_vmax, 5),
        )
        delta_bar = fig.colorbar(
            meshes[2],
            cax=cax_delta,
            orientation="horizontal",
            extend="both",
            ticks=np.linspace(-delta_vmax, delta_vmax, 5),
        )
        unit_label = display_units(variable).replace("Trend", "trend", 1)
        state_bar.set_label(f"OL and DA {unit_label}", labelpad=2)
        delta_bar.set_label(f"DA - OL {unit_label}", labelpad=2)
        state_bar.ax.tick_params(length=2.5, pad=1)
        delta_bar.ax.tick_params(length=2.5, pad=1)

    return save_figure(fig, basename)


def validate_trend_products() -> dict[str, object]:
    tile = pd.read_csv(RESULT_DIR / "phase1_state_trend_tile_summary.csv")
    domain = pd.read_csv(RESULT_DIR / "phase1_state_trend_domain_summary.csv")
    pairs = pd.read_csv(RESULT_DIR / "phase1_state_trend_pair_comparison.csv")

    expected_counts = {
        "PRECTOTCORRLAND": {"ol": 3719, "da": 3726, "delta": 0},
        "SFMC": {"ol": 6992, "da": 8966, "delta": 1412},
        "RZMC": {"ol": 4909, "da": 10329, "delta": 7892},
        "SNOMASLAND": {"ol": 12, "da": 5, "delta": 0},
        "SNODPLAND": {"ol": 7, "da": 6, "delta": 0},
    }
    for variable, series_counts in expected_counts.items():
        frame = tile[tile["variable"] == variable].set_index("series")
        actual = frame["n_significant_fdr"].to_dict()
        if any(int(actual[name]) != count for name, count in series_counts.items()):
            raise AssertionError(f"Trend FDR counts changed for {variable}: {actual}")

    rzmc_delta = tile[(tile["variable"] == "RZMC") & (tile["series"] == "delta")].iloc[0]
    if (int(rzmc_delta["n_significant_positive"]), int(rzmc_delta["n_significant_negative"])) != (7267, 625):
        raise AssertionError("RZMC DA-OL significant sign counts changed")

    scf = domain[domain["variable"] == "FRLANDSNO"].set_index("series")
    scf_slopes = scf["slope_per_year"].to_dict()
    expected_scf = {"ol": -0.000554, "da": -0.000549, "delta": -0.000006}
    for name, expected in expected_scf.items():
        if not np.isclose(scf_slopes[name], expected, atol=5.0e-7):
            raise AssertionError(f"SCF domain slope changed for {name}: {scf_slopes[name]}")

    precip = pairs[pairs["variable"] == "PRECTOTCORRLAND"].iloc[0]
    if (
        int(precip["n_ol_da_significant_overlap"]) != 3603
        or int(precip["n_overlap_same_direction"]) != 3603
        or int(precip["n_overlap_opposite_direction"]) != 0
        or not np.isclose(precip["ol_da_slope_correlation"], 0.9998, atol=5.0e-5)
    ):
        raise AssertionError("Precipitation trend-control cross-check changed")
    return {"tile_summary": tile, "domain_summary": domain, "pair_summary": pairs}


def shade_periods(ax: plt.Axes, periods: pd.DataFrame) -> None:
    for row in periods.itertuples():
        ax.axvspan(row.start, row.end + pd.Timedelta(days=1), color=PERIOD_COLORS[row.period_id], zorder=0)
        center = row.start + (row.end - row.start) / 2
        period_label = "P6 (SMAP Tb)" if row.period_id == "P6" else row.period_id
        ax.text(
            center,
            0.975,
            period_label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8.2,
            fontweight="bold",
            zorder=5,
        )
    for row in periods.iloc[1:].itertuples():
        is_p6 = row.period_id == "P6"
        ax.axvline(
            row.start,
            color="#8b1e1e" if is_p6 else "0.55",
            linewidth=1.35 if is_p6 else 0.65,
            linestyle="--",
            zorder=2,
        )


def build_assets() -> list[FigureAsset]:
    periods = load_periods()
    validate_trend_products()

    fig_x = plot_trend_rows(
        ["RZMC", "FRLANDSNO"],
        "fig16_longterm_rzmc_scf_trends",
        delta_scales=MAIN_TREND_DELTA_SCALES,
    )
    precip = plot_trend_rows(["PRECTOTCORRLAND"], "figS05_precipitation_trends")
    sfmc = plot_trend_rows(["SFMC"], "figS06_sfmc_trends")
    snow = plot_trend_rows(
        ["SNOMASLAND", "SNODPLAND"],
        "figS07_snow_mass_depth_trends",
        delta_scales=SNOW_TREND_DELTA_SCALES,
    )

    assets = [
        FigureAsset(
            "Figure 16: long-term state trends",
            *fig_x,
            caption=(
                "Figure 16. Long-term June 2000-May 2024 trends in (a-c) root-zone soil moisture (RZMC) and "
                "(d-f) snow-cover fraction (SCF) for the open-loop (OL), data-assimilation (DA), "
                "and paired DA-OL series. Trends are exact Theil-Sen slopes after trend-preserving "
                "removal of the calendar-month climatology. Black stippling denotes trends significant "
                "after Benjamini-Hochberg false-discovery-rate control at 0.05. RZMC uses the valid-land "
                "domain and SCF the Northern Hemisphere seasonal-snow domain. The DA-OL panels show the trend of the paired "
                "DA-OL series rather than the difference between independently estimated OL and DA trends. "
                "Regional RZMC trends are modified substantially by DA, whereas OL and DA exhibit nearly "
                "identical long-term SCF trends."
            ),
            sources=(
                "RZMC_{ol,da,delta}_valid_land_trend_statistics.nc",
                "FRLANDSNO_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
            ),
            provenance="Exact production slope; mapped significance is significant_fdr; Robinson projection; 60 S cutoff.",
        ),
        FigureAsset(
            "Supplemental Figure S5: precipitation trends",
            *precip,
            caption=(
                "Figure S5. Long-term precipitation trends for OL, DA, and the paired DA-OL series on valid land. "
                "Black stippling denotes production FDR significance. The common OL/DA pattern and null "
                "DA-OL result provide a forcing-control check."
            ),
            sources=("PRECTOTCORRLAND_{ol,da,delta}_valid_land_trend_statistics.nc",),
            provenance="Exact production slope and significant_fdr on valid land.",
        ),
        FigureAsset(
            "Supplemental Figure S6: SFMC trends",
            *sfmc,
            caption=(
                "Figure S6. Long-term surface soil-moisture trends for OL, DA, and the paired DA-OL series on valid land. "
                "Black stippling denotes production FDR significance."
            ),
            sources=("SFMC_{ol,da,delta}_valid_land_trend_statistics.nc",),
            provenance="Exact production slope and significant_fdr on valid land.",
        ),
        FigureAsset(
            "Supplemental Figure S7: snow mass and depth trends",
            *snow,
            caption=(
                "Figure S7. Long-term (a-c) snow-mass and (d-f) snow-depth trends for OL, DA, and the paired DA-OL "
                "series on the production Northern Hemisphere seasonal-snow mask. Black stippling denotes production FDR significance."
            ),
            sources=(
                "SNOMASLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
                "SNODPLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
            ),
            provenance="Exact production slope and significant_fdr on the Northern Hemisphere seasonal-snow mask.",
        ),
    ]
    write_report(assets, periods)
    return assets


def write_report(assets: list[FigureAsset], periods: pd.DataFrame) -> None:
    lines = [
        "# Trends And Observing-System Manuscript Figure Report",
        "",
        "These figures are plotting-only products generated from the accepted Phase 1 trend outputs. No analysis was recomputed.",
        "",
        "## Outputs And Captions",
        "",
    ]
    for asset in assets:
        lines.extend(
            [
                f"### {asset.name}",
                "",
                f"- PNG: `{asset.png.relative_to(REPO)}`",
                f"- PDF: `{asset.pdf.relative_to(REPO)}`",
                f"- Tracked review PNG: `{(DOC_FIG_DIR / asset.png.name).relative_to(REPO)}`",
                f"- Sources: {', '.join(f'`{source}`' for source in asset.sources)}",
                f"- Provenance: {asset.provenance}",
                "",
                f"**Draft caption.** {asset.caption}",
                "",
            ]
        )
    p9_start = periods.loc[periods["period_id"] == "P9", "start"].item().date()
    lines.extend(
        [
            "## Validation",
            "",
            "- RZMC FDR-significant counts reproduce 4,909 OL, 10,329 DA, and 7,892 paired DA-OL tiles; significant DA-OL signs reproduce 7,267 positive and 625 negative.",
            "- SCF area-weighted slopes reproduce -0.000554 yr-1 (OL), -0.000549 yr-1 (DA), and -0.000006 yr-1 (DA-OL).",
            "- Precipitation reproduces 3,719 OL and 3,726 DA significant tiles, 3,603 same-sign overlaps, slope correlation 0.9998, and zero DA-OL FDR tiles.",
            "- SFMC reproduces 6,992 OL, 8,966 DA, and 1,412 DA-OL FDR-significant tiles.",
            "- Snow mass reproduces 12 OL, 5 DA, and zero DA-OL FDR tiles; snow depth reproduces 7 OL, 6 DA, and zero DA-OL FDR tiles.",
            "- All maps use exact production Theil-Sen slopes and `significant_fdr`; pointwise confidence-interval exclusion is not used for mapped inference.",
            "",
            "## Plotting Choices",
            "",
            "- Maps follow the existing report convention: Robinson projection, 60 S cutoff, grey land, thin coastlines, segmented RdBu_r scales centered on a white zero bin, and black stippling.",
            "- OL and DA share a symmetric color scale within each row. Figure 16 retains a separate RZMC DA-OL scale but places SCF DA-OL on the OL/DA scale; snow mass/depth likewise retain the OL/DA scale so negligible differences are not visually exaggerated.",
            "",
            "## Discrepancies",
            "",
            "- No regenerated numerical result disagreed with `m21c_trends_breakpoints_report.md`.",
            f"- The task summary says P9 begins in November 2021, but the authoritative registry starts P9 on {p9_start} after P8 ends 2021-11-30. The figures use the registry date.",
            "",
            "PNGs are exported at 300 DPI. PDFs retain vector text, coastlines, markers, and intervals; dense map fields and stippling are rasterized within the vector PDF.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true", help="Validate accepted products without rendering")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validate_trend_products()
    if args.validate_only:
        print("Trend manuscript figure contracts: PASS")
        return
    assets = build_assets()
    print(f"Generated {len(assets)} manuscript figures")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Figure report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
