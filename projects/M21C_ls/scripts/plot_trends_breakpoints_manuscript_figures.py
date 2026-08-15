#!/usr/bin/env python3
"""Render final M21C trend and observing-system-transition figures.

This is a plotting-only workflow. It reads the accepted Phase 1 trend,
interrupted-series, and changepoint products and does not refit any analysis.
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
    delta_scale: str = "separate",
) -> tuple[Path, Path]:
    if delta_scale not in {"separate", "state"}:
        raise ValueError("delta_scale must be 'separate' or 'state'")
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
        delta_vmax = state_vmax if delta_scale == "state" else nice_limit(delta_values)
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


def validate_transition_products(periods: pd.DataFrame) -> dict[str, object]:
    coefficients = pd.read_csv(RESULT_DIR / "phase1_interrupted_series_coefficients.csv")
    selected = coefficients[
        (coefficients["coefficient"] == "level_change_P6")
        & coefficients["series_id"].isin(P6_SERIES)
    ].copy()
    if set(selected["series_id"]) != set(P6_SERIES) or len(selected) != len(P6_SERIES):
        raise AssertionError("P6 coefficient selection is incomplete or duplicated")
    selected = selected.set_index("series_id")
    for series_id, cfg in P6_SERIES.items():
        row = selected.loc[series_id]
        actual = row[["estimate", "ci_low_bootstrap", "ci_high_bootstrap"]].to_numpy(dtype=float)
        tolerance = 5.0e-4 if row["estimate_units"] == "kg m-2" else 5.0e-7
        if not np.allclose(actual, cfg["expected"], atol=tolerance):
            raise AssertionError(f"P6 coefficient changed for {series_id}: {actual}")
        if not bool(row["significant_fdr"]):
            raise AssertionError(f"P6 coefficient no longer survives boundary-family FDR: {series_id}")

    with xr.open_dataset(RESULT_DIR / "phase1_changepoint_monthly.nc") as source:
        monthly = source.load()
    if pd.DatetimeIndex(monthly["time"].values)[0] != periods.iloc[0]["start"]:
        raise AssertionError("Monthly transition series does not begin at P1")
    if pd.DatetimeIndex(monthly["time"].values)[-1] != periods.iloc[-1]["end"].replace(day=1):
        raise AssertionError("Monthly transition series does not end in the final P9 month")
    if not set(TIMESERIES).issubset(set(map(str, monthly["series_id"].values))):
        raise AssertionError("One or more Figure 17 production monthly series are missing")
    return {"coefficients": coefficients, "selected": selected, "monthly": monthly}


def validate_breakpoint_products(coefficients: pd.DataFrame) -> dict[str, object]:
    detections = pd.read_csv(RESULT_DIR / "phase1_changepoint_detections.csv", parse_dates=["break_date"])
    comparison = pd.read_csv(
        RESULT_DIR / "phase1_changepoint_boundary_comparison.csv",
        parse_dates=["boundary_date", "detected_date"],
    )
    accepted = detections[detections["accepted_detection"].astype(bool)].copy()
    if len(accepted) != 37:
        raise AssertionError(f"Expected 37 accepted breaks, found {len(accepted)}")
    role_counts = accepted["series_role"].value_counts().to_dict()
    if role_counts.get("primary_estimand", 0) != 37 or (accepted["source_series"].isin(["ol", "da"])).any():
        raise AssertionError(f"Accepted changepoint roles changed: {role_counts}")
    family_counts = accepted["source_series"].value_counts().to_dict()
    if family_counts.get("delta", 0) != 20 or family_counts.get("value", 0) != 17:
        raise AssertionError(f"Accepted changepoint source families changed: {family_counts}")

    april = accepted[accepted["break_date"] == pd.Timestamp("2015-04-01")]
    if april["series_id"].nunique() != 10:
        raise AssertionError(f"Expected ten primary April 2015 breaks, found {april['series_id'].nunique()}")
    significant_p6 = coefficients[
        (coefficients["coefficient"] == "level_change_P6")
        & coefficients["series_id"].isin(april["series_id"])
        & coefficients["significant_fdr"].astype(bool)
    ]
    if significant_p6["series_id"].nunique() != 9:
        raise AssertionError("April 2015 accepted-break/P6 known-date agreement is no longer 9 of 10")

    primary = comparison[comparison["series_role"] == "primary_estimand"]
    if int(primary["matched_within_primary_tolerance"].sum()) != 20:
        raise AssertionError("Expected 20 accepted matches within three months")
    if int(primary["matched_within_sensitivity_tolerance"].sum()) != 22:
        raise AssertionError("Expected 22 accepted matches within six months")
    if len(accepted) - int(primary["matched_within_sensitivity_tolerance"].sum()) != 15:
        raise AssertionError("Expected 15 accepted unmatched breaks")
    return {"detections": detections, "comparison": comparison, "accepted": accepted}


def standardized_monthly(monthly: xr.Dataset, series_id: str) -> np.ndarray:
    values = monthly["seasonal_adjusted"].sel(series=series_id).values.astype(float)
    mean = np.nanmean(values)
    std = np.nanstd(values, ddof=1)
    if not np.isfinite(std) or std <= 0:
        raise AssertionError(f"Cannot standardize production series {series_id}")
    return (values - mean) / std


def shade_periods(ax: plt.Axes, periods: pd.DataFrame) -> None:
    for row in periods.itertuples():
        ax.axvspan(row.start, row.end + pd.Timedelta(days=1), color=PERIOD_COLORS[row.period_id], zorder=0)
        center = row.start + (row.end - row.start) / 2
        ax.text(
            center,
            0.975,
            row.period_id,
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


def forest_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    series_ids: list[str],
    title: str,
    xlabel: str,
    multiplier: float,
    label: str,
) -> None:
    rows = frame.loc[series_ids]
    y = np.arange(len(rows))[::-1]
    estimates = rows["estimate"].to_numpy(dtype=float) * multiplier
    low = rows["ci_low_bootstrap"].to_numpy(dtype=float) * multiplier
    high = rows["ci_high_bootstrap"].to_numpy(dtype=float) * multiplier
    significant = rows["significant_fdr"].astype(bool).to_numpy()
    colors = np.where(significant, "#b73a35", "0.65")
    for index in range(len(rows)):
        ax.hlines(y[index], low[index], high[index], color=colors[index], linewidth=2.0, zorder=2)
        ax.scatter(
            estimates[index],
            y[index],
            s=42,
            facecolor=colors[index] if significant[index] else "white",
            edgecolor=colors[index],
            linewidth=1.2,
            zorder=3,
        )
    ax.axvline(0, color="0.25", linewidth=0.8, zorder=1)
    ax.set_yticks(y, [P6_SERIES[series_id]["label"] for series_id in series_ids])
    ax.set_ylim(-0.65, len(rows) - 0.35)
    xmax = max(float(np.nanmax(high)), 0.0)
    xmin = min(float(np.nanmin(low)), 0.0)
    span = xmax - xmin
    ax.set_xlim(xmin - 0.06 * span, xmax + 0.10 * span)
    ax.grid(axis="x", color="0.88", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_title(title, pad=7, fontweight="bold")
    ax.set_xlabel(xlabel)
    panel_label(ax, label, x=0.012, y=0.98)


def plot_observing_system_transitions(
    periods: pd.DataFrame,
    monthly: xr.Dataset,
    selected: pd.DataFrame,
) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(13.5, 7.6))
    grid = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.65, 1.0],
        left=0.07,
        right=0.985,
        bottom=0.085,
        top=0.975,
        hspace=0.30,
        wspace=0.38,
    )
    ax = fig.add_subplot(grid[0, :])
    shade_periods(ax, periods)
    time = pd.DatetimeIndex(monthly["time"].values)
    for series_id, (display, color) in TIMESERIES.items():
        ax.plot(time, standardized_monthly(monthly, series_id), color=color, linewidth=1.05, label=display, zorder=3)
    ax.axhline(0, color="0.35", linewidth=0.7, zorder=2)
    ax.set_xlim(periods.iloc[0]["start"], periods.iloc[-1]["end"] + pd.Timedelta(days=1))
    ax.set_ylabel("Standardized monthly value")
    ax.xaxis.set_major_locator(mdates.YearLocator(3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(axis="y", color="0.86", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", ncols=2, frameon=False, bbox_to_anchor=(0.01, 0.01))
    p6_start = periods.loc[periods["period_id"] == "P6", "start"].item()
    ax.text(
        p6_start + pd.Timedelta(days=50),
        0.06,
        "P6: SMAP Tb begins",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8.2,
        color="#7a1c1c",
        bbox={"boxstyle": "square,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        zorder=5,
    )
    panel_label(ax, "(a)")

    forest_panel(
        fig.add_subplot(grid[1, 0]),
        selected,
        ["RZMC_delta_valid_land__delta", "RZMC_delta_warm_snowfree_monthly__delta"],
        "RZMC state response",
        r"P6 level change ($\times 10^{-3}$ m$^3$ m$^{-3}$)",
        1000.0,
        "(b)",
    )
    forest_panel(
        fig.add_subplot(grid[1, 1]),
        selected,
        ["SFMC_INC_RMS_value_valid_land__value", "RZMC_INC_RMS_value_valid_land__value"],
        "Soil-moisture corrections",
        r"P6 level change ($\times 10^{-3}$ m$^3$ m$^{-3}$)",
        1000.0,
        "(c)",
    )
    forest_panel(
        fig.add_subplot(grid[1, 2]),
        selected,
        ["soil_water_abs_activity_value_valid_land__value", "soil_water_net_approx_value_valid_land__value"],
        "Soil-water activity",
        r"P6 level change (kg m$^{-2}$)",
        1.0,
        "(d)",
    )
    return save_figure(fig, "fig17_observing_system_transitions")


def breakpoint_display_label(variable: str, mask: str) -> str:
    variable_labels = {
        "PRECTOTCORRLAND": "Precipitation",
        "SNOMASLAND": "Snow mass",
        "SNODPLAND": "Snow depth",
        "FRLANDSNO": "SCF",
        "soil_water_net_approx": "Soil-water net",
        "soil_water_abs_activity": "Soil-water activity",
        "SFMC_INC_MEAN": "SFMC increment mean",
        "SFMC_INC_ABS_MEAN": "SFMC increment absolute mean",
        "SFMC_INC_RMS": "SFMC increment RMS",
        "RZMC_INC_MEAN": "RZMC increment mean",
        "RZMC_INC_ABS_MEAN": "RZMC increment absolute mean",
        "RZMC_INC_RMS": "RZMC increment RMS",
    }
    mask_labels = {
        "valid_land": "valid land",
        "seasonal_snow": "seasonal snow",
        "warm_snowfree_monthly": "warm snow-free",
        "locally_snowy_monthly": "locally snowy",
    }
    return f"{variable_labels.get(variable, variable)} | {mask_labels.get(mask, mask)}"


def plot_breakpoint_agreement(comparison: pd.DataFrame) -> tuple[Path, Path]:
    primary = comparison[comparison["series_role"] == "primary_estimand"].copy()
    series_meta = primary[["series_id", "variable", "mask"]].drop_duplicates()
    series_order = series_meta["series_id"].tolist()
    period_order = [f"P{i}" for i in range(2, 10)]
    matrix = np.full((len(series_order), len(period_order)), np.nan)
    for i, series_id in enumerate(series_order):
        for j, period_id in enumerate(period_order):
            row = primary[(primary["series_id"] == series_id) & (primary["period_id"] == period_id)].iloc[0]
            if bool(row["matched_within_sensitivity_tolerance"]):
                matrix[i, j] = row["signed_offset_months"]

    labels = [
        breakpoint_display_label(row.variable, row.mask)
        for row in series_meta.itertuples(index=False)
    ]
    colors = ["#31688e", "#86bdd0", "#ffffff", "#ee9b79", "#b73a35"]
    cmap = ListedColormap(colors)
    cmap.set_bad("#d9d9d9")
    norm = BoundaryNorm([-6.5, -3.5, -0.5, 0.5, 3.5, 6.5], cmap.N)
    fig, ax = plt.subplots(figsize=(8.4, 8.8), constrained_layout=True)
    mesh = ax.pcolormesh(
        np.ma.masked_invalid(matrix),
        cmap=cmap,
        norm=norm,
        edgecolors="white",
        linewidth=0.65,
    )
    xticklabels = ["P7\nexempt" if value == "P7" else value for value in period_order]
    ax.set_xticks(np.arange(len(period_order)) + 0.5, xticklabels)
    ax.set_yticks(np.arange(len(series_order)) + 0.5, labels)
    ax.invert_yaxis()
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                ax.text(j + 0.5, i + 0.5, f"{int(matrix[i, j]):+d}", ha="center", va="center", fontsize=7.4)
    p7_index = period_order.index("P7")
    ax.add_patch(
        Rectangle(
            (p7_index, 0),
            1,
            len(series_order),
            fill=False,
            hatch="////",
            edgecolor="#555555",
            linewidth=1.0,
        )
    )
    colorbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.035, fraction=0.042, ticks=[-5, -2, 0, 2, 5])
    colorbar.set_label("Accepted break minus known boundary (months)")
    ax.set_xlabel("Observing-system boundary")
    return save_figure(fig, "figSXX_breakpoint_boundary_agreement")


def build_assets() -> list[FigureAsset]:
    periods = load_periods()
    validate_trend_products()
    transition = validate_transition_products(periods)
    breakpoints = validate_breakpoint_products(transition["coefficients"])

    fig_x = plot_trend_rows(["RZMC", "FRLANDSNO"], "fig16_longterm_rzmc_scf_trends")
    fig_y = plot_observing_system_transitions(periods, transition["monthly"], transition["selected"])
    precip = plot_trend_rows(["PRECTOTCORRLAND"], "figSXX_precipitation_trends")
    sfmc = plot_trend_rows(["SFMC"], "figSXX_sfmc_trends")
    snow = plot_trend_rows(
        ["SNOMASLAND", "SNODPLAND"],
        "figSXX_snow_mass_depth_trends",
        delta_scale="state",
    )
    breaks = plot_breakpoint_agreement(breakpoints["comparison"])

    assets = [
        FigureAsset(
            "Main Figure 16: long-term state trends",
            *fig_x,
            caption=(
                "Long-term June 2000-May 2024 trends in (a-c) root-zone soil moisture (RZMC) and "
                "(d-f) snow-cover fraction (SCF) for the open-loop (OL), data-assimilation (DA), "
                "and paired DA-OL series. Trends are exact Theil-Sen slopes after trend-preserving "
                "removal of the calendar-month climatology. Black stippling denotes trends significant "
                "after Benjamini-Hochberg false-discovery-rate control at 0.05. RZMC uses the valid-land "
                "domain and SCF the Northern Hemisphere seasonal-snow domain. The DA-OL panels show the trend of the paired "
                "DA-OL series rather than the difference between independently estimated OL and DA trends."
            ),
            sources=(
                "RZMC_{ol,da,delta}_valid_land_trend_statistics.nc",
                "FRLANDSNO_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
            ),
            provenance="Exact production slope; mapped significance is significant_fdr; Robinson projection; 60 S cutoff.",
        ),
        FigureAsset(
            "Main Figure 17: observing-system transitions",
            *fig_y,
            caption=(
                "Changes in soil-water data-assimilation behavior across the P1-P9 observing-system "
                "periods. (a) Standardized area-weighted monthly RZMC DA-OL and soil-water analysis-"
                "correction diagnostics during June 2000-May 2024. Background shading denotes the P1-P9 "
                "periods defined in Fig. 1; the P5-P6 boundary in April 2015 marks the introduction of "
                "SMAP brightness-temperature assimilation. Each seasonally adjusted series is standardized "
                "by its full-record mean and sample standard deviation for visual comparison only. "
                "(b-d) Estimated P5-P6 level changes from the interrupted time-series analysis for RZMC "
                "DA-OL, soil-moisture analysis-correction RMS, and prognostic soil-water correction activity. "
                "Symbols show the estimate and horizontal bars the 95% fitted-AR(1) bootstrap interval. "
                "Statistical significance uses boundary-family false-discovery-rate control at 0.05."
            ),
            sources=(
                "phase1_changepoint_monthly.nc",
                "phase1_interrupted_series_coefficients.csv",
                "observing_system_registry.json",
            ),
            provenance="Monthly production seasonal_adjusted fields; display-only full-record z scores; native-unit P6 inference.",
        ),
        FigureAsset(
            "Supporting Figure: precipitation trends",
            *precip,
            caption=(
                "Long-term precipitation trends for OL, DA, and the paired DA-OL series on valid land. "
                "Black stippling denotes production FDR significance. The common OL/DA pattern and null "
                "DA-OL result provide a forcing-control check."
            ),
            sources=("PRECTOTCORRLAND_{ol,da,delta}_valid_land_trend_statistics.nc",),
            provenance="Exact production slope and significant_fdr on valid land.",
        ),
        FigureAsset(
            "Supporting Figure: SFMC trends",
            *sfmc,
            caption=(
                "Long-term surface soil-moisture trends for OL, DA, and the paired DA-OL series on valid land. "
                "Black stippling denotes production FDR significance."
            ),
            sources=("SFMC_{ol,da,delta}_valid_land_trend_statistics.nc",),
            provenance="Exact production slope and significant_fdr on valid land.",
        ),
        FigureAsset(
            "Supporting Figure: snow mass and depth trends",
            *snow,
            caption=(
                "Long-term (a-c) snow-mass and (d-f) snow-depth trends for OL, DA, and the paired DA-OL "
                "series on the production Northern Hemisphere seasonal-snow mask. Black stippling denotes production FDR significance."
            ),
            sources=(
                "SNOMASLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
                "SNODPLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc",
            ),
            provenance="Exact production slope and significant_fdr on the Northern Hemisphere seasonal-snow mask.",
        ),
        FigureAsset(
            "Supporting Figure: breakpoint-boundary agreement",
            *breaks,
            caption=(
                "Accepted independent changepoints relative to known P2-P9 dates for the primary Phase 1 "
                "estimands. Values are detected-minus-known months; blue is early, red late, white exact, "
                "and grey indicates no accepted match within six months. P7 is hatched because its 15-month "
                "duration is detection-exempt under the predeclared minimum-segment rule."
            ),
            sources=(
                "phase1_changepoint_boundary_comparison.csv",
                "phase1_changepoint_detections.csv",
            ),
            provenance="Accepted consensus breaks only; +/-3-month primary and +/-6-month sensitivity definitions retained.",
        ),
    ]
    write_report(assets, periods)
    return assets


def write_report(assets: list[FigureAsset], periods: pd.DataFrame) -> None:
    lines = [
        "# Trends And Observing-System Manuscript Figure Report",
        "",
        "These figures are plotting-only products generated from the accepted Phase 1 trend, interrupted-series, and changepoint outputs. No scientific analysis was refitted.",
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
            "- P6 begins 2015-04-01. All six plotted native-unit estimates, fitted-AR(1) bootstrap intervals, and boundary-family FDR flags reproduce the production coefficient table.",
            "- Ten primary series have accepted breaks exactly in April 2015; nine also have significant known-date P6 level changes; no accepted break occurs in paired OL or DA state controls.",
            "- The accepted-break inventory remains 37 total: 20 paired DA-OL and 17 correction diagnostics; 20 match within +/-3 months, two additional within +/-6 months, and 15 remain unmatched.",
            "- All maps use exact production Theil-Sen slopes and `significant_fdr`; pointwise confidence-interval exclusion is not used for mapped inference.",
            "",
            "## Plotting Choices",
            "",
            "- Maps follow the existing report convention: Robinson projection, 60 S cutoff, grey land, thin coastlines, segmented RdBu_r scales centered on a white zero bin, and black stippling.",
            "- OL and DA share a symmetric color scale within each row. DA-OL uses a separately labeled symmetric scale where needed; snow mass/depth retain the OL/DA scale so negligible differences are not visually exaggerated.",
            "- Figure 17 panel (a) shows unsmoothed monthly production series. Full-record z scoring is display-only; interrupted-series inference remains in native units.",
            "- The requested optional all-boundary decorative summary was not produced; the accepted breakpoint-agreement matrix already provides the auditable all-boundary view.",
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
    periods = load_periods()
    validate_trend_products()
    transition = validate_transition_products(periods)
    validate_breakpoint_products(transition["coefficients"])
    if args.validate_only:
        print("Trend and transition manuscript figure contracts: PASS")
        return
    assets = build_assets()
    print(f"Generated {len(assets)} manuscript figures")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Figure report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
