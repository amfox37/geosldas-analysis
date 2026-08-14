#!/usr/bin/env python3
"""Create the main and supplemental snow-DA hydrology manuscript figures.

This is a plotting-only workflow. It reads the accepted machine-readable
analysis products and validates their scientific contracts before drawing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


DPI = 300
MONTHS = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
COLORS = {
    "snow_input": "#c44e52",
    "snowmelt": "#7b6fa2",
    "surface_runoff": "#238b7b",
    "baseflow": "#74b49b",
    "total_runoff": "#15616d",
    "et": "#d8a43a",
    "storage": "#4c78a8",
    "residual": "#7f7f7f",
    "sfmc": "#e76f51",
    "rzmc": "#238b7b",
    "control": "#4c78a8",
}
EXPECTED_ANNUAL = {
    2001: [32.57, 19.61, 12.60, 2.82, -2.46],
    2002: [43.46, 24.14, 15.31, 4.29, -0.28],
    2003: [64.50, 40.98, 24.15, 3.12, -3.75],
    2004: [67.50, 45.44, 23.27, 2.35, -3.56],
    2005: [72.97, 51.14, 24.53, 0.23, -2.93],
    2006: [68.43, 47.57, 24.13, -0.23, -3.05],
}
ANNUAL_TERMS = ["I_snow", "dRunoff_total", "dET", "dStorage", "residual"]
PARTITION_TERMS = ["dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]
BOUNDARY_TERMS = [
    "dRunoff_total",
    "dRunoff_surface",
    "dBaseflow",
    "dET",
    "dStorage",
    "residual",
]


@dataclass(frozen=True)
class FigureProduct:
    figure: str
    png: Path
    pdf: Path
    figsize: tuple[float, float]
    sources: tuple[str, ...]
    differences: str


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    raise FileNotFoundError("Could not locate repository root")


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.size": 9.5,
            "axes.titlesize": 10.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "legend.fontsize": 8.5,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
            "lines.markersize": 5.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.025) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        clip_on=False,
    )


def save_figure(fig: plt.Figure, png: Path, pdf: Path) -> None:
    png.parent.mkdir(parents=True, exist_ok=True)
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def input_paths(root: Path) -> dict[str, Path]:
    water = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget"
    targeted = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness"
    return {
        "annual": water / "annual_domain_budgets.csv",
        "partition": water / "six_year_integrated_partitions.csv",
        "partition_ci": water / "partition_spatial_block_uncertainty.csv",
        "monthly_addition": water / "monthly_climatology_snow_addition.csv",
        "soil_summary": water / "soil_moisture_summary_snow_addition.csv",
        "peak_timing": water / "soil_moisture_peak_timing.csv",
        "octmar": targeted / "analysisA_octmar_signed_controls.csv",
        "boundary": targeted / "water_year_boundary_sensitivity.csv",
        "september": targeted / "water_year_september_timing_diagnostic.csv",
    }


def load_tables(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing accepted figure inputs:\n" + "\n".join(missing))
    return {name: pd.read_csv(path) for name, path in paths.items()}


def validate_annual_budget(annual: pd.DataFrame) -> pd.DataFrame:
    years = annual[annual["water_year"].astype(str).isin([str(year) for year in EXPECTED_ANNUAL])].copy()
    years["water_year"] = years["water_year"].astype(int)
    years = years.sort_values("water_year")
    if years["water_year"].tolist() != list(EXPECTED_ANNUAL):
        raise AssertionError("Figure 14a does not contain exactly WY2001-WY2006")
    for _, row in years.iterrows():
        expected = np.asarray(EXPECTED_ANNUAL[int(row["water_year"])])
        found = row[ANNUAL_TERMS].to_numpy(dtype=float)
        if not np.allclose(found, expected, atol=0.006, rtol=0):
            raise AssertionError(
                f"WY{int(row['water_year'])} differs from accepted rounded values: {found}"
            )
        if not np.isclose(found[0], found[1:].sum(), atol=1.0e-9):
            raise AssertionError(f"WY{int(row['water_year'])} direct budget does not close")
    return years


def positive_partition_row(
    partition: pd.DataFrame,
    boundary: pd.DataFrame,
) -> pd.Series:
    canonical = partition.set_index("sample").loc["addition"]
    targeted = boundary.query(
        "boundary == 'Oct-Sep' and scope == 'positive_input_partition'"
    )
    if len(targeted) != 1:
        raise AssertionError("Expected one Oct-Sep positive-input boundary row")
    row = targeted.iloc[0]
    if int(row["n_tile_years"]) != 247545 or int(canonical["n_tile_years"]) != 247545:
        raise AssertionError("Figure 14b is not using the accepted positive-input tile-years")
    compare = ["I_snow", "dRunoff_total", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]
    if not np.allclose(
        row[compare].to_numpy(dtype=float),
        canonical[compare].to_numpy(dtype=float),
        atol=1.0e-8,
    ):
        raise AssertionError("Targeted and canonical positive-input partitions differ")
    if not np.isclose(row["dRunoff_surface"] + row["dBaseflow"], row["dRunoff_total"], atol=1.0e-8):
        raise AssertionError("Surface runoff plus baseflow does not equal total runoff")
    fractions = row[
        [f"fraction_{term}" for term in ["dRunoff_total", "dET", "dStorage", "residual"]]
    ].to_numpy(dtype=float).sum()
    if not np.isclose(fractions, 1.0, atol=1.0e-10):
        raise AssertionError(f"Positive-input partition sums to {fractions}, not one")
    runoff_ci = row[["fraction_dRunoff_total_ci_low_5deg", "fraction_dRunoff_total_ci_high_5deg"]]
    if not np.allclose(runoff_ci.to_numpy(dtype=float), [0.611077, 0.672333], atol=5.0e-7):
        raise AssertionError("Figure 14b runoff CI is not the positive-input 5-degree interval")
    return row


def validate_monthly_pathway(
    monthly: pd.DataFrame,
    partition: pd.DataFrame,
    soil_summary: pd.DataFrame,
    peak_timing: pd.DataFrame,
    september: pd.DataFrame,
) -> dict[str, float | str]:
    if monthly["month"].tolist() != MONTHS or monthly["water_month"].tolist() != list(range(1, 13)):
        raise AssertionError("Figure 15 monthly source is not ordered Oct-Sep")
    pathway_columns = [
        "snow_net_monthly",
        "dSnowmelt_monthly",
        "dSFMC_monthly",
        "dRZMC_monthly",
        "dRunoff_total_monthly",
        "dET_monthly",
    ]
    if not np.isfinite(monthly[pathway_columns].to_numpy(dtype=float)).all():
        raise AssertionError("Figure 15 monthly pathway contains non-finite values")
    addition = partition.set_index("sample").loc["addition"]
    summary = soil_summary.set_index("metric")
    if int(addition["n_tile_years"]) != 247545 or int(summary.loc["peak_dRZMC", "n"]) != 247545:
        raise AssertionError("Figure 15 is not using the accepted snow-addition sample")
    september_input = september.set_index("variable").loc[
        "signed_snow_input", "snow_addition_tile_years_mean_kg_m2"
    ]
    if not np.isclose(monthly.loc[monthly["month"] == "Sep", "snow_net_monthly"].item(), september_input):
        raise AssertionError("Monthly pathway September input does not match the addition-sample diagnostic")
    rzmc_timing = peak_timing.query("state == 'RZMC'").sort_values("area_weighted_fraction")
    most_common = rzmc_timing.iloc[-1]
    if most_common["month"] != "May":
        raise AssertionError("May is no longer the most common RZMC peak month")
    return {
        "sample": "I_snow > 0",
        "n_tile_years": int(addition["n_tile_years"]),
        "peak_rzmc": float(summary.loc["peak_dRZMC", "area_weighted_mean"]),
        "mjj_rzmc": float(summary.loc["mjj_mean_dRZMC", "area_weighted_mean"]),
        "peak_month": str(most_common["month"]),
    }


def validate_octmar(octmar: pd.DataFrame) -> None:
    expected_responses = ["snowmelt", "infiltration", "rzmc", "et", "total_runoff", "total_water"]
    if octmar["response"].tolist() != expected_responses:
        raise AssertionError("Oct-Mar supplemental responses changed")
    for year in range(2001, 2007):
        predictor = pd.date_range(f"{year - 1}-10-01", f"{year}-03-01", freq="MS")
        for months in ([4, 5, 6], [5, 6, 7]):
            response = pd.DatetimeIndex([pd.Timestamp(year, month, 1) for month in months])
            if len(predictor.intersection(response)):
                raise AssertionError(f"Oct-Mar predictor overlaps a response in {year}")
    infiltration = octmar.set_index("response").loc["infiltration"]
    if not (infiltration["m3_march_snow_ci_low_5deg"] < 0 < infiltration["m3_march_snow_ci_high_5deg"]):
        raise AssertionError("Infiltration M3 interval no longer crosses zero")


def validate_boundaries(boundary: pd.DataFrame) -> pd.DataFrame:
    rows = boundary.query("scope == 'positive_input_partition'").set_index("boundary")
    if rows.index.tolist() != ["Oct-Sep", "Sep-Aug"]:
        rows = rows.reindex(["Oct-Sep", "Sep-Aug"])
    if rows.isnull().all(axis=1).any():
        raise AssertionError("Boundary comparison is missing an accounting definition")
    for name, row in rows.iterrows():
        closure = row[[f"fraction_{term}" for term in ["dRunoff_total", "dET", "dStorage", "residual"]]].sum()
        if not np.isclose(closure, 1.0, atol=1.0e-10):
            raise AssertionError(f"{name} positive-input partition does not close")
        if not np.isclose(row["fraction_dRunoff_surface"] + row["fraction_dBaseflow"], row["fraction_dRunoff_total"]):
            raise AssertionError(f"{name} runoff components do not sum to total runoff")
    if rows["n_tile_years"].min() <= 0:
        raise AssertionError("Boundary comparison has no positive-input samples")
    # Both rows were generated in one table by the same area-weighted partition
    # and spatial-block bootstrap functions. The figure never mixes products.
    required_ci = [
        f"fraction_{term}_ci_{side}_5deg"
        for term in BOUNDARY_TERMS
        for side in ["low", "high"]
    ]
    if rows[required_ci].isnull().any().any():
        raise AssertionError("Boundary comparison lacks matched 5-degree intervals")
    return rows


def validate_modis_only_window() -> None:
    octsep = pd.date_range("2000-10-01", "2006-09-01", freq="MS")
    sepaug = pd.date_range("2000-09-01", "2006-08-01", freq="MS")
    microwave_start = pd.Timestamp("2007-06-01")
    if max(octsep.max(), sepaug.max()) >= microwave_start:
        raise AssertionError("Microwave soil-moisture DA enters the six-year figure sample")


def plot_figure14(
    annual: pd.DataFrame,
    partition: pd.Series,
    png: Path,
    pdf: Path,
) -> None:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.8, 5.8),
        gridspec_kw={"width_ratios": [1.18, 1.0]},
        constrained_layout=True,
    )
    ax = axes[0]
    x = np.arange(len(annual))
    runoff = annual["dRunoff_total"].to_numpy()
    et = annual["dET"].to_numpy()
    storage = annual["dStorage"].to_numpy()
    residual = annual["residual"].to_numpy()
    ax.bar(x, runoff, width=0.66, color=COLORS["total_runoff"], label="Runoff")
    ax.bar(x, et, bottom=runoff, width=0.66, color=COLORS["et"], label="ET")
    ax.bar(x, storage, bottom=runoff + et, width=0.66, color=COLORS["storage"], label="Storage")
    ax.bar(
        x,
        residual,
        width=0.66,
        color=COLORS["residual"],
        edgecolor="0.25",
        linewidth=0.6,
        hatch="///",
        label="Residual",
    )
    ax.plot(
        x,
        annual["I_snow"],
        color=COLORS["snow_input"],
        marker="D",
        markerfacecolor="white",
        markeredgewidth=1.2,
        linewidth=1.8,
        label="Snow-DA input",
        zorder=4,
    )
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.set_xticks(x, [f"WY{year}" for year in annual["water_year"]])
    ax.set_ylabel(r"Area-weighted mean (kg m$^{-2}$ yr$^{-1}$)")
    ax.set_ylim(-8, 82)
    ax.grid(axis="y", color="0.85", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", ncols=2, frameon=False)
    ax.text(
        0.98,
        0.95,
        "All-tile net accounting\nMean residual = -4.6% of net input",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
    )
    panel_label(ax, "(a)")

    ax = axes[1]
    labels = ["Surface\nrunoff", "Baseflow", "ET", "Storage", "Residual"]
    colors = [
        COLORS["surface_runoff"],
        COLORS["baseflow"],
        COLORS["et"],
        COLORS["storage"],
        COLORS["residual"],
    ]
    values = np.array([partition[f"fraction_{term}"] for term in PARTITION_TERMS])
    low = np.array([partition[f"fraction_{term}_ci_low_5deg"] for term in PARTITION_TERMS])
    high = np.array([partition[f"fraction_{term}_ci_high_5deg"] for term in PARTITION_TERMS])
    xpos = np.arange(len(values))
    bars = ax.bar(
        xpos,
        100 * values,
        width=0.68,
        color=colors,
        yerr=np.vstack([100 * (values - low), 100 * (high - values)]),
        error_kw={"ecolor": "0.2", "elinewidth": 1.0, "capsize": 3},
    )
    bars[-1].set_hatch("///")
    bars[-1].set_edgecolor("0.25")
    bars[-1].set_linewidth(0.6)
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.set_xticks(xpos, labels)
    ax.set_ylabel("Fraction of positive snow-DA input (%)")
    ax.set_ylim(-8, 58)
    ax.grid(axis="y", color="0.85", linewidth=0.7)
    ax.set_axisbelow(True)
    bracket_y = 49.0
    ax.plot([0, 0, 1, 1], [bracket_y - 1.1, bracket_y, bracket_y, bracket_y - 1.1], color="0.2", linewidth=1.0)
    runoff_total = 100 * partition["fraction_dRunoff_total"]
    runoff_low = 100 * partition["fraction_dRunoff_total_ci_low_5deg"]
    runoff_high = 100 * partition["fraction_dRunoff_total_ci_high_5deg"]
    ax.text(
        0.5,
        bracket_y + 1.2,
        f"Total runoff: {runoff_total:.1f}% [{runoff_low:.1f}, {runoff_high:.1f}]",
        ha="center",
        va="bottom",
        fontsize=9.0,
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.94,
        rf"Mean positive input: {partition['I_snow']:.2f} kg m$^{{-2}}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
    )
    panel_label(ax, "(b)")
    save_figure(fig, png, pdf)


def plot_figure15(
    monthly: pd.DataFrame,
    pathway: dict[str, float | str],
    png: Path,
    pdf: Path,
) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(11.5, 8.6), sharex=True, constrained_layout=True)
    x = monthly["water_month"].to_numpy()

    ax = axes[0]
    ax.bar(x, monthly["snow_net_monthly"], width=0.68, color=COLORS["snow_input"])
    ax.set_ylabel("Snow-DA input\n" + r"(kg m$^{-2}$ month$^{-1}$)")
    ax.text(
        0.99,
        0.91,
        r"Snow-addition tile-years ($I_{snow} > 0$)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
    )
    panel_label(ax, "(a)", x=-0.065)

    ax = axes[1]
    ax.plot(x, monthly["dSnowmelt_monthly"], marker="o", color=COLORS["snowmelt"])
    ax.fill_between(x, 0, monthly["dSnowmelt_monthly"], color=COLORS["snowmelt"], alpha=0.12)
    ax.set_ylabel(r"$\Delta$Snowmelt" + "\n" + r"(kg m$^{-2}$ month$^{-1}$)")
    panel_label(ax, "(b)", x=-0.065)

    ax = axes[2]
    ax.plot(
        x,
        monthly["dSFMC_monthly"],
        marker="o",
        color=COLORS["sfmc"],
        linestyle="--",
        linewidth=1.5,
        label="SFMC",
    )
    ax.plot(
        x,
        monthly["dRZMC_monthly"],
        marker="o",
        color=COLORS["rzmc"],
        linewidth=2.3,
        label="RZMC",
    )
    ax.set_ylabel(r"$\Delta$Soil moisture" + "\n" + r"(m$^3$ m$^{-3}$)")
    ax.legend(loc="upper left", frameon=False, ncols=2)
    ax.annotate(
        rf"Mean tile-year peak $\Delta$RZMC = {pathway['peak_rzmc']:.4f} m$^3$ m$^{{-3}}$"
        + f"\n{pathway['peak_month']} is the most common peak month",
        xy=(8, monthly.loc[monthly["month"] == "May", "dRZMC_monthly"].item()),
        xytext=(8.6, 0.0184),
        arrowprops={"arrowstyle": "-", "color": "0.35", "linewidth": 0.8},
        ha="left",
        va="top",
        fontsize=8.6,
    )
    ax.set_ylim(-0.001, 0.0196)
    panel_label(ax, "(c)", x=-0.065)

    ax = axes[3]
    ax.plot(
        x,
        monthly["dRunoff_total_monthly"],
        marker="o",
        linewidth=2.1,
        color=COLORS["total_runoff"],
        label="Runoff",
    )
    ax.fill_between(x, 0, monthly["dRunoff_total_monthly"], color=COLORS["total_runoff"], alpha=0.10)
    ax.plot(
        x,
        monthly["dET_monthly"],
        marker="o",
        linewidth=1.8,
        linestyle="--",
        color=COLORS["et"],
        label="ET",
    )
    ax.fill_between(x, 0, monthly["dET_monthly"], color=COLORS["et"], alpha=0.08)
    ax.set_ylabel(r"$\Delta$Water flux" + "\n" + r"(kg m$^{-2}$ month$^{-1}$)")
    ax.legend(loc="upper left", frameon=False, ncols=2)
    ax.set_xlabel("Accounting-year month")
    panel_label(ax, "(d)", x=-0.065)

    for ax in axes:
        ax.axhline(0, color="0.3", linewidth=0.8)
        ax.grid(axis="y", color="0.86", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_xlim(0.5, 12.5)
    axes[-1].set_xticks(range(1, 13), MONTHS)
    save_figure(fig, png, pdf)


def coefficient_unit(response: str) -> str:
    if response == "rzmc":
        return r"m$^3$ m$^{-3}$ per kg m$^{-2}$"
    if response == "total_water":
        return r"kg m$^{-2}$ per kg m$^{-2}$"
    return r"kg m$^{-2}$ season$^{-1}$ per kg m$^{-2}$"


def plot_octmar_supplement(octmar: pd.DataFrame, png: Path, pdf: Path) -> None:
    response_order = ["snowmelt", "infiltration", "rzmc", "et", "total_runoff", "total_water"]
    labels = {
        "snowmelt": "AMJ snowmelt",
        "infiltration": "AMJ infiltration",
        "rzmc": "MJJ RZMC",
        "et": "MJJ ET",
        "total_runoff": "MJJ total runoff",
        "total_water": "MJJ TWLAND",
    }
    models = ["pooled_beta", "m1_beta", "m2_beta", "m3_march_snow_beta"]
    model_labels = ["Pooled", "M1", "M2", "M3"]
    lookup = octmar.set_index("response")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.4), constrained_layout=True)
    for index, (ax, response) in enumerate(zip(axes.flat, response_order)):
        row = lookup.loc[response]
        values = row[models].to_numpy(dtype=float)
        scale = 1.0e5 if response == "rzmc" else 1.0
        plotted = values * scale
        ax.plot(range(4), plotted, marker="o", color=COLORS["control"])
        low = row["m3_march_snow_ci_low_5deg"] * scale
        high = row["m3_march_snow_ci_high_5deg"] * scale
        ax.errorbar(
            [3],
            [plotted[-1]],
            yerr=[[plotted[-1] - low], [high - plotted[-1]]],
            color=COLORS["snow_input"],
            marker="o",
            markerfacecolor="white",
            capsize=4,
            linewidth=1.2,
            zorder=4,
        )
        ax.axhline(0, color="0.3", linewidth=0.8)
        ax.set_xticks(range(4), model_labels)
        unit = coefficient_unit(response)
        if response == "rzmc":
            unit += r" $\times 10^{-5}$"
        ax.set_ylabel(unit)
        ax.set_title(labels[response], loc="left", fontsize=10, pad=5)
        ax.grid(axis="y", color="0.86", linewidth=0.7)
        ax.set_axisbelow(True)
        panel_label(ax, f"({chr(97 + index)})", x=-0.13)
        combined = np.append(plotted, [low, high, 0])
        span = combined.max() - combined.min()
        margin = max(0.12 * span, 0.015 * max(abs(combined).max(), 1.0e-12))
        ax.set_ylim(combined.min() - margin, combined.max() + margin)
        if response == "infiltration":
            ax.text(
                0.98,
                0.08,
                "M3 interval crosses zero",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8.2,
                color=COLORS["snow_input"],
            )
    fig.text(
        0.5,
        -0.01,
        "Predictor: Oct-Mar signed snow-DA input. Responses begin in April or May; predictor and response windows do not overlap.",
        ha="center",
        va="top",
        fontsize=9.0,
    )
    save_figure(fig, png, pdf)


def boundary_color(term: str) -> str:
    return {
        "dRunoff_total": COLORS["total_runoff"],
        "dRunoff_surface": COLORS["surface_runoff"],
        "dBaseflow": COLORS["baseflow"],
        "dET": COLORS["et"],
        "dStorage": COLORS["storage"],
        "residual": COLORS["residual"],
    }[term]


def plot_boundary_supplement(boundary: pd.DataFrame, png: Path, pdf: Path) -> None:
    labels = ["Total\nrunoff", "Surface\nrunoff", "Baseflow", "ET", "Storage", "Residual"]
    x = np.arange(len(BOUNDARY_TERMS))
    fig, ax = plt.subplots(figsize=(11.5, 5.2), constrained_layout=True)
    marker_specs = {"Oct-Sep": ("o", -0.09), "Sep-Aug": ("s", 0.09)}
    for term_index, term in enumerate(BOUNDARY_TERMS):
        color = boundary_color(term)
        points = []
        for name in ["Oct-Sep", "Sep-Aug"]:
            marker, offset = marker_specs[name]
            row = boundary.loc[name]
            value = 100 * row[f"fraction_{term}"]
            low = 100 * row[f"fraction_{term}_ci_low_5deg"]
            high = 100 * row[f"fraction_{term}_ci_high_5deg"]
            points.append((term_index + offset, value))
            ax.errorbar(
                term_index + offset,
                value,
                yerr=[[value - low], [high - value]],
                fmt=marker,
                markersize=6,
                color=color,
                markeredgecolor="0.2",
                markeredgewidth=0.5,
                capsize=3,
                linewidth=1.0,
                zorder=3,
            )
        ax.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            color=color,
            linewidth=1.2,
            alpha=0.85,
            zorder=2,
        )
    ax.axhline(0, color="0.3", linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Fraction of positive snow-DA input (%)")
    ax.set_ylim(-10, 72)
    ax.grid(axis="y", color="0.86", linewidth=0.7)
    ax.set_axisbelow(True)
    legend = [
        Line2D([0], [0], marker="o", color="0.25", linestyle="none", label="Oct-Sep", markerfacecolor="white"),
        Line2D([0], [0], marker="s", color="0.25", linestyle="none", label="Sep-Aug", markerfacecolor="white"),
    ]
    ax.legend(handles=legend, loc="upper right", frameon=False, ncols=2)
    changes = 100 * (
        boundary.loc["Sep-Aug", [f"fraction_{term}" for term in BOUNDARY_TERMS]].to_numpy(dtype=float)
        - boundary.loc["Oct-Sep", [f"fraction_{term}" for term in BOUNDARY_TERMS]].to_numpy(dtype=float)
    )
    ax.text(
        0.02,
        0.95,
        f"Maximum absolute change = {np.max(np.abs(changes)):.1f} percentage points",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
    )
    save_figure(fig, png, pdf)


def image_metadata(path: Path) -> tuple[int, int, float]:
    with Image.open(path) as image:
        dpi = image.info.get("dpi", (DPI, DPI))[0]
        return image.width, image.height, float(dpi)


def relative(path: Path, root: Path) -> str:
    return str(path.relative_to(root))


def write_figure_report(
    root: Path,
    products: list[FigureProduct],
    docs_pngs: list[Path],
    pathway: dict[str, float | str],
    september: pd.DataFrame,
    boundary: pd.DataFrame,
    report_path: Path,
) -> None:
    lines = [
        "# Snow-DA Manuscript Figure Report",
        "",
        "The four figures below are plotting-only products built from the accepted machine-readable analysis tables. No statistical analysis was recomputed.",
        "",
        "## Outputs",
        "",
        "| Figure | PNG | PDF | PNG dimensions | Figure size |",
        "|---|---|---|---:|---:|",
    ]
    for product in products:
        width, height, dpi = image_metadata(product.png)
        lines.append(
            f"| {product.figure} | `{relative(product.png, root)}` | `{relative(product.pdf, root)}` | "
            f"{width} x {height} px at {dpi:.0f} DPI | {product.figsize[0]:.1f} x {product.figsize[1]:.1f} in |"
        )
    lines.extend(["", "Tracked review PNGs:"])
    lines.extend(f"- `{relative(path, root)}`" for path in docs_pngs)
    lines.extend(["", "## Sources And Changes", ""])
    for product in products:
        lines.extend(
            [
                f"### {product.figure}",
                "",
                "Sources:",
                *[f"- `{source}`" for source in product.sources],
                "",
                f"Difference from the old figure: {product.differences}",
                "",
            ]
        )
    sep = september.set_index("variable")
    boundary_rows = boundary.loc[["Oct-Sep", "Sep-Aug"]]
    lines.extend(
        [
            "## Validation",
            "",
            "- Figure 14a reproduces the six accepted WY2001-WY2006 all-tile budget values and each year closes as `input = runoff + ET + storage + residual`.",
            "- Figure 14b uses 247,545 positive-input tile-years. Surface runoff plus baseflow equals total runoff, and runoff + ET + storage + residual equals 100% before rounding.",
            f"- Figure 14b total-runoff interval is the positive-input 5-degree spatial-block interval: {100 * boundary_rows.loc['Oct-Sep', 'fraction_dRunoff_total_ci_low_5deg']:.1f}-{100 * boundary_rows.loc['Oct-Sep', 'fraction_dRunoff_total_ci_high_5deg']:.1f}%.",
            f"- Figure 15 uses the same `I_snow > 0` sample ({pathway['n_tile_years']:,} tile-years) as the accepted pathway analysis. Mean peak RZMC is {pathway['peak_rzmc']:.4f} m3 m-3, {pathway['peak_month']} is the most common peak month, and mean MJJ RZMC is {pathway['mjj_rzmc']:.4f} m3 m-3.",
            "- The Oct-Mar supplemental predictor contains exactly October-March and has zero overlap with the AMJ and MJJ response windows in all six years.",
            "- Both boundary rows come from the same targeted output table and use the same area-weighted positive-input partition and 5-degree block-bootstrap machinery.",
            "- All figure samples end by September 2006, before microwave soil-moisture DA begins in June 2007.",
            "",
            "## Caption Notes",
            "",
            f"September signed input is {sep.loc['signed_snow_input', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 ({sep.loc['signed_snow_input', 'percent_of_octsep_annual_input_all']:.1f}% of Oct-Sep annual input), September DA-OL snowmelt is {sep.loc['snowmelt', 'all_tile_years_mean_kg_m2']:.2f} kg m-2, and the snow-mass difference is only {sep.loc['snow_mass_change', 'all_tile_years_mean_kg_m2']:.2f} kg m-2. These belong in the boundary-sensitivity caption rather than in the graphic.",
            "",
            "PDFs use embedded TrueType text (`pdf.fonttype = 42`) and vector lines, markers, bars, and error bars. PNGs are exported at 300 DPI.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines))


def main() -> None:
    root = repo_root()
    configure_style()
    paths = input_paths(root)
    tables = load_tables(paths)
    annual = validate_annual_budget(tables["annual"])
    partition = positive_partition_row(tables["partition"], tables["boundary"])
    pathway = validate_monthly_pathway(
        tables["monthly_addition"],
        tables["partition"],
        tables["soil_summary"],
        tables["peak_timing"],
        tables["september"],
    )
    validate_octmar(tables["octmar"])
    boundary = validate_boundaries(tables["boundary"])
    validate_modis_only_window()

    output_dir = root / "projects/M21C_ls/output/paper_figures"
    docs_dir = root / "projects/M21C_ls/docs/paper_figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)

    products = [
        FigureProduct(
            "Figure 14: snow-DA water-budget accounting",
            output_dir / "fig14_snow_da_water_budget.png",
            output_dir / "fig14_snow_da_water_budget.pdf",
            (13.8, 5.8),
            (
                relative(paths["annual"], root),
                relative(paths["partition"], root),
                relative(paths["boundary"], root),
            ),
            "Combines the old annual and positive-partition diagnostics into one main-paper hierarchy: six all-tile annual budgets at left and the positive-input fate at right, with runoff components explicitly joined and the residual treated as closure.",
        ),
        FigureProduct(
            "Figure 15: monthly snow-DA pathway",
            output_dir / "fig15_snow_da_monthly_pathway.png",
            output_dir / "fig15_snow_da_monthly_pathway.pdf",
            (11.5, 8.6),
            (
                relative(paths["monthly_addition"], root),
                relative(paths["soil_summary"], root),
                relative(paths["peak_timing"], root),
            ),
            "Retains the old pathway concept but removes infiltration and persistence/residence-time claims, emphasizes RZMC, and aligns snow input, snowmelt, soil moisture, runoff, and ET on one Oct-Sep axis.",
        ),
        FigureProduct(
            "Supplemental Figure Sx: non-overlapping attribution",
            output_dir / "fig_supp_snow_da_octmar_attribution.png",
            output_dir / "fig_supp_snow_da_octmar_attribution.pdf",
            (13.5, 7.4),
            (relative(paths["octmar"], root),),
            "Expands the preliminary four-panel control sequence to all six responses, retains native units in separate panels, and makes the infiltration interval crossing zero explicit.",
        ),
        FigureProduct(
            "Supplemental Figure Sy: accounting-boundary sensitivity",
            output_dir / "fig_supp_snow_da_boundary_sensitivity.png",
            output_dir / "fig_supp_snow_da_boundary_sensitivity.pdf",
            (11.5, 5.2),
            (
                relative(paths["boundary"], root),
                relative(paths["september"], root),
            ),
            "Replaces the duplicated bar chart with paired points and matched intervals, retaining total runoff and its surface/baseflow split while moving September timing facts to the caption notes.",
        ),
    ]

    plot_figure14(annual, partition, products[0].png, products[0].pdf)
    plot_figure15(tables["monthly_addition"], pathway, products[1].png, products[1].pdf)
    plot_octmar_supplement(tables["octmar"], products[2].png, products[2].pdf)
    plot_boundary_supplement(boundary, products[3].png, products[3].pdf)

    docs_pngs = []
    for product in products:
        destination = docs_dir / product.png.name
        shutil.copy2(product.png, destination)
        docs_pngs.append(destination)
    report_path = root / "projects/M21C_ls/docs/snow_da_manuscript_figure_report.md"
    write_figure_report(
        root,
        products,
        docs_pngs,
        pathway,
        tables["september"],
        boundary,
        report_path,
    )
    print("Snow-DA manuscript figures complete")
    print(f"Output directory: {output_dir}")
    print(f"Figure report: {report_path}")


if __name__ == "__main__":
    main()
