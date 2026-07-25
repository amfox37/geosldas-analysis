#!/usr/bin/env python3
"""Plot AZ-domain ObsFcstAna OmF and observation-count summaries.

This follows the OFA plotting convention used in the CYGNSS and M21C notebooks:
use ``N_data`` to mask noisy OmF diagnostics, but plot ``N_data`` itself
without that mask. The default inputs are the local OLv8_M36_all_sensors_AZ
summary products staged under ``data/omf_compare_sums``.
"""

from __future__ import annotations

import argparse
import calendar
import csv
import math
import struct
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from netCDF4 import Dataset

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:  # pragma: no cover - depends on local plotting env
    ccrs = None
    cfeature = None


REPO_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = REPO_ROOT / "projects" / "CYGNSS_L1_AZ"
DESCRIBE_DIR = PROJECT_ROOT / "OLv8_M36_all_sensors_AZ_describe"
DEFAULT_STATS_DIR = REPO_ROOT / "data" / "omf_compare_sums" / "OL_AZ_all_sensors"
DEFAULT_TEMPORAL = DEFAULT_STATS_DIR / "OL_AZ_all_sensors_temporal_stats.nc4"
DEFAULT_MONTHLY = DEFAULT_STATS_DIR / "OL_AZ_all_sensors_monthly_stats.nc4"
DEFAULT_TILECOORD = DESCRIBE_DIR / "OLv8_M36_all_sensors_AZ.ldas_tilecoord.bin"
DEFAULT_OBSPARAM = DESCRIBE_DIR / "OLv8_M36_all_sensors_AZ.ldas_obsparam.txt"
DEFAULT_OUT_DIR = PROJECT_ROOT / "output" / "omf_summary_figures"

SPECIES_FAMILIES = {
    "SMOS": range(0, 4),
    "SMAP": range(4, 8),
    "ASCAT": range(8, 11),
    "CYGNSS": range(11, 13),
}

CYGNSS_L1_INDEX0 = 12
PANEL_LABELS = tuple(f"({chr(ord('a') + i)})" for i in range(26))


@dataclass(frozen=True)
class Species:
    index0: int
    name: str
    units: str
    fcst_units: str

    @property
    def label(self) -> str:
        return f"{self.index0 + 1}. {self.name}"


@dataclass(frozen=True)
class SpeciesGroup:
    name: str
    indices: tuple[int, ...]
    color: str
    linestyle: str = "-"
    linewidth: float = 2.5


PAPER_GROUPS = (
    SpeciesGroup("ASCAT", (8, 9, 10), "#ff7f0e", "--", 2.6),
    SpeciesGroup("SMOS", (0, 1, 2, 3), "#4daf4a", "-.", 2.5),
    SpeciesGroup("SMAP", (4, 5, 6, 7), "#9c755f", "-", 2.7),
    SpeciesGroup("CYGNSS L3", (11,), "#ff9da7", "-", 2.6),
    SpeciesGroup("CYGNSS L1", (12,), "#7b3294", ":", 2.8),
)

FIG01_OBS_DAY_GROUPS = (
    PAPER_GROUPS[0],
    PAPER_GROUPS[1],
    PAPER_GROUPS[2],
    SpeciesGroup("SMAP H", (4, 5), "#9c755f", "-", 2.5),
    PAPER_GROUPS[3],
    PAPER_GROUPS[4],
)


def read_exact(fp, nbytes: int) -> bytes:
    data = fp.read(nbytes)
    if len(data) != nbytes:
        raise EOFError(f"Expected {nbytes} bytes, got {len(data)}")
    return data


def read_record_tag(fp) -> int:
    return struct.unpack("<i", read_exact(fp, 4))[0]


def read_tilecoord(path: Path) -> dict[str, np.ndarray]:
    int_fields = {"tile_id", "typ", "pfaf", "i_indg", "j_indg"}
    fields = [
        "tile_id",
        "typ",
        "pfaf",
        "com_lon",
        "com_lat",
        "min_lon",
        "max_lon",
        "min_lat",
        "max_lat",
        "i_indg",
        "j_indg",
        "frac_cell",
        "frac_pfaf",
        "area",
        "elev",
    ]
    out: dict[str, np.ndarray] = {}
    with path.open("rb") as fp:
        tag = read_record_tag(fp)
        if tag != 4:
            raise ValueError(f"{path} N_tile record tag is {tag}, expected 4")
        n_tile = struct.unpack("<i", read_exact(fp, 4))[0]
        end_tag = read_record_tag(fp)
        if end_tag != tag:
            raise ValueError(f"{path} N_tile record tags do not match")
        out["N_tile"] = np.asarray(n_tile, dtype=np.int32)

        for field in fields:
            dtype = np.dtype("<i4") if field in int_fields else np.dtype("<f4")
            expected = n_tile * dtype.itemsize
            tag = read_record_tag(fp)
            if tag != expected:
                raise ValueError(f"{path} field {field} tag is {tag}, expected {expected}")
            out[field] = np.frombuffer(read_exact(fp, expected), dtype=dtype).copy()
            end_tag = read_record_tag(fp)
            if end_tag != tag:
                raise ValueError(f"{path} field {field} record tags do not match")
    return out


def clean_value(value: str) -> str:
    value = " ".join(value.strip().split())
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1]
    return value


def parse_obsparam(path: Path) -> list[Species]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    n_species = int(lines[0])
    block_len = (len(lines) - 1) // n_species
    if block_len < 30:
        raise ValueError(f"Unexpected obsparam block length {block_len} in {path}")

    species = []
    for index0 in range(n_species):
        block = [clean_value(v) for v in lines[1 + index0 * block_len : 1 + (index0 + 1) * block_len]]
        species.append(
            Species(
                index0=index0,
                name=block[0],
                units=block[18],
                fcst_units=block[20],
            )
        )
    return species


def read_nc_variables(path: Path, names: list[str]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    with Dataset(path) as ds:
        for name in names:
            values = ds.variables[name][:]
            if hasattr(values, "filled"):
                values = values.filled(np.nan)
            out[name] = np.asarray(values, dtype=float)
    return out


def read_months(path: Path) -> np.ndarray:
    with Dataset(path) as ds:
        yyyymm = np.asarray(ds.variables["yyyymm"][:], dtype=int)
    return np.asarray([datetime.strptime(str(v), "%Y%m") for v in yyyymm])


def family_for(index0: int) -> str:
    for name, indices in SPECIES_FAMILIES.items():
        if index0 in indices:
            return name
    return "Other"


def finite_percentile(values: np.ndarray, pct: float, default: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return default
    return float(np.nanpercentile(finite, pct))


def positive_norm(values: np.ndarray, percentile: float = 98.0) -> Normalize:
    vmax = finite_percentile(values, percentile, 1.0)
    return Normalize(vmin=0.0, vmax=max(vmax, 1.0e-12))


def range_norm(values: np.ndarray, lower: float = 2.0, upper: float = 98.0) -> Normalize:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    vmin = float(np.nanpercentile(finite, lower))
    vmax = float(np.nanpercentile(finite, upper))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return Normalize(vmin=0.0, vmax=1.0)
    if math.isclose(vmin, vmax):
        delta = max(abs(vmin) * 0.05, 1.0e-6)
        vmin -= delta
        vmax += delta
    return Normalize(vmin=vmin, vmax=vmax)


def signed_norm(values: np.ndarray, percentile: float = 98.0) -> TwoSlopeNorm:
    vmax = finite_percentile(np.abs(values), percentile, 1.0)
    vmax = max(vmax, 1.0e-12)
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def nobs_norm(values: np.ndarray) -> Normalize | LogNorm:
    finite = values[np.isfinite(values) & (values > 0)]
    if finite.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    if float(np.nanmax(finite)) / max(float(np.nanmin(finite)), 1.0) >= 100.0:
        return LogNorm(vmin=max(float(np.nanmin(finite)), 1.0), vmax=float(np.nanmax(finite)))
    return Normalize(vmin=0.0, vmax=float(np.nanmax(finite)))


def save_figure(fig, out_path: Path, tight: bool = False) -> None:
    kwargs = {"dpi": 180}
    if tight:
        kwargs["bbox_inches"] = "tight"
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in create_collection",
            category=RuntimeWarning,
            module="shapely.creation",
        )
        fig.savefig(out_path, **kwargs)


def axes_grid(n_items: int, ncols: int = 4, panel_size: tuple[float, float] = (3.0, 2.6)):
    nrows = int(math.ceil(n_items / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_size[0] * ncols, panel_size[1] * nrows))
    return fig, np.atleast_1d(axes).ravel()


def set_az_axes(ax, extent: tuple[float, float, float, float]) -> None:
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.85", linewidth=0.5)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("")


def cartopy_kwargs() -> dict:
    if ccrs is None:
        return {}
    return {"subplot_kw": {"projection": ccrs.PlateCarree()}}


def map_axes_grid(n_items: int, ncols: int = 2, panel_size: tuple[float, float] = (5.5, 3.7)):
    nrows = int(math.ceil(n_items / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_size[0] * ncols, panel_size[1] * nrows),
        **cartopy_kwargs(),
    )
    return fig, np.atleast_1d(axes).ravel()


def set_map_axes(ax, extent: tuple[float, float, float, float]) -> None:
    if ccrs is None or not hasattr(ax, "set_extent"):
        set_az_axes(ax, extent)
        return

    transform = ccrs.PlateCarree()
    ax.set_extent(extent, crs=transform)
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="0.94", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.55, edgecolor="0.25", zorder=2)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.45, edgecolor="0.35", zorder=2)
    states = cfeature.NaturalEarthFeature(
        "cultural",
        "admin_1_states_provinces_lakes",
        "50m",
        facecolor="none",
    )
    ax.add_feature(states, linewidth=0.45, edgecolor="0.35", zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.35, color="0.80", alpha=0.8)


def scatter_kwargs() -> dict:
    if ccrs is None:
        return {}
    return {"transform": ccrs.PlateCarree()}


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def value_stats(values: np.ndarray, mask: np.ndarray | None = None) -> dict[str, float | int]:
    valid = np.isfinite(values)
    if mask is not None:
        valid &= mask
    selected = values[valid]
    if selected.size == 0:
        return {"n": 0, "mean": np.nan, "std": np.nan, "min": np.nan, "max": np.nan}
    return {
        "n": int(selected.size),
        "mean": float(np.nanmean(selected)),
        "std": float(np.nanstd(selected)),
        "min": float(np.nanmin(selected)),
        "max": float(np.nanmax(selected)),
    }


def format_panel_stat(prefix: str, stats: dict[str, float | int], units: str = "") -> str:
    if stats["n"] == 0:
        return f"{prefix}: no data"
    suffix = f" {units}" if units else ""
    return f"{prefix}: {stats['mean']:.3g} +/- {stats['std']:.3g}{suffix} (n={stats['n']})"


def panel_mean_text(
    ax,
    values: np.ndarray,
    units: str = "",
    native_mask: np.ndarray | None = None,
    support_mask: np.ndarray | None = None,
    support_label: str = "CYG L1",
) -> None:
    lines = [format_panel_stat("all", value_stats(values, native_mask), units)]
    if support_mask is not None:
        mask = support_mask if native_mask is None else support_mask & native_mask
        lines.append(format_panel_stat(support_label, value_stats(values, mask), units))
    text = "\n".join(lines)
    ax.text(
        0.02,
        0.04,
        text,
        transform=ax.transAxes,
        fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="none", alpha=0.82),
        ha="left",
        va="bottom",
    )


def group_units(group: SpeciesGroup, species: list[Species], source: str = "obs") -> str:
    attr = "fcst_units" if source == "fcst" else "units"
    units = {getattr(species[index0], attr) for index0 in group.indices}
    return units.pop() if len(units) == 1 else "native units"


def days_covered(months: np.ndarray) -> int:
    return sum(calendar.monthrange(int(stamp.year), int(stamp.month))[1] for stamp in months)


def group_counts(n_data: np.ndarray, group: SpeciesGroup) -> np.ndarray:
    return np.nansum(n_data[:, group.indices], axis=1)


def cygnss_l1_support_mask(n_data: np.ndarray, nmin: int) -> np.ndarray:
    return np.isfinite(n_data[:, CYGNSS_L1_INDEX0]) & (n_data[:, CYGNSS_L1_INDEX0] >= nmin)


def weighted_group_mean(values: np.ndarray, n_data: np.ndarray, group: SpeciesGroup, nmin: int) -> np.ndarray:
    vals = values[:, group.indices].astype(float, copy=True)
    weights = n_data[:, group.indices].astype(float, copy=True)
    valid = np.isfinite(vals) & (weights >= nmin)
    weighted_sum = np.where(valid, vals * weights, 0.0).sum(axis=1)
    total_weight = np.where(valid, weights, 0.0).sum(axis=1)
    return np.divide(
        weighted_sum,
        total_weight,
        out=np.full(total_weight.shape, np.nan, dtype=float),
        where=total_weight > 0,
    )


def group_shared_range_norms(
    first_values: np.ndarray,
    second_values: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    groups: tuple[SpeciesGroup, ...],
    nmin: int,
    first_unit_source: str = "obs",
    second_unit_source: str = "fcst",
) -> dict[str, Normalize]:
    norms: dict[str, Normalize] = {}
    for group in groups:
        first_units = group_units(group, species, source=first_unit_source)
        second_units = group_units(group, species, source=second_unit_source)
        if first_units != second_units:
            continue
        first_group_values = weighted_group_mean(first_values, n_data, group, nmin)
        second_group_values = weighted_group_mean(second_values, n_data, group, nmin)
        norms[group.name] = range_norm(np.concatenate([first_group_values, second_group_values]))
    return norms


def month_edges(months: np.ndarray) -> tuple[datetime, datetime]:
    first = months[0]
    last = months[-1]
    if last.month == 12:
        end = datetime(last.year + 1, 1, 1)
    else:
        end = datetime(last.year, last.month + 1, 1)
    return datetime(first.year, first.month, 1), end


def shade_years(ax, months: np.ndarray) -> None:
    colors = ("#f4efe5", "#e8f2f8")
    start, end = month_edges(months)
    for i, year in enumerate(range(start.year, end.year + 1)):
        left = datetime(year, 1, 1)
        right = datetime(year + 1, 1, 1)
        if right <= start or left >= end:
            continue
        ax.axvspan(max(left, start), min(right, end), color=colors[i % len(colors)], alpha=0.35, zorder=0)


def format_date_axis(ax, months: np.ndarray) -> None:
    ax.set_xlim(*month_edges(months))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("right")


def padded_limits(values: np.ndarray, pad_fraction: float = 0.08) -> tuple[float, float] | None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None
    ymin = float(np.nanmin(finite))
    ymax = float(np.nanmax(finite))
    if math.isclose(ymin, ymax):
        pad = max(abs(ymin) * pad_fraction, 1.0e-6)
    else:
        pad = (ymax - ymin) * pad_fraction
    return ymin - pad, ymax + pad


def set_tight_xlim(ax, values: np.ndarray, pad_fraction: float = 0.08) -> None:
    limits = padded_limits(values, pad_fraction)
    if limits is not None:
        ax.set_xlim(*limits)


def set_tight_ylim(ax, values: np.ndarray, pad_fraction: float = 0.08) -> None:
    limits = padded_limits(values, pad_fraction)
    if limits is not None:
        ax.set_ylim(*limits)


def plot_group_obs_day_maps(
    n_data: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    months: np.ndarray,
    groups: tuple[SpeciesGroup, ...],
    extent: tuple[float, float, float, float],
    support_mask: np.ndarray,
    out_path: Path,
) -> None:
    n_days = days_covered(months)
    values_by_group = {group.name: group_counts(n_data, group) / n_days for group in groups}
    stacked = np.concatenate([values[np.isfinite(values)] for values in values_by_group.values()])
    norm = Normalize(vmin=0.0, vmax=max(float(np.nanpercentile(stacked, 98)), 1.0e-12))
    cmap = "viridis"

    fig, axes = map_axes_grid(len(groups), ncols=2)
    for i, (group, ax) in enumerate(zip(groups, axes)):
        vals = values_by_group[group.name]
        valid = np.isfinite(vals) & (vals > 0)
        ax.scatter(
            lon[valid],
            lat[valid],
            c=vals[valid],
            s=22,
            marker="s",
            cmap=cmap,
            norm=norm,
            linewidths=0,
            **scatter_kwargs(),
        )
        set_map_axes(ax, extent)
        ax.set_title(group.name, fontsize=12)
        panel_label(ax, PANEL_LABELS[i])
        panel_mean_text(ax, vals, "d$^{-1}$", support_mask=support_mask)

    for ax in axes[len(groups) :]:
        ax.axis("off")

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[: len(groups)], orientation="horizontal", fraction=0.045, pad=0.07)
    cbar.set_label("Monitored observations per day")
    fig.suptitle("Mean monitored observations per day by observing system", fontsize=14)
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_group_metric_maps(
    values: np.ndarray,
    n_data: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    species: list[Species],
    groups: tuple[SpeciesGroup, ...],
    extent: tuple[float, float, float, float],
    nmin: int,
    title_label: str,
    color_mode: str,
    support_mask: np.ndarray,
    out_path: Path,
    unit_source: str = "obs",
    norms_by_group: dict[str, Normalize] | None = None,
) -> None:
    fig, axes = map_axes_grid(len(groups), ncols=2)
    for i, (group, ax) in enumerate(zip(groups, axes)):
        vals = weighted_group_mean(values, n_data, group, nmin)
        valid = np.isfinite(vals)
        if norms_by_group is not None and group.name in norms_by_group:
            norm = norms_by_group[group.name]
            cmap = "viridis"
        elif color_mode == "signed":
            norm = signed_norm(vals)
            cmap = "RdBu_r"
        elif color_mode == "range":
            norm = range_norm(vals)
            cmap = "viridis"
        else:
            norm = positive_norm(vals)
            cmap = "magma"
        sc = ax.scatter(
            lon[valid],
            lat[valid],
            c=vals[valid],
            s=22,
            marker="s",
            cmap=cmap,
            norm=norm,
            linewidths=0,
            **scatter_kwargs(),
        )
        set_map_axes(ax, extent)
        units = group_units(group, species, source=unit_source)
        ax.set_title(group.name, fontsize=12)
        panel_label(ax, PANEL_LABELS[i])
        panel_mean_text(ax, vals, units, support_mask=support_mask)
        cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.060, pad=0.055)
        cbar.set_label(units, fontsize=8.5)
        cbar.ax.tick_params(labelsize=7)

    for ax in axes[len(groups) :]:
        ax.axis("off")

    fig.suptitle(f"Full-period {title_label} by observing system (N_data >= {nmin})", fontsize=14)
    fig.subplots_adjust(hspace=0.42, wspace=0.20)
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_group_count_time_series(
    months: np.ndarray,
    n_data: np.ndarray,
    groups: tuple[SpeciesGroup, ...],
    out_path: Path,
) -> None:
    group_series = {group.name: group_counts(n_data, group) for group in groups}
    total_obs = np.nansum(n_data, axis=1)

    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(11, 6.8), sharex=True, height_ratios=[1.05, 1.0])
    for ax in (ax_top, ax_bottom):
        shade_years(ax, months)
        ax.grid(True, color="0.85", linewidth=0.6, zorder=0)

    ax_top.plot(months, total_obs, color="0.15", linewidth=2.2)
    ax_top.set_ylabel("Total obs month$^{-1}$")
    ax_top.set_title("Monthly monitored observation counts")
    panel_label(ax_top, PANEL_LABELS[0])

    for group in groups:
        vals = np.where(group_series[group.name] > 0, group_series[group.name], np.nan)
        mean_val = np.nanmean(vals)
        ax_bottom.plot(
            months,
            vals,
            label=f"{group.name} ({mean_val:.1e})",
            color=group.color,
            linestyle=group.linestyle,
            linewidth=group.linewidth,
        )
    ax_bottom.set_ylabel("Obs month$^{-1}$")
    ax_bottom.set_xlabel("Year")
    panel_label(ax_bottom, PANEL_LABELS[1])
    ax_bottom.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.42), frameon=True, fontsize=9)
    format_date_axis(ax_bottom, months)

    fig.subplots_adjust(hspace=0.18, bottom=0.25)
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_group_metric_time_series(
    months: np.ndarray,
    values: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    groups: tuple[SpeciesGroup, ...],
    nmin: int,
    title_label: str,
    zero_line: bool,
    out_path: Path,
    unit_source: str = "obs",
) -> None:
    fig, axes = plt.subplots(len(groups), 1, figsize=(11, 1.95 * len(groups) + 1.2), sharex=True)
    axes = np.atleast_1d(axes)
    for i, (group, ax) in enumerate(zip(groups, axes)):
        shade_years(ax, months)
        vals = weighted_group_mean(values, n_data, group, nmin)
        mean_val = np.nanmean(vals)
        units = group_units(group, species, source=unit_source)
        ax.plot(
            months,
            vals,
            color=group.color,
            linestyle=group.linestyle,
            linewidth=group.linewidth,
            label=f"{group.name} ({mean_val:.3g} {units})",
        )
        ax.scatter(months, vals, color=group.color, s=16, zorder=3)
        if zero_line:
            ax.axhline(0.0, color="0.15", linestyle=":", linewidth=1.0)
        ax.grid(True, color="0.85", linewidth=0.6)
        ax.set_ylabel(units)
        ax.legend(loc="upper right", frameon=True, fontsize=9)
        panel_label(ax, PANEL_LABELS[i])

    axes[0].set_title(f"Monthly {title_label} by observing system (N_data >= {nmin})")
    axes[-1].set_xlabel("Year")
    format_date_axis(axes[-1], months)
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_group_observed_forecast_time_series(
    months: np.ndarray,
    obs_values: np.ndarray,
    fcst_values: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    groups: tuple[SpeciesGroup, ...],
    nmin: int,
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(len(groups), 1, figsize=(11, 1.95 * len(groups) + 1.2), sharex=True)
    axes = np.atleast_1d(axes)

    for i, (group, ax) in enumerate(zip(groups, axes)):
        shade_years(ax, months)
        obs_vals = weighted_group_mean(obs_values, n_data, group, nmin)
        fcst_vals = weighted_group_mean(fcst_values, n_data, group, nmin)
        obs_units = group_units(group, species, source="obs")
        fcst_units = group_units(group, species, source="fcst")

        obs_line = ax.plot(
            months,
            obs_vals,
            color=group.color,
            linestyle="-",
            linewidth=group.linewidth,
            marker="o",
            markersize=3.6,
            label=f"O ({np.nanmean(obs_vals):.3g} {obs_units})",
        )
        ax.set_ylabel(obs_units)

        if fcst_units == obs_units:
            fcst_ax = ax
        else:
            fcst_ax = ax.twinx()
            fcst_ax.set_ylabel(fcst_units)
            fcst_ax.tick_params(labelsize=8)

        fcst_line = fcst_ax.plot(
            months,
            fcst_vals,
            color="0.18",
            linestyle="--",
            linewidth=2.0,
            marker="s",
            markersize=3.2,
            label=f"F ({np.nanmean(fcst_vals):.3g} {fcst_units})",
        )

        ax.grid(True, color="0.85", linewidth=0.6)
        ax.set_title(group.name, fontsize=10, loc="left")
        panel_label(ax, PANEL_LABELS[i])
        handles = obs_line + fcst_line
        labels = [handle.get_label() for handle in handles]
        ax.legend(handles, labels, loc="upper right", frameon=True, fontsize=8.5)

    axes[0].set_title(f"Monthly observed and forecast mean by observing system (N_data >= {nmin})", fontsize=12)
    axes[-1].set_xlabel("Year")
    format_date_axis(axes[-1], months)
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_tile_elevation_map(
    elev: np.ndarray,
    min_lon: np.ndarray,
    max_lon: np.ndarray,
    min_lat: np.ndarray,
    max_lat: np.ndarray,
    i_indg: np.ndarray,
    j_indg: np.ndarray,
    extent: tuple[float, float, float, float],
    support_mask: np.ndarray,
    out_path: Path,
) -> None:
    vals = elev.astype(float, copy=True)
    native_mask = np.isfinite(vals)
    norm = range_norm(vals, lower=0.0, upper=100.0)

    cols = sorted(np.unique(i_indg), key=lambda col: float(np.nanmedian(min_lon[i_indg == col])))
    rows = sorted(np.unique(j_indg), key=lambda row: float(np.nanmedian(min_lat[j_indg == row])))
    col_lookup = {col: n for n, col in enumerate(cols)}
    row_lookup = {row: n for n, row in enumerate(rows)}

    grid = np.full((len(rows), len(cols)), np.nan, dtype=float)
    for tile_value, col, row in zip(vals, i_indg, j_indg):
        grid[row_lookup[row], col_lookup[col]] = tile_value

    col_min = np.asarray([np.nanmedian(min_lon[i_indg == col]) for col in cols], dtype=float)
    col_max = np.asarray([np.nanmedian(max_lon[i_indg == col]) for col in cols], dtype=float)
    row_min = np.asarray([np.nanmedian(min_lat[j_indg == row]) for row in rows], dtype=float)
    row_max = np.asarray([np.nanmedian(max_lat[j_indg == row]) for row in rows], dtype=float)

    lon_edges = np.empty(len(cols) + 1, dtype=float)
    lat_edges = np.empty(len(rows) + 1, dtype=float)
    lon_edges[0] = col_min[0]
    lon_edges[-1] = col_max[-1]
    lat_edges[0] = row_min[0]
    lat_edges[-1] = row_max[-1]
    lon_edges[1:-1] = 0.5 * (col_max[:-1] + col_min[1:])
    lat_edges[1:-1] = 0.5 * (row_max[:-1] + row_min[1:])

    fig, ax = plt.subplots(figsize=(7.2, 5.3), **cartopy_kwargs())
    sc = ax.pcolormesh(
        lon_edges,
        lat_edges,
        grid,
        cmap="terrain",
        norm=norm,
        shading="flat",
        **scatter_kwargs(),
    )
    set_map_axes(ax, extent)
    ax.set_title("GEOSldas tile elevation over AZ", fontsize=13)
    panel_mean_text(ax, vals, "m", native_mask=native_mask, support_mask=support_mask)
    cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.060, pad=0.065)
    cbar.set_label("Elevation (m)")
    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_cygnss_forecast_dual_axis(
    months: np.ndarray,
    f_mean: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    nmin: int,
    out_path: Path,
) -> None:
    l3_index0 = 11
    l1_index0 = CYGNSS_L1_INDEX0
    l3 = f_mean[:, l3_index0].astype(float, copy=True)
    l1 = f_mean[:, l1_index0].astype(float, copy=True)
    l3[n_data[:, l3_index0] < nmin] = np.nan
    l1[n_data[:, l1_index0] < nmin] = np.nan

    paired = np.isfinite(l3) & np.isfinite(l1)
    corr = np.nan
    if np.count_nonzero(paired) >= 2:
        corr = float(np.corrcoef(l3[paired], l1[paired])[0, 1])

    l3_color = "#1f78b4"
    l1_color = "#7b3294"

    fig = plt.figure(figsize=(11, 9.2))
    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.05, 1.10],
        width_ratios=[0.65, 1.90, 0.65],
        hspace=0.44,
    )
    ax_l3 = fig.add_subplot(gs[0, :])
    ax_scatter = fig.add_subplot(gs[1, 1])
    shade_years(ax_l3, months)
    ax_l3.grid(True, color="0.85", linewidth=0.6, zorder=0)

    l3_line = ax_l3.plot(
        months,
        l3,
        color=l3_color,
        linewidth=2.6,
        marker="o",
        markersize=4.8,
        label=f"{species[l3_index0].name} F",
    )
    ax_l3.set_ylabel(f"CYGNSS L3 forecast ({species[l3_index0].fcst_units})", color=l3_color)
    ax_l3.tick_params(axis="y", labelcolor=l3_color)
    set_tight_ylim(ax_l3, l3)

    ax_l1 = ax_l3.twinx()
    l1_line = ax_l1.plot(
        months,
        l1,
        color=l1_color,
        linewidth=2.6,
        linestyle="--",
        marker="s",
        markersize=4.4,
        label=f"{species[l1_index0].name} F",
    )
    ax_l1.set_ylabel(f"CYGNSS L1 forecast ({species[l1_index0].fcst_units})", color=l1_color)
    ax_l1.tick_params(axis="y", labelcolor=l1_color)
    set_tight_ylim(ax_l1, l1)

    handles = l3_line + l1_line
    labels = [handle.get_label() for handle in handles]
    ax_l3.legend(handles, labels, loc="upper left", frameon=True, fontsize=9)
    panel_label(ax_l3, PANEL_LABELS[0])

    title = f"Monthly CYGNSS forecast means over AZ (N_data >= {nmin})"
    if np.isfinite(corr):
        title += f"; r = {corr:.2f}"
    ax_l3.set_title(title)
    format_date_axis(ax_l3, months)

    ax_scatter.scatter(
        l3[paired],
        l1[paired],
        color="#2a9d8f",
        s=58,
        edgecolors="white",
        linewidths=0.7,
        zorder=3,
    )
    if np.count_nonzero(paired) >= 2:
        slope, intercept = np.polyfit(l3[paired], l1[paired], 1)
        fit_x = np.linspace(float(np.nanmin(l3[paired])), float(np.nanmax(l3[paired])), 100)
        ax_scatter.plot(
            fit_x,
            slope * fit_x + intercept,
            color="0.20",
            linestyle=":",
            linewidth=1.4,
            label=f"linear fit, slope={slope:.1f} dB/(m3 m-3)",
            zorder=2,
        )
        ax_scatter.legend(loc="upper left", frameon=True, fontsize=8.5)
    set_tight_xlim(ax_scatter, l3)
    set_tight_ylim(ax_scatter, l1)
    ax_scatter.grid(True, color="0.85", linewidth=0.6, zorder=0)
    ax_scatter.set_xlabel(f"CYGNSS L3 forecast ({species[l3_index0].fcst_units})")
    ax_scatter.set_ylabel(f"CYGNSS L1 forecast ({species[l1_index0].fcst_units})")
    ax_scatter.set_title("Forecast-mean scatter", fontsize=11, pad=10)
    ax_scatter.set_box_aspect(1)
    panel_label(ax_scatter, PANEL_LABELS[1])

    save_figure(fig, out_path, tight=True)
    plt.close(fig)


def plot_species_maps(
    values: np.ndarray,
    n_data: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    species: list[Species],
    extent: tuple[float, float, float, float],
    nmin: int,
    variable: str,
    support_mask: np.ndarray,
    out_path: Path,
) -> None:
    is_nobs = variable == "N_data"
    is_omf_mean = variable == "OmF_mean"
    is_o_mean = variable == "O_mean"
    fig, axes = axes_grid(len(species), ncols=4)
    for sp, ax in zip(species, axes):
        vals = values[:, sp.index0].copy()
        if not is_nobs:
            vals[n_data[:, sp.index0] < nmin] = np.nan
        if is_nobs:
            norm = nobs_norm(vals)
            cmap = "viridis"
        elif is_omf_mean:
            norm = signed_norm(vals)
            cmap = "RdBu_r"
        elif is_o_mean:
            norm = range_norm(vals)
            cmap = "viridis"
        else:
            norm = positive_norm(vals)
            cmap = "magma"
        valid = np.isfinite(vals)
        if is_nobs:
            valid &= vals > 0
        sc = ax.scatter(lon[valid], lat[valid], c=vals[valid], s=18, marker="s", cmap=cmap, norm=norm, linewidths=0)
        set_az_axes(ax, extent)
        ax.set_title(sp.label, fontsize=8)
        panel_mean_text(ax, vals, sp.units if not is_nobs else "count", native_mask=valid, support_mask=support_mask)
        cb = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.06, pad=0.08)
        cb.ax.tick_params(labelsize=6)
    for ax in axes[len(species) :]:
        ax.axis("off")
    label = "N_data" if is_nobs else f"{variable} masked where N_data < {nmin}"
    fig.suptitle(f"OLv8_M36_all_sensors_AZ {label}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, out_path)
    plt.close(fig)


def plot_species_time_series(
    months: np.ndarray,
    values: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    nmin: int,
    variable: str,
    out_path: Path,
) -> None:
    is_nobs = variable == "N_data"
    is_omf_mean = variable == "OmF_mean"
    is_o_mean = variable == "O_mean"
    fig, axes = axes_grid(len(species), ncols=4, panel_size=(3.2, 2.2))
    for sp, ax in zip(species, axes):
        vals = values[:, sp.index0].copy()
        if not is_nobs:
            vals[n_data[:, sp.index0] < nmin] = np.nan
        color = "#2563eb" if is_nobs else "#b2182b" if is_omf_mean else "#0f766e" if is_o_mean else "#b45309"
        ax.plot(months, vals, color=color, linewidth=1.6)
        ax.scatter(months, vals, color=color, s=9)
        if is_omf_mean:
            ax.axhline(0.0, color="0.2", linestyle=":", linewidth=0.8)
        mean_text = "no data" if not np.isfinite(vals).any() else f"mean {np.nanmean(vals):.3g}"
        unit = "" if is_nobs else f" ({sp.units})"
        ax.set_title(f"{sp.label}\n{mean_text}{unit}", fontsize=8)
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_ha("right")
    for ax in axes[len(species) :]:
        ax.axis("off")
    ylabel = "N_data" if is_nobs else f"{variable}, native units; masked where N_data < {nmin}"
    fig.suptitle(f"OLv8_M36_all_sensors_AZ monthly {ylabel}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, out_path)
    plt.close(fig)


def write_monthly_csv(
    path: Path,
    months: np.ndarray,
    o_mean: np.ndarray,
    f_mean: np.ndarray,
    omf_mean: np.ndarray,
    omf_stdv: np.ndarray,
    n_data: np.ndarray,
    species: list[Species],
    nmin: int,
) -> None:
    with path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(
            [
                "yyyymm",
                "species_index",
                "species",
                "family",
                "units",
                "fcst_units",
                "O_mean",
                "O_mean_masked",
                "F_mean",
                "F_mean_masked",
                "OmF_mean",
                "OmF_mean_masked",
                "OmF_stdv",
                "OmF_stdv_masked",
                "N_data",
            ]
        )
        for t, stamp in enumerate(months):
            yyyymm = stamp.strftime("%Y%m")
            for sp in species:
                count = n_data[t, sp.index0]
                obs_value = o_mean[t, sp.index0]
                fcst_value = f_mean[t, sp.index0]
                mean_value = omf_mean[t, sp.index0]
                stdv_value = omf_stdv[t, sp.index0]
                obs_masked = obs_value if count >= nmin else np.nan
                fcst_masked = fcst_value if count >= nmin else np.nan
                mean_masked = mean_value if count >= nmin else np.nan
                stdv_masked = stdv_value if count >= nmin else np.nan
                writer.writerow(
                    [
                        yyyymm,
                        sp.index0 + 1,
                        sp.name,
                        family_for(sp.index0),
                        sp.units,
                        sp.fcst_units,
                        obs_value,
                        obs_masked,
                        fcst_value,
                        fcst_masked,
                        mean_value,
                        mean_masked,
                        stdv_value,
                        stdv_masked,
                        count,
                    ]
                )


def append_map_stats(
    rows: list[list[object]],
    figure: str,
    panel: str,
    metric: str,
    units: str,
    values: np.ndarray,
    support_mask: np.ndarray,
    native_mask: np.ndarray | None = None,
) -> None:
    supports = (
        ("all_valid", native_mask),
        ("cygnss_l1_valid", support_mask if native_mask is None else native_mask & support_mask),
    )
    for support_name, mask in supports:
        stats = value_stats(values, mask)
        rows.append(
            [
                figure,
                panel,
                metric,
                units,
                support_name,
                stats["n"],
                stats["mean"],
                stats["std"],
                stats["min"],
                stats["max"],
            ]
        )


def write_map_support_csv(
    path: Path,
    temporal: dict[str, np.ndarray],
    months: np.ndarray,
    species: list[Species],
    groups: tuple[SpeciesGroup, ...],
    obs_day_groups: tuple[SpeciesGroup, ...],
    elevation: np.ndarray,
    nmin: int,
    support_mask: np.ndarray,
) -> None:
    rows: list[list[object]] = []
    n_days = days_covered(months)

    for group in obs_day_groups:
        append_map_stats(
            rows,
            "az_fig01_mean_monitored_observations_per_day_by_system.png",
            group.name,
            "N_data_per_day",
            "d-1",
            group_counts(temporal["N_data"], group) / n_days,
            support_mask,
        )

    append_map_stats(
        rows,
        "az_fig11_tile_elevation.png",
        "AZ tiles",
        "elev",
        "m",
        elevation,
        support_mask,
        native_mask=np.isfinite(elevation),
    )

    group_metrics = (
        ("az_fig03_full_period_omf_stddev_by_system.png", "OmF_stdv", "obs"),
        ("az_fig05_full_period_omf_mean_by_system.png", "OmF_mean", "obs"),
        ("az_fig07_full_period_o_mean_by_system.png", "O_mean", "obs"),
        ("az_fig09_full_period_f_mean_by_system.png", "F_mean", "fcst"),
    )
    for figure, metric, unit_source in group_metrics:
        for group in groups:
            append_map_stats(
                rows,
                figure,
                group.name,
                metric,
                group_units(group, species, source=unit_source),
                weighted_group_mean(temporal[metric], temporal["N_data"], group, nmin),
                support_mask,
            )

    species_metrics = (
        ("az_omf_stdv_maps_by_species.png", "OmF_stdv"),
        ("az_omf_mean_maps_by_species.png", "OmF_mean"),
        ("az_o_mean_maps_by_species.png", "O_mean"),
    )
    for figure, metric in species_metrics:
        for sp in species:
            vals = temporal[metric][:, sp.index0].copy()
            vals[temporal["N_data"][:, sp.index0] < nmin] = np.nan
            append_map_stats(rows, figure, sp.name, metric, sp.units, vals, support_mask)

    for sp in species:
        vals = temporal["N_data"][:, sp.index0].copy()
        native_mask = np.isfinite(vals) & (vals > 0)
        append_map_stats(
            rows,
            "az_nobs_maps_by_species.png",
            sp.name,
            "N_data",
            "count",
            vals,
            support_mask,
            native_mask=native_mask,
        )

    with path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["figure", "panel", "metric", "units", "support", "n_tile", "mean", "std", "min", "max"])
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temporal-stats", type=Path, default=DEFAULT_TEMPORAL)
    parser.add_argument("--monthly-stats", type=Path, default=DEFAULT_MONTHLY)
    parser.add_argument("--tilecoord", type=Path, default=DEFAULT_TILECOORD)
    parser.add_argument("--obsparam", type=Path, default=DEFAULT_OBSPARAM)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--nmin", type=int, default=20, help="Minimum N_data for plotting OmF diagnostics.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species = parse_obsparam(args.obsparam)
    tilecoord = read_tilecoord(args.tilecoord)
    lon = np.asarray(tilecoord["com_lon"], dtype=float)
    lat = np.asarray(tilecoord["com_lat"], dtype=float)
    elev = np.asarray(tilecoord["elev"], dtype=float)
    extent = (
        float(np.nanmin(tilecoord["min_lon"])) - 0.4,
        float(np.nanmax(tilecoord["max_lon"])) + 0.4,
        float(np.nanmin(tilecoord["min_lat"])) - 0.4,
        float(np.nanmax(tilecoord["max_lat"])) + 0.4,
    )

    temporal = read_nc_variables(args.temporal_stats, ["O_mean", "F_mean", "OmF_mean", "OmF_stdv", "N_data"])
    monthly = read_nc_variables(args.monthly_stats, ["O_mean", "F_mean", "OmF_mean", "OmF_stdv", "N_data"])
    months = read_months(args.monthly_stats)

    if temporal["N_data"].shape != (lon.size, len(species)):
        raise ValueError(
            "Temporal stats shape does not match tilecoord/species: "
            f"{temporal['N_data'].shape} vs ({lon.size}, {len(species)})"
        )

    support_mask = cygnss_l1_support_mask(temporal["N_data"], args.nmin)
    observed_forecast_range_norms = group_shared_range_norms(
        temporal["O_mean"],
        temporal["F_mean"],
        temporal["N_data"],
        species,
        PAPER_GROUPS,
        args.nmin,
    )

    plot_group_obs_day_maps(
        temporal["N_data"],
        lon,
        lat,
        months,
        FIG01_OBS_DAY_GROUPS,
        extent,
        support_mask,
        args.out_dir / "az_fig01_mean_monitored_observations_per_day_by_system.png",
    )
    plot_group_count_time_series(
        months,
        monthly["N_data"],
        PAPER_GROUPS,
        args.out_dir / "az_fig02_monthly_monitored_observation_counts.png",
    )
    plot_group_metric_maps(
        temporal["OmF_stdv"],
        temporal["N_data"],
        lon,
        lat,
        species,
        PAPER_GROUPS,
        extent,
        args.nmin,
        "O-F standard deviation",
        "positive",
        support_mask,
        args.out_dir / "az_fig03_full_period_omf_stddev_by_system.png",
    )
    plot_group_metric_time_series(
        months,
        monthly["OmF_stdv"],
        monthly["N_data"],
        species,
        PAPER_GROUPS,
        args.nmin,
        "O-F standard deviation",
        False,
        args.out_dir / "az_fig04_monthly_omf_stddev_by_system.png",
    )
    plot_group_metric_maps(
        temporal["OmF_mean"],
        temporal["N_data"],
        lon,
        lat,
        species,
        PAPER_GROUPS,
        extent,
        args.nmin,
        "O-F mean",
        "signed",
        support_mask,
        args.out_dir / "az_fig05_full_period_omf_mean_by_system.png",
    )
    plot_group_metric_time_series(
        months,
        monthly["OmF_mean"],
        monthly["N_data"],
        species,
        PAPER_GROUPS,
        args.nmin,
        "O-F mean",
        True,
        args.out_dir / "az_fig06_monthly_omf_mean_by_system.png",
    )
    plot_group_metric_maps(
        temporal["O_mean"],
        temporal["N_data"],
        lon,
        lat,
        species,
        PAPER_GROUPS,
        extent,
        args.nmin,
        "observed mean",
        "range",
        support_mask,
        args.out_dir / "az_fig07_full_period_o_mean_by_system.png",
        norms_by_group=observed_forecast_range_norms,
    )
    plot_group_observed_forecast_time_series(
        months,
        monthly["O_mean"],
        monthly["F_mean"],
        monthly["N_data"],
        species,
        PAPER_GROUPS,
        args.nmin,
        args.out_dir / "az_fig08_monthly_o_mean_by_system.png",
    )
    plot_group_metric_maps(
        temporal["F_mean"],
        temporal["N_data"],
        lon,
        lat,
        species,
        PAPER_GROUPS,
        extent,
        args.nmin,
        "forecast mean",
        "range",
        support_mask,
        args.out_dir / "az_fig09_full_period_f_mean_by_system.png",
        unit_source="fcst",
        norms_by_group=observed_forecast_range_norms,
    )
    plot_cygnss_forecast_dual_axis(
        months,
        monthly["F_mean"],
        monthly["N_data"],
        species,
        args.nmin,
        args.out_dir / "az_fig10_monthly_cygnss_l3_l1_forecast_mean.png",
    )
    plot_tile_elevation_map(
        elev,
        np.asarray(tilecoord["min_lon"], dtype=float),
        np.asarray(tilecoord["max_lon"], dtype=float),
        np.asarray(tilecoord["min_lat"], dtype=float),
        np.asarray(tilecoord["max_lat"], dtype=float),
        np.asarray(tilecoord["i_indg"], dtype=int),
        np.asarray(tilecoord["j_indg"], dtype=int),
        extent,
        support_mask,
        args.out_dir / "az_fig11_tile_elevation.png",
    )

    plot_species_maps(
        temporal["OmF_stdv"],
        temporal["N_data"],
        lon,
        lat,
        species,
        extent,
        args.nmin,
        "OmF_stdv",
        support_mask,
        args.out_dir / "az_omf_stdv_maps_by_species.png",
    )
    plot_species_maps(
        temporal["OmF_mean"],
        temporal["N_data"],
        lon,
        lat,
        species,
        extent,
        args.nmin,
        "OmF_mean",
        support_mask,
        args.out_dir / "az_omf_mean_maps_by_species.png",
    )
    plot_species_maps(
        temporal["O_mean"],
        temporal["N_data"],
        lon,
        lat,
        species,
        extent,
        args.nmin,
        "O_mean",
        support_mask,
        args.out_dir / "az_o_mean_maps_by_species.png",
    )
    plot_species_maps(
        temporal["N_data"],
        temporal["N_data"],
        lon,
        lat,
        species,
        extent,
        args.nmin,
        "N_data",
        support_mask,
        args.out_dir / "az_nobs_maps_by_species.png",
    )
    plot_species_time_series(
        months,
        monthly["OmF_stdv"],
        monthly["N_data"],
        species,
        args.nmin,
        "OmF_stdv",
        args.out_dir / "az_monthly_omf_stdv_by_species.png",
    )
    plot_species_time_series(
        months,
        monthly["OmF_mean"],
        monthly["N_data"],
        species,
        args.nmin,
        "OmF_mean",
        args.out_dir / "az_monthly_omf_mean_by_species.png",
    )
    plot_species_time_series(
        months,
        monthly["O_mean"],
        monthly["N_data"],
        species,
        args.nmin,
        "O_mean",
        args.out_dir / "az_monthly_o_mean_by_species.png",
    )
    plot_species_time_series(
        months,
        monthly["N_data"],
        monthly["N_data"],
        species,
        args.nmin,
        "N_data",
        args.out_dir / "az_monthly_nobs_by_species.png",
    )
    write_monthly_csv(
        args.out_dir / "az_monthly_species_summary.csv",
        months,
        monthly["O_mean"],
        monthly["F_mean"],
        monthly["OmF_mean"],
        monthly["OmF_stdv"],
        monthly["N_data"],
        species,
        args.nmin,
    )
    write_map_support_csv(
        args.out_dir / "az_map_support_summary.csv",
        temporal,
        months,
        species,
        PAPER_GROUPS,
        FIG01_OBS_DAY_GROUPS,
        elev,
        args.nmin,
        support_mask,
    )

    print(f"Wrote figures and CSV to {args.out_dir}")
    print(
        f"Species: {len(species)}; tiles: {lon.size}; months: {len(months)}; "
        f"Nmin: {args.nmin}; CYGNSS L1 support tiles: {int(np.sum(support_mask))}"
    )


if __name__ == "__main__":
    main()
