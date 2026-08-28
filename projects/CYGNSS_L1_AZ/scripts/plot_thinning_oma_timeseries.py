#!/usr/bin/env python3
"""Plot monthly OmF and OmA spread for a CYGNSS L1 density experiment."""
from __future__ import annotations

import csv
import argparse
import calendar
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np


PROJECT = Path(__file__).resolve().parent.parent
STATS = PROJECT / "output/thinning_expts"
OUT = STATS / "figures"
OL_MONTHLY = STATS / "spatial_stats_OL_paired_monitor_202001_202010.pkl"
OL_TEMPORAL = STATS / "temporal_stats_OL_paired_monitor_20200101_20201031.nc4"
NMIN = 10


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
        help="CYGNSS L1 density experiment to summarize",
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


def load_pickle(path: Path) -> dict[str, np.ndarray]:
    with path.open("rb") as stream:
        return pickle.load(stream)


def load_temporal(path: Path) -> dict[str, np.ndarray]:
    with Dataset(path) as nc:
        return {
            key: np.ma.filled(nc[key][:], np.nan).astype(float)
            for key in ("OmF_stdv", "OmA_stdv", "N_data")
        }


def monthly_group_series(data: dict[str, np.ndarray], key: str, group: Group) -> np.ndarray:
    values = np.asarray(data[key], dtype=float)[:, list(group.indices)]
    counts = np.asarray(data["N_data"], dtype=float)[:, list(group.indices)]
    valid = np.isfinite(values) & np.isfinite(counts) & (counts > 0)
    denominator = np.where(valid, counts, 0.0).sum(axis=1)
    return np.divide(
        np.where(valid, values * counts, 0.0).sum(axis=1),
        denominator,
        out=np.full(denominator.shape, np.nan),
        where=denominator > 0,
    )


def aggregate_group_value(
    data: dict[str, np.ndarray], key: str, group: Group, common: np.ndarray
) -> float:
    values = data[key][:, list(group.indices)]
    counts = data["N_data"][:, list(group.indices)]
    return float(
        np.sum(np.where(common, values * counts, 0.0))
        / np.sum(np.where(common, counts, 0.0))
    )


def main() -> None:
    args = parse_args()
    experiment = args.experiment
    label = experiment.capitalize()
    file_period, month_period, period_label = period_settings(args.end_month)
    if args.ol_support == "xmask" and experiment != "intermediate":
        raise ValueError("xmask OL support is available only for the intermediate experiment")
    ol_monthly_path = (
        STATS / f"spatial_stats_OL_paired_monitor_xmask_intermediate_{month_period}.pkl"
        if args.ol_support == "xmask"
        else STATS / f"spatial_stats_OL_paired_monitor_{month_period}.pkl"
    )
    ol_temporal_path = (
        STATS / f"temporal_stats_OL_paired_monitor_xmask_intermediate_{file_period}.nc4"
        if args.ol_support == "xmask"
        else STATS / f"temporal_stats_OL_paired_monitor_{file_period}.nc4"
    )
    ol_label = "OL O-F (xmask)" if args.ol_support == "xmask" else "OL O-F"
    output_tag = experiment + ("_xmask" if args.ol_support == "xmask" else "")
    if args.end_month != "202010":
        output_tag += f"_{month_period}"
    da_monthly_path = STATS / f"spatial_stats_DA_paired_{experiment}_{month_period}.pkl"
    da_temporal_path = STATS / f"temporal_stats_DA_paired_{experiment}_{file_period}.nc4"
    OUT.mkdir(parents=True, exist_ok=True)
    ol_monthly = load_pickle(ol_monthly_path)
    da_monthly = load_pickle(da_monthly_path)
    ol_temporal = load_temporal(ol_temporal_path)
    da_temporal = load_temporal(da_temporal_path)

    months = [str(value) for value in ol_monthly["date_vec"]]
    x = np.arange(len(months))
    tick_step = max(1, int(np.ceil(len(months) / 8)))
    tick_positions = np.arange(0, len(months), tick_step)
    if tick_positions[-1] != len(months) - 1:
        tick_positions = np.append(tick_positions, len(months) - 1)
    tick_labels = [f"{months[index][:4]}-{months[index][4:]}" for index in tick_positions]
    summary_rows = []

    fig, axes = plt.subplots(
        len(GROUPS),
        2,
        figsize=(13.5, 16.5),
        gridspec_kw={"width_ratios": (1.25, 1.0), "hspace": 0.72, "wspace": 0.22},
    )

    for row_index, group in enumerate(GROUPS):
        ol_f = monthly_group_series(ol_monthly, "OmF_stdv", group)
        da_f = monthly_group_series(da_monthly, "OmF_stdv", group)
        da_a = monthly_group_series(da_monthly, "OmA_stdv", group)
        background_reduction = 100.0 * (ol_f - da_f) / ol_f
        analysis_reduction = 100.0 * (da_f - da_a) / da_f

        ax_values, ax_percent = axes[row_index]
        ax_values.plot(x, ol_f, "--o", color="0.35", lw=1.7, ms=3.5, label=ol_label)
        ax_values.plot(x, da_f, "-o", color="#2b6cb0", lw=1.9, ms=3.7, label=f"{label} O-F")
        ax_values.plot(x, da_a, "-s", color="#e07b39", lw=1.9, ms=3.5, label=f"{label} O-A")
        ax_values.set_ylabel(f"standard deviation ({group.units})")
        ax_values.set_title(f"({chr(97 + row_index)}) {group.name}", loc="left", fontweight="bold")

        ax_percent.plot(
            x,
            background_reduction,
            "-o",
            color="#2b6cb0",
            lw=1.8,
            ms=3.6,
            label="O-F reduction vs OL",
        )
        ax_percent.plot(
            x,
            analysis_reduction,
            "-s",
            color="#e07b39",
            lw=1.8,
            ms=3.4,
            label="O-F to O-A reduction",
        )
        ax_percent.axhline(0.0, color="0.3", lw=0.9)
        ax_percent.set_ylabel("reduction (%)")

        for axis in (ax_values, ax_percent):
            axis.set_xticks(tick_positions)
            axis.set_xticklabels(tick_labels, rotation=35, ha="right")
            axis.grid(color="0.90", lw=0.7)
            axis.set_axisbelow(True)
        if row_index == 0:
            ax_values.legend(frameon=False, ncols=3, fontsize=8.5)
            ax_percent.legend(frameon=False, fontsize=8.5)
        if row_index == len(GROUPS) - 1:
            ax_values.set_xlabel("month")
            ax_percent.set_xlabel("month")

        idx = list(group.indices)
        common = (
            np.isfinite(ol_temporal["OmF_stdv"][:, idx])
            & np.isfinite(da_temporal["OmF_stdv"][:, idx])
            & np.isfinite(da_temporal["OmA_stdv"][:, idx])
            & (ol_temporal["N_data"][:, idx] >= NMIN)
            & (da_temporal["N_data"][:, idx] >= NMIN)
        )
        ol_f_all = aggregate_group_value(ol_temporal, "OmF_stdv", group, common)
        da_f_all = aggregate_group_value(da_temporal, "OmF_stdv", group, common)
        da_a_all = aggregate_group_value(da_temporal, "OmA_stdv", group, common)
        summary_rows.append(
            {
                "system": group.name,
                "common_species_tile_cells": int(common.sum()),
                "ol_omf_stdv": ol_f_all,
                "experiment_omf_stdv": da_f_all,
                "experiment_oma_stdv": da_a_all,
                "omf_reduction_vs_ol_percent": 100.0 * (ol_f_all - da_f_all) / ol_f_all,
                "oma_reduction_from_omf_percent": 100.0 * (da_f_all - da_a_all) / da_f_all,
                "months_oma_lower_than_omf": int(np.sum(analysis_reduction > 0)),
            }
        )

    fig.suptitle(
        f"{label} CYGNSS L1 experiment: monthly O-F and O-A spread, {period_label}\n"
        "positive percentages indicate reduced residual variability",
        fontsize=14,
        y=0.995,
    )
    fig.text(
        0.5,
        0.005,
        (
            "CYGNSS L1 OL uses the intermediate retained-event xmask (0.24% fewer events); "
            "intermediate O-F and O-A use the same retained events."
            if args.ol_support == "xmask"
            else f"CYGNSS L1 OL and {experiment} O-F use different event samples when thinning changes "
            f"support; {experiment} O-F and O-A use the same retained events."
        ),
        ha="center",
        fontsize=8.5,
        color="0.35",
    )
    fig.subplots_adjust(top=0.94, bottom=0.06)
    figure_path = OUT / f"{output_tag}_omf_oma_monthly_timeseries.png"
    fig.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    summary_path = OUT / f"{output_tag}_omf_oma_summary.csv"
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"wrote {figure_path.relative_to(PROJECT)}")
    print(f"wrote {summary_path.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
