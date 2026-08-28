#!/usr/bin/env python3
"""
Step 1 figure: actual_incr vs pred_incr hexbin density, one panel per
CYGNSS L1 DA experiment, with the 1:1 line and the through-origin best fit.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

PKL_PATH = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/_cygl1_gain_consistency.pkl"
OUT_PNG = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/figures/cygl1_gain_consistency.png"

EXPERIMENTS = ["cygl1assim", "cygl1assim_halfR", "cygl1assim_quarterR"]
TITLES = {
    "cygl1assim": "cygl1assim\n(errstd 4.4 dB)",
    "cygl1assim_halfR": "cygl1assim_halfR\n(errstd 2.2 dB)",
    "cygl1assim_quarterR": "cygl1assim_quarterR\n(errstd 1.1 dB)",
}

BLUE_RAMP = ["#eaf2fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_CMAP = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP, N=256)
INK_PRIMARY = "#1a1a1a"
INK_MUTED = "#6b6b6b"
FIT_COLOR = "#d9730d"  # warm accent, distinct from the sequential blue density fill
ONE_TO_ONE_COLOR = "#9a9a9a"


def through_origin_fit(x, y):
    denom = np.sum(x * x)
    slope = np.sum(x * y) / denom if denom > 0 else np.nan
    resid = y - slope * x
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return slope, r2


def main():
    df = pd.read_pickle(PKL_PATH)

    lim = np.nanpercentile(np.abs(np.concatenate([df["pred_incr"].values, df["actual_incr"].values])), 99.5)
    lim = float(np.ceil(lim * 10) / 10)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharex=True, sharey=True)

    hb_last = None
    for ax, exp in zip(axes, EXPERIMENTS):
        sub = df[df["experiment"] == exp]
        x = sub["pred_incr"].values
        y = sub["actual_incr"].values
        slope, r2 = through_origin_fit(x, y)
        n = len(x)

        hb_last = ax.hexbin(
            x, y, gridsize=70, extent=[-lim, lim, -lim, lim],
            cmap=SEQ_CMAP, bins="log", mincnt=1, linewidths=0,
        )

        line_x = np.array([-lim, lim])
        ax.plot(line_x, line_x, color=ONE_TO_ONE_COLOR, lw=1.5, ls="--", zorder=3)
        ax.plot(line_x, slope * line_x, color=FIT_COLOR, lw=2, zorder=4)

        ax.axhline(0, color=ONE_TO_ONE_COLOR, lw=0.6, zorder=1)
        ax.axvline(0, color=ONE_TO_ONE_COLOR, lw=0.6, zorder=1)

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(TITLES[exp], fontsize=11, color=INK_PRIMARY)
        ax.set_xlabel("pred_incr = K·d", fontsize=10, color=INK_PRIMARY)

        ax.text(
            0.04, 0.96,
            f"slope = {slope:.3f}\nR² = {r2:.3f}\nN = {n:,}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
            color=INK_PRIMARY,
        )

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(INK_MUTED)
        ax.tick_params(colors=INK_MUTED, labelsize=9)

    axes[0].set_ylabel("actual_incr = ana − fcst", fontsize=10, color=INK_PRIMARY)

    legend_elems = [
        Line2D([0], [0], color=ONE_TO_ONE_COLOR, lw=1.5, ls="--", label="1:1 (gain matches exactly)"),
        Line2D([0], [0], color=FIT_COLOR, lw=2, label="through-origin best fit"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=2, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, -0.02))

    cbar = fig.colorbar(hb_last, ax=axes, orientation="vertical", fraction=0.025, pad=0.02)
    cbar.set_label("observation count (log scale)", fontsize=9.5, color=INK_PRIMARY)
    cbar.ax.tick_params(colors=INK_MUTED, labelsize=8.5)

    fig.suptitle(
        "CYGNSS L1 gain-consistency check: does the analysis update apply its own saved K·d?",
        fontsize=12.5, color=INK_PRIMARY, y=1.03,
    )

    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
