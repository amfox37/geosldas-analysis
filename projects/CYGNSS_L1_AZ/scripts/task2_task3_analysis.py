#!/usr/bin/env python3
"""
Task 2 (operator transfer function) and Task 3 (observable coherence) analysis,
per the 2026-08-20 brief on CYGNSS L1 vs L3.

Reads the paired extraction produced by extract_cygl1_l3_pairs.py and computes:

Task 2 (forecast-only, unaffected by obs quality):
  - Pooled scatter/hexbin of F(L1) [dB] vs F(L3) [m3/m3] at matched tile-times.
  - A handful of individual-tile scatters.
  - Per-tile Spearman rho between F(L1) and F(L3); distribution summary.
  - Pooled R^2 (F(L1) explained by F(L3), linear fit) and median per-tile R^2.

Task 3 (matched observations, same assimilation window):
  - Count of matched (tile,time) pairs with both O(L1) and O(L3) present, and
    number of distinct tiles.
  - Per-tile temporal correlation (Pearson) between O(L1) and O(L3).
  - On the identical matched samples, per-tile correlation of each species'
    O against its own F (for comparison).

Outputs: printed summary + PNG figures under output/figures/.
"""
import gzip
from io import StringIO

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_PATH = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygnss_l1_l3_paired_ofa_2020.csv.gz"
FIG_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/figures"

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GRID_GRAY = "#c9c8c2"
TEXT_SECONDARY = "#52514e"

import os
os.makedirs(FIG_DIR, exist_ok=True)


def read_csv_skip_header(path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    n_hdr = sum(1 for line in lines if line.startswith("#"))
    return pd.read_csv(StringIO("".join(lines[n_hdr:])))


def style_axes(ax):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID_GRAY)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9)
    ax.xaxis.label.set_color(TEXT_SECONDARY)
    ax.yaxis.label.set_color(TEXT_SECONDARY)


def main():
    df = read_csv_skip_header(CSV_PATH)
    df["datetime"] = pd.to_datetime(df["datetime"])

    l1 = df[df["species_name"] == "CYGNSS_L1_DDM3X5_CROP_SCALAR"]
    l3 = df[df["species_name"] == "CYGNSS_SM_6hr"]

    # ---------------- TASK 2: forecast-only merge ----------------
    l1_f = l1.dropna(subset=["fcst"])[["tile_id", "datetime", "fcst"]].rename(columns={"fcst": "F_L1"})
    l3_f = l3.dropna(subset=["fcst"])[["tile_id", "datetime", "fcst"]].rename(columns={"fcst": "F_L3"})
    merged_f = pd.merge(l1_f, l3_f, on=["tile_id", "datetime"], how="inner")
    print(f"[Task2] matched tile-times with both forecasts: {len(merged_f)}  across {merged_f['tile_id'].nunique()} tiles")

    # Pooled stats
    r_pearson, p_pearson = stats.pearsonr(merged_f["F_L1"], merged_f["F_L3"])
    r_spearman_pooled, _ = stats.spearmanr(merged_f["F_L1"], merged_f["F_L3"])
    slope, intercept, r_lin, p_lin, se = stats.linregress(merged_f["F_L3"], merged_f["F_L1"])
    print(f"[Task2] pooled Pearson r = {r_pearson:.4f}  (R^2 = {r_pearson**2:.4f})")
    print(f"[Task2] pooled Spearman rho = {r_spearman_pooled:.4f}")
    print(f"[Task2] pooled linear fit F_L1 ~ F_L3: slope={slope:.4f} dB per m3/m3, intercept={intercept:.4f}, R^2={r_lin**2:.4f}")

    # Per-tile Spearman
    tile_rhos = []
    for tid, g in merged_f.groupby("tile_id"):
        if len(g) < 10:
            continue
        rho, _ = stats.spearmanr(g["F_L1"], g["F_L3"])
        r2_tile = stats.pearsonr(g["F_L1"], g["F_L3"])[0] ** 2 if len(g) > 2 else np.nan
        tile_rhos.append((tid, len(g), rho, r2_tile))
    tile_rho_df = pd.DataFrame(tile_rhos, columns=["tile_id", "N", "spearman_rho", "pearson_r2"])
    print(f"[Task2] per-tile Spearman rho (N tiles with >=10 matched obs: {len(tile_rho_df)}):")
    print(tile_rho_df["spearman_rho"].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]))
    print(f"[Task2] median per-tile Pearson R^2 (F_L1 vs F_L3): {tile_rho_df['pearson_r2'].median():.4f}")
    print(f"[Task2] fraction of tiles with spearman_rho > 0.5: {(tile_rho_df['spearman_rho'] > 0.5).mean():.3f}")
    print(f"[Task2] fraction of tiles with spearman_rho < 0:   {(tile_rho_df['spearman_rho'] < 0).mean():.3f}")

    tile_rho_df.to_csv(f"{FIG_DIR}/../task2_per_tile_spearman.csv", index=False)

    # ---- Figure 1: pooled hexbin F(L1) vs F(L3) ----
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    hb = ax.hexbin(merged_f["F_L3"], merged_f["F_L1"], gridsize=60, cmap="Blues", mincnt=1, bins="log")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("log10(count)", color=TEXT_SECONDARY, fontsize=9)
    ax.set_xlabel("F(CYGNSS L3) forecast soil moisture [m$^3$ m$^{-3}$]")
    ax.set_ylabel("F(CYGNSS L1) forecast DDM scalar [dB]")
    ax.set_title(f"Operator transfer function: F(L1) vs F(L3)\nOLv8_M36_all_sensors_AZ_scaled, 2020, N={len(merged_f)}, pooled r={r_pearson:.2f}",
                 fontsize=10, color="#0b0b0b")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/task2_fig1_pooled_hexbin_F_L1_vs_F_L3.png")
    plt.close(fig)

    # ---- Figure 2: small multiples, handful of individual tiles ----
    top_tiles = tile_rho_df.sort_values("N", ascending=False)["tile_id"].head(6).tolist()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), dpi=150, sharex=False, sharey=False)
    for ax, tid in zip(axes.ravel(), top_tiles):
        g = merged_f[merged_f["tile_id"] == tid]
        rho = tile_rho_df.loc[tile_rho_df["tile_id"] == tid, "spearman_rho"].values[0]
        ax.scatter(g["F_L3"], g["F_L1"], s=10, color=BLUE, alpha=0.6, edgecolors="none")
        ax.set_title(f"tile {tid}  (N={len(g)}, rho={rho:.2f})", fontsize=9, color="#0b0b0b")
        style_axes(ax)
    fig.supxlabel("F(CYGNSS L3) [m$^3$ m$^{-3}$]", fontsize=10, color=TEXT_SECONDARY)
    fig.supylabel("F(CYGNSS L1) [dB]", fontsize=10, color=TEXT_SECONDARY)
    fig.suptitle("Operator transfer function at individual tiles (highest-N tiles)", fontsize=11, color="#0b0b0b")
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/task2_fig2_individual_tiles_scatter.png")
    plt.close(fig)

    # ---- Figure 3: distribution of per-tile Spearman rho ----
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=150)
    ax.hist(tile_rho_df["spearman_rho"], bins=40, color=BLUE, edgecolor="white", linewidth=0.3)
    med = tile_rho_df["spearman_rho"].median()
    ax.axvline(med, color=ORANGE, linewidth=2, linestyle="--")
    ax.text(med, ax.get_ylim()[1]*0.92, f" median={med:.2f}", color=ORANGE, fontsize=9)
    ax.set_xlabel("Per-tile Spearman rho, F(L1) vs F(L3)")
    ax.set_ylabel("Number of tiles")
    ax.set_title("Task 2: per-tile forecast-forecast correlation distribution", fontsize=10, color="#0b0b0b")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/task2_fig3_spearman_distribution.png")
    plt.close(fig)

    # ---------------- TASK 3: matched observations ----------------
    l1_o = l1.dropna(subset=["obs"])[["tile_id", "datetime", "obs", "fcst"]].rename(columns={"obs": "O_L1", "fcst": "F_L1"})
    l3_o = l3.dropna(subset=["obs"])[["tile_id", "datetime", "obs", "fcst"]].rename(columns={"obs": "O_L3", "fcst": "F_L3"})
    merged_o = pd.merge(l1_o, l3_o, on=["tile_id", "datetime"], how="inner")
    n_pairs = len(merged_o)
    n_tiles = merged_o["tile_id"].nunique()
    print(f"\n[Task3] matched (tile,time) pairs with BOTH O(L1) and O(L3): {n_pairs}  across {n_tiles} tiles")

    task3_rows = []
    for tid, g in merged_o.groupby("tile_id"):
        if len(g) < 10:
            continue
        r_oo, _ = stats.pearsonr(g["O_L1"], g["O_L3"])
        # own-forecast correlation on the SAME identical samples
        g_l1_valid = g.dropna(subset=["O_L1", "F_L1"])
        g_l3_valid = g.dropna(subset=["O_L3", "F_L3"])
        r_l1_of = stats.pearsonr(g_l1_valid["O_L1"], g_l1_valid["F_L1"])[0] if len(g_l1_valid) >= 10 else np.nan
        r_l3_of = stats.pearsonr(g_l3_valid["O_L3"], g_l3_valid["F_L3"])[0] if len(g_l3_valid) >= 10 else np.nan
        task3_rows.append((tid, len(g), r_oo, r_l1_of, r_l3_of))

    task3_df = pd.DataFrame(task3_rows, columns=["tile_id", "N", "r_O_L1_O_L3", "r_O_L1_F_L1", "r_O_L3_F_L3"])
    print(f"[Task3] tiles with >=10 matched obs-obs pairs: {len(task3_df)}")
    print("[Task3] per-tile O(L1) vs O(L3) correlation distribution:")
    print(task3_df["r_O_L1_O_L3"].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]))
    print("[Task3] per-tile O(L1) vs F(L1) correlation distribution (same samples):")
    print(task3_df["r_O_L1_F_L1"].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]))
    print("[Task3] per-tile O(L3) vs F(L3) correlation distribution (same samples):")
    print(task3_df["r_O_L3_F_L3"].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]))

    task3_df.to_csv(f"{FIG_DIR}/../task3_per_tile_correlations.csv", index=False)

    # ---- Figure 4: Task3 comparison boxplots ----
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=150)
    data = [
        task3_df["r_O_L1_O_L3"].dropna().values,
        task3_df["r_O_L1_F_L1"].dropna().values,
        task3_df["r_O_L3_F_L3"].dropna().values,
    ]
    labels = ["O(L1) vs O(L3)\n(cross-product)", "O(L1) vs F(L1)\n(own fcst)", "O(L3) vs F(L3)\n(own fcst)"]
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False, widths=0.5)
    colors = [AQUA, BLUE, ORANGE]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
        patch.set_edgecolor(c)
    for med_line in bp["medians"]:
        med_line.set_color("#0b0b0b")
    ax.axhline(0, color=GRID_GRAY, linewidth=1)
    ax.set_ylabel("Per-tile Pearson correlation")
    ax.set_title(f"Task 3: does the observable itself carry the signal?\nN={n_pairs} matched obs pairs, {len(task3_df)} tiles (N>=10)", fontsize=10, color="#0b0b0b")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/task3_fig1_correlation_comparison.png")
    plt.close(fig)

    return dict(
        n_forecast_pairs=len(merged_f),
        n_forecast_tiles=merged_f["tile_id"].nunique(),
        pooled_pearson_r=r_pearson,
        pooled_r2=r_pearson**2,
        pooled_spearman=r_spearman_pooled,
        median_tile_spearman=tile_rho_df["spearman_rho"].median(),
        median_tile_r2=tile_rho_df["pearson_r2"].median(),
        n_obs_pairs=n_pairs,
        n_obs_tiles=n_tiles,
        median_r_oo=task3_df["r_O_L1_O_L3"].median(),
        median_r_l1_of=task3_df["r_O_L1_F_L1"].median(),
        median_r_l3_of=task3_df["r_O_L3_F_L3"].median(),
    )


if __name__ == "__main__":
    result = main()
    print("\n=== SUMMARY ===")
    for k, v in result.items():
        print(f"{k}: {v}")
