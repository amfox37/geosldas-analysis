#!/usr/bin/env python3
"""Compare SFMC and RZMC adjacent-period differences by region and boundary."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
ROOT = HERE.parent
OUT = ROOT / "output" / "regional_rzmc_transitions"
FIGDIR = ROOT / "docs" / "regional_rzmc_diagnostic_figures"
REG = json.loads((ROOT/"config"/"regional_rzmc_regions.json").read_text())["regions"]

def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    order = [r["region_id"] for r in REG]
    lab = {r["region_id"]: r["label"] for r in REG}
    rows = [lab[r] for r in order]
    trs = [f"P{i+1}−P{i}" for i in range(1, 9)]
    S = 1000.0
    d = {v: pd.read_csv(OUT/f"regional_{v}_period_diffs_fdr.csv") for v in ("rzmc", "sfmc")}

    fig, axes = plt.subplots(1, 2, figsize=(14.2, 4.3), sharey=True)
    vmax = 9.0
    for ax, v, title in zip(axes, ("rzmc", "sfmc"),
                            ("(a) RZMC  (0–100 cm)", "(b) SFMC  (0–5 cm)")):
        piv = d[v].pivot(index="region", columns="transition", values="diff").reindex(rows)[trs]*S
        sig = d[v].pivot(index="region", columns="transition", values="sig_fdr").reindex(rows)[trs]
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(trs)):
                val = piv.values[i, j]
                is_sig = bool(sig.values[i, j])
                if is_sig:
                    ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                                               edgecolor="black", lw=2.0, zorder=3))
                ax.text(j, i, f"{val:+.1f}", ha="center", va="center", zorder=4,
                        fontsize=7.6, fontweight="bold" if is_sig else "normal",
                        color="black" if abs(val) < 5 else "white")
        ax.set_xticks(range(len(trs)), trs, fontsize=8.5)
        ax.set_yticks(range(len(rows)), rows, fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.set_xticks(np.arange(-0.5, len(trs), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
        ax.grid(which="minor", color="white", lw=1.2); ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=axes, orientation="vertical", pad=0.015, shrink=0.85, extend="both")
    cb.set_label(r"adjacent-period difference in DA $-$ OL ($\times10^{-3}$ m$^3$ m$^{-3}$)")
    fig.suptitle("Root-zone and surface soil moisture: DA $-$ OL change between successive "
                 "observing-system periods\n(black outline = FDR $q<0.05$ within that transition)",
                 fontsize=12, y=1.06)
    fig.savefig(FIGDIR/"sfmc_rzmc_transition_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig); print("wrote sfmc_rzmc_transition_comparison.png")

if __name__ == "__main__":
    main()
