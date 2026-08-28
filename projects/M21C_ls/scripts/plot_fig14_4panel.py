#!/usr/bin/env python3
"""Figure 14 (4-panel): absolute OL/DA annual water balance, the DA-OL difference
budget, and the partition of the positive snow-DA input.

Seasonal-snow domain, WY2001-WY2006. Reads only production output from
water_year_snow_da_budget.py; run that first. Storage endpoints come from
catch_internal_rst restarts and peat free-water change from PEATCLSM_FSWCHANGE,
both wired into that pipeline.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
BUDGET = ROOT / "output" / "monthly_synthesis_diagnostics" / "water_year_snow_da_budget"
YEAR_ROWS = lambda f: f[f["water_year"].astype(str).str.fullmatch(r"\d+")]

# (a), (b): per-run absolute balances
_abs = YEAR_ROWS(pd.read_csv(BUDGET / "annual_absolute_budgets.csv"))
W = {
    run: _abs[_abs["run"] == run].sort_values("water_year")[
        ["precipitation", "I_snow", "ET", "runoff_surface", "baseflow",
         "storage", "peat_free_standing_water", "residual"]
    ].to_numpy()
    for run in ("OL", "DA")
}
# (c): DA - OL differential budget
Cd = YEAR_ROWS(pd.read_csv(BUDGET / "annual_domain_budgets.csv")).sort_values("water_year")[
    ["I_snow", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "dPeatFreeStandingWater", "residual"]
].to_numpy()
# (d): snow-addition partition with 5-degree spatial-block intervals
_part = pd.read_csv(BUDGET / "six_year_integrated_partitions.csv").set_index("sample").loc["addition"]
_unc = pd.read_csv(BUDGET / "partition_spatial_block_uncertainty.csv")
_unc = _unc[(_unc["sample"] == "addition") & (_unc["block_degrees"] == 5.0)].set_index("variable")
B = np.array([
    [_part[f"fraction_{term}"], _unc.loc[term, "ci_low"], _unc.loc[term, "ci_high"]]
    for term in ["dRunoff_surface", "dBaseflow", "dET", "dStorage",
                 "dPeatFreeStandingWater", "residual", "dRunoff_total"]
])
POSITIVE_INPUT = float(_part["I_snow"])
YEARS = [f"WY{y}" for y in range(2001, 2007)]

C = dict(et="#e0a02c", rs="#2a9d8f", rb="#83c5be", ds="#4a72b0",
         fsw="#8d6bb1", res="#8a8a8a", P="#c0392b", I="#2b6cb0")
TERMS = (("et", "Evapotranspiration", None), ("rs", "Surface runoff", None),
         ("rb", "Baseflow", None), ("ds", "Storage change", None),
         ("fsw", "Peatland free-standing water", None), ("res", "Residual", "//"))

fig = plt.figure(figsize=(14.0, 10.4))
gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.13,
                      height_ratios=[1.0, 0.92], left=0.06, right=0.985,
                      top=0.925, bottom=0.115)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1], sharey=axA)
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
x = np.arange(6)


def stack(ax, series, width=0.62):
    pos = np.zeros(6); neg = np.zeros(6)
    for vals, (key, lab, h) in zip(series, TERMS):
        up = np.where(vals > 0, vals, 0.0); dn = np.where(vals < 0, vals, 0.0)
        ax.bar(x, up, width, bottom=pos, color=C[key], edgecolor="0.25",
               lw=0.5, hatch=h, zorder=2)
        ax.bar(x, dn, width, bottom=neg, color=C[key], edgecolor="0.25",
               lw=0.5, hatch=h, zorder=2)
        pos += up; neg += dn
    return pos, neg


# ---- (a), (b) absolute balances ----
for ax, run, ttl in ((axA, "OL", "(a) Open loop"), (axB, "DA", "(b) Data assimilation")):
    P, I, ET, RS, RB, DS, FSW, RES = W[run].T
    stack(ax, (ET, RS, RB, DS, FSW, RES))
    ax.plot(x, P, "-o", color=C["P"], mfc="white", mec=C["P"], ms=6.5, lw=1.7, zorder=5)
    if run == "DA":
        ax.plot(x, P + I, "-^", color=C["I"], mfc="white", mec=C["I"], ms=6.5, lw=1.7, zorder=5)
    ax.axhline(0, color="0.3", lw=0.8, zorder=3)
    ax.set_xticks(x, YEARS, fontsize=9.5)
    ax.grid(axis="y", color="0.9", lw=0.7); ax.set_axisbelow(True)
    ax.set_title(ttl, loc="left", fontweight="bold", fontsize=12)
    mres, min_ = RES.mean(), (P + I).mean()
    line2 = f"Mean residual {mres:+.2f} ({100*mres/min_:+.2f}%)".replace("-", "\u2212")
    ax.text(0.985, 0.975,
            f"Mean input {min_:.0f} kg m$^{{-2}}$ yr$^{{-1}}$\n" + line2,
            transform=ax.transAxes, ha="right", va="top", fontsize=9.2)
axA.set_ylabel(r"Area-weighted water flux (kg m$^{-2}$ yr$^{-1}$)")
axA.set_ylim(-30, 830)
plt.setp(axB.get_yticklabels(), visible=False)

# ---- (c) DA - OL difference ----
Isnow = Cd[:, 0]
pos, neg = stack(axC, (Cd[:, 3], Cd[:, 1], Cd[:, 2], Cd[:, 4], Cd[:, 5], Cd[:, 6]))
axC.plot(x, Isnow, "-D", color=C["P"], mfc="white", mec=C["P"], ms=6.5, lw=1.7, zorder=5)
axC.axhline(0, color="0.3", lw=0.8, zorder=3)
axC.set_xticks(x, YEARS, fontsize=9.5)
axC.grid(axis="y", color="0.9", lw=0.7); axC.set_axisbelow(True)
axC.set_ylabel(r"DA $-$ OL (kg m$^{-2}$ yr$^{-1}$)")
axC.set_ylim(min(neg.min() * 1.7, -7), max(pos.max(), Isnow.max()) * 1.30)
axC.set_title("(c) DA $-$ OL difference budget", loc="left", fontweight="bold", fontsize=12)
mres = Cd[:, 6].mean()
axC.text(0.985, 0.975,
         f"Mean net snow-DA input {Isnow.mean():.1f} kg m$^{{-2}}$ yr$^{{-1}}$\n"
         + f"Mean residual {mres:+.2f} ({100*mres/Isnow.mean():+.1f}% of net input)"
           .replace("-", "\u2212"),
         transform=axC.transAxes, ha="right", va="top", fontsize=9.2)

# ---- (d) partition ----
labs = ["Surface\nrunoff", "Baseflow", "ET", "Storage", "Free-standing\nwater", "Residual"]
cols = [C["rs"], C["rb"], C["et"], C["ds"], C["fsw"], C["res"]]
v, lo, hi = B[:6, 0] * 100, B[:6, 1] * 100, B[:6, 2] * 100
xb = np.arange(6)
bars = axD.bar(xb, v, 0.62, color=cols, edgecolor="0.25", lw=0.5, zorder=2)
bars[-1].set_hatch("//")
axD.errorbar(xb, v, yerr=[v - lo, hi - v], fmt="none", ecolor="0.15",
             capsize=4, lw=1.3, zorder=3)
axD.axhline(0, color="0.3", lw=0.8, zorder=3)
axD.set_xticks(xb, labs, fontsize=9)
axD.set_ylabel("Fraction of input (%)")
axD.grid(axis="y", color="0.9", lw=0.7); axD.set_axisbelow(True)
rt, rtl, rth = B[6] * 100
top = max(hi[0], hi[1]) + 3.0
axD.plot([0, 0, 1, 1], [top - 1.2, top, top, top - 1.2], color="0.2", lw=1.1, zorder=4)
axD.text(0.5, top + 0.8, f"Total runoff\n{rt:.1f}% [{rtl:.1f}, {rth:.1f}]",
         ha="center", va="bottom", fontsize=9.8, fontweight="bold", linespacing=1.25)
axD.text(0.985, 0.985, f"Mean positive input: {POSITIVE_INPUT:.1f} kg m$^{{-2}}$ yr$^{{-1}}$",
         transform=axD.transAxes, ha="right", va="top", fontsize=9.2)
axD.set_title("(d) Partition of the added water", loc="left",
              fontweight="bold", fontsize=12)
axD.set_ylim(min(lo.min() - 3, -7), top + 11)

# ---- shared legend ----
handles = [Patch(facecolor=C[k], edgecolor="0.25", lw=0.5, hatch=h, label=l)
           for k, l, h in TERMS]
handles += [Line2D([], [], color=C["P"], marker="o", mfc="white", ms=6.5, lw=1.7,
                   label="Precipitation (a, b)"),
            Line2D([], [], color=C["I"], marker="^", mfc="white", ms=6.5, lw=1.7,
                   label="Precipitation + snow-DA input (b)"),
            Line2D([], [], color=C["P"], marker="D", mfc="white", ms=6.5, lw=1.7,
                   label="Net snow-DA input (c)")]
fig.legend(handles=handles, loc="lower center", ncols=5, frameon=False,
           fontsize=9.6, bbox_to_anchor=(0.5, 0.005), columnspacing=1.6)
fig.suptitle("Land water budget over the Northern Hemisphere seasonal-snow domain, "
             "WY2001–WY2006", fontsize=13, y=0.972)

stem = "fig14_snow_da_water_budget_4panel"
for directory, suffix in (("docs", ".png"), ("output", ".png"), ("output", ".pdf")):
    out = ROOT / directory / "paper_figures" / f"{stem}{suffix}"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print("wrote", out.relative_to(ROOT))
