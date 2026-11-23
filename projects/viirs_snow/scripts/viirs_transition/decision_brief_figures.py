"""Figures and numbers for docs/viirs_decision_brief.md (python -I; set MPLCONFIGDIR to a writable dir under land_da).

Plain-language version of analyses B, C, D for the operational decision: MODIS vs VIIRS, and Terra+Aqua (two sensors) vs one sensor.
Writes into docs/viirs_transition_report_figures/:
  brief_1_water_added.png     net snow water added by assimilation (kg m-2 per snow year), split into ET / runoff / storage, Aug 2023 - Jul 2026
  brief_2_update_structure.png  why: updates per day and total adjustment per tile
  brief_3_skill.png           skill benefit of VIIRS over MODIS, and of a second satellite (benefit > 0 = better)
and prints the numbers quoted in the brief (python stdout) so the text can be checked against them.
Inputs: D_budget_tables.csv, B_summary.md, C_skill_contrasts.csv, IMS_skill_contrasts_fixed.csv (flip-corrected IMS).
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc

FIG = sc.FIG
YEARS = {"SY4": "2023-24", "SY5": "2024-25", "SY6": "2025-26"}
CFG = [("MODIS\nTerra only", "T_control", "#8c6bb1"), ("MODIS\nTerra+Aqua", "TA_control", "#3f007d"),
       ("VIIRS\nJ2 only", "T_to_J2", "#fc8d59"), ("VIIRS\nJ1+J2", "TA_to_J1J2", "#b30000")]

# ---------------- Figure 1: water added ----------------
d = pd.read_csv(FIG / "D_budget_tables.csv")
d = d[d.kind == "C0"]
def val(arm, var, sy):
    return float(d[(d.contrast == f"{arm}-OL") & (d["var"] == var) & (d.sy == sy)]["mean"].iloc[0])
W = {}
for lab, arm, _ in CFG:
    W[arm] = {v: np.mean([val(arm, v, s) for s in YEARS]) for v in ["I_snow", "dET", "dRunoff", "dStorage", "RZMC_mean", "SWE_max"]}
    W[arm]["years"] = [val(arm, "I_snow", s) for s in YEARS]
ref = W["TA_control"]["I_snow"]
fig, ax = plt.subplots(figsize=(8.6, 5.2))
for i, (lab, arm, col) in enumerate(CFG):
    b = 0
    for v, hatch, a in [("dET", "", 1.0), ("dRunoff", "//", 0.65), ("dStorage", "..", 0.4)]:
        x = max(W[arm][v], 0)
        ax.bar(i, W[arm][v], bottom=b, color=col, alpha=a, hatch=hatch, edgecolor="white", width=0.62)
        b += W[arm][v]
    ax.plot([i] * 3, W[arm]["years"], "k.", ms=7, zorder=5)
    ax.text(i, max(W[arm]["years"]) + 4, f"{W[arm]['I_snow']:.0f}\n({100 * W[arm]['I_snow'] / ref:.0f}% of MODIS Terra+Aqua)", ha="center", fontsize=8.5)
ax.set_xticks(range(len(CFG))); ax.set_xticklabels([c[0] for c in CFG])
ax.set_ylabel("Net snow water added by assimilation\n(kg m-2 = mm of water, per snow year)")
ax.set_ylim(0, 170)
h = [plt.Rectangle((0, 0), 1, 1, fc="0.3", alpha=1.0), plt.Rectangle((0, 0), 1, 1, fc="0.3", alpha=0.65, hatch="//", ec="white"),
     plt.Rectangle((0, 0), 1, 1, fc="0.3", alpha=0.4, hatch="..", ec="white"), plt.Line2D([], [], marker=".", ls="", color="k")]
ax.legend(h, ["goes to evapotranspiration", "goes to runoff", "stays in storage", "one dot per snow year"], fontsize=8.5, loc="upper right", frameon=False)
ax.set_title("How much extra water does snow assimilation add?\n(mean of Aug 2023 - Jul 2026; difference from a run without snow assimilation)", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "brief_1_water_added.png", dpi=150); plt.close(fig)

# ---------------- Figure 2: update structure (from B_summary.md) ----------------
txt = (FIG / "B_summary.md").read_text()
def parse(section, ncols):
    s = txt.split(section, 1)[1].split("\n## ", 1)[0]
    rows = [[c.strip() for c in l.strip().strip("|").split("|")] for l in s.splitlines() if l.startswith("| ") and not l.startswith("| arm") and not l.startswith("|---")]
    return rows
b1 = parse("## B1.", 7); b3 = parse("## B3.", 7)
two = {a: np.mean([float(r[5]) for r in b1 if r[0] == a and r[1] in YEARS]) for _, a, _ in CFG}
absinc = {a: np.mean([float(r[3]) for r in b3 if r[0] == a and r[1] in YEARS]) for _, a, _ in CFG}
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
for i, (lab, arm, col) in enumerate(CFG):
    axs[0].bar(i, two[arm], color=col, width=0.62); axs[0].text(i, two[arm] + 1.2, f"{two[arm]:.0f}%", ha="center", fontsize=9)
    axs[1].bar(i, absinc[arm], color=col, width=0.62); axs[1].text(i, absinc[arm] + 2, f"{absinc[arm]:.0f}", ha="center", fontsize=9)
for ax in axs:
    ax.set_xticks(range(len(CFG))); ax.set_xticklabels([c[0] for c in CFG], fontsize=9)
axs[0].set_ylabel("% of tile-days with an update\nthat get two updates in the same day"); axs[0].set_ylim(0, 70)
axs[1].set_ylabel("Total adjustment to snow water per tile\n(kg m-2 per year, sum of absolute increments)"); axs[1].set_ylim(0, 165)
axs[0].set_title("Two updates in one day", fontsize=10); axs[1].set_title("Total size of the adjustments", fontsize=10)
fig.suptitle("Why: only MODIS Terra+Aqua updates twice a day (north of 30N, mean Aug 2023 - Jul 2026)", fontsize=10.5)
fig.tight_layout(); fig.savefig(FIG / "brief_2_update_structure.png", dpi=150); plt.close(fig)

# ---------------- Figure 3: skill benefit ----------------
ct = pd.read_csv(FIG / "C_skill_contrasts.csv")
fx = pd.read_csv(FIG / "IMS_skill_contrasts_fixed.csv").assign(dataset="IMS")
ct = pd.concat([ct[ct.dataset != "IMS"], fx[[c for c in ct.columns if c in fx.columns]]], ignore_index=True)
ct = ct[(ct.season == "ALL") | ct.season.isna()]
PAN = [("IMS", "accuracy", "IMS snow-cover accuracy\n(percentage points)", -100.0, "pp"),
       ("snotel_SWE", "rmse", "SNOTEL snow water equivalent error\n(RMSE, mm)", 1.0, "mm"),
       ("ghcn_SNWD", "rmse", "GHCN snow depth error\n(RMSE, mm)", 1.0, "mm")]
LBL = {"SY1": "2020-21", "SY2": "2021-22", "SY3": "2022-23", "SY4": "2023-24", "SY5": "2024-25", "SY6": "2025-26"}
SYC = {"SY1": "#9ecae1", "SY2": "#6baed6", "SY3": "#4292c6", "SY4": "#fdae6b", "SY5": "#fd8d3c", "SY6": "#e6550d"}
SW = [("Terra only:\nMODIS vs VIIRS J2", "C1 T_control-T_to_J2", "SY4-6", ["SY4", "SY5", "SY6"]),
      ("Terra+Aqua:\nMODIS vs VIIRS J1+J2", "C2 TA_control-TA_to_J1J2", "SY4-6", ["SY4", "SY5", "SY6"])]
CT = [("MODIS: Terra+Aqua\nvs Terra only (6 years)", "C3 MODIS_T-MODIS_TA", "SY1-6", [f"SY{i}" for i in range(1, 7)]),
      ("VIIRS: J1+J2\nvs J2 only (3 years)", "C4a T_to_J2-TA_to_J1J2 (post)", "SY4-6", ["SY4", "SY5", "SY6"])]
def row(ds, m, name, win):
    r = ct[(ct.dataset == ds) & (ct.metric == m) & (ct.contrast == name) & (ct.window == win)]
    return None if r.empty else r.iloc[0]
fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.8))
NUM = []
for c, (ds, m, lab, sg, unit) in enumerate(PAN):
    ref0 = row(ds, m, "C0 TA_to_J1J2-OL", "SY1-6")
    gain = (-ref0.delta if sg > 0 else ref0.delta * 100.0)
    for r_i, series in enumerate([SW, CT]):
        ax = axes[r_i, c]
        for k, (sl, name, pw, yrs) in enumerate(series):
            y = len(series) - 1 - k
            for yr in yrs:
                rr = row(ds, m, name, yr)
                if rr is not None:
                    ax.plot(rr.delta * sg, y + (int(yr[2:]) - 3.5) * 0.04, "o", color=SYC[yr], ms=6, zorder=3)
            p = row(ds, m, name, pw)
            lo, hi = sorted([p.ci_lo * sg, p.ci_hi * sg])
            ax.errorbar(p.delta * sg, y, xerr=[[p.delta * sg - lo], [hi - p.delta * sg]], fmt="s", color="k", ms=7, capsize=4, zorder=4)
            NUM.append((ds, m, name, p.delta * sg, lo, hi, gain))
        ax.axvline(0, color="0.4", lw=1); ax.grid(axis="x", color="0.9")
        ax.set_yticks(range(len(series))); ax.set_yticklabels([s[0] for s in series][::-1] if c == 0 else [], fontsize=8.5)
        ax.set_ylim(-0.6, len(series) - 0.4)
        ax.set_xlabel(f"benefit ({unit}); right = better", fontsize=8.5)
        if r_i == 0:
            ax.set_title(f"{lab}\nfor scale: assimilating snow cover gains {gain:.2f} {unit}", fontsize=9)
fig.tight_layout(rect=[0, 0.08, 1, 0.93], h_pad=4.0)
for r_i, ttl in enumerate(["Does switching from MODIS to VIIRS change skill?   Right of zero = VIIRS better (Aug 2023 - Jul 2026)",
                           "Does a second satellite improve skill?   Right of zero = two satellites better"]):
    y1 = axes[r_i, 0].get_position().y1
    fig.text(0.01, y1 + (0.085 if r_i == 0 else 0.03), ttl, fontsize=11, fontweight="bold", ha="left")
hd = [plt.Line2D([], [], marker="o", ls="", color=SYC[k], label=LBL[k]) for k in SYC] + [plt.Line2D([], [], marker="s", ls="", color="k", label="all years pooled, 95% range")]
fig.legend(handles=hd, loc="lower center", ncol=7, fontsize=8.5, frameon=False, title="dots = single snow years (Aug-Jul)", title_fontsize=8.5)
fig.savefig(FIG / "brief_3_skill.png", dpi=150); plt.close(fig)

# ---------------- numbers for the brief ----------------
print("WATER (mean Aug 2023-Jul 2026, kg m-2 per snow year; % of MODIS T+A)")
for lab, arm, _ in CFG:
    w = W[arm]
    print(f"  {arm:11s} I_snow {w['I_snow']:.1f} ({100*w['I_snow']/ref:.0f}%)  ET {w['dET']:.1f} runoff {w['dRunoff']:.1f} storage {w['dStorage']:.1f}  RZMC {w['RZMC_mean']:.4f}  peakSWE {w['SWE_max']:.2f}  years {[round(x,1) for x in w['years']]}")
print("UPDATES % two-per-day:", {k: round(v, 1) for k, v in two.items()}, " abs increment kg m-2/yr:", {k: round(v, 1) for k, v in absinc.items()})
print("SKILL benefit (pooled, CI lo, hi, DA-OL gain)")
for n in NUM: print("  ", n[0], n[1], n[2], [round(x, 3) for x in n[3:]])
