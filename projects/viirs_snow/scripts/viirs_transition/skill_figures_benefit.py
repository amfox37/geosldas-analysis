"""Reader-friendly skill figure for analysis C (python -I; set MPLCONFIGDIR to a writable dir under land_da).

Two questions, one sign convention. Every dot/square is the BENEFIT of the newer / larger configuration, so right of zero always
means "better":
  row 1  Does the MODIS -> VIIRS switch change skill?   benefit = VIIRS arm minus its never-switched MODIS control (SY4-6)
  row 2  Does a second satellite help?                  benefit = (Terra+Aqua or J1+J2) minus (Terra-only or J2-only)
Benefit = RMSE(first) - RMSE(second) for RMSE metrics (mm), and -(accuracy delta) x 100 for IMS accuracy (percentage points),
because the contrasts in C_skill_contrasts.csv are 'first minus second' and positive accuracy delta = first term better.
The DA - OL gain of TA_to_J1J2 over SY1-6 is quoted in each column title for scale (not drawn, it would squash the axes).
Reads C_skill_contrasts.csv and IMS_skill_contrasts_fixed.csv (flip-corrected IMS); writes figC_5_benefit_switch_and_constellation.png.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc

ct = pd.read_csv(sc.FIG / "C_skill_contrasts.csv")
fx = pd.read_csv(sc.FIG / "IMS_skill_contrasts_fixed.csv").assign(dataset="IMS")
ct = pd.concat([ct[ct.dataset != "IMS"], fx[[c for c in ct.columns if c in fx.columns]]], ignore_index=True)
ct = ct[(ct.season == "ALL") | ct.season.isna()]  # station data sets carry no season label (all non-JJA days)

PANELS = [("IMS", "accuracy", "IMS accuracy\n(percentage points; IMS not independent of imagery)", -100.0),
          ("snotel_SWE", "rmse", "SNOTEL SWE RMSE\n(mm; western US/Alaska)", 1.0),
          ("snotel_SNWD", "rmse", "SNOTEL snow depth RMSE\n(mm; western US/Alaska)", 1.0),
          ("ghcn_SNWD", "rmse", "GHCN snow depth RMSE\n(mm; sparse stations)", 1.0)]
# (label, contrast name, pooled window, per-year windows)
SWITCH = [("Terra: MODIS control vs T_to_J2\n(C1, SY4-6)", "C1 T_control-T_to_J2", "SY4-6", ["SY4", "SY5", "SY6"]),
          ("Terra+Aqua: MODIS control vs\nTA_to_J1J2 (C2, SY4-6)", "C2 TA_control-TA_to_J1J2", "SY4-6", ["SY4", "SY5", "SY6"])]
COUNT = [("MODIS: Terra+Aqua vs Terra-only\n(C3, SY1-6)", "C3 MODIS_T-MODIS_TA", "SY1-6", [f"SY{i}" for i in range(1, 7)]),
         ("VIIRS: J1+J2 vs J2-only\n(C4a, SY4-6)", "C4a T_to_J2-TA_to_J1J2 (post)", "SY4-6", ["SY4", "SY5", "SY6"])]
SYC = {"SY1": "#9ecae1", "SY2": "#6baed6", "SY3": "#4292c6", "SY4": "#fdae6b", "SY5": "#fd8d3c", "SY6": "#e6550d"}

def row(ds, m, name, win):
    r = ct[(ct.dataset == ds) & (ct.metric == m) & (ct.contrast == name) & (ct.window == win)]
    return None if r.empty else r.iloc[0]

fig, axes = plt.subplots(2, 4, figsize=(17, 7.4))
for c, (ds, m, lab, sg) in enumerate(PANELS):
    ref = row(ds, m, "C0 TA_to_J1J2-OL", "SY1-6")
    # DA - OL gain: RMSE delta is negative when DA is better, IMS accuracy delta positive -> sign handled below
    gain = (-ref.delta if sg > 0 else ref.delta * 100.0)
    for r_i, (series, ttl) in enumerate([(SWITCH, "Does the MODIS -> VIIRS switch change skill?\nright of zero = VIIRS better"),
                                         (COUNT, "Does a second satellite help?\nright of zero = more satellites better")]):
        ax = axes[r_i, c]
        for k, (sl, name, pw, years) in enumerate(series):
            y = len(series) - 1 - k
            for yr in years:
                rr = row(ds, m, name, yr)
                if rr is not None:
                    ax.plot(rr.delta * sg, y + (int(yr[2:]) - 3.5) * 0.04, "o", color=SYC[yr], ms=6, zorder=3)
            p = row(ds, m, name, pw)
            lo, hi = sorted([p.ci_lo * sg, p.ci_hi * sg])
            ax.errorbar(p.delta * sg, y, xerr=[[p.delta * sg - lo], [hi - p.delta * sg]], fmt="s", color="k", ms=7, capsize=4, zorder=4)
        ax.axvline(0, color="0.4", lw=1)
        ax.set_yticks(range(len(series)))
        ax.set_yticklabels([s[0] for s in series][::-1] if c == 0 else [], fontsize=8.5)
        ax.set_ylim(-0.6, len(series) - 0.4)
        ax.grid(axis="x", color="0.9")
        ax.set_xlabel(lab.split("\n")[0].replace("RMSE", "RMSE benefit").replace("accuracy", "accuracy benefit") + (" (pp)" if ds == "IMS" else " (mm)"), fontsize=8.5)
        if r_i == 0:
            unit = "pp" if ds == "IMS" else "mm"
            ax.set_title(f"{lab}\nfor scale, DA - OL gain = {gain:.2f} {unit}", fontsize=9)

handles = [plt.Line2D([], [], marker="o", ls="", color=SYC[k], label=k) for k in SYC]
handles += [plt.Line2D([], [], marker="s", ls="", color="k", label="pooled window, 95% CI (5-deg block bootstrap)")]
fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=8.5, frameon=False)
fig.suptitle("Skill benefit of VIIRS over MODIS (top) and of a second satellite (bottom); dots = single snow years", fontsize=12, y=0.995)
fig.text(0.5, 0.045, "Bottom row C4a caveat: T_to_J2 and TA_to_J1J2 differ in history before 2023-08-01, so C3 (MODIS, six years) is the cleaner one-vs-two-satellites estimate. "
         "Differences are tenths of a mm up to about 1 mm, against station RMSE of ~114-582 mm.", ha="center", fontsize=8)
fig.tight_layout(rect=[0, 0.08, 1, 0.93], h_pad=4.0)
for r_i, ttl in enumerate(["Does the MODIS -> VIIRS switch change skill?   Right of zero = VIIRS better",
                           "Does a second satellite help?   Right of zero = more satellites better"]):
    y1 = axes[r_i, 0].get_position().y1
    fig.text(0.01, y1 + (0.085 if r_i == 0 else 0.03), ttl, fontsize=11, fontweight="bold", ha="left")
out = sc.FIG / "figC_5_benefit_switch_and_constellation.png"
fig.savefig(out, dpi=150)
print(out)
