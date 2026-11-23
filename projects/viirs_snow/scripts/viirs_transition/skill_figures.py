"""Forest plots and a short results table for analysis C from C_skill_contrasts.csv / C_skill_by_arm.csv (python -I)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
ct = pd.read_csv(sc.FIG / "C_skill_contrasts.csv"); ab = pd.read_csv(sc.FIG / "C_skill_by_arm.csv")
# IMS rows: use the flip-corrected scoring (IMS_skill_contrasts_fixed.csv; 2024+ IMS arrays are vertically flipped in the original C scoring)
_fx = pd.read_csv(sc.FIG / "IMS_skill_contrasts_fixed.csv").assign(dataset="IMS")
ct = pd.concat([ct[ct.dataset != "IMS"], _fx[[c for c in ct.columns if c in _fx.columns]]], ignore_index=True)
ORDER = ["C0 T_to_J2-OL", "C0 TA_to_J1J2-OL", "C0 MODIS_T-OL", "C0 MODIS_TA-OL", "C0 T_control-OL", "C0 TA_control-OL",
         "C1 T_control-T_to_J2", "C2 TA_control-TA_to_J1J2", "C3 MODIS_T-MODIS_TA", "C4 DoD (T_to_J2-TA_to_J1J2)",
         "C4a T_to_J2-TA_to_J1J2 (post)", "C4b T_to_J2-TA_to_J1J2 (pre)"]
POOL = {"C0": None, "C1": "SY4-6", "C2": "SY4-6", "C3": "SY1-6", "C4": None}
def pooled_win(name):
    if name.startswith("C4 DoD"): return "DoD"
    if name.startswith("C0"): return "SY4-6" if "control" in name else "SY1-6"
    if name.startswith("C4b"): return "SY1-3"
    if name.startswith("C4a"): return "SY4-6"
    return POOL[name[:2]]
SYC = ["#9ecae1", "#6baed6", "#4292c6", "#fdae6b", "#fd8d3c", "#e6550d"]
CAV = {"IMS": "IMS is NOT independent of MODIS/VIIRS imagery", "snotel_SWE": "SNOTEL: western US/Alaska only", "snotel_SNWD": "SNOTEL: western US/Alaska only",
       "ghcn_SNWD": "GHCN: sparse, mostly mid-latitude"}
UNIT = {"IMS": {"accuracy": "accuracy", "csi": "CSI"}, }
def forest(ds, metrics, fname, title, units):
    fig, axes = plt.subplots(1, len(metrics), figsize=(5.2 * len(metrics) + 2.5, 5.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, m in zip(axes, metrics):
        for yi, name in enumerate(ORDER):
            y = len(ORDER) - 1 - yi
            sub = ct[(ct.dataset == ds) & (ct.contrast == name) & (ct.metric == m) & (ct.season.fillna("ALL") == "ALL" if "season" in ct else True)]
            if sub.empty: continue
            pw = pooled_win(name); p = sub[sub.window == pw]
            for k in range(6):
                s = sub[sub.window == f"SY{k+1}"]
                if len(s): ax.plot(s.delta.values[0], y + 0.18, "o", ms=3.5, color=SYC[k], zorder=3, label=f"SY{k+1}" if (yi == 0 or name == ORDER[6]) else None)
            if len(p):
                r = p.iloc[0]; ax.errorbar(r.delta, y, xerr=[[r.delta - r.ci_lo], [r.ci_hi - r.delta]], fmt="s", color="k", ms=5, capsize=3, zorder=4)
        ax.axvline(0, color="grey", lw=0.8); ax.set_title(m + f" ({units.get(m, '')})", fontsize=10); ax.grid(axis="x", alpha=0.3)
    axes[0].set_yticks(range(len(ORDER))); axes[0].set_yticklabels(ORDER[::-1], fontsize=8)
    h, l = axes[-1].get_legend_handles_labels(); d = dict(zip(l, h)); axes[-1].legend(d.values(), d.keys(), fontsize=7, loc="lower right", title="per snow year", title_fontsize=7)
    fig.suptitle(f"{title}\nBlack square = pooled window with 95% CI (paired 5-degree block bootstrap); dots = single snow years; delta = first minus second term\n{CAV[ds]}",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.86)); fig.savefig(sc.FIG / fname, dpi=150); plt.close(fig)
ct["season"] = ct["season"] if "season" in ct else "ALL"
ct["season"] = ct["season"].fillna("ALL")
forest("IMS", ["accuracy", "csi"], "figC_1_ims_contrasts.png", "IMS snow cover: contrasts C0-C4 (positive accuracy/CSI delta = first term better)", {"accuracy": "-", "csi": "-"})
forest("snotel_SWE", ["rmse", "ubrmse"], "figC_2_snotel_swe_contrasts.png", "SNOTEL SWE (non-JJA): contrasts (negative delta = first term better)", {"rmse": "mm", "ubrmse": "mm"})
forest("snotel_SNWD", ["rmse", "ubrmse"], "figC_3_snotel_depth_contrasts.png", "SNOTEL snow depth (non-JJA): contrasts (negative delta = first term better)", {"rmse": "mm", "ubrmse": "mm"})
forest("ghcn_SNWD", ["rmse", "ubrmse"], "figC_4_ghcn_depth_contrasts.png", "GHCN snow depth (non-JJA): contrasts (negative delta = first term better)", {"rmse": "mm", "ubrmse": "mm"})
# compact pooled table
rows = []
primary = {"IMS": "accuracy", "snotel_SWE": "rmse", "snotel_SNWD": "rmse", "ghcn_SNWD": "rmse"}
for name in ORDER:
    r = {"contrast": name, "window": pooled_win(name)}
    for ds, m in primary.items():
        s = ct[(ct.dataset == ds) & (ct.contrast == name) & (ct.metric == m) & (ct.window == r["window"]) & (ct.season == "ALL")]
        r[f"{ds} {m}"] = "" if s.empty else f"{s.delta.iloc[0]:+.4g} [{s.ci_lo.iloc[0]:+.3g},{s.ci_hi.iloc[0]:+.3g}]"
    rows.append(r)
pd.DataFrame(rows).to_csv(sc.FIG / "C_pooled_contrasts_table.csv", index=False)
print(pd.DataFrame(rows).to_string())
