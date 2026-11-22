"""Monthly IMS accuracy time series (NH >= 30N and all cells) for the arms, with the original (as scored in C) and the
flip-corrected IMS; DA-OL by month; monthly hit/false-alarm/miss shares. Reads <OUT>/imsfix/ims_<arm>.npz."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
cfg = sc.load_arms(); arms = list(cfg["arms"])
rows = []
for a in arms:
    z = np.load(sc.OUT / "imsfix" / f"ims_{a}.npz"); d = pd.DatetimeIndex(z["days"]); dst = z["dst"].astype(float); f = z["found"]
    for v, vn in enumerate(("as_scored", "fixed")):
        for r, rn in enumerate(("NH30", "all")):
            x = dst[:, v, r, :]
            df = pd.DataFrame(x[f], columns=list("ABCD"), index=d[f]); df["ym"] = df.index.to_period("M")
            m = df.groupby("ym")[list("ABCD")].sum()
            m["arm"], m["ims"], m["region"] = a, vn, rn
            rows.append(m.reset_index())
M = pd.concat(rows, ignore_index=True); N = M[list("ABCD")].sum(axis=1)
M["accuracy"] = (M.A + M.D) / N; M["hit_rate"] = M.A / (M.A + M.C); M["far"] = M.B / (M.A + M.B); M["csi"] = M.A / (M.A + M.B + M.C)
M["ym"] = M.ym.astype(str)
M.to_csv(sc.FIG / "IMS_monthly_accuracy.csv", index=False)
# DA - OL monthly accuracy
P = M[M.region == "NH30"].pivot_table(index=["ims", "ym"], columns="arm", values="accuracy")
for a in arms:
    if a != "OL": P[f"{a}-OL"] = P[a] - P["OL"]
P.reset_index().to_csv(sc.FIG / "IMS_monthly_DA_minus_OL.csv", index=False)
fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
cols = {"OL": "#666666", "T_to_J2": "#1f77b4", "TA_to_J1J2": "#d62728", "T_control": "#17becf", "TA_control": "#ff7f0e"}
for k, vn in enumerate(("as_scored", "fixed")):
    for a in arms:
        s = M[(M.arm == a) & (M.ims == vn) & (M.region == "NH30")]
        ax[k].plot(pd.PeriodIndex(s.ym, freq="M").to_timestamp(), s.accuracy, color=cols[a], label=a, lw=1.4)
    ax[k].set_ylabel("NH>=30N IMS accuracy"); ax[k].set_title("IMS as scored in analysis C (2024+ not flipped)" if k == 0 else "IMS with 2024+ vertical flip corrected")
    ax[k].axvline(pd.Timestamp("2024-01-01"), color="k", ls=":", lw=1)
ax[0].legend(ncol=5, fontsize=8)
for a in arms:
    if a == "OL": continue
    s = P.loc["fixed", f"{a}-OL"].dropna()
    ax[2].plot(pd.PeriodIndex(s.index, freq="M").to_timestamp(), s.values, color=cols[a], label=f"{a} - OL")
ax[2].axhline(0, color="k", lw=.7); ax[2].set_ylabel("accuracy DA - OL (fixed IMS)"); ax[2].legend(ncol=4, fontsize=8)
fig.tight_layout(); fig.savefig(sc.FIG / "IMS_3_monthly_accuracy.png", dpi=110)
pd.set_option("display.width", 200)
print(P.loc["fixed"].round(4).tail(40).to_string())
