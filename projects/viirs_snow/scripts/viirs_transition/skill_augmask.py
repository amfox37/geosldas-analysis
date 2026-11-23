"""E2: re-score C1 and C2 (IMS, SNOTEL, GHCN) with the conservative 1-28 Aug 2023 mask (analysis E recommendation).
Reuses skill_score.py unchanged (imported as a module) and the existing skill_data npz.
 - stations: the 28 days are set to NaN in the observations (station selection is left as in C, i.e. made before the mask);
 - IMS: the per-cell counts are stored by snow year x season only, so the counts of the 28 days (SY4 x JJA), re-extracted by
   skill_augmask_ims.py for the four arms of C1/C2, are subtracted from skill_<arm>.npz.
Both passes (nomask = reproduces C, augmask) use the same fresh RNG seed so the block-bootstrap draws are paired.
Flags: same logic as rule2.py (same sign in every snow year AND pooled bootstrap CI excludes 0); rule2.py is not imported
because it writes rule2_*.csv. Outputs: augmask_skill_contrasts.csv, augmask_skill_by_arm.csv, augmask_flags.csv,
augmask_delta.csv in docs/viirs_transition_report_figures/. Run with python -I."""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
import skill_score as ss

KEEP = ["C1 T_control-T_to_J2", "C2 TA_control-TA_to_J1J2"]
ss.CONTRASTS = {k: v for k, v in ss.CONTRASTS.items() if k in KEEP}
ARMS4 = ["T_to_J2", "TA_to_J1J2", "T_control", "TA_control"]
MASK_DAYS = (ss.DAYS >= "2023-08-01") & (ss.DAYS <= "2023-08-28")
Z0 = ss.Z
Z0d = {a: {k: z[k] for k in z.files} for a, z in Z0.items()}
Z1 = {a: dict(d) for a, d in Z0d.items()}
tot = 0
for a in ARMS4:
    aug = np.load(sc.OUT / f"augmask_ims_cnt_{a}.npz")["cnt"].astype(np.int32)       # 4, cell
    c = Z1[a]["cnt"].astype(np.int32); c[:, 3, 2, :] -= aug
    assert c.min() >= 0, a
    Z1[a]["cnt"] = c; tot += int(aug.sum())
print("IMS cell-days removed (4 arms):", tot, flush=True)

abs_rows, ct_rows = [], []
station_data = {k: ss.build_dataset(k) for k in ("snotel", "ghcn")}
for label in ("nomask", "augmask"):
    ss.RNG = np.random.default_rng(42)
    ss.Z = Z0d if label == "nomask" else Z1
    a, c = ss.run_ims(); a.insert(0, "dataset", "IMS"); c.insert(0, "dataset", "IMS")
    A, C = [a], [c]
    for kind, data in station_data.items():
        d2 = {v: dict(D, o=D["o"].copy()) for v, D in data.items()}
        if label == "augmask":
            for v, D in d2.items():
                D["o"][MASK_DAYS] = np.nan
                print(f"{kind} {v}: stations {D['lat'].size}, S-hemisphere stations {(D['lat'] < 0).sum()}; obs-days masked {int(np.isfinite(data[v]['o'][MASK_DAYS]).sum())} "
                      f"of which in S hemisphere {int(np.isfinite(data[v]['o'][MASK_DAYS][:, D['lat'] < 0]).sum())}", flush=True)
        a, c = ss.run_stations(kind, d2); A.append(a); C.append(c)
    for df in A: df.insert(0, "pass", label); abs_rows.append(df)
    for df in C: df.insert(0, "pass", label); ct_rows.append(df)
abs_all = pd.concat(abs_rows, ignore_index=True); ct_all = pd.concat(ct_rows, ignore_index=True)
abs_all.to_csv(sc.FIG / "augmask_skill_by_arm.csv", index=False); ct_all.to_csv(sc.FIG / "augmask_skill_contrasts.csv", index=False)

# ---- flags (rule2 logic) and change table
PRIM = {"IMS": "accuracy", "snotel_SWE": "rmse", "snotel_SNWD": "rmse", "ghcn_SNWD": "rmse"}
YRS = ["SY4", "SY5", "SY6"]
ct_all["season"] = ct_all["season"].fillna("ALL")
cta = ct_all[ct_all.season == "ALL"]
def g(p, ds, c, w, m):
    s = cta[(cta["pass"] == p) & (cta.dataset == ds) & (cta.contrast == c) & (cta.window == w) & (cta.metric == m)]
    return s.iloc[0]
rows = []
for p in ("nomask", "augmask"):
    for c in KEEP:
        for ds, m in PRIM.items():
            r = g(p, ds, c, "SY4-6", m); ci = (r.ci_lo > 0) or (r.ci_hi < 0)
            peryear = {y: float(g(p, ds, c, y, m).delta) for y in YRS}
            same = all(np.sign(v) == np.sign(r.delta) for v in peryear.values())
            rows.append(dict(**{"pass": p}, contrast=c, dataset=ds, metric=m, delta=r.delta, ci_lo=r.ci_lo, ci_hi=r.ci_hi,
                             **{f"delta_{y}": v for y, v in peryear.items()}, same_sign_all_years=same, ci_excludes_0=ci,
                             verdict="difference" if (same and ci) else "no resolvable difference"))
fl = pd.DataFrame(rows); fl.to_csv(sc.FIG / "augmask_flags.csv", index=False)
w = fl.pivot_table(index=["contrast", "dataset", "metric"], columns="pass", values=["delta", "ci_lo", "ci_hi", "delta_SY4"], aggfunc="first")
w.columns = [f"{a}_{b}" for a, b in w.columns]; w = w.reset_index()
v = fl.pivot_table(index=["contrast", "dataset", "metric"], columns="pass", values="verdict", aggfunc="first").reset_index()
w = w.merge(v.rename(columns={"nomask": "verdict_nomask", "augmask": "verdict_augmask"}), on=["contrast", "dataset", "metric"])
w["d_delta"] = w["delta_augmask"] - w["delta_nomask"]; w["verdict_changed"] = w.verdict_nomask != w.verdict_augmask
w.to_csv(sc.FIG / "augmask_delta.csv", index=False)
pd.set_option("display.width", 250)
print(w[["contrast", "dataset", "metric", "delta_nomask", "delta_augmask", "d_delta", "verdict_nomask", "verdict_augmask"]].round(5).to_string())
# IMS season/abs checks
x = abs_all[(abs_all.dataset == "IMS") & (abs_all.window == "SY4") & (abs_all.season == "ALL") & abs_all.unit.isin(ARMS4)]
print(x.pivot_table(index="unit", columns="pass", values=["N", "accuracy"]).round(5).to_string())
