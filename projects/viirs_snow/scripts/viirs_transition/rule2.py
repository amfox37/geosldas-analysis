"""CONFIRMED decision rule (amfox, 2026-11-21; see plan 'How results will be read'). Python -I.
Skill (IMS, SNOTEL, GHCN): difference iff same sign in every snow year of the window (3; 6 for C3) AND the pooled paired
  block-bootstrap CI excludes 0. C4 DoD has no per-year values: CI-only, labelled.
Water balance: flux terms I_snow, dET, dRunoff, SWE_max, RZMC_mean keep the old rule (|A-B| >= 10% of DA-OL reference in
  every year AND same sign in every year). dStorage, residual, dFSW, dP, dSnowmelt, SWE_end, SFMC, TWLAND: descriptive, no flag.
Reads existing C_/D_ csv; writes rule2_skill.csv, rule2_water_balance.csv, rule2_changes.csv in docs/viirs_transition_report_figures/."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import skill_common as sc
FIG = sc.FIG
ct = pd.read_csv(FIG / "C_skill_contrasts.csv"); ct["season"] = ct["season"].fillna("ALL"); ct = ct[ct.season == "ALL"]
pre = pd.read_csv(FIG / "C_decision_rule_preliminary.csv")
PRIM = {"IMS": "accuracy", "snotel_SWE": "rmse", "snotel_SNWD": "rmse", "ghcn_SNWD": "rmse"}
SPEC = {"C1 T_control-T_to_J2": ("SY4-6", ["SY4", "SY5", "SY6"]), "C2 TA_control-TA_to_J1J2": ("SY4-6", ["SY4", "SY5", "SY6"]),
        "C3 MODIS_T-MODIS_TA": ("SY1-6", [f"SY{i}" for i in range(1, 7)]), "C4 DoD (T_to_J2-TA_to_J1J2)": ("DoD", None),
        "C4a T_to_J2-TA_to_J1J2 (post)": ("SY4-6", ["SY4", "SY5", "SY6"])}
def g(ds, c, w, m):
    s = ct[(ct.dataset == ds) & (ct.contrast == c) & (ct.window == w) & (ct.metric == m)]
    return None if s.empty else s.iloc[0]
rows = []
for c, (pw, yrs) in SPEC.items():
    for ds, m in PRIM.items():
        r = g(ds, c, pw, m)
        if r is None: continue
        ci = (r.ci_lo > 0) or (r.ci_hi < 0)
        if yrs is None: same, verdict = "n/a", ("difference (CI only)" if ci else "no resolvable difference")
        else:
            same = all(np.sign(g(ds, c, y, m).delta) == np.sign(r.delta) for y in yrs)
            verdict = "difference" if (same and ci) else "no resolvable difference"
        rows.append(dict(contrast=c, dataset=ds, metric=m, delta=r.delta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, same_sign_all_years=same, ci_excludes_0=ci, verdict=verdict))
sk = pd.DataFrame(rows); sk.to_csv(FIG / "rule2_skill.csv", index=False)
chg = []
m = sk.merge(pre[["contrast", "dataset", "metric", "verdict"]], on=["contrast", "dataset", "metric"], suffixes=("", "_prelim"))
for _, r in m.iterrows():
    a = r.verdict_prelim.startswith("difference"); b = r.verdict.startswith("difference")
    if a != b: chg.append(dict(analysis="C skill", contrast=r.contrast, metric=f"{r.dataset} {r.metric}", prelim=r.verdict_prelim, confirmed=r.verdict))
d = pd.read_csv(FIG / "D_decision_rule.csv")
FLAG = ["I_snow", "dET", "dRunoff", "SWE_max", "RZMC_mean"]
d["flagged_by_rule"] = d["var"].isin(FLAG)
d["verdict"] = np.where(d.flagged_by_rule, np.where(d.resolvable, "difference", "no resolvable difference"), "descriptive (no flag)")
d.to_csv(FIG / "rule2_water_balance.csv", index=False)
for _, r in d.iterrows():
    if r.resolvable and not r.flagged_by_rule:
        chg.append(dict(analysis="D water", contrast=r.contrast, metric=r["var"], prelim="resolvable", confirmed="descriptive (no flag)"))
pd.DataFrame(chg).to_csv(FIG / "rule2_changes.csv", index=False)
pd.set_option("display.width", 250)
print(sk.round(4).to_string()); print(d[d.flagged_by_rule][["contrast", "var", "mean_AminusB", "ratio", "same_sign_all_years", "ge10pct_all_years", "verdict"]].round(3).to_string()); print(pd.DataFrame(chg).to_string())
