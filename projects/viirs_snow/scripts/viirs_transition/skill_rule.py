"""PRELIMINARY application of the proposed decision rule (plan 'How results will be read'; NOT yet confirmed by amfox):
a contrast is a difference only if |delta| >= 10% of the mean |DA-OL| of the two compared arms over the same pooled window
and delta has the same sign in all snow years of the window. Writes C_decision_rule_preliminary.csv (python -I)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import skill_common as sc
ct = pd.read_csv(sc.FIG / "C_skill_contrasts.csv"); ct["season"] = ct["season"].fillna("ALL"); ct = ct[ct.season == "ALL"]
PRIM = {"IMS": "accuracy", "snotel_SWE": "rmse", "snotel_SNWD": "rmse", "ghcn_SNWD": "rmse"}
SPEC = {  # contrast -> (pooled window, snow years, units for the reference DA-OL effect and window used for it)
 "C1 T_control-T_to_J2": ("SY4-6", ["SY4", "SY5", "SY6"], [("C0 T_control-OL", "SY4-6"), ("C0 T_to_J2-OL", "SY4-6")]),
 "C2 TA_control-TA_to_J1J2": ("SY4-6", ["SY4", "SY5", "SY6"], [("C0 TA_control-OL", "SY4-6"), ("C0 TA_to_J1J2-OL", "SY4-6")]),
 "C3 MODIS_T-MODIS_TA": ("SY1-6", [f"SY{i}" for i in range(1, 7)], [("C0 MODIS_T-OL", "SY1-6"), ("C0 MODIS_TA-OL", "SY1-6")]),
 "C4 DoD (T_to_J2-TA_to_J1J2)": ("DoD", None, [("C0 T_to_J2-OL", "SY1-6"), ("C0 TA_to_J1J2-OL", "SY1-6")]),
 "C4a T_to_J2-TA_to_J1J2 (post)": ("SY4-6", ["SY4", "SY5", "SY6"], [("C0 T_to_J2-OL", "SY4-6"), ("C0 TA_to_J1J2-OL", "SY4-6")])}
def g(ds, c, w, m):
    s = ct[(ct.dataset == ds) & (ct.contrast == c) & (ct.window == w) & (ct.metric == m)]
    return None if s.empty else s.iloc[0]
rows = []
for c, (pw, sys_, refs) in SPEC.items():
    for ds, m in PRIM.items():
        r = g(ds, c, pw, m)
        if r is None: continue
        def refval(rc, rw):
            x = g(ds, rc, rw, m)
            if x is not None: return abs(x.delta)
            return abs(np.mean([g(ds, rc, f"SY{i}", m).delta for i in (4, 5, 6)]))   # SY4-6 pooled not evaluated: mean of per-year deltas
        ref = np.mean([refval(rc, rw) for rc, rw in refs])
        if sys_ is None:
            sign_ok = None
            # DoD has no per-year values: use sign of the pooled CI instead
            ci_excl0 = (r.ci_lo > 0) or (r.ci_hi < 0); same = "n/a (CI excludes 0)" if ci_excl0 else "n/a (CI includes 0)"
        else:
            ds_ = [g(ds, c, y, m).delta for y in sys_]; sign_ok = bool(all(np.sign(d) == np.sign(r.delta) for d in ds_)); same = str(sign_ok)
        big = abs(r.delta) >= 0.10 * ref
        rows.append(dict(contrast=c, dataset=ds, metric=m, delta=r.delta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, ref_DA_minus_OL=ref, pct_of_ref=100 * abs(r.delta) / ref if ref else np.nan,
                         ge_10pct=big, same_sign_all_years=same,
                         verdict="difference" if (big and sign_ok is True) else ("no resolvable difference" if sign_ok is not None else ("difference-candidate (CI excl. 0, >=10%)" if big and ci_excl0 else "no resolvable difference"))))
out = pd.DataFrame(rows); out.to_csv(sc.FIG / "C_decision_rule_preliminary.csv", index=False)
pd.set_option("display.width", 250); print(out.round(4).to_string())
