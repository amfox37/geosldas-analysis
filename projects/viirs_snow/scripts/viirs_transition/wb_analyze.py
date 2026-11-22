#!/usr/bin/env python3
"""Analysis D, step 2: water-balance differences for contrasts C0-C4 (docs/viirs_transition_plan.md).

Input: output/water_balance/wb_<arm>.nc from wb_extract.py. Output: tables and figures with prefix D_ in
docs/viirs_transition_report_figures/.

Budget definitions follow projects/M21C_ls/scripts/water_year_snow_da_budget.py (differential budget A - B, per tile,
per snow year):  residual = dI_snow - dET - dRunoff - dStorage - dFSW, where
  I_snow = summed WESNN1-3_INCR, Runoff = RUNSURFLAND + BASEFLOWLAND, dStorage = change of the A-B TWLAND difference
  between 1 Aug restarts, dFSW = PEATCLSM_FSWCHANGE total. For C0 (DA - OL) B = OL, whose I_snow is 0.
Domain: Northern-Hemisphere seasonal-snow mask as in M21C_ls (lat > 20N, snow possible, mean JJA SCF < 0.2), rebuilt
from the arms' own SWE/SCF (the M21C tile count 48067 is for the M21C record, not expected here).
Domain means are area weighted. The 5-degree block bootstrap (M21C_ls block_codes/confidence_interval) captures
spatial uncertainty only; with 3 (C1, C2, C4) or 6 (C3) snow years the year-to-year spread is the larger uncertainty,
so every snow year is reported separately.
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset

REPO = "/discover/nobackup/projects/land_da/geosldas-analysis"
PROJ = f"{REPO}/projects/viirs_snow"
sys.path.insert(0, f"{REPO}/projects/M21C_ls/scripts")
from analysis_a_robustness import block_codes, confidence_interval  # noqa: E402

CATCHMENT_DEF = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def"
SY_IDS = ["SY1", "SY2", "SY3", "SY4", "SY5", "SY6"]
# variable -> (label, kg m-2 unless stated)
VARS = ["I_snow", "dET", "dRunoff", "dStorage", "dFSW", "residual", "dP", "dSnowmelt", "SWE_max", "SWE_end",
        "RZMC_mean", "SFMC_mean", "TWLAND_mean"]
KEY = ["I_snow", "dET", "dRunoff", "dStorage", "SWE_max", "SWE_end", "RZMC_mean"]
REPL = 1000


def load_arm(arm, d):
    with Dataset(f"{d}/wb_{arm}.nc") as o:
        out = {"sy": list(o.variables["sy_id"][:]), "missing": json.loads(o.missing),
               "lat": np.array(o["lat"][:], float), "lon": np.array(o["lon"][:], float),
               "rst_year": list(o["endpoint_year"][:])}
        for k, v in o.variables.items():
            if k in ("sy_id", "lat", "lon", "endpoint_year", "month_order"):
                continue
            out[k] = np.ma.filled(v[:].astype(np.float64), np.nan)
    return out


def series(arms, spec):
    """spec: list of (arm, sy list) -> dict sy_id -> per-tile record (dict of fields)."""
    rec = {}
    for arm, sys_ in spec:
        a = arms[arm]
        for sid in sys_:
            i = a["sy"].index(sid)
            r = {k: a[k][i] for k in ["P", "ET", "RUNSURF", "BASEFLOW", "SNOWMELT", "FSW", "WCHANGE", "I_snow",
                                      "I_snow_abs", "n_incr", "SWE_max", "SWE_end"]}
            r["RZMC_mean"] = np.nanmean(a["m_RZMC"][i], 0)
            r["SFMC_mean"] = np.nanmean(a["m_SFMC"][i], 0)
            r["TWLAND_mean"] = np.nanmean(a["m_TWLAND"][i], 0)
            r["TW0"], r["TW1"] = a["TWLAND_rst"][i], a["TWLAND_rst"][i + 1]
            for k in ["SWE", "RZMC", "SFMC", "TWLAND", "SCF", "ET", "RUNOFF", "SNOWMELT", "P", "I_snow"]:
                r["m_" + k] = a["m_" + k][i]
            rec[sid] = r
    return rec


def diff(A, B, sid):
    a, b = A[sid], B[sid]
    d = {"I_snow": a["I_snow"] - b["I_snow"],
         "dET": a["ET"] - b["ET"],
         "dRunoff": (a["RUNSURF"] + a["BASEFLOW"]) - (b["RUNSURF"] + b["BASEFLOW"]),
         "dStorage": (a["TW1"] - b["TW1"]) - (a["TW0"] - b["TW0"]),
         "dFSW": a["FSW"] - b["FSW"], "dP": a["P"] - b["P"], "dSnowmelt": a["SNOWMELT"] - b["SNOWMELT"]}
    d["residual"] = d["I_snow"] - d["dET"] - d["dRunoff"] - d["dStorage"] - d["dFSW"]
    for k in ["SWE_max", "SWE_end", "RZMC_mean", "SFMC_mean", "TWLAND_mean"]:
        d[k] = a[k] - b[k]
    d["monthly"] = {k: a["m_" + k] - b["m_" + k] for k in ["SWE", "RZMC", "ET", "RUNOFF", "SNOWMELT", "I_snow"]}
    return d


def wmean(x, w, keep):
    v = keep & np.isfinite(x)
    return float(np.sum(x[v] * w[v]) / np.sum(w[v])) if v.any() else np.nan


def boot(x, w, keep, codes, rng):
    v = keep & np.isfinite(x)
    nb = codes.max() + 1
    S = np.bincount(codes[v], weights=(x * w)[v], minlength=nb)
    W = np.bincount(codes[v], weights=w[v], minlength=nb)
    idx = rng.integers(0, nb, size=(REPL, nb))
    m = S[idx].sum(1) / np.maximum(W[idx].sum(1), 1e-30)
    return confidence_interval(m)


def mdtable(df):
    cols = list(df.columns)
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
                     + ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=f"{PROJ}/output/water_balance")
    ap.add_argument("--figdir", default=f"{PROJ}/docs/viirs_transition_report_figures")
    a = ap.parse_args()
    cfg = json.load(open(f"{PROJ}/config/arms.json"))
    names = list(cfg["arms"])
    arms = {n: load_arm(n, a.data) for n in names}
    lat, lon = arms["OL"]["lat"], arms["OL"]["lon"]
    area = np.loadtxt(CATCHMENT_DEF, usecols=(6,), skiprows=1)
    assert area.size == lat.size
    # seasonal-snow mask (M21C_ls definition) from the 6-year OL record
    ol = arms["OL"]
    swe_max = np.nanmax(ol["SWE_max"], 0)
    scf = ol["m_SCF"]  # (sy, month[Aug..Jul], tile)
    scf_max = np.nanmax(scf.reshape(-1, lat.size), 0)
    jja = np.nanmean(scf[:, [10, 11, 0], :].reshape(-1, lat.size), 0)  # Jun, Jul, Aug
    mask = (lat > 20) & ((scf_max > 0.05) | (swe_max > 5)) & (jja < 0.2) & np.isfinite(lat)
    print("mask tiles:", int(mask.sum()))
    codes = block_codes(pd.DataFrame({"lat": lat, "lon": lon}), 5.0)
    rng = np.random.default_rng(21022)

    S1, S3, S4, S6 = ["SY1", "SY2", "SY3"], ["SY4", "SY5", "SY6"], None, SY_IDS
    ser = {
        "OL": series(arms, [("OL", SY_IDS)]),
        "T_to_J2": series(arms, [("T_to_J2", SY_IDS)]),
        "TA_to_J1J2": series(arms, [("TA_to_J1J2", SY_IDS)]),
        "T_control": series(arms, [("T_control", S3)]),
        "TA_control": series(arms, [("TA_control", S3)]),
        "MODIS_T": series(arms, [("T_to_J2", S1), ("T_control", S3)]),
        "MODIS_TA": series(arms, [("TA_to_J1J2", S1), ("TA_control", S3)]),
    }
    rows, mrows = [], []

    def record(label, d, sid, kind):
        for k in VARS:
            x = d[k]
            lo, hi = boot(x, area, mask, codes, rng)
            rows.append(dict(kind=kind, contrast=label, sy=sid, var=k, mean=wmean(x, area, mask), lo=lo, hi=hi))

    # C0: DA - OL
    c0 = {}
    for arm, s in [("T_to_J2", ser["T_to_J2"]), ("TA_to_J1J2", ser["TA_to_J1J2"]),
                   ("MODIS_T", ser["MODIS_T"]), ("MODIS_TA", ser["MODIS_TA"]),
                   ("T_control", ser["T_control"]), ("TA_control", ser["TA_control"])]:
        for sid in s:
            d = diff(s, ser["OL"], sid)
            c0[(arm, sid)] = d
            record(f"{arm}-OL", d, sid, "C0")
            for k, v in d["monthly"].items():
                for m in range(12):
                    mrows.append(dict(contrast=f"{arm}-OL", sy=sid, var=k, month_idx=m, mean=wmean(v[m], area, mask)))
    # contrasts
    cdef = {"C1": ("T_control", "T_to_J2", S3), "C2": ("TA_control", "TA_to_J1J2", S3),
            "C3": ("MODIS_T", "MODIS_TA", SY_IDS), "C4": ("T_to_J2", "TA_to_J1J2", S3)}
    cd = {}
    for c, (A, B, sy) in cdef.items():
        for sid in sy:
            d = diff(ser[A], ser[B], sid)
            cd[(c, sid)] = d
            record(f"{c}:{A}-{B}", d, sid, c)
            for k, v in d["monthly"].items():
                for m in range(12):
                    mrows.append(dict(contrast=f"{c}:{A}-{B}", sy=sid, var=k, month_idx=m, mean=wmean(v[m], area, mask)))
    # C4 difference of differences: (J2 vs J1J2 in SY4-6) - (same T vs TA difference in SY1-3, mean)
    base = {k: np.nanmean([diff(ser["T_to_J2"], ser["TA_to_J1J2"], s)[k] for s in S1], 0) for k in VARS}
    for sid in S3:
        d = {k: cd[("C4", sid)][k] - base[k] for k in VARS}
        record("C4dd:(T_to_J2-TA_to_J1J2)[SY4-6] - mean[SY1-3]", d, sid, "C4dd")
    df = pd.DataFrame(rows)
    os.makedirs(a.figdir, exist_ok=True)
    df.to_csv(f"{a.figdir}/D_budget_tables.csv", index=False)
    pd.DataFrame(mrows).to_csv(f"{a.figdir}/D_monthly_domain_means.csv", index=False)

    # decision rule: |A-B| >= 10% of the DA-OL effect (mean of the two arms' |domain mean DA-OL| in that SY), same sign
    ref_arm = {"C1": ("T_control", "T_to_J2"), "C2": ("TA_control", "TA_to_J1J2"),
               "C3": ("MODIS_T", "MODIS_TA"), "C4": ("T_to_J2", "TA_to_J1J2")}
    dec = []
    for c, (A, B, sy) in cdef.items():
        for k in VARS:
            vals, refs = [], []
            for sid in sy:
                vals.append(df[(df.contrast == f"{c}:{A}-{B}") & (df.sy == sid) & (df["var"] == k)]["mean"].iloc[0])
                r = [abs(df[(df.contrast == f"{x}-OL") & (df.sy == sid) & (df["var"] == k)]["mean"].iloc[0]) for x in ref_arm[c]]
                refs.append(np.mean(r))
            vals, refs = np.array(vals), np.array(refs)
            same = bool(np.all(vals > 0) or np.all(vals < 0))
            big = bool(np.all(np.abs(vals) >= 0.10 * refs))
            dec.append(dict(contrast=c, var=k, mean_AminusB=vals.mean(), min=vals.min(), max=vals.max(),
                            mean_ref_DAminusOL=refs.mean(), ratio=vals.mean() / refs.mean() if refs.mean() else np.nan,
                            same_sign_all_years=same, ge10pct_all_years=big, resolvable=same and big))
    dec = pd.DataFrame(dec)
    dec.to_csv(f"{a.figdir}/D_decision_rule.csv", index=False)

    # ---- figures
    sub = lambda kind, lab, k: df[(df.kind == kind) & (df.contrast == lab) & (df["var"] == k)]
    # D_1: C0 budget per snow year, TA and T series
    terms = ["I_snow", "dET", "dRunoff", "dStorage", "dFSW", "residual"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, arm in zip(axs, ["MODIS_T", "MODIS_TA"]):
        w = 0.13
        for j, k in enumerate(terms):
            s = sub("C0", f"{arm}-OL", k)
            ax.bar(np.arange(len(s)) + (j - 2.5) * w, s["mean"], w, label=k)
        ax.set_xticks(range(6)); ax.set_xticklabels(SY_IDS); ax.axhline(0, c="k", lw=.5)
        ax.set_title(f"{arm} minus OL (domain mean, seasonal-snow mask)")
    axs[0].set_ylabel("kg m-2 per snow year"); axs[0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(f"{a.figdir}/D_1_C0_budget_modis.png", dpi=130); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, arm in zip(axs, ["T_to_J2", "TA_to_J1J2"]):
        for j, k in enumerate(terms):
            s = sub("C0", f"{arm}-OL", k)
            ax.bar(np.arange(len(s)) + (j - 2.5) * 0.13, s["mean"], 0.13, label=k)
        ax.set_xticks(range(6)); ax.set_xticklabels(SY_IDS); ax.axhline(0, c="k", lw=.5)
        ax.axvline(2.5, c="gray", ls=":"); ax.set_title(f"{arm} minus OL (switch to VIIRS after SY3)")
    axs[0].set_ylabel("kg m-2 per snow year"); axs[0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(f"{a.figdir}/D_2_C0_budget_viirs_arms.png", dpi=130); plt.close(fig)
    # D_3: contrasts C1-C4 for key variables per SY with bootstrap CI
    fig, axs = plt.subplots(2, 4, figsize=(15, 6))
    for ax, k in zip(axs.flat, KEY):
        for ci, (c, (A, B, sy)) in enumerate(cdef.items()):
            s = sub(c, f"{c}:{A}-{B}", k)
            x = np.arange(len(s)) * 1.0 + ci * 0.0
            xs = [SY_IDS.index(i) + (ci - 1.5) * 0.18 for i in s["sy"]]
            ax.errorbar(xs, s["mean"], yerr=[s["mean"] - s["lo"], s["hi"] - s["mean"]], fmt="o", ms=3, label=c)
        ax.axhline(0, c="k", lw=.5); ax.set_xticks(range(6)); ax.set_xticklabels(SY_IDS, fontsize=7)
        ax.set_title(k);
    axs.flat[0].legend(fontsize=7)
    axs.flat[-1].axis("off")
    fig.suptitle("A - B domain means by snow year (bars: spatial 5-degree block bootstrap 95% CI only)\n"
                 "C1 T_control-T_to_J2; C2 TA_control-TA_to_J1J2; C3 MODIS_T-MODIS_TA; C4 T_to_J2-TA_to_J1J2", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{a.figdir}/D_3_contrasts_by_snow_year.png", dpi=130); plt.close(fig)
    # D_4: seasonal cycle of DA-OL SWE and runoff (monthly domain means), TA vs T MODIS and VIIRS (SY4-6 mean)
    md = pd.DataFrame(mrows)
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.8))
    lab = ["Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"]
    for ax, k in zip(axs, ["SWE", "RUNOFF", "RZMC"]):
        for contrast, col in [("T_control-OL", "C0"), ("TA_control-OL", "C1"), ("T_to_J2-OL", "C2"), ("TA_to_J1J2-OL", "C3")]:
            s = md[(md.contrast == contrast) & (md["var"] == k) & (md.sy.isin(S3))].groupby("month_idx")["mean"].mean()
            ax.plot(range(12), s.values, label=contrast.replace("-OL", ""))
        ax.axhline(0, c="k", lw=.5); ax.set_xticks(range(12)); ax.set_xticklabels(lab, fontsize=7)
        ax.set_title(f"DA - OL, {k}, mean SY4-SY6")
    axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(f"{a.figdir}/D_4_seasonal_cycle_SY4_6.png", dpi=130); plt.close(fig)
    # D_5: maps of pooled differences (C1, C2, C3: runoff, SWE_end) and C0 TA_to_J1J2 I_snow
    def pooled(A, B, sy, k):
        return np.nanmean([diff(ser[A], ser[B], s)[k] for s in sy], 0)
    panels = [("TA_to_J1J2-OL, I_snow", pooled("TA_to_J1J2", "OL", S3, "I_snow")),
              ("TA_to_J1J2-OL, dRunoff", pooled("TA_to_J1J2", "OL", S3, "dRunoff")),
              ("C2 TA_control-TA_to_J1J2, I_snow", pooled("TA_control", "TA_to_J1J2", S3, "I_snow")),
              ("C2, dRunoff", pooled("TA_control", "TA_to_J1J2", S3, "dRunoff")),
              ("C1 T_control-T_to_J2, I_snow", pooled("T_control", "T_to_J2", S3, "I_snow")),
              ("C3 MODIS_T-MODIS_TA, I_snow", pooled("MODIS_T", "MODIS_TA", SY_IDS, "I_snow"))]
    fig, axs = plt.subplots(2, 3, figsize=(15, 7))
    for ax, (t, f) in zip(axs.flat, panels):
        v = np.where(mask, f, np.nan)
        lim = np.nanpercentile(np.abs(v), 98) or 1
        sc = ax.scatter(lon, lat, c=v, s=1.5, cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_xlim(-180, 180); ax.set_ylim(20, 85); ax.set_title(t + " (kg m-2 / snow year, mean SY4-6 or SY1-6)", fontsize=8)
        fig.colorbar(sc, ax=ax, shrink=.8)
    fig.tight_layout(); fig.savefig(f"{a.figdir}/D_5_maps_pooled.png", dpi=120); plt.close(fig)

    # ---- markdown summary
    def table(kind, lab, ks):
        out = ["| SY | " + " | ".join(ks) + " |", "|---|" + "---|" * len(ks)]
        for sid in SY_IDS:
            r = []
            for k in ks:
                s = sub(kind, lab, k); s = s[s.sy == sid]
                r.append(f"{s['mean'].iloc[0]:.2f}" if len(s) else "")
            if any(r):
                out.append(f"| {sid} | " + " | ".join(r) + " |")
        return "\n".join(out)
    md_out = [f"# Analysis D tables (auto-generated by wb_analyze.py)\n\nMask tiles: {int(mask.sum())}. Units kg m-2 per snow year "
              f"(RZMC in m3 m-3). Domain means area weighted.\n"]
    ks = ["I_snow", "dET", "dRunoff", "dStorage", "dFSW", "residual", "dP", "SWE_max", "SWE_end", "RZMC_mean"]
    for arm in ["MODIS_T", "MODIS_TA", "T_to_J2", "TA_to_J1J2"]:
        md_out += [f"\n## C0: {arm} - OL\n", table("C0", f"{arm}-OL", ks)]
    for c, (A, B, sy) in cdef.items():
        md_out += [f"\n## {c}: {A} - {B}\n", table(c, f"{c}:{A}-{B}", ks)]
    md_out += ["\n## C4 difference of differences\n", table("C4dd", "C4dd:(T_to_J2-TA_to_J1J2)[SY4-6] - mean[SY1-3]", ks)]
    md_out += ["\n## Decision rule (proposed: |A-B| >= 10% of mean |DA-OL| in each year, same sign in all years)\n",
               mdtable(dec.round(3))]
    miss = {n: [m for m in arms[n]["missing"]] for n in names}
    md_out += ["\n## Input gaps recorded by wb_extract.py\n", "```", json.dumps(miss), "```"]
    open(f"{a.figdir}/D_water_balance_tables.md", "w").write("\n".join(md_out))
    print("done")


if __name__ == "__main__":
    main()
