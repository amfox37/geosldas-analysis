#!/usr/bin/env python3
"""Analysis D additions (step 2): partition fractions, soil-moisture response, significance maps, SY1 and gap-day
sensitivities. Input: output/water_balance/wb_<arm>.nc (wb_extract.py) and wb2_gap_<arm>.nc (wb2_extract.py).
Output: docs/viirs_transition_report_figures/D2_*. Reuses wb_analyze.py (budget definitions, mask, block bootstrap).
Partition fraction = (domain-mean term difference) / (domain-mean I_snow difference), ratio of sums, resampled by 5-degree
blocks (the numerator and denominator resampled together); it is only meaningful where the I_snow difference is clearly
non-zero (flagged 'unstable' when the 95% CI of the I_snow difference includes 0).
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wb_analyze as W  # noqa: E402
from wb_analyze import PROJ, CATCHMENT_DEF, SY_IDS, REPL, block_codes, confidence_interval, load_arm, series, diff, wmean, mdtable  # noqa: E402

FIG = f"{PROJ}/docs/viirs_transition_report_figures"
S1, S3 = ["SY1", "SY2", "SY3"], ["SY4", "SY5", "SY6"]
TERMS = ["dET", "dRunoff", "dStorage", "dFSW", "residual"]
SEAS = {"SON": [1, 2, 3], "DJF": [4, 5, 6], "MAM": [7, 8, 9], "JJA": [10, 11, 0], "AMJ": [8, 9, 10]}  # month_idx: Aug=0
rng = np.random.default_rng(7)


def build(arms):
    ser = {"OL": series(arms, [("OL", SY_IDS)]), "T_to_J2": series(arms, [("T_to_J2", SY_IDS)]),
           "TA_to_J1J2": series(arms, [("TA_to_J1J2", SY_IDS)]),
           "T_control": series(arms, [("T_control", S3)]), "TA_control": series(arms, [("TA_control", S3)]),
           "MODIS_T": series(arms, [("T_to_J2", S1), ("T_control", S3)]),
           "MODIS_TA": series(arms, [("TA_to_J1J2", S1), ("TA_control", S3)])}
    c0 = [("T_to_J2", "OL", SY_IDS), ("TA_to_J1J2", "OL", SY_IDS), ("MODIS_T", "OL", SY_IDS), ("MODIS_TA", "OL", SY_IDS),
          ("T_control", "OL", S3), ("TA_control", "OL", S3)]
    cc = [("C1", "T_control", "T_to_J2", S3), ("C2", "TA_control", "TA_to_J1J2", S3),
          ("C3", "MODIS_T", "MODIS_TA", SY_IDS), ("C4", "T_to_J2", "TA_to_J1J2", S3)]
    items = [(f"C0:{a}-{b}", a, b, sy) for a, b, sy in c0] + [(f"{c}:{a}-{b}", a, b, sy) for c, a, b, sy in cc]
    return ser, items


def boot_ratio(num, den, w, keep, codes):
    v = keep & np.isfinite(num) & np.isfinite(den)
    nb = codes.max() + 1
    Sn = np.bincount(codes[v], weights=(num * w)[v], minlength=nb); Sd = np.bincount(codes[v], weights=(den * w)[v], minlength=nb)
    idx = rng.integers(0, nb, size=(REPL, nb))
    r = Sn[idx].sum(1) / Sd[idx].sum(1)
    d = Sd[idx].sum(1) / np.bincount(codes[v], weights=w[v], minlength=nb)[idx].sum(1)
    lo, hi = np.percentile(r, [2.5, 97.5])
    dlo, dhi = np.percentile(d, [2.5, 97.5])
    return Sn.sum() / Sd.sum(), lo, hi, (dlo <= 0 <= dhi)


def boot_mean(x, w, keep, codes):
    v = keep & np.isfinite(x)
    nb = codes.max() + 1
    S = np.bincount(codes[v], weights=(x * w)[v], minlength=nb); Wt = np.bincount(codes[v], weights=w[v], minlength=nb)
    idx = rng.integers(0, nb, size=(REPL, nb))
    m = S[idx].sum(1) / Wt[idx].sum(1)
    return np.percentile(m, [2.5, 97.5])


def pooled_diff(ser, a, b, sy, k):
    return np.nanmean([diff(ser[a], ser[b], s)[k] for s in sy], 0)


def main():
    cfg = json.load(open(f"{PROJ}/config/arms.json"))
    d0 = f"{PROJ}/output/water_balance"
    arms = {n: load_arm(n, d0) for n in cfg["arms"]}
    lat, lon = arms["OL"]["lat"], arms["OL"]["lon"]
    area = np.loadtxt(CATCHMENT_DEF, usecols=(6,), skiprows=1)
    ol = arms["OL"]
    scf = ol["m_SCF"]
    mask = ((lat > 20) & ((np.nanmax(scf.reshape(-1, lat.size), 0) > 0.05) | (np.nanmax(ol["SWE_max"], 0) > 5))
            & (np.nanmean(scf[:, [10, 11, 0], :].reshape(-1, lat.size), 0) < 0.2) & np.isfinite(lat))
    print("mask tiles", int(mask.sum()), flush=True)
    codes = block_codes(pd.DataFrame({"lat": lat, "lon": lon}), 5.0)
    ser, items = build(arms)
    out = []

    # ---------------- (1) partition fractions and (2) soil-moisture response, per contrast and SY + pooled
    prow, srow = [], []
    for label, a, b, sy in items:
        ds = {s: diff(ser[a], ser[b], s) for s in sy}
        ds["pooled"] = {k: np.nanmean([ds[s][k] for s in sy], 0) for k in ds[sy[0]] if k != "monthly"}
        ds["pooled"]["monthly"] = {k: np.nanmean([ds[s]["monthly"][k] for s in sy], 0) for k in ds[sy[0]]["monthly"]}
        for sid, d in ds.items():
            for t in TERMS:
                r, lo, hi, unst = boot_ratio(d[t], d["I_snow"], area, mask, codes)
                prow.append(dict(contrast=label, sy=sid, term=t, fraction=r, lo=lo, hi=hi, unstable=unst,
                                 I_snow_mean=wmean(d["I_snow"], area, mask)))
            for var in ["RZMC", "SFMC"]:
                for sname, mi in list(SEAS.items()) + [("year", list(range(12)))]:
                    x = np.nanmean(d["monthly"][var][mi], 0) if var in d["monthly"] else None
                    if x is None:
                        continue
                    lo, hi = boot_mean(x, area, mask, codes)
                    srow.append(dict(contrast=label, sy=sid, var=var, season=sname, mean=wmean(x, area, mask), lo=lo, hi=hi))
    # SFMC is not in diff()['monthly']; compute directly from the monthly means
    srow = [r for r in srow if r["var"] == "RZMC"]
    for label, a, b, sy in items:
        for sid in sy + ["pooled"]:
            use = sy if sid == "pooled" else [sid]
            for sname, mi in list(SEAS.items()) + [("year", list(range(12)))]:
                x = np.nanmean([np.nanmean(ser[a][s]["m_SFMC"][mi] - ser[b][s]["m_SFMC"][mi], 0) for s in use], 0)
                lo, hi = boot_mean(x, area, mask, codes)
                srow.append(dict(contrast=label, sy=sid, var="SFMC", season=sname, mean=wmean(x, area, mask), lo=lo, hi=hi))
    pdf, sdf = pd.DataFrame(prow), pd.DataFrame(srow)
    pdf.to_csv(f"{FIG}/D2_partition_fractions.csv", index=False)
    sdf.to_csv(f"{FIG}/D2_soil_moisture_response.csv", index=False)

    # figure: partition fractions, C0 pooled (SY1-6 per arm) + per SY for MODIS_TA
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    ax = axs[0]
    labs = [l for l, *_ in items if l.startswith("C0:") and ("MODIS" in l or "to_" in l)]
    for j, t in enumerate(TERMS):
        for i, l in enumerate(labs):
            s = pdf[(pdf.contrast == l) & (pdf.sy == "pooled") & (pdf.term == t)].iloc[0]
            ax.errorbar(i + (j - 2) * 0.15, s.fraction, yerr=[[s.fraction - s.lo], [s.hi - s.fraction]], fmt="o", ms=3,
                        color=f"C{j}", label=t if i == 0 else None)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels([l.split(":")[1] for l in labs], rotation=25, fontsize=7)
    ax.axhline(0, c="k", lw=.5); ax.set_ylabel("fraction of I_snow"); ax.legend(fontsize=7); ax.set_title("C0 partition (pooled over years; block-bootstrap CI)")
    ax = axs[1]
    for j, t in enumerate(TERMS):
        s = pdf[(pdf.contrast == "C0:MODIS_TA-OL") & (pdf.sy != "pooled") & (pdf.term == t)]
        ax.plot(range(len(s)), s.fraction, "o-", color=f"C{j}", label=t)
    ax.set_xticks(range(6)); ax.set_xticklabels(SY_IDS); ax.axhline(0, c="k", lw=.5); ax.set_title("MODIS_TA - OL partition by snow year")
    fig.tight_layout(); fig.savefig(f"{FIG}/D2_1_partition.png", dpi=130); plt.close(fig)

    # figure: soil moisture response
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for ax, var in zip(axs, ["RZMC", "SFMC"]):
        labs2 = [l for l, *_ in items if l.split(":")[0] in ("C0", "C1", "C2", "C3", "C4") and ("MODIS" not in l or l.startswith(("C0", "C3")))]
        labs2 = [l for l in labs2 if not (l.startswith("C0") and l.split(":")[1].startswith(("T_control", "TA_control")))]
        for i, l in enumerate(labs2):
            for j, sn in enumerate(["SON", "DJF", "MAM", "JJA"]):
                s = sdf[(sdf.contrast == l) & (sdf.sy == "pooled") & (sdf["var"] == var) & (sdf.season == sn)].iloc[0]
                ax.errorbar(i + (j - 1.5) * 0.15, s["mean"], yerr=[[s["mean"] - s.lo], [s.hi - s["mean"]]], fmt="o", ms=3,
                            color=f"C{j}", label=sn if i == 0 else None)
        ax.set_xticks(range(len(labs2))); ax.set_xticklabels([l.replace("C0:", "DA-OL:") for l in labs2], rotation=40, ha="right", fontsize=6)
        ax.axhline(0, c="k", lw=.5); ax.set_title(f"{var} difference (m3 m-3), seasonal means, pooled years"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/D2_2_soil_moisture.png", dpi=130); plt.close(fig)

    # ---------------- (3) significance maps for C3 and C2
    def sigmap(label, a, b, sy, var, test):
        per = np.array([diff(ser[a], ser[b], s)[var] for s in sy])  # (n, tile)
        m = np.nanmean(per, 0)
        same = (np.all(per > 0, 0) | np.all(per < 0, 0))
        sig = same.copy()
        if test:
            t, p = stats.ttest_1samp(per, 0, axis=0, nan_policy="omit")
            p = np.where(np.isfinite(p), p, 1.0)
            pm = p[mask]; o = np.argsort(pm); q = np.empty_like(pm); n = len(pm)
            q[o] = np.minimum.accumulate((pm[o] * n / (np.arange(n) + 1))[::-1])[::-1]
            pbh = np.ones_like(p); pbh[mask] = q
            sig = same & (pbh < 0.05)
        return m, sig

    srows = []
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    for r, (c, a, b, sy, test) in enumerate([("C3", "MODIS_T", "MODIS_TA", SY_IDS, True), ("C2", "TA_control", "TA_to_J1J2", S3, False)]):
        for col, var in enumerate(["I_snow", "dRunoff", "SWE_end"]):
            m, sig = sigmap(c, a, b, sy, var, test)
            v = np.where(mask, m, np.nan); lim = np.nanpercentile(np.abs(v), 98) or 1
            ax = axs[r, col]
            sc = ax.scatter(lon, lat, c=v, s=1.5, cmap="RdBu_r", vmin=-lim, vmax=lim)
            k = mask & sig
            ax.scatter(lon[k], lat[k], s=0.4, c="k", alpha=.35)
            ax.set_xlim(-180, 180); ax.set_ylim(20, 85)
            crit = "same sign in all %d years + BH-FDR t-test p<0.05" % len(sy) if test else "same sign in all %d years (no t-test, n=3)" % len(sy)
            ax.set_title(f"{c} {a}-{b}: {var}\nblack = {crit}", fontsize=8); fig.colorbar(sc, ax=ax, shrink=.8)
            srows.append(dict(contrast=c, var=var, frac_tiles_marked=float(k.sum() / mask.sum()),
                              frac_same_sign=float((mask & (np.array([diff(ser[a], ser[b], s)[var] for s in sy]) > 0).all(0)).sum() / mask.sum()
                                                   + (mask & (np.array([diff(ser[a], ser[b], s)[var] for s in sy]) < 0).all(0)).sum() / mask.sum())))
    fig.tight_layout(); fig.savefig(f"{FIG}/D2_3_significance_maps_C3_C2.png", dpi=110); plt.close(fig)
    sg = pd.DataFrame(srows); sg.to_csv(f"{FIG}/D2_significance_summary.csv", index=False)

    # ---------------- (4a) SY1 sensitivity and spin-up check
    sens = []
    for label, a, b, sy in items:
        if not label.startswith(("C0:MODIS", "C0:T_to", "C0:TA_to", "C3")):
            continue
        for tag, use in [("all", sy), ("excl_SY1", [s for s in sy if s != "SY1"])]:
            for t in ["I_snow", "dET", "dRunoff", "dStorage", "dFSW", "residual", "RZMC_mean"]:
                x = np.nanmean([diff(ser[a], ser[b], s)[t] for s in use], 0)
                lo, hi = boot_mean(x, area, mask, codes)
                sens.append(dict(contrast=label, years=tag, var=t, mean=wmean(x, area, mask), lo=lo, hi=hi))
    sens_df = pd.DataFrame(sens)
    # spin-up check: DA-OL TWLAND difference at 1 Aug endpoints and monthly trajectory
    traj, endp = {}, []
    for arm in ["MODIS_TA", "MODIS_T", "TA_to_J1J2", "T_to_J2"]:
        ts = []
        for s in SY_IDS:
            if s not in ser[arm]:
                continue
            ts.append([wmean(ser[arm][s]["m_TWLAND"][m] - ser["OL"][s]["m_TWLAND"][m], area, mask) for m in range(12)])
        traj[arm] = np.array(ts)
        sids = [s for s in SY_IDS if s in ser[arm]]
        for s in sids:
            endp.append(dict(arm=arm, sy=s, TW_diff_1Aug_start=wmean(ser[arm][s]["TW0"] - ser["OL"][s]["TW0"], area, mask),
                             TW_diff_1Aug_end=wmean(ser[arm][s]["TW1"] - ser["OL"][s]["TW1"], area, mask),
                             dStorage=wmean(diff(ser[arm], ser["OL"], s)["dStorage"], area, mask)))
    endp_df = pd.DataFrame(endp)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for arm, c in [("MODIS_TA", "C0"), ("MODIS_T", "C1")]:
        axs[0].plot(np.arange(traj[arm].size), traj[arm].ravel(), c=c, label=f"{arm} - OL")
    axs[0].set_xticks(np.arange(0, 72, 12)); axs[0].set_xticklabels(["Aug " + str(2020 + i) for i in range(6)], fontsize=7)
    axs[0].set_ylabel("TWLAND difference (kg m-2, monthly mean)"); axs[0].legend(); axs[0].set_title("Accumulated storage difference vs OL")
    for arm in ["MODIS_TA", "MODIS_T"]:
        s = endp_df[endp_df.arm == arm]; axs[1].plot(range(len(s)), s.dStorage, "o-", label=f"{arm} dStorage")
    axs[1].set_xticks(range(6)); axs[1].set_xticklabels(SY_IDS); axs[1].axhline(0, c="k", lw=.5); axs[1].legend(); axs[1].set_title("dStorage per snow year (DA-OL)")
    fig.tight_layout(); fig.savefig(f"{FIG}/D2_4_SY1_storage_spinup.png", dpi=130); plt.close(fig)

    # ---------------- (4b) gap-day sensitivity: remove gap/boundary-day fluxes from every arm (flux terms only)
    gapn = {}
    for n in arms:
        with Dataset(f"{d0}/wb2_gap_{n}.nc") as g:
            gapn[n] = {k: np.ma.filled(g[k][:].astype(np.float64), 0.0) for k in ["P", "ET", "RUNOFF", "SNOWMELT", "FSW", "I_snow"]}
            gapn[n]["n"] = list(g["n_gap_days"][:])
    base_vals = {}
    def flux_table(tag):
        res = []
        ser2, items2 = build(arms)
        for label, a, b, sy in items2:
            for t in ["I_snow", "dET", "dRunoff"]:
                x = np.nanmean([diff(ser2[a], ser2[b], s)[t] for s in sy], 0)
                lo, hi = boot_mean(x, area, mask, codes)
                res.append(dict(contrast=label, case=tag, var=t, mean=wmean(x, area, mask), lo=lo, hi=hi))
        return res
    gres = flux_table("original")
    for n, a in arms.items():  # in-place correction, all arms, own snow years
        for i, sid in enumerate(a["sy"]):
            a["I_snow"][i] -= gapn[n]["I_snow"][i]; a["ET"][i] -= gapn[n]["ET"][i]
            a["RUNSURF"][i] -= gapn[n]["RUNOFF"][i]  # RUNOFF = RUNSURF + BASEFLOW combined, removed from RUNSURF
    gres += flux_table("gap_days_removed")
    gdf = pd.DataFrame(gres)
    gw = gdf.pivot_table(index=["contrast", "var"], columns="case", values="mean").reset_index()
    gw["change"] = gw["gap_days_removed"] - gw["original"]
    gw.to_csv(f"{FIG}/D2_gapday_sensitivity.csv", index=False)
    sens_df.to_csv(f"{FIG}/D2_SY1_sensitivity.csv", index=False); endp_df.to_csv(f"{FIG}/D2_storage_endpoints.csv", index=False)

    # ---------------- markdown
    f = lambda x: f"{x:.2f}"
    md = ["# Analysis D additions (auto-generated by wb2_analyze.py)\n",
          f"Mask tiles: {int(mask.sum())}. Fractions are domain-mean term difference / domain-mean I_snow difference.\n",
          "## Partition of added snow water (C0, pooled over years)\n"]
    pt = pdf[(pdf.sy == "pooled") & pdf.contrast.str.startswith("C0")].copy()
    pt["txt"] = pt.apply(lambda r: f"{r.fraction:.2f} [{r.lo:.2f}, {r.hi:.2f}]", axis=1)
    md.append(mdtable(pt.pivot(index="contrast", columns="term", values="txt").reset_index()[["contrast"] + TERMS]))
    md.append("\n## Partition for contrasts C1-C4 (pooled; 'unstable' where the I_snow difference CI includes 0)\n")
    pt = pdf[(pdf.sy == "pooled") & ~pdf.contrast.str.startswith("C0")].copy()
    pt["txt"] = pt.apply(lambda r: f"{r.fraction:.2f} [{r.lo:.2f}, {r.hi:.2f}]" + (" unstable" if r.unstable else ""), axis=1)
    md.append(mdtable(pt.pivot(index="contrast", columns="term", values="txt").reset_index()[["contrast"] + TERMS]))
    md.append("\n## Root-zone soil moisture difference (m3 m-3), pooled, by season\n")
    st = sdf[(sdf.sy == "pooled") & (sdf["var"] == "RZMC")].copy()
    st["txt"] = st.apply(lambda r: f"{r['mean']:.4f} [{r.lo:.4f}, {r.hi:.4f}]", axis=1)
    md.append(mdtable(st.pivot(index="contrast", columns="season", values="txt").reset_index()[["contrast", "SON", "DJF", "MAM", "AMJ", "JJA", "year"]]))
    md.append("\n## Significance maps: fraction of mask tiles marked\n"); md.append(mdtable(sg.round(3)))
    md.append("\n## SY1 sensitivity (domain mean, kg m-2 per snow year; RZMC m3 m-3)\n")
    sp = sens_df.copy(); sp["txt"] = sp.apply(lambda r: f"{r['mean']:.3f}", axis=1)
    md.append(mdtable(sp.pivot_table(index=["contrast", "var"], columns="years", values="mean", aggfunc="first").reset_index().round(3)))
    md.append("\n## Storage difference vs OL at the 1 Aug endpoints (spin-up check)\n"); md.append(mdtable(endp_df.round(2)))
    md.append("\n## Gap-day sensitivity (flux terms; gap + boundary days removed from every arm)\n"); md.append(mdtable(gw.round(3)))
    open(f"{FIG}/D2_water_balance_additions.md", "w").write("\n".join(md))
    print("done")


if __name__ == "__main__":
    main()
