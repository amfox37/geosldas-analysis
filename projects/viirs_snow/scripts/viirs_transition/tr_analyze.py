#!/usr/bin/env python3
"""Analysis E step 2: transition behaviour after 2023-08-01 (docs/viirs_transition_plan.md, section E).

Inputs: output/transition/tr_<year>.nc (tr_extract.py), output/incr_stats/{incr_stats,innov}_<arm>.nc (analysis B),
        output/water_balance/wb_<arm>.nc (analysis D, monthly SWE).
Outputs: docs/viirs_transition_report_figures/E_summary.md, E_weekly_tables.csv (git-ignored), E_1..E_4 *.png.
Weeks are 7-day bins counted from 1 Aug. NH = bands >= 30N (bands N_30_45, N_45_60, N_60_90) as in analysis B.
python3 -I tr_analyze.py
"""
import csv, datetime as dt, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
FIG = PROJECT / "docs" / "viirs_transition_report_figures"
OUT = PROJECT / "output"
COL = {"OL": "#898781", "T_to_J2": "#2a78d6", "TA_to_J1J2": "#eb6834", "T_control": "#1baf7a", "TA_control": "#6250d6"}
NHB = [2, 3, 4]
NB, NL = 5, 8
WEND = 20   # weeks from 1 Aug to ~15 Dec


def dd(x):
    return np.array([dt.date(int(v) // 10000, int(v) // 100 % 100, int(v) % 100) for v in x])


def week_of(days, year):
    return np.array([(d - dt.date(year, 8, 1)).days // 7 for d in days])


def wsum(vals, wk, n=WEND):
    return np.array([np.nansum(vals[wk == w]) if (wk == w).any() else np.nan for w in range(n)])


def wmean(vals, wk, n=WEND):
    return np.array([np.nanmean(vals[wk == w]) if (wk == w).any() else np.nan for w in range(n)])


def grp_nh(a):
    a = a.reshape(a.shape[:-1] + (NB, NL))
    return a[..., NHB, :].sum((-1, -2))


def main():
    cfg = json.load(open(PROJECT / "config" / "arms.json"))
    rows = []
    md = ["# Analysis E: transition behaviour after 2023-08-01", "",
          "Weeks are 7-day bins from 1 Aug (week 0 = 1-7 Aug; week 8 starts 26 Sep). NH = bands >= 30N. One member.", ""]
    # ---------- A. when does snow exist? (domain means, 2023 window) ----------
    T = {y: netCDF4.Dataset(OUT / "transition" / f"tr_{y}.nc") for y in (2023, 2024, 2025)}
    arms = json.loads(T[2023].arms); pairs = [tuple(p) for p in json.loads(T[2023].pairs)]
    lat = T[2023]["lat"][:]
    bidx = np.digitize(lat, [-90, 0, 30, 45, 60, 90.01]) - 1
    ntile_b = np.array([(bidx == b).sum() for b in range(5)], float)
    nh_w = ntile_b[NHB] / ntile_b[NHB].sum()
    days23 = dd(T[2023]["date"][:]); wk23 = week_of(days23, 2023)
    swe = T[2023]["swe"][:]; nsn = T[2023]["nsnow"][:]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    x = np.array([(d - dt.date(2023, 8, 1)).days for d in days23])
    for j, a in enumerate(arms):
        ax[0].plot(x, (swe[:, j, NHB] * nh_w).sum(1), color=COL[a], label=a)
        ax[1].plot(x, nsn[:, j, NHB].sum(1), color=COL[a])
        ax[2].plot(x, nsn[:, j, 0], color=COL[a])
    for k, t in enumerate(["NH (>=30N) mean SWE, kg m-2", "NH tiles with SWE > 1 kg m-2", "Southern Hemisphere tiles with SWE > 1 kg m-2"]):
        ax[k].set_title(t, fontsize=9); ax[k].axvline(0, color="k", lw=0.6); ax[k].axvline(56, color="k", lw=0.6, ls=":")
        ax[k].set_xlabel("days since 2023-08-01")
    ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(FIG / "E_1_snow_state_2023.png", dpi=130); plt.close(fig)
    # monthly share of seasonal SWE from analysis D (NH tiles)
    md += ["## 1. How much snow is there in the first 4-8 weeks?", "",
           "Share of the snow-year integral of NH monthly-mean SWE (analysis D) that falls in Aug+Sep, and of absolute snow increments (analysis B) in the first 8 weeks (1 Aug - 25 Sep), NH. Snow years SY4-SY6.", "",
           "| arm | SY | Aug+Sep share of SWE integral | first-8-week share of NH abs increment | first-8-week NH abs increment (kg m-2 summed over tiles) |", "|---|---|---|---|---|"]
    SYR = cfg["snow_years"]
    for a in ["OL", "T_to_J2", "TA_to_J1J2", "T_control", "TA_control"]:
        w = netCDF4.Dataset(OUT / "water_balance" / f"wb_{a}.nc"); w.set_auto_mask(False)
        wl = w["lat"][:]; nhm = wl >= 30; sids = list(w["sy_id"][:])
        ms = w["m_SWE"][:]
        sh = {}
        for i, s in enumerate(sids):
            if s not in ("SY4", "SY5", "SY6"):
                continue
            tot = np.nansum(np.nanmean(ms[i][:, nhm], 1)); sh[s] = np.nansum(np.nanmean(ms[i][:2][:, nhm], 1)) / tot
        ab = {}
        if a != "OL":
            S = netCDF4.Dataset(OUT / "incr_stats" / f"incr_stats_{a}.nc"); S.set_auto_mask(False)
            dts = dd(S["date"][:]); sa = grp_nh(S["step_abs"][:]).sum(1)
            for s in SYR:
                if s["id"] not in ("SY4", "SY5", "SY6"):
                    continue
                b0 = dt.date.fromisoformat(s["beg"]); m = np.array([(d >= b0) and (d <= dt.date.fromisoformat(s["end"])) for d in dts])
                m8 = np.array([(d >= b0) and ((d - b0).days < 56) for d in dts])
                ab[s["id"]] = (sa[m8].sum() / sa[m].sum(), sa[m8].sum())
        for s in ("SY4", "SY5", "SY6"):
            if s in sh:
                md.append(f"| {a} | {s} | {sh[s]*100:.2f}% | " + (f"{ab[s][0]*100:.2f}% | {ab[s][1]:.3g} |" if s in ab else "- | - |"))
    # ---------- B. divergence of C1/C2 from the common start state ----------
    def nh_rms(rmsb, nb=ntile_b):
        return np.sqrt((rmsb[..., NHB] ** 2 * ntile_b[NHB]).sum(-1) / ntile_b[NHB].sum())
    pidx = {p: i for i, p in enumerate(pairs)}
    ref = {"C1": ("T_control", "OL"), "C2": ("TA_control", "OL")}
    c = {"C1": ("T_control", "T_to_J2"), "C2": ("TA_control", "TA_to_J1J2")}
    rms = T[2023]["dswe_rms"][:]; scf = T[2023]["dscf_rms"][:]
    md += ["", "## 2. Divergence of the MODIS control from the VIIRS arm after the shared 2023-08-01 restart (SWE, kg m-2, rms over tiles)", "",
           "Both arms start from the identical state, so the step at t = 0 is zero by construction and the question is how fast the difference grows and when it is no longer limited by the absence of snow. Ratio = rms(control - VIIRS arm) / rms(control - OL), the reference DA effect.", "",
           "| week (from 1 Aug) | C1 NH rms | C1 NH ratio | C2 NH rms | C2 NH ratio | C1 SH rms | C2 SH rms | NH tiles with SWE>1 (T_to_J2) |", "|---|---|---|---|---|---|---|---|"]
    tabs = {}
    for k, p in c.items():
        r = nh_rms(rms[:, pidx[p], :]); rr = nh_rms(rms[:, pidx[ref[k]], :])
        tabs[k] = (wmean(r, wk23), wmean(rr, wk23), wmean(rms[:, pidx[p], 0], wk23), wmean(rms[:, pidx[ref[k]], 0], wk23))
    nsw = wmean(nsn[:, arms.index("T_to_J2"), NHB].sum(1), wk23)
    for w in range(WEND):
        md.append(f"| {w} | {tabs['C1'][0][w]:.3f} | {tabs['C1'][0][w]/tabs['C1'][1][w]:.2f} | {tabs['C2'][0][w]:.3f} | {tabs['C2'][0][w]/tabs['C2'][1][w]:.2f} | {tabs['C1'][2][w]:.3f} | {tabs['C2'][2][w]:.3f} | {nsw[w]:.0f} |")
        rows.append(["divergence", w, tabs["C1"][0][w], tabs["C1"][1][w], tabs["C2"][0][w], tabs["C2"][1][w], tabs["C1"][2][w], tabs["C2"][2][w]])
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for k, col in (("C1", "#2a78d6"), ("C2", "#eb6834")):
        ax[0].plot(range(WEND), tabs[k][0], color=col, label=f"{k}: control - VIIRS arm"); ax[0].plot(range(WEND), tabs[k][1], color=col, ls="--", label=f"{k}: control - OL (DA effect)")
        ax[1].plot(range(WEND), tabs[k][0] / tabs[k][1], color=col, label=k)
        ax[2].plot(range(WEND), tabs[k][2], color=col); ax[2].plot(range(WEND), tabs[k][3], color=col, ls="--")
    ax[0].set_title("NH rms SWE difference, 2023", fontsize=9); ax[1].set_title("ratio (contrast / DA effect), NH", fontsize=9)
    ax[2].set_title("SH rms SWE difference (dashed: DA effect)", fontsize=9)
    for a_ in ax:
        a_.set_xlabel("week since 2023-08-01"); a_.axvline(8, color="k", lw=0.6, ls=":")
    ax[0].legend(fontsize=6); ax[1].set_ylim(0, 3); ax[1].legend(fontsize=7); fig.tight_layout()
    fig.savefig(FIG / "E_2_divergence_2023.png", dpi=130); plt.close(fig)
    # SCF too
    md += ["", "SCF rms difference (NH), weeks 0-19: C1 " + ", ".join(f"{v:.3f}" for v in wmean(nh_rms(scf[:, pidx[c['C1']], :]), wk23)) +
           "; C2 " + ", ".join(f"{v:.3f}" for v in wmean(nh_rms(scf[:, pidx[c['C2']], :]), wk23)), ""]
    # same calendar window in 2024, 2025 (arms have diverged by then: shows what the differences look like later)
    md += ["Same calendar weeks in 2024 and 2025 (no longer a shared start state; the difference then includes everything accumulated in the previous year). NH rms SWE difference C1 / C2:", "",
           "| week | C1 2023 | C1 2024 | C1 2025 | C2 2023 | C2 2024 | C2 2025 |", "|---|---|---|---|---|---|---|"]
    ser = {}
    for y in (2023, 2024, 2025):
        dsy = dd(T[y]["date"][:]); wk = week_of(dsy, y); r = T[y]["dswe_rms"][:]
        for k, p in c.items():
            ser[(k, y)] = wmean(nh_rms(r[:, pidx[p], :]), wk)
    for w in range(WEND):
        md.append(f"| {w} | " + " | ".join(f"{ser[(k,y)][w]:.3f}" for k in ("C1", "C2") for y in (2023, 2024, 2025)) + " |")
    # ---------- C. increments around the switch ----------
    md += ["", "## 3. Increments: step at the switch relative to the same calendar weeks in other years", "",
           "NH sum of absolute snow increments per week (kg m-2 summed over NH tiles, millions), first 8 weeks (1 Aug - 25 Sep), per snow year. 2020-2022 are the pre-switch (MODIS) years of the T_to_J2 / TA_to_J1J2 arms; 2023-2025 are the post-switch years (controls stay MODIS).", "",
           "| arm | " + " | ".join(str(y) for y in range(2020, 2026)) + " | 2023 / mean(2024,2025) | 2023 / mean(2020-22) |", "|---|---|---|---|---|---|---|---|---|"]
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6), sharey=False)
    inc = {}
    for ai, a in enumerate(["T_to_J2", "TA_to_J1J2", "T_control", "TA_control"]):
        S = netCDF4.Dataset(OUT / "incr_stats" / f"incr_stats_{a}.nc"); S.set_auto_mask(False)
        dts = dd(S["date"][:]); sa = grp_nh(S["step_abs"][:]).sum(1); sn = grp_nh(S["step_n"][:]).sum(1)
        vals = {}
        for y in range(2020, 2026):
            m = np.array([(d.year == y) and (dt.date(y, 8, 1) <= d <= dt.date(y, 12, 15)) and not (y == 2023 and d == dt.date(2023, 8, 1)) for d in dts])
            if m.sum() < 100:
                continue
            wkk = week_of(dts[m], y); w_abs = wsum(sa[m], wkk); w_n = wsum(sn[m], wkk)
            vals[y] = (w_abs, w_n); axs[ai].plot(range(WEND), w_abs / 1e6, color=("#d03a2f" if y == 2023 else None), lw=2 if y == 2023 else 1, label=str(y))
        axs[ai].set_title(a, fontsize=9); axs[ai].axvline(8, color="k", lw=0.6, ls=":"); axs[ai].set_xlabel("week since 1 Aug")
        inc[a] = vals
        e8 = {y: np.nansum(v[0][:8]) / 1e6 for y, v in vals.items()}
        late = [e8[y] for y in (2024, 2025) if y in e8]; early = [e8[y] for y in (2020, 2021, 2022) if y in e8]
        md.append(f"| {a} | " + " | ".join(f"{e8[y]:.2f}" if y in e8 else "-" for y in range(2020, 2026)) +
                  f" | {e8[2023]/np.mean(late):.2f} | " + (f"{e8[2023]/np.mean(early):.2f} |" if early else "- |"))
        for y, v in vals.items():
            for w in range(WEND):
                rows.append(["incr_abs_nh", a, y, w, v[0][w], v[1][w]])
    axs[0].legend(fontsize=6, ncol=2); axs[0].set_ylabel("NH sum |incr| per week, 1e6 kg m-2")
    fig.tight_layout(); fig.savefig(FIG / "E_3_increments_by_week.png", dpi=130); plt.close(fig)
    # ---------- C2. same-state ratio VIIRS arm / MODIS control, weekly ----------
    md += ["", "Same-state ratio of NH weekly absolute increments, VIIRS arm / MODIS control (T_to_J2 / T_control; TA_to_J1J2 / TA_control), by week from 1 Aug 2023, and for the same weeks in 2024 and 2025 (state no longer shared):", "",
           "| week | T 2023 | T 2024 | T 2025 | TA 2023 | TA 2024 | TA 2025 |", "|---|---|---|---|---|---|---|"]
    rat = {}
    for (v_, c_, tag) in (("T_to_J2", "T_control", "T"), ("TA_to_J1J2", "TA_control", "TA")):
        for y in (2023, 2024, 2025):
            rat[(tag, y)] = inc[v_][y][0] / inc[c_][y][0]
    for w in range(WEND):
        md.append(f"| {w} | " + " | ".join(f"{rat[(t, y)][w]:.2f}" for t in ("T", "TA") for y in (2023, 2024, 2025)) + " |")
        rows.append(["incr_ratio", w] + [rat[(t, y)][w] for t in ("T", "TA") for y in (2023, 2024, 2025)])
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for (tag, y), v in rat.items():
        ax.plot(range(WEND), v, color=("#2a78d6" if tag == "T" else "#eb6834"), ls={2023: "-", 2024: "--", 2025: ":"}[y], label=f"{tag} {y}")
    ax.axvline(8, color="k", lw=0.6, ls=":"); ax.set_xlabel("week since 1 Aug"); ax.set_title("NH |incr|: VIIRS arm / MODIS control", fontsize=9)
    ax.legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig(FIG / "E_5_increment_ratio.png", dpi=130); plt.close(fig)
    # ---------- D. new sensor spin-up: obs count and O-F of J1/J2 ----------
    md += ["", "## 4. New sensor: J1 / J2 obs counts and O-F after the switch (TA_to_J1J2 arm, NH tiles)", "",
           "Weekly mean daily tile count and rms(obs - forecast) in SCF units. If the new sensors needed a spin-up, 2023 would differ from 2024-25 in the same weeks.", "",
           "| species | week | n/day 2023 | n/day 2024 | n/day 2025 | rms O-F 2023 | rms O-F 2024 | rms O-F 2025 |", "|---|---|---|---|---|---|---|---|"]
    I = netCDF4.Dataset(OUT / "incr_stats" / "innov_TA_to_J1J2.nc"); I.set_auto_mask(False)
    idts = dd(I["date"][:]); st = I["stats"][:]   # (day, species, cls, band, var) species order Aqua, Terra, J1, J2
    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    for si, (sp, ix) in enumerate([("J1", 2), ("J2", 3)]):
        n = {}; r = {}
        for y in (2023, 2024, 2025):
            m = np.array([(d.year == y) and (dt.date(y, 8, 1) <= d <= dt.date(y, 12, 15)) and d != dt.date(2023, 8, 1) for d in idts])
            wkk = week_of(idts[m], y)
            nn = st[m][:, ix, 0][:, NHB, 0].sum(1); s2 = st[m][:, ix, 0][:, NHB, 2].sum(1)
            nd_ = np.array([ (wkk == w).sum() for w in range(WEND)])
            n[y] = wsum(nn, wkk) / np.maximum(nd_, 1); r[y] = np.sqrt(wsum(s2, wkk) / np.maximum(wsum(nn, wkk), 1))
            axs[si, 0].plot(range(WEND), n[y], label=str(y), lw=2 if y == 2023 else 1, color=("#d03a2f" if y == 2023 else None))
            axs[si, 1].plot(range(WEND), r[y], label=str(y), lw=2 if y == 2023 else 1, color=("#d03a2f" if y == 2023 else None))
        axs[si, 0].set_title(f"{sp}: NH obs per day", fontsize=9); axs[si, 1].set_title(f"{sp}: NH rms O-F (SCF)", fontsize=9)
        for w in range(WEND):
            md.append(f"| {sp} | {w} | {n[2023][w]:.0f} | {n[2024][w]:.0f} | {n[2025][w]:.0f} | {r[2023][w]:.3f} | {r[2024][w]:.3f} | {r[2025][w]:.3f} |")
            rows.append(["innov", sp, w, n[2023][w], n[2024][w], n[2025][w], r[2023][w], r[2024][w], r[2025][w]])
    for a_ in axs.ravel():
        a_.set_xlabel("week since 1 Aug"); a_.axvline(8, color="k", lw=0.6, ls=":")
    axs[0, 0].legend(fontsize=7); fig.tight_layout(); fig.savefig(FIG / "E_4_new_sensor_innov.png", dpi=130); plt.close(fig)
    with open(FIG / "E_weekly_tables.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    (FIG / "E_summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
