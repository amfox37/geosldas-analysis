#!/usr/bin/env python3
"""
Three-arm (CYGNSS L1 / CYGNSS L3 / SMAP Tb assimilation) mechanism comparison.

Motivation: CYGNSS L1's own OmF_stdv improves by ~13% vs OL in all three arms, even though L1
is only assimilated in one of them (verified 2026-09-25 via ObsFcstAna assim_flag counts). This
script tests whether the arms actually reach that number by the same route:

  (1) daily state:   per-arm daily SFMC/RZMC minus OL -- how big is each arm's state change,
                     and do the arms change the state in the same places/times (per-tile
                     correlation of the DA-minus-OL series, raw and high-pass filtered)?
  (2) increments:    per-arm 3-hourly analysis-minus-forecast SFMC/RZMC (inst3_1d_lndfcstana_Nt)
                     -- frequency, magnitude, spatial pattern, cross-arm agreement at the same
                     tile/cycle, and what the L3/SMAP arms do at the L1 arm's own update tiles.
  (3) L1 O-F split:  every CYGNSS L1 obs matched obs-for-obs across OL and all arms, OmF_stdv %
                     change vs OL binned by time since the L1 arm's last update at that tile and
                     by time since the arm's OWN last update; var(O-F) decomposition into
                     std(O), std(F), corr(O,F); per-tile skill agreement across arms.

Usage:
  cygl1_three_arm_mechanism.py \\
      --arm L1=DAv8_M36_AZ_paired_cygl1_dense075_coh05 \\
      --arm L3=DAv8_M36_AZ_paired_cygl1_dense075_coh05_L3assim \\
      --arm SMAP=DAv8_M36_all_sensors_AZ_scaled_smapassim_6mo \\
      --ol-expid OLv8_M36_AZ_paired_monitor --ol-root <dir holding the extracted OL exp dir> \\
      --start 20200101 --end 20220101 --out-tag three_arm_24mo

The first --arm is treated as the reference arm for the "time since L1-arm update" lag.
Loaded arrays are cached under <out-dir>/cache/ (reuse with the same args; --no-cache to rebuild).
"""
import argparse
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
import netCDF4 as nc  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
OUT_ROOT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"

CYGL1_DESCR = "CYGNSS_L1_DDM3X5_CROP_SCALAR"
FILL_THRESH = 1e14
INCR_EPS = 1e-6  # |ANA-FCST| above this (m3/m3) counts as an update at that tile/cycle
HIGHPASS_DAYS = 31

LAG_BINS = [
    (0, 1, "<1 d"),
    (1, 3, "1-3 d"),
    (3, 7, "3-7 d"),
    (7, 14, "7-14 d"),
    (14, np.inf, ">14 d"),
    (np.nan, np.nan, "never"),
]

# dataviz reference palette: categorical slots 1-3 (validated all-pairs), diverging blue<->red
ARM_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK_MUTED, GRID = "#1f1f1e", "#6b6b66", "#e4e3df"
DIVERGING = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f1a3a2", "#b8302f"])
SEQUENTIAL = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#2a78d6", "#1c5cab", "#0d366b"])


# ----------------------------------------------------------------------------- loading

def exp_dir(root, expid):
    return os.path.join(root, expid, "output", DOMAIN)


def month_dirs(base, start, end):
    out = []
    y, m = start.year, start.month
    while (y, m) < (end.year, end.month) or (y, m) == (end.year, end.month):
        out.append(os.path.join(base, f"Y{y}", f"M{m:02d}"))
        m += 1
        if m == 13:
            m, y = 1, y + 1
    return out


def stamp_to_dt(fname):
    stamp = os.path.basename(fname).split(".")[-2].rstrip("z")
    d, hm = stamp.split("_")
    return pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:8]} {hm[:2]}:{hm[2:]}")


def files_in_range(root, expid, sub, pattern, start, end):
    base = os.path.join(exp_dir(root, expid), sub, "ens_avg")
    fs = []
    for d in month_dirs(base, start, end):
        fs.extend(glob.glob(os.path.join(d, f"{expid}.*{pattern}.*.nc4")))
    fs = sorted(f for f in fs if start <= stamp_to_dt(f) < end)
    if not fs:
        sys.exit(f"ERROR: no {pattern} files for {expid} under {base}")
    return fs


def load_daily(args):
    root, expid, start, end = args
    fs = files_in_range(root, expid, "cat", "tavg24_1d_lnd_Nt", start, end)
    dates, sf, rz = [], [], []
    for f in fs:
        with nc.Dataset(f) as d:
            sf.append(np.asarray(d["SFMC"][0, :], dtype="f8"))
            rz.append(np.asarray(d["RZMC"][0, :], dtype="f8"))
            if not dates:
                lon, lat = np.asarray(d["lon"][:]), np.asarray(d["lat"][:])
        dates.append(stamp_to_dt(f).normalize())
    return dict(dates=np.array(dates, dtype="datetime64[ns]"), SFMC=np.array(sf), RZMC=np.array(rz),
                lon=lon, lat=lat)


def load_incr(args):
    root, expid, start, end = args
    fs = files_in_range(root, expid, "cat", "inst3_1d_lndfcstana_Nt", start, end)
    times, sf, rz = [], [], []
    for f in fs:
        with nc.Dataset(f) as d:
            sf.append(np.asarray(d["SFMC_ANA"][0, :], "f8") - np.asarray(d["SFMC_FCST"][0, :], "f8"))
            rz.append(np.asarray(d["RZMC_ANA"][0, :], "f8") - np.asarray(d["RZMC_FCST"][0, :], "f8"))
        times.append(stamp_to_dt(f))
    return dict(times=np.array(times, dtype="datetime64[ns]"), dSFMC=np.array(sf), dRZMC=np.array(rz))


def load_l1_ofa(args):
    root, expid, start, end = args
    fs = files_in_range(root, expid, "ana", "ldas_ObsFcstAna", start, end)
    rows = []
    for f in fs:
        dt = stamp_to_dt(f)
        with nc.Dataset(f) as d:
            if d.dimensions["n_obs"].size == 0:
                continue
            names = {str(s): int(i) for s, i in zip(d["obsparam_descr"][:], d["obsparam_species_id"][:])}
            if CYGL1_DESCR not in names:
                continue
            sp = np.asarray(d["species"][:])
            m = sp == names[CYGL1_DESCR]
            if not m.any():
                continue
            obs = np.asarray(d["obs"][:])[m]
            fc = np.asarray(d["fcst"][:])[m]
            ok = (np.abs(obs) < FILL_THRESH) & (np.abs(fc) < FILL_THRESH)
            tile = np.asarray(d["tilenum"][:])[m] - 1  # 0-based row into the tile-space arrays
            af = np.asarray(d["assim_flag"][:])[m]
            for k in np.nonzero(ok)[0]:
                rows.append((int(tile[k]), dt, float(obs[k]), float(fc[k]), int(af[k])))
    return pd.DataFrame(rows, columns=["tile", "time", "obs", "fcst", "assim"])


def cached(path, fn, argtuple, use_cache):
    if use_cache and os.path.exists(path):
        if path.endswith(".npz"):
            z = np.load(path, allow_pickle=False)
            return {k: z[k] for k in z.files}
        return pd.read_parquet(path)
    out = fn(argtuple)
    if path.endswith(".npz"):
        np.savez_compressed(path, **out)
    else:
        out.to_parquet(path)
    return out


# ----------------------------------------------------------------------------- helpers

def md(df, floatfmt=".3f", index=False):
    """Minimal DataFrame -> markdown table (tabulate isn't installed in GEOSpyD)."""
    if index:
        df = df.reset_index()
    fmt = lambda v: format(v, floatfmt[1:]) if isinstance(v, (float, np.floating)) and not np.isnan(v) else ("" if isinstance(v, (float, np.floating)) else str(v))
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, sep] + body)


def pct(a, b):
    return 100.0 * (a - b) / b


def rowcorr(a, b):
    """Per-column Pearson correlation of two (time, tile) arrays."""
    a = a - a.mean(0)
    b = b - b.mean(0)
    den = np.sqrt((a * a).sum(0) * (b * b).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (a * b).sum(0) / den, np.nan)


def highpass(x, days):
    return x - pd.DataFrame(x).rolling(days, center=True, min_periods=days // 2).mean().to_numpy()


def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def map_panels(fname, lon, lat, fields, titles, cmap, vmin, vmax, cbar_label, suptitle, marks=None):
    fig, axes = plt.subplots(1, len(fields), figsize=(4.2 * len(fields), 4.2), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, fld, t in zip(axes, fields, titles):
        sc = ax.scatter(lon, lat, c=fld, cmap=cmap, vmin=vmin, vmax=vmax, s=34, marker="s", linewidths=0)
        if marks is not None:
            ax.scatter(lon[marks], lat[marks], s=4, c=INK, linewidths=0)
        ax.set_title(t, fontsize=10, color=INK, loc="left")
        ax.set_aspect("equal")
        style_axes(ax)
        ax.grid(False)
    cb = fig.colorbar(sc, ax=axes, shrink=0.85)
    cb.set_label(cbar_label, color=INK_MUTED, fontsize=9)
    cb.ax.tick_params(labelsize=8, colors=INK_MUTED)
    fig.suptitle(suptitle, fontsize=11, color=INK, x=0.01, ha="left")
    fig.savefig(fname, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------- analyses

def analysis_state(tags, ol, daily, l1_tiles, out, rep, figdir):
    rep.append("## (1) Daily state: arm minus OL\n")
    for var in ("SFMC", "RZMC"):
        delta = {t: daily[t][var] - ol[var] for t in tags}
        hp = {t: highpass(delta[t], HIGHPASS_DAYS) for t in tags}
        valid_hp = ~np.isnan(hp[tags[0]]).any(1)
        rows = []
        for t in tags:
            rms = np.sqrt(np.nanmean(delta[t] ** 2, 0))
            rows.append(dict(arm=t, var=var,
                             rms_delta_median_tile=np.median(rms),
                             rms_delta_median_L1tiles=np.median(rms[l1_tiles]),
                             mean_delta_domain=np.nanmean(delta[t]),
                             rms_highpass_median_tile=np.median(np.sqrt(np.nanmean(hp[t][valid_hp] ** 2, 0)))))
        df = pd.DataFrame(rows)
        rep.append(f"### {var}: size of each arm's change vs OL (m3/m3)\n")
        rep.append(md(df, ".5f") + "\n")
        df.to_csv(os.path.join(out, f"state_{var}_magnitude.csv"), index=False)

        prow = []
        cmaps = {}
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                a, b = tags[i], tags[j]
                r_raw = rowcorr(delta[a], delta[b])
                r_hp = rowcorr(hp[a][valid_hp], hp[b][valid_hp])
                r_dom = np.corrcoef(delta[a].mean(1), delta[b].mean(1))[0, 1]
                cmaps[f"{a} vs {b}"] = r_hp
                prow.append(dict(pair=f"{a}-{b}", var=var,
                                 r_raw_median_tile=np.nanmedian(r_raw),
                                 r_highpass_median_tile=np.nanmedian(r_hp),
                                 r_highpass_median_L1tiles=np.nanmedian(r_hp[l1_tiles]),
                                 r_domain_mean_series=r_dom))
        pdf = pd.DataFrame(prow)
        rep.append(f"### {var}: do the arms change the state the same way? (correlation of DA-minus-OL series)\n")
        rep.append(md(pdf, ".3f") + "\n")
        pdf.to_csv(os.path.join(out, f"state_{var}_pair_corr.csv"), index=False)

        lon, lat = ol["lon"], ol["lat"]
        rmsf = [np.sqrt(np.nanmean(delta[t] ** 2, 0)) for t in tags]
        map_panels(os.path.join(figdir, f"state_{var}_rms_delta_maps.png"), lon, lat, rmsf,
                   [f"{t} assim" for t in tags], SEQUENTIAL, 0, np.nanpercentile(np.concatenate(rmsf), 98),
                   f"RMS daily {var} change vs OL (m3/m3)",
                   f"Daily {var}: RMS of arm minus OL (dots = tiles with CYGNSS L1 obs)", marks=l1_tiles)
        map_panels(os.path.join(figdir, f"state_{var}_pair_corr_maps.png"), lon, lat, list(cmaps.values()),
                   list(cmaps.keys()), DIVERGING, -1, 1, "correlation",
                   f"Daily {var}: per-tile correlation of high-pass ({HIGHPASS_DAYS} d) arm-minus-OL series")


def last_update_lag(upd_mask, upd_times, tiles, times):
    """Days since the last update strictly before `times` at `tiles` (NaN if none)."""
    lag = np.full(len(tiles), np.nan)
    tnum = upd_times.astype("datetime64[m]").astype("int64")
    q = np.asarray(times, dtype="datetime64[m]").astype("int64")
    for tile in np.unique(tiles):
        sel = np.nonzero(tiles == tile)[0]
        ut = tnum[upd_mask[:, tile]]
        if ut.size == 0:
            continue
        k = np.searchsorted(ut, q[sel], side="left") - 1
        ok = k >= 0
        lag[sel[ok]] = (q[sel[ok]] - ut[k[ok]]) / 1440.0
    return lag


def analysis_incr(tags, incr, l1_events, out, rep, figdir, lon, lat, l1_tiles):
    rep.append("## (2) Increments: 3-hourly analysis minus forecast\n")
    times = incr[tags[0]]["times"]
    for t in tags[1:]:
        if not np.array_equal(incr[t]["times"], times):
            sys.exit("ERROR: increment time axes differ between arms")
    rows = []
    for t in tags:
        for var in ("dSFMC", "dRZMC"):
            x = incr[t][var]
            upd = np.abs(incr[t]["dSFMC"]) > INCR_EPS
            rows.append(dict(arm=t, var=var[1:],
                             frac_tile_cycles_updated=upd.mean(),
                             tiles_ever_updated=int(upd.any(0).sum()),
                             mean_abs_incr_when_updated=np.abs(x[upd]).mean() if upd.any() else np.nan,
                             net_incr_per_tile_per_yr_median=np.median(x.sum(0)) / 2.0,
                             frac_wetting=(x[upd] > 0).mean() if upd.any() else np.nan))
    df = pd.DataFrame(rows)
    rep.append("### Frequency and size of updates (update = |dSFMC| > 1e-6 m3/m3)\n")
    rep.append(md(df, ".5f") + "\n")
    df.to_csv(os.path.join(out, "incr_summary.csv"), index=False)

    prow = []
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            a, b = tags[i], tags[j]
            xa, xb = incr[a]["dSFMC"], incr[b]["dSFMC"]
            ua, ub = np.abs(xa) > INCR_EPS, np.abs(xb) > INCR_EPS
            both = ua & ub
            netA, netB = xa.sum(0), xb.sum(0)
            # daily-summed increments, correlated per tile over time
            day = times.astype("datetime64[D]")
            _, inv = np.unique(day, return_inverse=True)
            da = np.zeros((inv.max() + 1, xa.shape[1])); np.add.at(da, inv, xa)
            db = np.zeros_like(da); np.add.at(db, inv, xb)
            prow.append(dict(pair=f"{a}-{b}",
                             overlap_frac_of_A_updates=both.sum() / max(ua.sum(), 1),
                             overlap_frac_of_B_updates=both.sum() / max(ub.sum(), 1),
                             same_cycle_corr=np.corrcoef(xa[both], xb[both])[0, 1] if both.sum() > 2 else np.nan,
                             same_cycle_sign_agree=(np.sign(xa[both]) == np.sign(xb[both])).mean() if both.any() else np.nan,
                             daily_sum_corr_median_tile=np.nanmedian(rowcorr(da, db)),
                             net_incr_spatial_corr=np.corrcoef(netA, netB)[0, 1]))
    pdf = pd.DataFrame(prow)
    rep.append("### Cross-arm agreement of SFMC increments\n")
    rep.append(md(pdf, ".3f") + "\n")
    pdf.to_csv(os.path.join(out, "incr_pair_agreement.csv"), index=False)

    # what do the other arms do at the reference arm's own L1 update tile-cycles?
    ref = tags[0]
    tix = {v: i for i, v in enumerate(times.astype("datetime64[m]").astype("int64"))}
    ev = l1_events.drop_duplicates(["tile", "time"])
    ci = np.array([tix.get(v, -1) for v in ev["time"].values.astype("datetime64[m]").astype("int64")])
    ok = ci >= 0
    ci, ti = ci[ok], ev["tile"].values[ok]
    erow = []
    xref = incr[ref]["dSFMC"][ci, ti]
    for t in tags:
        x = incr[t]["dSFMC"]
        same = x[ci, ti]
        # +-1 day window sum at the same tile
        win = np.array([x[max(c - 8, 0):c + 9, tl].sum() for c, tl in zip(ci, ti)])
        erow.append(dict(arm=t, n_events=len(ci),
                         frac_updated_same_cycle=(np.abs(same) > INCR_EPS).mean(),
                         mean_abs_incr_same_cycle=np.abs(same).mean(),
                         mean_abs_incr_pm1day=np.abs(win).mean(),
                         corr_with_ref_same_cycle=np.corrcoef(xref, same)[0, 1] if t != ref and same.std() > 0 else np.nan,
                         corr_with_ref_pm1day=np.corrcoef(xref, win)[0, 1] if t != ref and win.std() > 0 else np.nan,
                         sign_agree_with_ref_pm1day=(np.sign(win[np.abs(win) > INCR_EPS]) == np.sign(xref[np.abs(win) > INCR_EPS])).mean() if t != ref else np.nan))
    edf = pd.DataFrame(erow)
    rep.append(f"### SFMC increments at the {ref} arm's own L1-assimilation tile-cycles\n")
    rep.append(md(edf, ".4f") + "\n")
    edf.to_csv(os.path.join(out, "incr_at_ref_L1_events.csv"), index=False)

    freq = [(np.abs(incr[t]["dSFMC"]) > INCR_EPS).mean(0) * 8 * 365.25 for t in tags]
    map_panels(os.path.join(figdir, "incr_update_frequency_maps.png"), lon, lat, freq,
               [f"{t} assim" for t in tags], SEQUENTIAL, 0, np.nanpercentile(np.concatenate(freq), 98),
               "SFMC updates per year", "How often each tile gets an SFMC update (dots = tiles with CYGNSS L1 obs)",
               marks=l1_tiles)
    net = [incr[t]["dSFMC"].sum(0) / 2.0 for t in tags]
    vmax = np.nanpercentile(np.abs(np.concatenate(net)), 98)
    map_panels(os.path.join(figdir, "incr_net_sfmc_maps.png"), lon, lat, net,
               [f"{t} assim" for t in tags], DIVERGING.reversed(), -vmax, vmax,
               "net SFMC increment per year (m3/m3; blue = wetting)",
               "Net SFMC increment per tile per year")


def analysis_l1_omf(tags, ol_ofa, arm_ofa, incr, out, rep, figdir, lon, lat, min_tile_obs):
    rep.append("## (3) CYGNSS L1 O-F, matched obs-for-obs across OL and all arms\n")
    key = ["tile", "time", "obs"]
    m = ol_ofa[key + ["fcst"]].rename(columns={"fcst": "fcst_OL"}).assign(obs=lambda d: d.obs.round(5))
    for t in tags:
        a = arm_ofa[t][key + ["fcst"]].rename(columns={"fcst": f"fcst_{t}"}).assign(obs=lambda d: d.obs.round(5))
        m = m.merge(a, on=key, how="inner")
    m = m.drop_duplicates(key)
    rep.append(f"Matched L1 obs in all {len(tags) + 1} runs: N = {len(m)} "
               f"(OL {len(ol_ofa)}, " + ", ".join(f"{t} {len(arm_ofa[t])}" for t in tags) + ")\n")

    o = m["obs"].to_numpy()
    rows = []
    for t in ["OL"] + tags:
        f = m[f"fcst_{t}"].to_numpy()
        rows.append(dict(run=t, omf_std=np.std(o - f), omf_mean=np.mean(o - f), std_obs=np.std(o),
                         std_fcst=np.std(f), corr_obs_fcst=np.corrcoef(o, f)[0, 1]))
    df = pd.DataFrame(rows)
    base = df.iloc[0]
    for c in ("omf_std", "std_fcst"):
        df[f"{c}_pct_vs_OL"] = pct(df[c], base[c])
    rep.append("### Overall, with var(O-F) = var(O) + var(F) - 2 cov(O,F) decomposition\n")
    rep.append(md(df, ".4f") + "\n")
    df.to_csv(os.path.join(out, "l1_omf_decomposition.csv"), index=False)

    # counterfactual: OmF_stdv if only corr changed (std_fcst held at OL) vs only std_fcst changed
    crow = []
    so, sfo, ro = base["std_obs"], base["std_fcst"], base["corr_obs_fcst"]
    for _, r in df.iloc[1:].iterrows():
        only_r = np.sqrt(so ** 2 + sfo ** 2 - 2 * r["corr_obs_fcst"] * so * sfo)
        only_s = np.sqrt(so ** 2 + r["std_fcst"] ** 2 - 2 * ro * so * r["std_fcst"])
        crow.append(dict(run=r["run"], actual_pct=pct(r["omf_std"], base["omf_std"]),
                         corr_change_only_pct=pct(only_r, base["omf_std"]),
                         fcst_std_change_only_pct=pct(only_s, base["omf_std"])))
    cdf = pd.DataFrame(crow)
    rep.append("### Which part of the O-F change comes from better correlation vs a different forecast spread?\n")
    rep.append(md(cdf, ".2f") + "\n")
    cdf.to_csv(os.path.join(out, "l1_omf_counterfactual.csv"), index=False)

    # lag since last update
    ref = tags[0]
    upd = {t: np.abs(incr[t]["dSFMC"]) > INCR_EPS for t in tags}
    itimes = incr[ref]["times"]
    tiles = m["tile"].to_numpy()
    times = m["time"].to_numpy()
    lag_ref = last_update_lag(upd[ref], itimes, tiles, times)
    lag_own = {t: last_update_lag(upd[t], itimes, tiles, times) for t in tags}

    def binned(lag, t):
        out_rows = []
        for lo, hi, lab in LAG_BINS:
            sel = np.isnan(lag) if lab == "never" else (lag >= lo) & (lag < hi)
            n = int(sel.sum())
            if n < 30:
                out_rows.append(dict(bin=lab, n=n, pct=np.nan))
                continue
            s_arm = np.std(o[sel] - m[f"fcst_{t}"].to_numpy()[sel])
            s_ol = np.std(o[sel] - m["fcst_OL"].to_numpy()[sel])
            out_rows.append(dict(bin=lab, n=n, pct=pct(s_arm, s_ol)))
        return pd.DataFrame(out_rows)

    brow = []
    for t in tags:
        for kind, lag in ((f"since last {ref}-arm update", lag_ref), ("since arm's own last update", lag_own[t])):
            b = binned(lag, t)
            b.insert(0, "lag_kind", kind)
            b.insert(0, "arm", t)
            brow.append(b)
    bdf = pd.concat(brow, ignore_index=True)
    rep.append("### L1 OmF_stdv % change vs OL, by time since last SFMC update at the obs tile\n")
    piv = bdf.pivot_table(index=["lag_kind", "bin"], columns="arm", values="pct", sort=False)
    npiv = bdf.pivot_table(index=["lag_kind", "bin"], columns="arm", values="n", sort=False)
    rep.append(md(piv, ".2f", index=True) + "\n\nN per bin:\n\n" + md(npiv, ".0f", index=True) + "\n")
    bdf.to_csv(os.path.join(out, "l1_omf_by_lag.csv"), index=False)

    # per-tile skill agreement
    g = m.assign(**{f"e_{t}": o - m[f"fcst_{t}"] for t in ["OL"] + tags}).groupby("tile")
    cnt = g.size()
    keep = cnt[cnt >= min_tile_obs].index
    if len(keep) < 10:
        rep.append(f"(per-tile comparison skipped: only {len(keep)} tiles with >= {min_tile_obs} matched obs)\n")
        return
    ts = pd.DataFrame({t: pct(g[f"e_{t}"].std().loc[keep], g["e_OL"].std().loc[keep]) for t in tags})
    trow = []
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            a, b = tags[i], tags[j]
            trow.append(dict(pair=f"{a}-{b}", n_tiles=len(ts), pearson=ts[a].corr(ts[b]),
                             spearman=ts[a].corr(ts[b], method="spearman")))
    tdf = pd.DataFrame(trow)
    rep.append(f"### Do the same tiles improve? Per-tile L1 OmF_stdv % change vs OL, correlated across arms "
               f"(tiles with >= {min_tile_obs} matched obs)\n")
    rep.append(md(ts.describe().T[["mean", "50%", "std"]].rename(columns={"50%": "median"}), ".2f", index=True) + "\n\n")
    rep.append(md(tdf, ".3f") + "\n")
    ts.to_csv(os.path.join(out, "l1_omf_per_tile_pct.csv"))

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    labs = [b[2] for b in LAG_BINS]
    for ax, kind in zip(axes, [f"since last {ref}-arm update", "since arm's own last update"]):
        for t, col in zip(tags, ARM_COLORS):
            b = bdf[(bdf.arm == t) & (bdf.lag_kind == kind)]
            ax.plot(range(len(labs)), b["pct"], color=col, lw=2, marker="o", ms=6, label=f"{t} assim")
        ax.axhline(0, color=INK_MUTED, lw=1)
        ax.set_xticks(range(len(labs)), labs)
        ax.set_title(f"Time {kind}", fontsize=10, color=INK, loc="left")
        style_axes(ax)
    axes[0].set_ylabel("L1 OmF std, % change vs OL", color=INK_MUTED, fontsize=9)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("CYGNSS L1 O-F improvement by time since the last update at the obs tile",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.savefig(os.path.join(figdir, "l1_omf_by_lag.png"), dpi=150)
    plt.close(fig)

    others = tags[1:]
    fig, axes = plt.subplots(1, len(others), figsize=(4.4 * len(others), 4.2), constrained_layout=True)
    axes = np.atleast_1d(axes)
    lim = np.nanpercentile(np.abs(ts.to_numpy()), 99)
    for ax, t, col in zip(axes, others, ARM_COLORS[1:]):
        ax.scatter(ts[ref], ts[t], s=14, color=col, alpha=0.75, linewidths=0)
        ax.plot([-lim, lim], [-lim, lim], color=INK_MUTED, lw=1)
        ax.axhline(0, color=GRID, lw=1); ax.axvline(0, color=GRID, lw=1)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_xlabel(f"{ref} assim: % change vs OL", color=INK_MUTED, fontsize=9)
        ax.set_ylabel(f"{t} assim: % change vs OL", color=INK_MUTED, fontsize=9)
        r = ts[ref].corr(ts[t])
        ax.set_title(f"{ref} vs {t} (r = {r:.2f})", fontsize=10, color=INK, loc="left")
        style_axes(ax)
    fig.suptitle("Per-tile L1 OmF std change vs OL: do the same tiles improve?",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.savefig(os.path.join(figdir, "l1_omf_per_tile_scatter.png"), dpi=150)
    plt.close(fig)

    tl = ts.index.to_numpy()
    map_panels(os.path.join(figdir, "l1_omf_per_tile_maps.png"), lon[tl], lat[tl], [ts[t].to_numpy() for t in tags],
               [f"{t} assim" for t in tags], DIVERGING, -lim, lim, "L1 OmF std % change vs OL (blue = better)",
               "Per-tile CYGNSS L1 OmF std change vs OL")


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", action="append", required=True, help="TAG=EXPID (repeat; first is the L1 reference arm)")
    ap.add_argument("--ol-expid", default="OLv8_M36_AZ_paired_monitor")
    ap.add_argument("--ol-root", default=EXP_ROOT, help="dir containing the OL exp dir (e.g. an archive extraction)")
    ap.add_argument("--exp-root", default=EXP_ROOT)
    ap.add_argument("--start", default="20200101")
    ap.add_argument("--end", default="20220101")
    ap.add_argument("--out-tag", default="three_arm")
    ap.add_argument("--min-tile-obs", type=int, default=30, help="min matched L1 obs for a tile in the per-tile comparison")
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    arms = dict(a.split("=", 1) for a in args.arm)
    tags = list(arms)
    start, end = pd.Timestamp(args.start), pd.Timestamp(args.end)
    out = os.path.join(OUT_ROOT, "three_arm_mechanism", args.out_tag)
    figdir = os.path.join(out, "figures")
    cache = os.path.join(out, "cache")
    for d in (out, figdir, cache):
        os.makedirs(d, exist_ok=True)
    uc = not args.no_cache

    jobs = {("daily", "OL"): (load_daily, (args.ol_root, args.ol_expid, start, end), "npz"),
            ("ofa", "OL"): (load_l1_ofa, (args.ol_root, args.ol_expid, start, end), "parquet")}
    for t, e in arms.items():
        jobs[("daily", t)] = (load_daily, (args.exp_root, e, start, end), "npz")
        jobs[("incr", t)] = (load_incr, (args.exp_root, e, start, end), "npz")
        jobs[("ofa", t)] = (load_l1_ofa, (args.exp_root, e, start, end), "parquet")
    res = {}
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {k: ex.submit(cached, os.path.join(cache, f"{k[0]}_{k[1]}_{args.start}_{args.end}.{ext}"), fn, a, uc)
                for k, (fn, a, ext) in jobs.items()}
        for k, f in futs.items():
            res[k] = f.result()
            print(f"loaded {k}", flush=True)

    ol = res[("daily", "OL")]
    daily = {t: res[("daily", t)] for t in tags}
    for t in tags:
        if not (np.array_equal(daily[t]["dates"], ol["dates"]) and np.array_equal(daily[t]["lon"], ol["lon"])):
            sys.exit(f"ERROR: daily axes of {t} differ from OL")
    lon, lat = ol["lon"], ol["lat"]
    ol_ofa = res[("ofa", "OL")]
    arm_ofa = {t: res[("ofa", t)] for t in tags}
    incr = {t: res[("incr", t)] for t in tags}
    l1_tiles = np.zeros(len(lon), bool)
    l1_tiles[ol_ofa["tile"].unique()] = True
    ref_events = arm_ofa[tags[0]].query("assim == 1")
    if ref_events.empty:
        print(f"WARNING: reference arm {tags[0]} has no assimilated L1 obs", flush=True)

    rep = [f"# Three-arm mechanism comparison ({args.start} to {args.end})\n",
           "Arms: " + ", ".join(f"**{t}** = `{e}`" for t, e in arms.items()) + f"; OL = `{args.ol_expid}`\n",
           f"{len(lon)} tiles, {l1_tiles.sum()} with CYGNSS L1 obs. "
           f"L1 assimilation events in {tags[0]} arm: {len(ref_events)}.\n"]
    analysis_state(tags, ol, daily, l1_tiles, out, rep, figdir)
    analysis_incr(tags, incr, ref_events, out, rep, figdir, lon, lat, l1_tiles)
    analysis_l1_omf(tags, ol_ofa, arm_ofa, incr, out, rep, figdir, lon, lat, args.min_tile_obs)
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(rep))
    print(f"wrote {out}/report.md", flush=True)


if __name__ == "__main__":
    main()
