"""Analysis C scoring: skill of the five VIIRS-transition arms (and the two MODIS continuous series) against
IMS snow cover (per-cell contingency), SNOTEL (SWE, depth) and GHCN (depth); contrasts C0-C4 (config/arms.json).
Reads the compact extractions from skill_extract.py. Paired 5-degree spatial-block bootstrap (cells / stations).
Outputs: tables C_*.csv and figures figC_*.png in docs/viirs_transition_report_figures/. Run with python -I.
"""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
warnings.filterwarnings("ignore")
import os
NBOOT = int(os.environ.get("NBOOT", 1000))
RNG = np.random.default_rng(42)
MIN_VALID_DAYS = 360      # 1500 d / 25 y record (9131 d) scaled to the 2191-d window
MIN_SNOW_DAYS_PER_YEAR = 5.0
MIN_PAIRS_WINDOW = 60     # paired non-JJA days for a station to count in a window
MIN_PAIRS_SY = 30
cfg = sc.load_arms()
ARMS = list(cfg["arms"])
DAYS = pd.date_range("2020-08-01", "2026-07-31")
SY = sc.snow_year_index(DAYS); SEAS = sc.season_index(DAYS)
NT = len(DAYS)

# ---------- units (arm or continuous series) and coverage by snow year
UNITS = {"OL": [("OL", range(6))], "T_to_J2": [("T_to_J2", range(6))], "TA_to_J1J2": [("TA_to_J1J2", range(6))],
         "T_control": [("T_control", range(3, 6))], "TA_control": [("TA_control", range(3, 6))],
         "MODIS_T": [("T_to_J2", range(3)), ("T_control", range(3, 6))],
         "MODIS_TA": [("TA_to_J1J2", range(3)), ("TA_control", range(3, 6))]}
UNIT_SY = {u: sorted(s for _, r in p for s in r) for u, p in UNITS.items()}
SYLAB = [f"SY{i+1}" for i in range(6)]
WIN = {"SY1-6": list(range(6)), "SY1-3": [0, 1, 2], "SY4-6": [3, 4, 5]} | {SYLAB[i]: [i] for i in range(6)}

CONTRASTS = {}   # name -> (terms [(sign, unit, windowkey|'pre'|'post')], windows to evaluate, description)
for u, w in [("T_to_J2", "SY1-6"), ("TA_to_J1J2", "SY1-6"), ("MODIS_T", "SY1-6"), ("MODIS_TA", "SY1-6"), ("T_control", "SY4-6"), ("TA_control", "SY4-6")]:
    CONTRASTS[f"C0 {u}-OL"] = ([(1, u, "*"), (-1, "OL", "*")], [w] + [s for s in SYLAB if int(s[2:]) - 1 in UNIT_SY[u]], "DA impact (DA - OL)")
CONTRASTS["C1 T_control-T_to_J2"] = ([(1, "T_control", "*"), (-1, "T_to_J2", "*")], ["SY4-6"] + SYLAB[3:], "MODIS Terra - VIIRS J2 (same start state 2023-08-01)")
CONTRASTS["C2 TA_control-TA_to_J1J2"] = ([(1, "TA_control", "*"), (-1, "TA_to_J1J2", "*")], ["SY4-6"] + SYLAB[3:], "MODIS Terra+Aqua - VIIRS J1+J2 (same start state)")
CONTRASTS["C3 MODIS_T-MODIS_TA"] = ([(1, "MODIS_T", "*"), (-1, "MODIS_TA", "*")], ["SY1-6"] + SYLAB, "Terra - Terra+Aqua, MODIS, 6 y, same start")
CONTRASTS["C4 DoD (T_to_J2-TA_to_J1J2)"] = ([(1, "T_to_J2", "SY4-6"), (-1, "TA_to_J1J2", "SY4-6"), (-1, "T_to_J2", "SY1-3"), (1, "TA_to_J1J2", "SY1-3")], ["DoD"],
                                            "(T-only minus T+A) in SY4-6 minus the same in SY1-3: J2 vs J1+J2 as a difference of differences")
CONTRASTS["C4a T_to_J2-TA_to_J1J2 (post)"] = ([(1, "T_to_J2", "*"), (-1, "TA_to_J1J2", "*")], ["SY4-6"] + SYLAB[3:], "J2-only minus J1+J2, no baseline correction (different histories)")
CONTRASTS["C4b T_to_J2-TA_to_J1J2 (pre)"] = ([(1, "T_to_J2", "*"), (-1, "TA_to_J1J2", "*")], ["SY1-3"] + SYLAB[:3], "Terra-only minus Terra+Aqua in SY1-3 (baseline for C4)")

def term_window(term_w, win):
    return WIN[win] if term_w == "*" else WIN[term_w]

def load_arm(arm):
    return np.load(sc.OUT / f"skill_{arm}.npz")
Z = {a: load_arm(a) for a in ARMS}
tc = sc.read_tilecoord_quiet()

# ---------- helpers: assemble unit series on the common day axis
def unit_series(u, key):
    """(T, S) array NaN outside the unit's coverage."""
    parts = UNITS[u]; z0 = Z[parts[0][0]][key]
    out = np.full((NT,) + z0.shape[1:], np.nan, np.float32)
    for arm, sys_ in parts:
        zd = pd.DatetimeIndex(Z[arm]["days"]); idx = DAYS.get_indexer(zd)
        a = Z[arm][key]
        sel = np.isin(SY[idx], list(sys_)) & (idx >= 0)
        out[idx[sel]] = a[sel]
    return out

def block_ids(lat, lon, size=5.0):
    return (np.floor((lat + 90) / size).astype(int) * 1000 + np.floor((lon + 180) / size).astype(int))

def ci_from_boot(b):
    return np.nanpercentile(b, [2.5, 97.5])

# =====================================================================================================================
# IMS
# =====================================================================================================================
def ims_part():
    lat = np.asarray(tc["com_lat"], float)[Z["OL"]["rep"]]; lon = np.asarray(tc["com_lon"], float)[Z["OL"]["rep"]]
    ncell = lat.size
    cnt = {}  # unit -> (4, SY, season, cell) int, zero where the unit has no coverage
    for u, parts in UNITS.items():
        c = np.zeros((4, 6, 4, ncell), np.int32)
        for arm, sys_ in parts:
            for s in sys_: c[:, s] = Z[arm]["cnt"][:, s]
        cnt[u] = c
    # NH seasonal-snow mask: >=10 IMS snow days (A+C of OL, all SY, all seasons), and not permanent (snow on <=90% of valid days)
    a, b, c_, d = [cnt["OL"][q].sum(axis=(0, 1)) for q in range(4)]
    snowd = a + c_; valid = a + b + c_ + d
    elig = (snowd >= 10) & (snowd <= 0.9 * valid)
    # season weights: NH exclude JJA, SH exclude DJF  (ALL); season-specific rows use NH cells only
    w_all = np.ones((4, ncell), bool); w_all[2, lat >= 0] = False; w_all[0, lat < 0] = False
    blk = block_ids(lat, lon); ub, binv = np.unique(blk, return_inverse=True)
    def agg(u, sy_list, season):
        c = cnt[u][:, sy_list]                       # 4,sy,season,cell
        if season == "ALL": w = w_all
        else:
            w = np.zeros((4, ncell), bool); w[sc.SEASONS.index(season), lat >= 0] = True
        x = (c * w[None, None]).sum(axis=(1, 2)) * elig[None]     # 4, cell
        out = np.zeros((4, ub.size)); 
        for q in range(4): out[q] = np.bincount(binv, weights=x[q], minlength=ub.size)
        return out                                     # 4, block
    return dict(agg=agg, nb=ub.size, elig=elig, ncell=int(elig.sum()), lat=lat, lon=lon)

def ims_scores(x):                                      # x (...,4) counts A,B,C,D
    A, B, C, D = [x[..., q] for q in range(4)]
    N = A + B + C + D
    with np.errstate(all="ignore"):
        return {"accuracy": (A + D) / N, "hit_rate": A / (A + C), "false_alarm_ratio": B / (A + B),
                "csi": A / (A + B + C)}

IMS_METRICS = ["accuracy", "csi", "hit_rate", "false_alarm_ratio"]

def run_ims():
    P = ims_part(); agg = P["agg"]
    print("IMS eligible cells:", P["ncell"], "blocks:", P["nb"])
    rows_abs, rows_ct = [], []
    # cache block counts
    cache = {}
    def blkcnt(u, win, season):
        k = (u, win, season)
        if k not in cache: cache[k] = agg(u, WIN[win], season)
        return cache[k]
    for u in UNITS:
        for win in ["SY1-6", "SY1-3", "SY4-6"] + SYLAB:
            if not set(WIN[win]) <= set(UNIT_SY[u]): continue
            for season in ["ALL"] + ["DJF", "MAM", "SON"]:
                x = blkcnt(u, win, season).sum(axis=1)
                sc_ = ims_scores(x)
                rows_abs.append(dict(unit=u, window=win, season=season, N=int(x.sum()), **{m: float(sc_[m]) for m in IMS_METRICS}))
    for name, (terms, wins, desc) in CONTRASTS.items():
        for win in wins:
            for season in ["ALL"] + (["DJF", "MAM", "SON"] if win in ("SY1-6", "SY4-6", "SY1-3", "DoD") else []):
                if win == "DoD" and season != "ALL": continue
                bc = []
                ok = True
                for sg, u, tw in terms:
                    w = win if tw == "*" else tw
                    if win == "DoD": w = tw
                    if not set(WIN[w]) <= set(UNIT_SY[u]): ok = False; break
                    bc.append((sg, blkcnt(u, w, season)))
                if not ok: continue
                def stat(idx):
                    out = {m: 0.0 for m in IMS_METRICS}
                    for sg, b in bc:
                        s = ims_scores(b[:, idx].sum(axis=1))
                        for m in IMS_METRICS: out[m] += sg * s[m]
                    return out
                pt = stat(slice(None))
                boots = {m: [] for m in IMS_METRICS}
                for _ in range(NBOOT):
                    idx = RNG.integers(0, P["nb"], P["nb"]); s = stat(idx)
                    for m in IMS_METRICS: boots[m].append(s[m])
                for m in IMS_METRICS:
                    lo, hi = ci_from_boot(np.array(boots[m]))
                    rows_ct.append(dict(contrast=name, window=win, season=season, metric=m, delta=float(pt[m]), ci_lo=float(lo), ci_hi=float(hi), n_blocks=P["nb"], n_cells=P["ncell"]))
    return pd.DataFrame(rows_abs), pd.DataFrame(rows_ct)

# =====================================================================================================================
# Station data (SNOTEL SWE, SNOTEL depth, GHCN depth)
# =====================================================================================================================
def station_metrics(o, m, mask):
    """per-station N, bias, rmse, ubrmse from (T,S) obs/model and boolean mask."""
    n = mask.sum(axis=0).astype(float)
    with np.errstate(all="ignore"):
        o = np.where(mask, o, 0.0); m = np.where(mask, m, 0.0)
        om = o.sum(0) / n; mm = m.sum(0) / n
        e = np.where(mask, m - o, 0.0)
        bias = e.sum(0) / n; rmse = np.sqrt((e ** 2).sum(0) / n)
        ub = np.sqrt((np.where(mask, (m - mm) - (o - om), 0.0) ** 2).sum(0) / n)
    return n, bias, rmse, ub

def season_mask(lat):
    ok = np.ones((NT, lat.size), bool)
    nh = lat >= 0
    ok[np.ix_(SEAS == 2, nh)] = False       # NH: exclude JJA
    ok[np.ix_(SEAS == 0, ~nh)] = False      # SH: exclude DJF
    return ok

def build_dataset(kind):
    if kind == "snotel":
        df, st = sc.load_snotel(); tiles = pd.read_csv(sc.OUT / "snotel_candidate_tiles.csv")
        assert (tiles.station.values == st.station.values).all()
        vars_ = {"SWE": ("swe_mm", "sn_swe", 1.0), "SNWD": ("snwd_mm", "sn_snod", 1000.0)}
        mkey = {"SWE": "sn_swe", "SNWD": "sn_snod"}
    else:
        df, st = sc.load_ghcn(); tiles = pd.read_csv(sc.OUT / "ghcn_candidate_tiles.csv")
        assert (tiles.station.values == st.station.values).all()
        vars_ = {"SNWD": ("snwd_mm", None, 1.0)}
    st = st.copy(); st["dist_km"] = tiles.dist_km.values; st["elev_diff_m"] = tiles.elev_diff_m.values
    keep = (st.dist_km <= sc.MAX_DIST_KM) & (st.elev_diff_m.abs() <= sc.MAX_ABS_ELEV_DIFF_M)
    pos = {s: i for i, s in enumerate(st.station)}
    out = {}
    for var, (col, _, scale) in vars_.items():
        d = df[["date", "station", col]].dropna()
        d = d[(d[col] >= 0) & (d[col] < 5000)]
        d = d[d.station.isin(set(st.station[keep]))]
        piv = d.pivot_table(index="date", columns="station", values=col, aggfunc="mean").reindex(DAYS)
        stns = list(piv.columns); ix = np.array([pos[s] for s in stns])
        o = piv.to_numpy(np.float32)
        # model on the unit axis
        mod = {}
        for u in UNITS:
            if kind == "snotel":
                a = unit_series(u, "sn_swe" if var == "SWE" else "sn_snod")[:, ix]
                if var == "SNWD": a = a * 1000.0                         # m -> mm
            else:
                a = unit_series(u, "gh_snod")[:, ix] * np.clip(unit_series(u, "gh_scf")[:, ix], 0, 1) * 1000.0   # cell-mean depth, mm
            mod[u] = a.astype(np.float32)
        lat = st.lat.values[ix].astype(float); lon = st.lon.values[ix].astype(float)
        # station selection (6-y window): >= MIN_VALID_DAYS valid obs days, >= 5 snow days per snow year on average, OL model finite
        okm = np.isfinite(o) & np.isfinite(mod["OL"])
        nvalid = np.isfinite(o).sum(0)
        snowdays = np.array([(np.where(okm, o, 0)[SY == y] > 0).sum(0) for y in range(6)]).mean(0)
        sm = season_mask(lat)
        npair = (okm & sm).sum(0)
        sel = (nvalid >= MIN_VALID_DAYS) & (snowdays >= MIN_SNOW_DAYS_PER_YEAR) & (npair >= MIN_VALID_DAYS)
        print(f"{kind} {var}: candidates {len(stns)} -> selected {int(sel.sum())} (valid>={MIN_VALID_DAYS}, snow d/y>={MIN_SNOW_DAYS_PER_YEAR})")
        out[var] = dict(o=o[:, sel], mod={u: a[:, sel] for u, a in mod.items()}, lat=lat[sel], lon=lon[sel],
                        names=np.array(stns)[sel], n_cand=len(stns))
    return out

def run_stations(kind, data):
    rows_abs, rows_ct = [], []
    for var, D in data.items():
        o, mod, lat, lon = D["o"], D["mod"], D["lat"], D["lon"]
        S = lat.size
        sm = season_mask(lat); blk = block_ids(lat, lon); ub, binv = np.unique(blk, return_inverse=True); nb = ub.size
        okobs = np.isfinite(o) & sm
        cache = {}
        def metr(u, win, wmask_units):
            k = (u, win)
            if k in cache: return cache[k]
            sy_ok = np.isin(SY, WIN[win])[:, None]
            mask = okobs & sy_ok
            for uu in wmask_units: mask = mask & np.isfinite(mod[uu])
            n, b, r, ubr = station_metrics(o, mod[u], mask)
            thr = MIN_PAIRS_SY if len(WIN[win]) == 1 else MIN_PAIRS_WINDOW
            good = n >= thr
            res = dict(n=n, absbias=np.abs(b), rmse=r, ubrmse=ubr, bias=b, good=good)
            cache[k] = res; return res
        def avail_units(win): return [u for u in UNITS if set(WIN[win]) <= set(UNIT_SY[u])]
        METR = ["rmse", "ubrmse", "absbias", "bias"]
        for u in UNITS:
            for win in ["SY1-6", "SY1-3", "SY4-6"] + SYLAB:
                if not set(WIN[win]) <= set(UNIT_SY[u]): continue
                m_ = metr(u, win, avail_units(win)); g = m_["good"]
                rows_abs.append(dict(dataset=f"{kind}_{var}", unit=u, window=win, n_stations=int(g.sum()),
                                     **{m: float(np.nanmean(m_[m][g])) for m in METR}))
        for name, (terms, wins, desc) in CONTRASTS.items():
            for win in wins:
                ds, good = [], np.ones(S, bool); ok = True
                for sg, u, tw in terms:
                    w = win if (tw == "*" and win != "DoD") else tw
                    if not set(WIN[w]) <= set(UNIT_SY[u]): ok = False; break
                    m_ = metr(u, w, avail_units(w)); good &= m_["good"]; ds.append((sg, m_))
                if not ok: continue
                for m in METR[:3]:
                    d = sum(sg * mm[m] for sg, mm in ds)
                    dg = d[good]; bg = binv[good]
                    pt = float(np.nanmean(dg))
                    # block bootstrap over 5-degree blocks of the selected stations
                    ublk = np.unique(bg); per = [np.where(bg == b)[0] for b in ublk]
                    boots = []
                    for _ in range(NBOOT):
                        pick = RNG.integers(0, len(ublk), len(ublk))
                        idx = np.concatenate([per[p] for p in pick]); boots.append(np.nanmean(dg[idx]))
                    lo, hi = ci_from_boot(np.array(boots))
                    rows_ct.append(dict(dataset=f"{kind}_{var}", contrast=name, window=win, metric=m, delta=pt, ci_lo=float(lo), ci_hi=float(hi),
                                        n_stations=int(good.sum()), n_blocks=len(ublk)))
    return pd.DataFrame(rows_abs), pd.DataFrame(rows_ct)

if __name__ == "__main__":
    sc.OUT.mkdir(exist_ok=True)
    tabs = {}
    a, c = run_ims(); a.insert(0, "dataset", "IMS"); c.insert(0, "dataset", "IMS")
    tabs["ims"] = (a, c)
    for kind in ("snotel", "ghcn"):
        data = build_dataset(kind)
        tabs[kind] = run_stations(kind, data)
    abs_all = pd.concat([t[0] for t in tabs.values()], ignore_index=True)
    ct_all = pd.concat([t[1] for t in tabs.values()], ignore_index=True)
    abs_all.to_csv(sc.FIG / "C_skill_by_arm.csv", index=False); ct_all.to_csv(sc.FIG / "C_skill_contrasts.csv", index=False)
    print(abs_all.shape, ct_all.shape)
