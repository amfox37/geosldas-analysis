#!/usr/bin/env python3
"""Tile-level 'benefit table': DA arm minus OL, sign so positive = DA better.

Arms: L1 = DA_L1_full_xc015_coh040216_fixedop (errstd 2.75, coherency filter), L3 = DA_L3_fixedop,
OL = OLv8_M36_AZ_fixedop.  Period 2020-2022.  Key: tile_axis_index (0..908).

Always written (local data):
  * tile keys: tile_id, idx0 (derived from com_lon/com_lat on the EASEv2 M36 global grid, 964 x 406,
    j*964+i, same definition as build_iv_tc_tile_map.py), com_lon, com_lat, n_l1_obs, in_l1_mask
    (n_l1_obs >= 100, as build_iv_tc_tile_map.py).
  * Monitor O-F: per-tile % change of O-F stdv (sign flipped) for SMOS Tb (species axis 0-3), SMAP Tb (4-7),
    ASCAT (8-10), from temporal_stats nc4 of the stats bundle; DA file vs the OL xmask file of the same tag.
    Species of a group are pooled per tile N-weighted:  var = sum N (s^2+m^2)/sum N - pooled_mean^2
    (same pooling as the bundle README / spatial pkl).  Also the arm's own L1 species (axis 12) and L3 (11).
    Whole period and per year (2020, 2021, 2022).  Sample sizes: N_data per tile (DA and OL).
Written only if --iv-tc-root is available (it lives on Discover, not local):
  * TC (SMAP L3 primary, ASCAT H121 secondary), IV, via iv_tc_az_maps.compare_tc/compare_iv
    (NMIN TC = 20 paired days, IV = 100 paired days, lag 2, in_l1_mask, baseline sigma2 floor 10th pct).
Per-tile dMSE / alpha_opt: NOT available locally (noise_gain csv is domain-wide per month only).
Seasonal: domain-level O-F and dMSE/alpha written to benefit_table_seasonal_domain.csv (per-tile seasonal
needs per-tile monthly stats or ens_avg OFA files, which are on Discover).
"""
from __future__ import annotations
import argparse, pickle, sys
from pathlib import Path
import numpy as np, pandas as pd
import netCDF4 as nc
from pyproj import Transformer

PROJ = Path(__file__).resolve().parents[2]
STATS = PROJ / "output/cygl1_fixedop_2020_2022_stats"
OUT = PROJ / "output/benefit_drivers"
TAGS = {"L1": "full_xc015_coh040216_fixedop", "L3": "l3_fixedop"}
GROUPS = {"smos": [0, 1, 2, 3], "smap": [4, 5, 6, 7], "ascat": [8, 9, 10], "l3obs": [11], "l1obs": [12]}
PERIODS = {"all": ("stats", "20200101_20221231"), "2020": ("stats/per_year", "20200101_20201231"),
           "2021": ("stats/per_year", "20210101_20211231"), "2022": ("stats/per_year", "20220101_20221231")}
M36_NX, M36_NY, CELL = 964, 406, 36032.220840584
NMIN_L1_OBS = 100


def tile_keys():
    t = pd.read_csv(STATS / "tile_coords.csv")
    x, y = Transformer.from_crs(4326, 6933, always_xy=True).transform(t.com_lon.values, t.com_lat.values)
    i = np.floor(x / CELL + M36_NX / 2).astype(int); j = np.floor(M36_NY / 2 - y / CELL).astype(int)
    t["idx0"] = j * M36_NX + i
    assert t.idx0.nunique() == len(t)
    c = pd.read_csv(PROJ / "output/ismn_fixedop_20200101_20221231/tile_l1_obs_count_2020_2022.csv")
    t = t.merge(c.rename(columns={"tile_index": "tile_axis_index"}), on="tile_axis_index", how="left")
    t["in_l1_mask"] = t.n_l1_obs >= NMIN_L1_OBS
    return t[["tile_axis_index", "tile_id", "idx0", "com_lon", "com_lat", "n_l1_obs", "in_l1_mask"]]


def read_nc(path):
    with nc.Dataset(path) as d:
        out = {}
        for k in ("N_data", "OmF_mean", "OmF_stdv"):
            a = np.ma.filled(d[k][:], np.nan).astype(float)
            a[a <= -9998] = np.nan
            out[k] = a
    return out


def pooled(d, sp):
    N, m, s = (d[k][:, sp] for k in ("N_data", "OmF_mean", "OmF_stdv"))
    ok = np.isfinite(N) & np.isfinite(m) & np.isfinite(s) & (N > 0)
    Nz = np.where(ok, N, 0.0)
    tot = Nz.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pm = (Nz * np.where(ok, m, 0)).sum(1) / tot
        v = (Nz * np.where(ok, s ** 2 + m ** 2, 0)).sum(1) / tot - pm ** 2
    sd = np.sqrt(np.clip(v, 0, None)); sd[tot == 0] = np.nan
    return sd, tot


def omf_columns(keys):
    cols = {}
    for arm, tag in TAGS.items():
        for per, (sub, span) in PERIODS.items():
            da = read_nc(STATS / sub / f"temporal_stats_DA_paired_{tag}_{span}.nc4")
            ol = read_nc(STATS / sub / f"temporal_stats_OL_paired_monitor_xmask_{tag}_{span}.nc4")
            for g, sp in GROUPS.items():
                sd_da, n_da = pooled(da, sp); sd_ol, n_ol = pooled(ol, sp)
                with np.errstate(invalid="ignore", divide="ignore"):
                    pct = -100.0 * (sd_da - sd_ol) / sd_ol
                sfx = "" if per == "all" else f"_{per}"
                cols[f"omf_{g}_pct_better_{arm}{sfx}"] = pct
                cols[f"omf_{g}_N_{arm}{sfx}"] = n_da
                cols[f"omf_{g}_N_OLxmask_{arm}{sfx}"] = n_ol
    df = pd.DataFrame(cols)
    # L1 minus L3 (whole period + per year), only where both arms have a value; positive = L1 better
    for g in ("smos", "smap", "ascat"):
        for per in PERIODS:
            sfx = "" if per == "all" else f"_{per}"
            df[f"omf_{g}_pct_better_L1minusL3{sfx}"] = df[f"omf_{g}_pct_better_L1{sfx}"] - df[f"omf_{g}_pct_better_L3{sfx}"]
    return df


def tc_iv_columns(root, keys):
    sys.path.insert(0, str(PROJ / "scripts"))
    import iv_tc_az_maps as m
    tm = keys.copy()
    tc = {r: m.load_tc(m.tc_file(root, r), tm) for r in ("OL", "L1_filt275", "L3")}
    sens = ("smap_l3", "ascat_h121")
    iv = {(s, r): m.load_iv(m.iv_file(root, s, r), tm) for s in sens for r in ("OL", "L1_filt275", "L3")}
    out = {}
    for arm, run in (("L1", "L1_filt275"), ("L3", "L3")):
        c = m.compare_tc(tc["OL"], tc[run], tm)
        out[f"tc_dR2mod_{arm}"] = c.d_r2
        out[f"tc_pct_sig2_better_{arm}"] = c.pct_sig2_inv
        out[f"tc_R2mod_OL_{arm}"] = c.r2_ref
        out[f"tc_N_{arm}"] = np.minimum(tc["OL"].N.to_numpy(), tc[run].N.to_numpy())
        for s in sens:
            ok, d = m.compare_iv(iv[(s, "OL")], iv[(s, run)], tm)
            out[f"iv_dAnomR_{s}_{arm}"] = d
            out[f"iv_N_{s}_{arm}"] = np.minimum(iv[(s, "OL")].N.to_numpy(), iv[(s, run)].N.to_numpy())
    df = pd.DataFrame(out)
    for c_ in [c for c in df.columns if c.endswith("_L1") and not c.startswith(("tc_N", "iv_N", "tc_R2"))]:
        base = c_[:-3]
        df[base + "_L1minusL3"] = df[c_] - df[base + "_L3"]
    return df


def seasonal_domain():
    season = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
              9: "SON", 10: "SON", 11: "SON"}
    gsp = {"SMOS": [0, 1, 2, 3], "SMAP": [4, 5, 6, 7], "ASCAT": [8, 9, 10], "L3obs": [11], "L1obs": [12]}
    rows = []
    for arm, tag in TAGS.items():
        da = pickle.load(open(STATS / f"stats/spatial_stats_DA_paired_{tag}_202001_202212.pkl", "rb"))
        ol = pickle.load(open(STATS / f"stats/spatial_stats_OL_paired_monitor_xmask_{tag}_202001_202212.pkl", "rb"))
        mon = np.array([int(s[4:]) for s in da["date_vec"]])
        def pst(d, rows_, sp):
            N, m, s = (np.asarray(d[k])[rows_][:, sp] for k in ("N_data", "OmF_mean", "OmF_stdv"))
            ok = (N > 0) & np.isfinite(m) & np.isfinite(s)
            N, m, s = N[ok], m[ok], s[ok]
            pm = (N * m).sum() / N.sum()
            return np.sqrt(max((N * (s ** 2 + m ** 2)).sum() / N.sum() - pm ** 2, 0)), N.sum()
        ng = pd.read_csv(STATS / f"summary/noise_gain_{tag}_202001_202212.csv")
        ng = ng[ng.month != "all"].copy(); ng["mo"] = ng.month.str[5:].astype(int)
        for sea in ("DJF", "MAM", "JJA", "SON"):
            sel = np.array([season[m_] == sea for m_ in mon])
            for g, sps in gsp.items():
                per_sp = [(pst(da, sel, s), pst(ol, sel, s)) for s in sps]
                pct = np.mean([-100 * (a[0] - b[0]) / b[0] for a, b in per_sp])  # mean of per-species %, as bundle summary
                r = dict(arm=arm, season=sea, group=g, omf_pct_better=pct, N_da=sum(a[1] for a, b in per_sp))
                gg = ng[(ng.grp == ("L3" if g == "L3obs" else "L1" if g == "L1obs" else g)) & ng.mo.map(season).eq(sea)]
                if len(gg):
                    w = gg.N.to_numpy()
                    dm, nz, gn = (np.average(gg[k], weights=w) for k in ("dMSE", "noise", "gain"))
                    r.update(dMSE_pct_better=-dm, noise=nz, gain=gn, alpha_opt_approx=(-gn / 2) / nz, N_ofa=int(w.sum()))
                rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iv-tc-root", type=Path,
                    default=Path("/discover/nobackup/projects/land_da/cygl1_operator_test/iv_tc_output"))
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    keys = tile_keys()
    df = pd.concat([keys, omf_columns(keys)], axis=1)
    if (a.iv_tc_root / "step4_tc").exists():
        df = pd.concat([df, tc_iv_columns(a.iv_tc_root, keys)], axis=1)
    else:
        print(f"WARNING: {a.iv_tc_root} not found; TC and IV columns NOT written (Discover only).")
    df.to_csv(a.out / "benefit_table.csv", index=False)
    seasonal_domain().to_csv(a.out / "benefit_table_seasonal_domain.csv", index=False)
    print(df.shape, "in mask:", int(df.in_l1_mask.sum()))


if __name__ == "__main__":
    main()
