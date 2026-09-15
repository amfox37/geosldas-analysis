#!/usr/bin/env python3
"""
Event-study test for local skill transfer from CYGNSS L1 assimilation to the
withheld monitor-only species (SMOS/SMAP Tb, ASCAT, CYGNSS_SM_6hr).

Motivation: the pooled OmF-vs-OL comparison (score_cygl1_arm.py) mixes every
monitor obs together regardless of how close in time/space it is to a CygL1
update -- if any real local transfer exists it could be diluted to the
noise-floor result seen so far. This script instead matches every monitor obs
(same tile, same DA vs OL cross-comparison) to the most recent CygL1 update at
that tile and bins DA-vs-OL skill by time-since-that-update, to see whether
skill is better right after a nearby update and decays away, vs. flat/zero at
every lag (which would argue against any local transfer channel at all).

Usage:
  cygl1_local_skill_decay.py --da-expid DAv8_M36_AZ_paired_cygl1_dense075_coh05 \\
                              --arm-tag dense075_coh05 \\
                              --start 20200101 --end 20220101
"""
import argparse
import glob
import gzip
import hashlib
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_cygl1_paired_gain_consistency import read_tilecoord  # noqa: E402

DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"

XCOMPACT_DEG = 1.25  # obs_param_nml xcompact=ycompact for this project's update_type=3d;
# a tile at (dlon,dlat) from an obs is inside that obs's admission ellipse iff
# (dlon/xcompact)^2 + (dlat/ycompact)^2 < 1 -- see thin_cygl1_nested_density_6mo.py's own
# docstring for the same ellipse, established there as the actual GEOSldas check_compact()
# admission window (raw lon/lat degree differences, no cos(lat) scaling -- matches this
# project's own established convention, not a literal-geometry correction).

CYGL1_DESCR = "CYGNSS_L1_DDM3X5_CROP_SCALAR"
TB_DESCRS = [
    "SMOS_fit_Tbh_A", "SMOS_fit_Tbh_D", "SMOS_fit_Tbv_A", "SMOS_fit_Tbv_D",
    "SMAP_L1C_Tbh_A", "SMAP_L1C_Tbh_D", "SMAP_L1C_Tbv_A", "SMAP_L1C_Tbv_D",
]
SM_DESCRS = ["ASCAT_HSAF_META_SM", "ASCAT_HSAF_METB_SM", "ASCAT_HSAF_METC_SM", "CYGNSS_SM_6hr"]
MONITOR_DESCRS = TB_DESCRS + SM_DESCRS
GROUPS = {"Tb": TB_DESCRS, "SM": SM_DESCRS}

FILL_THRESH = 1e14

LAG_BINS = [
    (0, 0.125, "same_cycle(<=3h)"),
    (0.125, 1, "<1day"),
    (1, 3, "1-3day"),
    (3, 7, "3-7day"),
    (7, 14, "7-14day"),
    (14, np.inf, ">14day"),
]


def month_range(start, end):
    ym = []
    y, m = start.year, start.month
    while (y, m) < (end.year, end.month):
        ym.append((y, m))
        m += 1
        if m == 13:
            m = 1
            y += 1
    return ym


def read_ofa_dir(exp_id, start, end, want_cygl1_events, want_monitor_obs):
    """One pass per file: optionally pull CygL1 assim events and/or monitor obs/fcst."""
    exp_dir = os.path.join(EXP_ROOT, exp_id)
    tilecoord_file = os.path.join(exp_dir, "output", DOMAIN, "rc_out", f"{exp_id}.ldas_tilecoord.bin")
    ana_root = os.path.join(exp_dir, "output", DOMAIN, "ana", "ens_avg")
    tc = read_tilecoord(tilecoord_file)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]
    tile_lonlat = {int(tile_id_by_row[i]): (float(lon_by_row[i]), float(lat_by_row[i])) for i in range(n_tile)}

    files = []
    for yy, mm in month_range(start, end):
        files.extend(glob.glob(os.path.join(ana_root, f"Y{yy}", f"M{mm:02d}", f"{exp_id}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files = sorted(files)
    print(f"[{exp_id}] {len(files)} OFA files", flush=True)
    if not files:
        print(f"ERROR: no OFA files for {exp_id}", file=sys.stderr)
        sys.exit(1)

    cygl1_rows = []
    monitor_rows = []

    for i, fpath in enumerate(files):
        if i % 500 == 0:
            print(f"    [{exp_id}] ...{i}/{len(files)}", flush=True)
        fname = os.path.basename(fpath)
        stamp = fname.split(".")[-2].rstrip("z")
        yyyymmdd, hhmm = stamp.split("_")
        dt = pd.Timestamp(year=int(yyyymmdd[0:4]), month=int(yyyymmdd[4:6]), day=int(yyyymmdd[6:8]),
                           hour=int(hhmm[0:2]), minute=int(hhmm[2:4]))

        with nc.Dataset(fpath) as f:
            descr = np.array(f.variables["obsparam_descr"][:])
            spid = np.array(f.variables["obsparam_species_id"][:])
            name_to_id = {str(d): int(s) for d, s in zip(descr, spid)}
            species = np.array(f.variables["species"][:])
            tilenum = np.array(f.variables["tilenum"][:])
            obs = np.array(f.variables["obs"][:])
            fcst = np.array(f.variables["fcst"][:])
            assim_flag = np.array(f.variables["assim_flag"][:])

            if want_cygl1_events and CYGL1_DESCR in name_to_id:
                wanted_id = name_to_id[CYGL1_DESCR]
                mask = (species == wanted_id) & (assim_flag == 1)
                if np.any(mask):
                    tn = tilenum[mask]
                    ov = obs[mask]
                    row_idx = np.clip(tn - 1, 0, n_tile - 1)
                    bad = (tn - 1 < 0) | (tn - 1 >= n_tile) | (np.abs(ov) > FILL_THRESH)
                    tid = tile_id_by_row[row_idx]
                    for k in range(len(tn)):
                        if bad[k]:
                            continue
                        cygl1_rows.append((int(tid[k]), dt))

            if want_monitor_obs:
                for descr_name in MONITOR_DESCRS:
                    if descr_name not in name_to_id:
                        continue
                    wanted_id = name_to_id[descr_name]
                    mask = species == wanted_id
                    if not np.any(mask):
                        continue
                    tn = tilenum[mask]
                    ov = obs[mask]
                    fv = fcst[mask]
                    row_idx = np.clip(tn - 1, 0, n_tile - 1)
                    bad = (tn - 1 < 0) | (tn - 1 >= n_tile) | (np.abs(ov) > FILL_THRESH) | (np.abs(fv) > FILL_THRESH)
                    tid = tile_id_by_row[row_idx]
                    for k in range(len(tn)):
                        if bad[k]:
                            continue
                        monitor_rows.append((descr_name, int(tid[k]), dt, float(ov[k]), float(fv[k])))

    cygl1_df = pd.DataFrame(cygl1_rows, columns=["tile_id", "datetime"])
    monitor_df = pd.DataFrame(monitor_rows, columns=["species", "tile_id", "datetime", "obs", "fcst"])
    return cygl1_df, monitor_df, tile_lonlat


def build_neighbor_expanded_updates(cygl1_df, tile_lonlat, monitor_tile_ids, radius_deg):
    """For every monitor tile, the sorted union of CygL1 update times from itself and every
    tile whose centroid falls inside its xcompact/ycompact admission ellipse -- i.e. the same
    ellipse a CygL1 obs at that tile would actually use to influence neighboring tiles in
    GEOSldas's own local update (see XCOMPACT_DEG comment above)."""
    cygl1_tile_ids = sorted(cygl1_df["tile_id"].unique().tolist())
    monitor_tile_ids = sorted(set(monitor_tile_ids))
    print(f"Building neighbor-expanded update lookup: {len(monitor_tile_ids)} monitor tiles x "
          f"{len(cygl1_tile_ids)} CygL1-update tiles, radius={radius_deg}deg", flush=True)

    cyg_lonlat = np.array([tile_lonlat[t] for t in cygl1_tile_ids])  # (Nc, 2)
    mon_lonlat = np.array([tile_lonlat[t] for t in monitor_tile_ids])  # (Nm, 2)

    dlon = mon_lonlat[:, 0:1] - cyg_lonlat[None, :, 0]  # (Nm, Nc)
    dlat = mon_lonlat[:, 1:2] - cyg_lonlat[None, :, 1]
    inside = (dlon / radius_deg) ** 2 + (dlat / radius_deg) ** 2 < 1.0

    times_by_cygtile = {tid: g.sort_values("datetime")["datetime"].values.astype("datetime64[ns]")
                         for tid, g in cygl1_df.groupby("tile_id")}

    updates_by_monitor_tile = {}
    for mi, mtile in enumerate(monitor_tile_ids):
        neighbor_idx = np.where(inside[mi])[0]
        if len(neighbor_idx) == 0:
            continue
        arrs = [times_by_cygtile[cygl1_tile_ids[ci]] for ci in neighbor_idx]
        merged = np.sort(np.concatenate(arrs))
        updates_by_monitor_tile[mtile] = merged
    return updates_by_monitor_tile


def assign_lag_bin(hours):
    if pd.isna(hours):
        return "never"
    days = hours / 24.0
    for lo, hi, label in LAG_BINS:
        if lo <= days < hi:
            return label
    return "never"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--da-expid", default="DAv8_M36_AZ_paired_cygl1_dense075_coh05")
    ap.add_argument("--ol-expid", default="OLv8_M36_AZ_paired_monitor")
    ap.add_argument("--arm-tag", default="dense075_coh05")
    ap.add_argument("--start", default="20200101")
    ap.add_argument("--end", default="20220101")
    ap.add_argument("--radius-deg", type=float, default=XCOMPACT_DEG,
                     help="admission-ellipse radius (xcompact=ycompact) for 'was there a "
                          "nearby CygL1 update', default = this project's actual localization setting")
    args = ap.parse_args()

    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)

    print(f"=== DA arm: {args.da_expid} ===")
    cygl1_df, da_monitor_df, tile_lonlat = read_ofa_dir(args.da_expid, start, end, want_cygl1_events=True, want_monitor_obs=True)
    print(f"CygL1 assim events: {len(cygl1_df)}; DA monitor obs rows: {len(da_monitor_df)}")

    print(f"=== OL arm: {args.ol_expid} ===")
    _, ol_monitor_df, _ = read_ofa_dir(args.ol_expid, start, end, want_cygl1_events=False, want_monitor_obs=True)
    print(f"OL monitor obs rows: {len(ol_monitor_df)}")

    # join DA and OL monitor obs on (species, tile_id, datetime)
    m = da_monitor_df.merge(ol_monitor_df, on=["species", "tile_id", "datetime"], suffixes=("_da", "_ol"))
    print(f"Matched DA/OL monitor events: {len(m)} (DA had {len(da_monitor_df)}, OL had {len(ol_monitor_df)})")
    m["of_da"] = m["obs_da"] - m["fcst_da"]
    m["of_ol"] = m["obs_ol"] - m["fcst_ol"]

    # neighbor-expanded (within radius_deg admission ellipse) sorted CygL1 update times per
    # monitor tile, for searchsorted lookup -- this is the actual localization-window test,
    # not just same-tile
    updates_by_tile = build_neighbor_expanded_updates(cygl1_df, tile_lonlat, m["tile_id"].unique(), args.radius_deg)

    m_dt = m["datetime"].values.astype("datetime64[ns]")
    lag_hours = np.full(len(m), np.nan)
    for tid, g_idx in m.groupby("tile_id").groups.items():
        upd_times = updates_by_tile.get(tid)
        if upd_times is None or len(upd_times) == 0:
            continue
        idx_arr = np.array(g_idx)
        ev_times = m_dt[m.index.get_indexer(idx_arr)]
        pos = np.searchsorted(upd_times, ev_times, side="right") - 1
        valid = pos >= 0
        diffs = np.full(len(ev_times), np.nan)
        diffs[valid] = (ev_times[valid] - upd_times[pos[valid]]) / np.timedelta64(1, "h")
        lag_hours[m.index.get_indexer(idx_arr)] = diffs

    m["lag_hours"] = lag_hours
    m["lag_bin"] = m["lag_hours"].apply(assign_lag_bin)
    m["group"] = m["species"].apply(lambda s: "Tb" if s in TB_DESCRS else "SM")

    bin_order = [b[2] for b in LAG_BINS] + ["never"]
    m["lag_bin"] = pd.Categorical(m["lag_bin"], categories=bin_order, ordered=True)

    lines = []
    lines.append(f"CYGNSS L1 local-skill-decay event study: {args.arm_tag} vs OL")
    lines.append(f"period: {args.start}-{args.end}, N matched monitor events = {len(m)}")
    lines.append("")
    header = f"{'lag_bin':20s}{'group':6s}{'N':>10s}{'OmF_stdv(DA)':>14s}{'OmF_stdv(OL)':>14s}{'%vsOL':>10s}"
    lines.append(header)
    for group in ["Tb", "SM"]:
        for lag_bin in bin_order:
            sub = m[(m["group"] == group) & (m["lag_bin"] == lag_bin)]
            if len(sub) < 20:
                lines.append(f"{lag_bin:20s}{group:6s}{len(sub):>10d}{'--':>14s}{'--':>14s}{'--':>10s}")
                continue
            s_da = sub["of_da"].std()
            s_ol = sub["of_ol"].std()
            pct = 100 * (s_da - s_ol) / s_ol
            lines.append(f"{lag_bin:20s}{group:6s}{len(sub):>10d}{s_da:>14.4f}{s_ol:>14.4f}{pct:>+10.2f}")
        lines.append("")

    summary_text = "\n".join(lines)
    print()
    print(summary_text)

    os.makedirs(OUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUT_DIR, f"cygl1_local_skill_decay_{args.arm_tag}.csv.gz")
    summary_path = os.path.join(OUT_DIR, f"cygl1_local_skill_decay_{args.arm_tag}_summary.txt")
    with gzip.open(csv_path, "wt") as fo:
        fo.write(f"# {args.da_expid} vs {args.ol_expid}, {args.start}-{args.end}\n")
        m.to_csv(fo, index=False)
    with open(summary_path, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {csv_path} and {summary_path}")


if __name__ == "__main__":
    main()
