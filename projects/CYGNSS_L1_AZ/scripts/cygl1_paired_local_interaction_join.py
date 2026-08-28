#!/usr/bin/env python3
"""
CYGNSS L1 paired thinning-density experiment: compute the per-observation
local-interaction-count covariate for all three DA arms (sparse/intermediate/
dense) and join it onto tile membership, producing a lookup table
(arm, datetime, tilenum, local_interaction_count) that
cygl1_paired_event_diagnostics.py left-joins onto the per-event OFA table
from extract_cygl1_paired_gain_consistency.py.

Reuses (imports, does not duplicate) load_all_obs / build_sparse /
build_intermediate_candidate / local_interaction_counts from
thin_cygl1_nested_density_6mo.py, which already runs over the full
BEG_DATE=20200101 / END_DATE=20201031 range (module constants) -- necessary
because a May obs's local-interaction count depends on its window-mates
regardless of when the window falls. We compute over the full Jan-Oct range
but only EMIT the May-Oct (Y2020 M05-M10) subset of kept obs, matching the
evaluation window used by the event-diagnostics script.

DA-sparse:       kept_mask = build_sparse(df) (min_sep=5.0deg) -- trivially 0
                  local-interaction by construction; recomputed anyway as a
                  cheap double-check, cross-checked against
                  cygl1_sparse_6mo_thinning_log.csv (3515 rows).
DA-intermediate:  kept_mask = build_intermediate_candidate(df, sparse_mask, 2.40)
                  -- cross-checked against cygl1_intermediate_6mo_thinning_log.csv
                  (8789 rows).
DA-dense:         kept_mask = all True (full unthinned stream, 181,630 obs over
                  Jan-Oct) -- the expensive O(n_per_window^2) case; this is the
                  single most important number for the "does DA-dense reproduce
                  the R-sweep-style failure" question.

Tile-membership join: for each retained obs (lon,lat,anchor-as-datetime), find
tilenum by bbox lookup (min_lon<=lon<=max_lon and min_lat<=lat<=max_lat) into
the tilecoord arrays, 0-based row index -> tilenum = row_index+1. Any
tilecoord file works (all 4 arms' tilecoord.bin are byte-identical size /
same AZ grid) -- uses the sparse arm's file.

anchor (unix epoch seconds) -> datetime via pd.Timestamp(anchor, unit='s'),
verified to match OFA filename cycle stamps exactly.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_cygl1_nested_density_6mo as thin  # noqa: E402
from extract_cygl1_paired_gain_consistency import read_tilecoord  # noqa: E402

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
LOOKUP_CSV = os.path.join(OUT_DIR, "cygl1_paired_local_interaction_lookup.csv")

TILECOORD_FILE = ("/discover/nobackup/projects/land_da/cygl1_operator_test/"
                   "DAv8_M36_AZ_paired_cygl1_sparse/output/SMAP_EASEv2_M36_GLOBAL/rc_out/"
                   "DAv8_M36_AZ_paired_cygl1_sparse.ldas_tilecoord.bin")

SPARSE_LOG_CSV = os.path.join(OUT_DIR, "cygl1_sparse_6mo_thinning_log.csv")
INTERMEDIATE_LOG_CSV = os.path.join(OUT_DIR, "cygl1_intermediate_6mo_thinning_log.csv")

EVAL_MONTHS = {5, 6, 7, 8, 9, 10}


def tilenum_for_points(lon, lat, tc):
    """Vectorized bbox lookup: for each (lon,lat), find the tilecoord row whose
    bounding box contains it. Returns 1-based tilenum, or -1 if no match /
    ambiguous match resolved to first hit."""
    min_lon, max_lon = tc["min_lon"], tc["max_lon"]
    min_lat, max_lat = tc["min_lat"], tc["max_lat"]
    n = len(lon)
    tilenum = np.full(n, -1, dtype=np.int64)
    n_multi = 0
    for i in range(n):
        hits = np.where((lon[i] >= min_lon) & (lon[i] <= max_lon) &
                         (lat[i] >= min_lat) & (lat[i] <= max_lat))[0]
        if len(hits) == 0:
            continue
        if len(hits) > 1:
            n_multi += 1
        tilenum[i] = hits[0] + 1
        if i % 20000 == 0:
            print(f"    tile lookup {i}/{n}", flush=True)
    print(f"  tile lookup done: {n} points, {int((tilenum < 0).sum())} unmatched, "
          f"{n_multi} matched >1 tile (used first hit)", flush=True)
    return tilenum


def emit_rows(df, kept_mask, counts, arm, tc):
    kept_df = df[kept_mask].copy()
    kept_df["local_interaction_count"] = counts[kept_mask]
    kept_df["datetime"] = pd.to_datetime(kept_df["anchor"], unit="s")
    eval_df = kept_df[kept_df["datetime"].dt.month.isin(EVAL_MONTHS)].copy()
    print(f"[{arm}] kept={kept_mask.sum()}, in May-Oct eval window={len(eval_df)}", flush=True)

    print(f"[{arm}] resolving tilenum for {len(eval_df)} May-Oct obs...", flush=True)
    eval_df["tilenum"] = tilenum_for_points(eval_df["lon"].values, eval_df["lat"].values, tc)
    eval_df["arm"] = arm
    return eval_df[["arm", "datetime", "tilenum", "local_interaction_count",
                     "file_idx", "obs_idx", "lon", "lat", "anchor"]]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Reading tilecoord for tile-membership lookup...", flush=True)
    tc = read_tilecoord(TILECOORD_FILE)
    print(f"  N_tile = {tc['N_tile']}", flush=True)

    print("Loading full Jan-Oct 2020 CYGNSS L1 obs stream (thin_cygl1_nested_density_6mo.load_all_obs)...",
          flush=True)
    df, src_paths, dates = thin.load_all_obs()

    print("Building DA-sparse mask (min_sep=5.0deg)...", flush=True)
    sparse_mask = thin.build_sparse(df)
    print(f"  DA-sparse kept: {sparse_mask.sum()} (expect 3515)", flush=True)

    print("Building DA-intermediate mask (min_sep=2.40deg, nested superset of sparse)...", flush=True)
    inter_mask = thin.build_intermediate_candidate(df, sparse_mask, 2.40)
    print(f"  DA-intermediate kept: {inter_mask.sum()} (expect 8789)", flush=True)

    dense_mask = np.ones(len(df), dtype=bool)
    print(f"  DA-dense kept: {dense_mask.sum()} (all obs, full stream)", flush=True)

    print("Computing local_interaction_counts for DA-sparse (expect trivially all 0)...", flush=True)
    sparse_counts = thin.local_interaction_counts(df, sparse_mask)
    sc = sparse_counts[sparse_mask]
    print(f"  DA-sparse local_interaction_count: min={sc.min()}, max={sc.max()}, "
          f"all_zero={bool(np.all(sc == 0))}", flush=True)

    print("Computing local_interaction_counts for DA-intermediate...", flush=True)
    inter_counts = thin.local_interaction_counts(df, inter_mask)
    ic = inter_counts[inter_mask]
    print(f"  DA-intermediate local_interaction_count: median={np.median(ic):.2f}, "
          f"p90={np.percentile(ic,90):.2f}, max={ic.max()}", flush=True)

    print("Computing local_interaction_counts for DA-dense (expensive, full unthinned stream "
          f"{dense_mask.sum()} obs -- this is the single most important number for the "
          "dose-response question; may take a while)...", flush=True)
    dense_counts = thin.local_interaction_counts(df, dense_mask)
    dc = dense_counts[dense_mask]
    print(f"  DA-dense local_interaction_count: median={np.median(dc):.2f}, "
          f"p90={np.percentile(dc,90):.2f}, max={dc.max()}", flush=True)

    # --- cross-checks against existing thinning logs ---
    print("\nCross-checking against existing thinning-log CSVs...", flush=True)
    sparse_log = pd.read_csv(SPARSE_LOG_CSV)
    sparse_kept_keys = set(zip(df.loc[sparse_mask, "file_idx"], df.loc[sparse_mask, "obs_idx"]))
    log_keys = set(zip(sparse_log["file_idx"], sparse_log["obs_idx"]))
    print(f"  sparse: recomputed kept N={len(sparse_kept_keys)} vs log N={len(log_keys)}, "
          f"identical_set={sparse_kept_keys == log_keys}", flush=True)

    inter_log = pd.read_csv(INTERMEDIATE_LOG_CSV)
    inter_kept_keys = set(zip(df.loc[inter_mask, "file_idx"], df.loc[inter_mask, "obs_idx"]))
    log_keys2 = set(zip(inter_log["file_idx"], inter_log["obs_idx"]))
    print(f"  intermediate: recomputed kept N={len(inter_kept_keys)} vs log N={len(log_keys2)}, "
          f"identical_set={inter_kept_keys == log_keys2}", flush=True)

    # --- emit May-Oct lookup rows per arm ---
    out_dfs = []
    out_dfs.append(emit_rows(df, sparse_mask, sparse_counts, "sparse", tc))
    out_dfs.append(emit_rows(df, inter_mask, inter_counts, "intermediate", tc))
    out_dfs.append(emit_rows(df, dense_mask, dense_counts, "dense", tc))

    lookup = pd.concat(out_dfs, ignore_index=True)
    n_unmatched = int((lookup["tilenum"] < 0).sum())
    print(f"\nTotal lookup rows: {len(lookup)}; unmatched tile lookups: {n_unmatched}", flush=True)

    lookup.to_csv(LOOKUP_CSV, index=False)
    print(f"Wrote {LOOKUP_CSV}")
    lookup.to_pickle(os.path.join(OUT_DIR, "_cygl1_paired_local_interaction_lookup.pkl"))
    return lookup


if __name__ == "__main__":
    main()
