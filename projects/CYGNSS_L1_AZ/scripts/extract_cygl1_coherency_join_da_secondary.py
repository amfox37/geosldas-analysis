#!/usr/bin/env python3
"""
Step 2b (secondary/comparison check) for the CYGNSS L1 coherency follow-up.

Repeats the exact same join as extract_cygl1_coherency_join.py (same staging-file
candidate index, same GEOSldas-reader-faithful window/owner-tile/nearest-distance
matching, same QC-pass-CSV coherency lookup by exact sample_id/ch_id key) but against
the DA run DAv8_M36_all_sensors_AZ_scaled_cygl1assim instead of the OL run
OLv8_M36_all_sensors_AZ_scaled, restricted to the same Jan-Feb 2020 window.

Species index re-verified directly from THIS experiment's own obsparam (do not assume
it matches the OL run's index): CYGNSS_L1_DDM3X5_CROP_SCALAR = species 13 here too
(same as the OL run, confirmed independently -- coincidence, not an assumption), with
obsparam_assim=1 (actually assimilated, unlike the OL run's assim=0) and
obsparam_errstd=4.4 dB (vs 3.0 dB in the OL run -- these are two different experiment
configs, not the same R value).

This is SECONDARY evidence only: O-F in a run that actually assimilates the species
reflects post-increment "damage" (the analysis nudges the model state, which the next
cycle's forecast partially retains), not pure operator fit. Treat any difference from
the OL-run (primary) result with that caveat foremost in mind.

Reuses read_tilecoord / epoch_seconds / load_staging_index / step1_duplicate_check /
find_best_candidate / load_qc_pass_csv unchanged from extract_cygl1_coherency_join.py.
"""
import glob
import os
import sys
import hashlib
import subprocess

import numpy as np
import netCDF4 as nc
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_cygl1_coherency_join import (
    read_tilecoord, load_staging_index, step1_duplicate_check,
    find_best_candidate, load_qc_pass_csv, DTSTEP_ASSIM,
    L1_SPECIES_NAME, L3_SPECIES_NAME,
)

EXP_DIR = "/discover/nobackup/projects/land_da/cygl1_operator_test/DAv8_M36_all_sensors_AZ_scaled_cygl1assim"
EXP_ID = "DAv8_M36_all_sensors_AZ_scaled_cygl1assim"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
TILECOORD_FILE = os.path.join(EXP_DIR, "output", DOMAIN, "rc_out", f"{EXP_ID}.ldas_tilecoord.bin")
ANA_ROOT = os.path.join(EXP_DIR, "output", DOMAIN, "ana", "ens_avg")

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygnss_l1_coherency_joined_202001_202002_DAv8_secondary.csv.gz")


def main():
    print(f"[1] Reading tilecoord: {TILECOORD_FILE}", flush=True)
    tc = read_tilecoord(TILECOORD_FILE)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]
    ig_by_row = tc["i_indg"]
    jg_by_row = tc["j_indg"]
    print(f"    N_tile = {n_tile}")

    dates = pd.date_range("2020-01-01", "2020-02-29", freq="D")
    print(f"[2] Loading staging index for {len(dates)} days (2020-01-01 .. 2020-02-29)", flush=True)
    index, n_files, n_rows_total, n_ok_total = load_staging_index(dates)
    step1_duplicate_check(index, DTSTEP_ASSIM)

    files = sorted(glob.glob(os.path.join(ANA_ROOT, "Y2020", "M01", f"{EXP_ID}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files += sorted(glob.glob(os.path.join(ANA_ROOT, "Y2020", "M02", f"{EXP_ID}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files = sorted(files)
    print(f"[3] Found {len(files)} OFA files for Jan-Feb 2020", flush=True)

    rows = []
    n_l1_obs_seen = 0
    n_no_owner_tile_match = 0
    n_no_candidate_in_window = 0
    n_matched = 0
    n_multi_candidate_matched = 0
    species_name_to_id_seen = {}

    for i, fpath in enumerate(files):
        if i % 100 == 0:
            print(f"    ...{i}/{len(files)}", flush=True)
        fname = os.path.basename(fpath)
        stamp = fname.split(".")[-2].rstrip("z")
        yyyymmdd, hhmm = stamp.split("_")
        dt = pd.Timestamp(
            year=int(yyyymmdd[0:4]), month=int(yyyymmdd[4:6]), day=int(yyyymmdd[6:8]),
            hour=int(hhmm[0:2]), minute=int(hhmm[2:4]),
        )
        t_center = int(dt.value // 10**9)
        t_low = t_center - DTSTEP_ASSIM // 2
        t_up = t_center + DTSTEP_ASSIM // 2

        with nc.Dataset(fpath) as f:
            descr = np.array(f.variables["obsparam_descr"][:])
            spid = np.array(f.variables["obsparam_species_id"][:])
            name_to_id = {str(d): int(s) for d, s in zip(descr, spid)}
            if L1_SPECIES_NAME not in name_to_id:
                continue
            species_name_to_id_seen[L1_SPECIES_NAME] = name_to_id[L1_SPECIES_NAME]
            l1_id = name_to_id[L1_SPECIES_NAME]
            l3_id = name_to_id.get(L3_SPECIES_NAME)

            species = np.array(f.variables["species"][:])
            mask_l1 = species == l1_id
            if not np.any(mask_l1):
                continue

            tilenum_all = np.array(f.variables["tilenum"][:])
            tilenum = tilenum_all[mask_l1]
            obs = np.array(f.variables["obs"][:])[mask_l1]
            obsvar = np.array(f.variables["obsvar"][:])[mask_l1]
            fcst = np.array(f.variables["fcst"][:])[mask_l1]
            ana = np.array(f.variables["ana"][:])[mask_l1]
            assim_flag = np.array(f.variables["assim_flag"][:])[mask_l1]

            l3_fcst_by_tilenum = {}
            if l3_id is not None:
                mask_l3 = species == l3_id
                if np.any(mask_l3):
                    tn_l3 = tilenum_all[mask_l3]
                    fcst_l3 = np.array(f.variables["fcst"][:])[mask_l3]
                    for tn, fv in zip(tn_l3, fcst_l3):
                        l3_fcst_by_tilenum[int(tn)] = float(fv)

            for arr in (obs, fcst, ana):
                arr[arr > 1e14] = np.nan

            for k in range(len(tilenum)):
                n_l1_obs_seen += 1
                tn = int(tilenum[k])
                row_idx = tn - 1
                if row_idx < 0 or row_idx >= n_tile:
                    n_no_owner_tile_match += 1
                    continue
                ig = int(ig_by_row[row_idx])
                jg = int(jg_by_row[row_idx])
                tile_id = int(tile_id_by_row[row_idx])

                cand = find_best_candidate(index, ig, jg, t_low, t_up)
                if cand is None:
                    n_no_candidate_in_window += 1
                    continue
                if cand["n_candidates_in_window"] > 1:
                    n_multi_candidate_matched += 1

                qc = load_qc_pass_csv(cand["year"], cand["day"], cand["sc_num"])
                coherency_state = np.nan
                coherency_ratio = np.nan
                qc_found = False
                if qc is not None:
                    key = (cand["sample_id"], cand["ch_id"])
                    if key in qc.index:
                        qc_row = qc.loc[key]
                        coherency_state = float(qc_row["coherency_state"])
                        coherency_ratio = float(qc_row["coherency_ratio"])
                        qc_found = True

                n_matched += 1
                rows.append((
                    dt, tile_id, float(lon_by_row[row_idx]), float(lat_by_row[row_idx]),
                    float(obs[k]), float(obsvar[k]), float(fcst[k]), float(ana[k]), int(assim_flag[k]),
                    l3_fcst_by_tilenum.get(tn, np.nan),
                    cand["year"], cand["day"], cand["sc_num"], cand["sample_id"], cand["ch_id"],
                    cand["distance_km"], cand["n_candidates_in_window"], cand["obs_db"],
                    coherency_state, coherency_ratio, qc_found,
                ))

    print(f"[4] L1 obs seen in OFA files: {n_l1_obs_seen}")
    print(f"    dropped, bad tilenum: {n_no_owner_tile_match}")
    print(f"    dropped, no matching raw candidate in window: {n_no_candidate_in_window}")
    print(f"    matched to a raw candidate: {n_matched}  (of which >1 candidate competed: {n_multi_candidate_matched})")

    df = pd.DataFrame(rows, columns=[
        "datetime", "tile_id", "lon", "lat",
        "obs_db", "obsvar", "fcst_db", "ana_db", "assim_flag",
        "l3_sm6hr_fcst",
        "raw_year", "raw_day", "raw_sc_num", "raw_sample_id", "raw_ch_id",
        "raw_owner_distance_km", "raw_n_candidates_in_window", "raw_obs_db_staging",
        "coherency_state", "coherency_ratio", "qc_pass_csv_found",
    ])
    n_qc_missing = int((~df["qc_pass_csv_found"]).sum())
    print(f"[5] Rows with no QC-pass-CSV coherency match: {n_qc_missing} / {len(df)}")
    print(f"[6] assim_flag value counts (DA run -- expect some ==1, unlike the OL run):")
    print(df["assim_flag"].value_counts())

    df = df.sort_values(["datetime", "tile_id"]).reset_index(drop=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    header_lines = [
        f"# experiment_path: {EXP_DIR}",
        f"# experiment_id: {EXP_ID}",
        "# SECONDARY/comparison check only -- see extract_cygl1_coherency_join.py for the PRIMARY OL-run join.",
        "# O-F here reflects post-increment state (assim=1 for this species in this run), not pure operator fit.",
        f"# n_l1_obs_seen_in_ofa: {n_l1_obs_seen}",
        f"# n_matched_to_raw_candidate: {n_matched}",
        f"# n_multi_candidate_matched: {n_multi_candidate_matched}",
        f"# n_qc_pass_csv_coherency_missing: {n_qc_missing}",
    ]
    import gzip
    with gzip.open(CSV_PATH, "wt") as fo:
        for line in header_lines:
            fo.write(line + "\n")
        df.to_csv(fo, index=False)
    sha = hashlib.sha256()
    with open(CSV_PATH, "rb") as fb:
        for chunk in iter(lambda: fb.read(1 << 20), b""):
            sha.update(chunk)
    with open(CSV_PATH + ".sha256", "w") as fo:
        fo.write(f"{sha.hexdigest()}  {os.path.basename(CSV_PATH)}\n")
    print(f"[7] Wrote {CSV_PATH} and .sha256")
    return df


if __name__ == "__main__":
    main()
