#!/usr/bin/env python3
"""
Join CYGNSS L1 (CYGNSS_L1_DDM3X5_CROP_SCALAR) OFA obs/fcst/ana records from the
OLv8_M36_all_sensors_AZ_scaled open-loop run to the raw-granule `coherency_state`/
`coherency_ratio` SDS fields, for Jan-Feb 2020. Follow-up to the L1-vs-L3 operator
diagnosis at projects/CYGNSS_L1_AZ/runs/cygl1_operator_diagnosis.md: does the native
CYGNSS L1 coherency flag predict where the crude flat-Fresnel-reflectivity operator
(mwRTM_get_lr_reflectivity, no roughness/vegetation correction) fits well vs. poorly?

--------------------------------------------------------------------------------
Step 0 finding (see also runs/cygl1_coherency_stratification.md):
    The GEOSldas-facing staging files at
    /discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/Y2020/M0{1,2}/
    cygnss_l1_ddm3x5_crop_scalar_m36_<date>_all_cyg.nc4
    do NOT carry coherency_state/coherency_ratio (confirmed via ncdump -h: only
    sample_id/ch_id/sc_num/year/day/observed_y_db/sp_nearest_tile_*/tile_ig/tile_jg/
    coefficient support arrays -- no coherency fields, no product_json payload either,
    that field is empty-string for all ok rows). So this is NOT a one-step join; we
    must go through the earlier-pipeline per-(day,spacecraft) QC-pass CSVs at
    /discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images/
    cygnss_qc_m36_window_counts_<YYYYMMDD>_cyg<NN>/cygnss_l1_qc_pass_<YYYYMMDD>_cyg<NN>.csv
    which DO carry coherency_state/coherency_ratio, keyed by (sample_id, ch_id) and
    unique on that key within each (day, spacecraft) file (checked; 0 duplicates in a
    spot check of 20200101_cyg01, N=51095 rows).

Step 1 finding: the staging "all_cyg" daily merge files are NOT one-obs-per-tile-window.
    Multiple raw candidate obs (different spacecraft, different sample_id/ch_id) CAN
    map to the same owner tile within the same 3-h assimilation window (this is why
    GEOSldas's own reader, read_obs_cygnss_l1_scalar() in clsm_ensupd_read_obs.F90,
    line ~3606-3623, has an explicit "keep the candidate with min sp_nearest_tile_
    distance_km, drop the rest" dedup step, with an N_duplicate counter). So a
    (tile_ig,tile_jg,3h-window) key is many-to-one on the raw-sample side, and there is
    NO single coherency_state that belongs to a staged M36 record "for free" -- we must
    reproduce GEOSldas's own tie-break (nearest specular point to tile center) to know
    WHICH raw sample (and therefore which coherency_state) is the one actually used by
    the forward operator for a given OFA (tile_id, window) row. See empirical duplicate
    count distribution printed by this script (search "[Step1]" in stdout) and recorded
    in the markdown report.

Join logic (exactly reproduces read_obs_cygnss_l1_scalar(), clsm_ensupd_read_obs.F90):
  1. OFA file timestamp T (e.g. ..._ObsFcstAna.20200101_0300z.nc4) is the assimilation
     window CENTER. dtstep_assim = 10800 s (3 h; confirmed by 3-hourly OFA cadence and
     by the dtstep_assim(4)=10800 entry in clsm_bias_routines.F90). Window is the
     half-open interval (T - 5400s, T + 5400s].
  2. Owner tile = the local tile whose (i_indg, j_indg) tilecoord entry matches a raw
     obs's (sp_nearest_tile_ig, sp_nearest_tile_jg). tilecoord read the same way as
     extract_cygl1_l3_pairs.py (read_tilecoord()); OFA tilenum is 1-based, so tilenum-1
     indexes the 0-based tilecoord arrays -- same convention as that script.
  3. Raw-obs absolute timestamp = Jan-1-00:00:00 of the raw obs's own `year` field, plus
     (day-1)*86400 + ddm_timestamp_utc_sec seconds -- using the per-row year/day fields
     from the staging file, NOT the filename's target date (the source code comment notes
     the daily merge folds in the previous day's 22:30-24:00 UTC tail, so a row's true
     origin day can differ from the file's target_date_utc).
  4. Among all staging rows (status==1) in a given day's "all_cyg" file whose
     (sp_nearest_tile_ig, sp_nearest_tile_jg) matches the OFA row's owner tile and whose
     absolute timestamp falls in the window, keep the one with the SMALLEST
     sp_nearest_tile_distance_km (ties broken by first-encountered, as in the Fortran
     do-loop). This reproduces GEOSldas's own duplicate-owner-tile resolution exactly.
  5. That row's (year, day, sc_num, sample_id, ch_id) is looked up directly (exact,
     non-fuzzy) against the (day, spacecraft) QC-pass CSV to fetch coherency_state /
     coherency_ratio. No value-matching / floating-point join is used anywhere in this
     script -- the entire chain is integer/exact-key joins (tile ig/jg, sample_id, ch_id,
     year, day, sc_num), specifically to avoid the ~10% repeated-obs-value ambiguity
     flagged for this dataset.
  6. Windows never straddle a day boundary for Jan 1 - Feb 29 2020 inclusive (checked:
     first window is 20200101_0300z, i.e. low bound 20200101_0130z; last window used is
     20200229_2100z, i.e. up bound 20200229_2230z) -- so only the single matching day's
     staging file is ever needed as the candidate pool for a given OFA row.

Species resolved by NAME from each OFA file's own obsparam_descr/obsparam_species_id
(never a hardcoded index), per project convention. CYGNSS_SM_6hr fcst is also pulled
(same tile/window) for soil-moisture-equivalent context (Step 4 "deliverable" ask),
reusing the same tile-id resolution as extract_cygl1_l3_pairs.py.

Output: gzip CSV with one row per joined CYGNSS_L1 OFA observation, plus a
diagnostics/dropped-record summary printed to stdout and saved as a header block.
"""
import glob
import os
import struct
import sys
import hashlib
import subprocess
from collections import defaultdict

import numpy as np
import netCDF4 as nc
import pandas as pd

EXP_DIR = "/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_all_sensors_AZ_scaled"
EXP_ID = "OLv8_M36_all_sensors_AZ_scaled"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
TILECOORD_FILE = os.path.join(EXP_DIR, "output", DOMAIN, "rc_out", f"{EXP_ID}.ldas_tilecoord.bin")
ANA_ROOT = os.path.join(EXP_DIR, "output", DOMAIN, "ana", "ens_avg")

STAGING_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1"
QC_PASS_ROOT = "/discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images"

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygnss_l1_coherency_joined_202001_202002.csv.gz")

DTSTEP_ASSIM = 10800  # seconds; 3-h assimilation cadence (see docstring)

L1_SPECIES_NAME = "CYGNSS_L1_DDM3X5_CROP_SCALAR"
L3_SPECIES_NAME = "CYGNSS_SM_6hr"


def read_tilecoord(fname):
    machfmt = "<"
    tile_coord = {}
    with open(fname, "rb") as ifp:
        _ = struct.unpack(f"{machfmt}i", ifp.read(4))[0]
        tile_coord["N_tile"] = struct.unpack(f"{machfmt}i", ifp.read(4))[0]
        _ = struct.unpack(f"{machfmt}i", ifp.read(4))[0]

        Nt = tile_coord["N_tile"]
        fields = [
            "tile_id", "typ", "pfaf", "com_lon", "com_lat", "min_lon", "max_lon",
            "min_lat", "max_lat", "i_indg", "j_indg", "frac_cell", "frac_pfaf",
            "area", "elev",
        ]
        for field in fields:
            this_dtype = "i" if field in ["tile_id", "typ", "pfaf", "i_indg", "j_indg"] else "f"
            _ = struct.unpack(f"{machfmt}i", ifp.read(4))[0]
            tile_coord[field] = np.frombuffer(ifp.read(Nt * 4), dtype=f"{machfmt}{this_dtype}").copy()
            _ = struct.unpack(f"{machfmt}i", ifp.read(4))[0]
    return tile_coord


def epoch_seconds(year, day, sec):
    """Seconds since <year>-01-01 00:00:00 UTC, given 1-based day-of-year and
    intra-day seconds -- matches the Fortran augment_date_time((day-1)*86400+sec,
    Jan-1-of-year) construction in read_obs_cygnss_l1_scalar()."""
    base = np.datetime64(f"{int(year):04d}-01-01T00:00:00")
    return (base + np.timedelta64(1, "s") * (int(round((day - 1) * 86400 + sec)))).astype("datetime64[s]").astype(np.int64)


def load_staging_index(dates):
    """Load all Jan-Feb 2020 staging 'all_cyg' files into one flat index of
    candidate raw obs, keyed for fast (ig,jg) lookup. Returns a dict
    (ig,jg) -> list of row-dicts."""
    index = defaultdict(list)
    n_files = 0
    n_rows_total = 0
    n_ok_total = 0
    for d in dates:
        ymd = d.strftime("%Y%m%d")
        yyyy, mm = d.strftime("%Y"), d.strftime("%m")
        fpath = os.path.join(STAGING_ROOT, f"Y{yyyy}", f"M{mm}", f"cygnss_l1_ddm3x5_crop_scalar_m36_{ymd}_all_cyg.nc4")
        if not os.path.exists(fpath):
            print(f"WARNING: missing staging file {fpath}", file=sys.stderr)
            continue
        n_files += 1
        with nc.Dataset(fpath) as f:
            status = np.array(f.variables["status"][:])
            n_rows_total += len(status)
            ok = status == 1
            n_ok_total += int(ok.sum())
            sp_ig = np.array(f.variables["sp_nearest_tile_ig"][:])[ok]
            sp_jg = np.array(f.variables["sp_nearest_tile_jg"][:])[ok]
            dist = np.array(f.variables["sp_nearest_tile_distance_km"][:])[ok]
            year = np.array(f.variables["year"][:])[ok]
            day = np.array(f.variables["day"][:])[ok]
            sc_num = np.array(f.variables["sc_num"][:])[ok]
            sample_id = np.array(f.variables["sample_id"][:])[ok]
            ch_id = np.array(f.variables["ch_id"][:])[ok]
            ddm_sec = np.array(f.variables["ddm_timestamp_utc_sec"][:])[ok]
            obs_db = np.array(f.variables["observed_y_db"][:])[ok]

        epoch = np.array([epoch_seconds(y, dd, s) for y, dd, s in zip(year, day, ddm_sec)])

        for i in range(len(sp_ig)):
            index[(int(sp_ig[i]), int(sp_jg[i]))].append({
                "epoch": int(epoch[i]),
                "distance_km": float(dist[i]),
                "year": int(year[i]),
                "day": int(day[i]),
                "sc_num": int(sc_num[i]),
                "sample_id": int(sample_id[i]),
                "ch_id": int(ch_id[i]),
                "obs_db": float(obs_db[i]),
                "staging_file_date": ymd,
            })

    print(f"[Staging] loaded {n_files} daily files, {n_rows_total} total rows, {n_ok_total} status==1 rows")
    return index, n_files, n_rows_total, n_ok_total


def step1_duplicate_check(index, dtstep_assim):
    """Empirical many-to-one check (Step 1): for each (ig,jg) group, bin candidate
    obs into non-overlapping dtstep_assim-wide windows anchored at the observed
    epoch (approximate binning: floor to dtstep_assim grid) and count how many
    candidates land in the same (ig,jg,window) bin."""
    bin_counts = defaultdict(int)
    for (ig, jg), rows in index.items():
        for r in rows:
            wbin = r["epoch"] // dtstep_assim
            bin_counts[(ig, jg, wbin)] += 1
    counts = np.array(list(bin_counts.values()))
    dist = pd.Series(counts).value_counts().sort_index()
    n_multi = int((counts > 1).sum())
    n_total_bins = len(counts)
    print(f"[Step1] (tile,window)-bins with >1 candidate raw obs: {n_multi} / {n_total_bins} "
          f"({100.0*n_multi/max(n_total_bins,1):.2f}%)")
    print(f"[Step1] distribution of candidates-per-(tile,window)-bin:\n{dist}")
    return dist, n_multi, n_total_bins


def find_best_candidate(index, ig, jg, t_low, t_up):
    rows = index.get((ig, jg))
    if not rows:
        return None
    best = None
    best_dist = None
    n_candidates = 0
    for r in rows:
        if t_low < r["epoch"] <= t_up:
            n_candidates += 1
            if best_dist is None or r["distance_km"] < best_dist:
                best = r
                best_dist = r["distance_km"]
    if best is not None:
        best = dict(best)
        best["n_candidates_in_window"] = n_candidates
    return best


_qc_cache = {}


def load_qc_pass_csv(year, day, sc_num):
    key = (year, day, sc_num)
    if key in _qc_cache:
        return _qc_cache[key]
    dt = np.datetime64(f"{year:04d}-01-01") + np.timedelta64(day - 1, "D")
    ymd = pd.Timestamp(dt).strftime("%Y%m%d")
    fdir = os.path.join(QC_PASS_ROOT, f"cygnss_qc_m36_window_counts_{ymd}_cyg{sc_num:02d}")
    fpath = os.path.join(fdir, f"cygnss_l1_qc_pass_{ymd}_cyg{sc_num:02d}.csv")
    if not os.path.exists(fpath):
        _qc_cache[key] = None
        return None
    df = pd.read_csv(fpath, usecols=["sample_id", "ch_id", "coherency_state", "coherency_ratio"])
    dup = df.duplicated(subset=["sample_id", "ch_id"]).sum()
    if dup > 0:
        print(f"WARNING: {dup} duplicate (sample_id,ch_id) rows in {fpath}; keeping first", file=sys.stderr)
        df = df.drop_duplicates(subset=["sample_id", "ch_id"], keep="first")
    df = df.set_index(["sample_id", "ch_id"])
    _qc_cache[key] = df
    return df


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

    step1_dist, n_multi, n_total_bins = step1_duplicate_check(index, DTSTEP_ASSIM)

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
        t_center = int(dt.value // 10**9)  # unix epoch seconds
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
            if l3_id is not None:
                species_name_to_id_seen[L3_SPECIES_NAME] = l3_id

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

            # co-located L3 fcst (Step 4 context), keyed by tilenum within this file
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

    # Sanity: staging obs_db (used to select the winning candidate) should equal
    # the OFA-side obs_db up to the fixed z-score scaling offset applied by
    # scale_obs_cygl1scal_zscore -- NOT expected to be numerically identical.
    # We do NOT rely on this for the join (see docstring); reported only as a
    # secondary consistency signal.
    df["obs_minus_raw_obs_db"] = df["obs_db"] - df["raw_obs_db_staging"]

    df = df.sort_values(["datetime", "tile_id"]).reset_index(drop=True)

    print("[6] Row counts by coherency_state:")
    print(df["coherency_state"].value_counts(dropna=False).sort_index())

    os.makedirs(OUT_DIR, exist_ok=True)
    header_lines = [
        f"# experiment_path: {EXP_DIR}",
        f"# experiment_id: {EXP_ID}",
        f"# window: Jan 1 - Feb 29 2020 (OFA file timestamps)",
        f"# ofa_species: {L1_SPECIES_NAME} (species_id resolved per-file, seen={species_name_to_id_seen.get(L1_SPECIES_NAME)})",
        f"# l3_context_species: {L3_SPECIES_NAME} (species_id seen={species_name_to_id_seen.get(L3_SPECIES_NAME)})",
        f"# staging_root: {STAGING_ROOT}",
        f"# qc_pass_root: {QC_PASS_ROOT}",
        f"# dtstep_assim_seconds: {DTSTEP_ASSIM}",
        f"# n_l1_obs_seen_in_ofa: {n_l1_obs_seen}",
        f"# n_dropped_bad_tilenum: {n_no_owner_tile_match}",
        f"# n_dropped_no_candidate_in_window: {n_no_candidate_in_window}",
        f"# n_matched_to_raw_candidate: {n_matched}",
        f"# n_multi_candidate_matched (owner-tile duplicate resolved by nearest distance): {n_multi_candidate_matched}",
        f"# n_qc_pass_csv_coherency_missing: {n_qc_missing}",
        "# join logic: see script docstring (exact tile ig/jg + sample_id/ch_id keys, no value-matching)",
        "# columns: datetime,tile_id,lon,lat,obs_db,obsvar,fcst_db,ana_db,assim_flag,l3_sm6hr_fcst,"
        "raw_year,raw_day,raw_sc_num,raw_sample_id,raw_ch_id,raw_owner_distance_km,"
        "raw_n_candidates_in_window,raw_obs_db_staging,coherency_state,coherency_ratio,"
        "qc_pass_csv_found,obs_minus_raw_obs_db",
    ]
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")

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
