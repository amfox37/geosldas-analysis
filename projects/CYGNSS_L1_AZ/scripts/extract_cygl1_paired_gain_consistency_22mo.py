#!/usr/bin/env python3
"""
22-month extension of extract_cygl1_paired_gain_consistency.py, restricted to
the ONE arm that was actually extended (DA-intermediate; sparse/dense stopped
at Oct 2020). Evaluation window widened from May-Oct 2020 to May 2020 - Dec
2021 (Jan-Apr 2020 remains excluded as spin-up).

    intermediate: DAv8_M36_AZ_paired_cygl1_intermediate (min_sep=2.40deg)

Same extraction logic as the original script (species resolved by name, same
fill/tilenum sanity checks), just widened month glob across 2 calendar years.
"""
import glob
import os
import sys
import hashlib
import subprocess
import gzip

import numpy as np
import netCDF4 as nc
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_cygl1_paired_gain_consistency import read_tilecoord  # noqa: E402

DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
SPECIES_NAME_WANTED = "CYGNSS_L1_DDM3X5_CROP_SCALAR"

ARMS = {
    "intermediate": "DAv8_M36_AZ_paired_cygl1_intermediate",
}
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"

# May 2020 - Dec 2021 (Jan-Apr 2020 spin-up excluded)
YEAR_MONTHS_WANTED = [(2020, m) for m in range(5, 13)] + [(2021, m) for m in range(1, 13)]

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygl1_paired_gain_consistency_22mo_ofa.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_paired_gain_consistency_22mo_summary.txt")

FILL_THRESH = 1e14


def extract_arm(arm, exp_id):
    exp_dir = os.path.join(EXP_ROOT, exp_id)
    tilecoord_file = os.path.join(exp_dir, "output", DOMAIN, "rc_out", f"{exp_id}.ldas_tilecoord.bin")
    ana_root = os.path.join(exp_dir, "output", DOMAIN, "ana", "ens_avg")

    print(f"[{arm}] Reading tilecoord: {tilecoord_file}", flush=True)
    tc = read_tilecoord(tilecoord_file)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]

    files = []
    for yy, mm in YEAR_MONTHS_WANTED:
        files.extend(glob.glob(os.path.join(ana_root, f"Y{yy}", f"M{mm:02d}",
                                             f"{exp_id}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files = sorted(files)
    print(f"[{arm}] Found {len(files)} OFA files (May2020-Dec2021)", flush=True)
    if len(files) == 0:
        print(f"ERROR: no OFA files found for {arm} -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_id_seen = None
    n_bad_tilenum = 0

    for i, fpath in enumerate(files):
        if i % 500 == 0:
            print(f"    [{arm}] ...{i}/{len(files)}", flush=True)
        fname = os.path.basename(fpath)
        stamp = fname.split(".")[-2].rstrip("z")
        yyyymmdd, hhmm = stamp.split("_")
        dt = pd.Timestamp(
            year=int(yyyymmdd[0:4]), month=int(yyyymmdd[4:6]), day=int(yyyymmdd[6:8]),
            hour=int(hhmm[0:2]), minute=int(hhmm[2:4]),
        )

        with nc.Dataset(fpath) as f:
            descr = np.array(f.variables["obsparam_descr"][:])
            spid = np.array(f.variables["obsparam_species_id"][:])
            name_to_id = {str(d): int(s) for d, s in zip(descr, spid)}
            if SPECIES_NAME_WANTED not in name_to_id:
                continue
            wanted_id = name_to_id[SPECIES_NAME_WANTED]
            species_id_seen = wanted_id

            species = np.array(f.variables["species"][:])
            assim_flag = np.array(f.variables["assim_flag"][:])
            mask = (species == wanted_id) & (assim_flag == 1)
            if not np.any(mask):
                continue

            tilenum = np.array(f.variables["tilenum"][:])[mask]
            obs = np.array(f.variables["obs"][:])[mask]
            fcst = np.array(f.variables["fcst"][:])[mask]
            ana = np.array(f.variables["ana"][:])[mask]
            obsvar = np.array(f.variables["obsvar"][:])[mask]
            fcstvar = np.array(f.variables["fcstvar"][:])[mask]
            anavar = np.array(f.variables["anavar"][:])[mask]

            stacked = np.stack([obs, fcst, ana, obsvar, fcstvar, anavar])
            bad_fill = np.any(stacked > FILL_THRESH, axis=0)

            row_idx = tilenum - 1
            bad_tile = (row_idx < 0) | (row_idx >= n_tile)
            n_bad_tilenum += int(bad_tile.sum())
            bad = bad_fill | bad_tile
            row_idx_c = np.clip(row_idx, 0, n_tile - 1)
            tile_id_vals = tile_id_by_row[row_idx_c]
            lon_vals = lon_by_row[row_idx_c]
            lat_vals = lat_by_row[row_idx_c]

            for k in range(len(obs)):
                if bad[k]:
                    continue
                rows.append((
                    arm, dt, int(tilenum[k]), int(tile_id_vals[k]), float(lon_vals[k]), float(lat_vals[k]),
                    float(obs[k]), float(fcst[k]), float(ana[k]),
                    float(obsvar[k]), float(fcstvar[k]), float(anavar[k]),
                ))

    print(f"[{arm}] Assembled {len(rows)} assimilated CYGNSS L1 rows "
          f"(species_id={species_id_seen}), n_bad_tilenum={n_bad_tilenum}", flush=True)
    df = pd.DataFrame(rows, columns=[
        "arm", "datetime", "tilenum", "tile_id", "lon", "lat",
        "obs", "fcst", "ana", "obsvar", "fcstvar", "anavar",
    ])
    return df, species_id_seen, exp_dir, tilecoord_file, ana_root


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_dfs = []
    provenance = []
    species_ids = {}

    for arm, exp_id in ARMS.items():
        df, sp_id, exp_dir, tilecoord_file, ana_root = extract_arm(arm, exp_id)
        all_dfs.append(df)
        species_ids[arm] = sp_id
        provenance.append((arm, exp_id, exp_dir, tilecoord_file, ana_root))

    df = pd.concat(all_dfs, ignore_index=True)
    df["d"] = df["obs"] - df["fcst"]
    df["K"] = df["fcstvar"] / (df["fcstvar"] + df["obsvar"])
    df["pred_incr"] = df["K"] * df["d"]
    df["actual_incr"] = df["ana"] - df["fcst"]
    df = df.sort_values(["arm", "datetime", "tilenum"]).reset_index(drop=True)

    summary_lines = []
    summary_lines.append("CYGNSS L1 paired thinning-density experiment: OFA extraction, "
                          "DA-intermediate only, May2020-Dec2021 (22-month extension)")
    summary_lines.append(f"species_name = {SPECIES_NAME_WANTED}")
    summary_lines.append("")
    summary_lines.append(f"{'arm':14s} {'species_id':>10s} {'N':>8s}")
    for arm in ARMS:
        n = int((df["arm"] == arm).sum())
        summary_lines.append(f"{arm:14s} {species_ids[arm]:>10d} {n:>8d}")

    # through-origin and OLS regression, actual_incr ~ pred_incr
    sub = df.dropna(subset=["pred_incr", "actual_incr"])
    x = sub["pred_incr"].values
    y = sub["actual_incr"].values
    slope_to = float(np.sum(x * y) / np.sum(x * x))
    yhat_to = slope_to * x
    ss_res_to = np.sum((y - yhat_to) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2_to = 1 - ss_res_to / ss_tot

    A = np.vstack([x, np.ones_like(x)]).T
    (slope_ols, intercept_ols), _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    yhat_ols = slope_ols * x + intercept_ols
    r2_ols = 1 - np.sum((y - yhat_ols) ** 2) / ss_tot

    summary_lines.append("")
    summary_lines.append(f"Through-origin regression actual_incr ~ pred_incr: slope={slope_to:.4f}, R2={r2_to:.4f}, N={len(sub)}")
    summary_lines.append(f"OLS regression (with intercept): slope={slope_ols:.4f}, intercept={intercept_ols:.4f}, R2={r2_ols:.4f}")

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    header_lines = ["# CYGNSS L1 paired thinning-density experiment: OFA extraction, "
                    "DA-intermediate, May2020-Dec2021 (22-month extension)"]
    header_lines.append(f"# species_name_wanted: {SPECIES_NAME_WANTED}")
    for arm, exp_id, exp_dir, tilecoord_file, ana_root in provenance:
        header_lines.append(f"# arm: {arm} = {exp_id}")
        header_lines.append(f"#   path: {exp_dir}")
        header_lines.append(f"#   tilecoord_file: {tilecoord_file}")
        header_lines.append(f"#   ofa_root: {ana_root}")
        header_lines.append(f"#   species_id_resolved: {species_ids[arm]}")
    header_lines.append("# months restricted to May2020-Dec2021 (evaluation window; Jan-Apr 2020 spin-up excluded)")
    header_lines.append("# formulas: d=obs-fcst; K=fcstvar/(fcstvar+obsvar); pred_incr=K*d; actual_incr=ana-fcst")
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")
    header_lines.append("# columns: arm,datetime,tilenum,tile_id,lon,lat,obs,fcst,ana,obsvar,fcstvar,anavar,d,K,pred_incr,actual_incr")

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

    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")

    print(f"\nWrote {CSV_PATH}, .sha256, and {SUMMARY_PATH}")
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_paired_gain_consistency_22mo.pkl"))
    return df


if __name__ == "__main__":
    main()
