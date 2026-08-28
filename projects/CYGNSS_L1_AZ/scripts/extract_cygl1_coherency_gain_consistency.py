#!/usr/bin/env python3
"""
CYGNSS L1 coherency-screening experiment: extract assimilated (assim_flag==1)
CYGNSS_L1_DDM3X5_CROP_SCALAR ObsFcstAna (OFA) rows from the two DA arms:

    coherency_screened: DAv8_M36_AZ_paired_cygl1_coherency_screened  (coherency_ratio>=0.523, one-sided)
    coherency_randmatch: DAv8_M36_AZ_paired_cygl1_coherency_randmatch (random-matched control, same per-window count)

Full calendar year 2020 (Jan-Dec, no spin-up excluded, per the experiment's
design spec). Same extraction pattern as extract_cygl1_paired_gain_consistency.py,
generalized to these two arms and the full-year period; local-interaction-count
join is not part of this experiment (that covariate was specific to the
thinning-density family) and is omitted here.

Species resolved by NAME (obsparam_descr/obsparam_species_id) per file, never
a hardcoded index, per project convention.
"""
import glob
import os
import struct
import sys
import hashlib
import subprocess
import gzip

import numpy as np
import netCDF4 as nc
import pandas as pd

DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
SPECIES_NAME_WANTED = "CYGNSS_L1_DDM3X5_CROP_SCALAR"

ARMS = {
    "coherency_screened": "DAv8_M36_AZ_paired_cygl1_coherency_screened",
    "coherency_randmatch": "DAv8_M36_AZ_paired_cygl1_coherency_randmatch",
}
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"

MONTHS_WANTED = [f"M{m:02d}" for m in range(1, 13)]

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygl1_coherency_gain_consistency_ofa.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_coherency_gain_consistency_summary.txt")

FILL_THRESH = 1e14


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
    for mm in MONTHS_WANTED:
        files.extend(glob.glob(os.path.join(ana_root, "Y2020", mm, f"{exp_id}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files = sorted(files)
    print(f"[{arm}] Found {len(files)} OFA files (Y2020 M01-M12)", flush=True)
    if len(files) == 0:
        print(f"ERROR: no OFA files found for {arm} -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_id_seen = None
    n_bad_tilenum = 0

    for i, fpath in enumerate(files):
        if i % 200 == 0:
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
    summary_lines.append("CYGNSS L1 coherency-screening experiment: OFA extraction, 2 DA arms, Jan-Dec 2020")
    summary_lines.append(f"species_name = {SPECIES_NAME_WANTED}")
    summary_lines.append("")
    summary_lines.append(f"{'arm':22s} {'species_id':>10s} {'N':>8s}")
    for arm in ARMS:
        n = int((df["arm"] == arm).sum())
        summary_lines.append(f"{arm:22s} {species_ids[arm]:>10d} {n:>8d}")
    summary_text = "\n".join(summary_lines)
    print(summary_text)

    header_lines = ["# CYGNSS L1 coherency-screening experiment: OFA extraction, 2 DA arms, Jan-Dec 2020"]
    header_lines.append(f"# species_name_wanted: {SPECIES_NAME_WANTED}")
    for arm, exp_id, exp_dir, tilecoord_file, ana_root in provenance:
        header_lines.append(f"# arm: {arm} = {exp_id}")
        header_lines.append(f"#   path: {exp_dir}")
        header_lines.append(f"#   tilecoord_file: {tilecoord_file}")
        header_lines.append(f"#   ofa_root: {ana_root}")
        header_lines.append(f"#   species_id_resolved: {species_ids[arm]}")
    header_lines.append("# months: Y2020 M01-M12 (full calendar year, no spin-up excluded per design spec)")
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
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_coherency_gain_consistency.pkl"))
    return df


if __name__ == "__main__":
    main()
