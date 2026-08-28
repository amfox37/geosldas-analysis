#!/usr/bin/env python3
"""
Step 1 of the revised (2026-08-21) CYGNSS L1 single-observation diagnosis: a
Kalman-gain-consistency check built entirely from existing ObsFcstAna (OFA)
output -- no rerun, no Fortran change.

For every assimilated (assim_flag==1) CYGNSS_L1_DDM3X5_CROP_SCALAR observation
in the three DA experiments below, OFA already saves obsvar (= R, post-scaling)
and fcstvar (= HPH^T through the real nonlinear operator, i.e. ensemble
variance in observation space) alongside obs/fcst/ana. That is enough to
compute, per observation:

    d          = obs - fcst                       (innovation)
    K          = fcstvar / (fcstvar + obsvar)      (scalar gain)
    pred_incr  = K * d                             (gain-predicted obs-space increment)
    actual_incr= ana - fcst                        (actual saved obs-space increment)

Regressing actual_incr against pred_incr tells us whether the analysis update
is applying its own saved gain correctly:
    slope ~= 1, no offset   -> update matches its own K,d; any problem is upstream
                                (R, ensemble spread, or the operator itself)
    slope ~= 1, with offset -> some additive term is present beyond K*d
    slope ~= 0               -> update ignores its own gain (a real defect)
    slope <  0               -> sign error

Caveat: this identity is exact only for a single isolated observation updating
a single tile. With multi-tile-footprint observations and multiple
simultaneous observations per cycle, the true update is a joint multi-obs
solve, so this is an approximation -- but the sign and slope remain
diagnostic. See memory `cygl1-single-obs-diagnosis-roadmap` for the full
staged plan this is step 1 of.

Experiments (all under /discover/nobackup/projects/land_da/cygl1_operator_test/):
  DAv8_M36_all_sensors_AZ_scaled_cygl1assim            (CYGNSS L1 errstd 4.4 dB, full R)
  DAv8_M36_all_sensors_AZ_scaled_cygl1assim_halfR      (errstd 2.2 dB)
  DAv8_M36_all_sensors_AZ_scaled_cygl1assim_quarterR   (errstd 1.1 dB)

Species resolved by NAME per file (obsparam_descr/obsparam_species_id), never
by a hardcoded index, per project convention.
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

EXPERIMENTS = {
    "cygl1assim": "DAv8_M36_all_sensors_AZ_scaled_cygl1assim",
    "cygl1assim_halfR": "DAv8_M36_all_sensors_AZ_scaled_cygl1assim_halfR",
    "cygl1assim_quarterR": "DAv8_M36_all_sensors_AZ_scaled_cygl1assim_quarterR",
}
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygl1_gain_consistency_ofa.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_gain_consistency_summary.txt")

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


def extract_experiment(exp_key, exp_id):
    exp_dir = os.path.join(EXP_ROOT, exp_id)
    tilecoord_file = os.path.join(exp_dir, "output", DOMAIN, "rc_out", f"{exp_id}.ldas_tilecoord.bin")
    ana_root = os.path.join(exp_dir, "output", DOMAIN, "ana", "ens_avg")

    print(f"[{exp_key}] Reading tilecoord: {tilecoord_file}", flush=True)
    tc = read_tilecoord(tilecoord_file)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]

    files = sorted(glob.glob(os.path.join(ana_root, "Y????", "M??", f"{exp_id}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    print(f"[{exp_key}] Found {len(files)} OFA files", flush=True)
    if len(files) == 0:
        print(f"ERROR: no OFA files found for {exp_key} -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_id_seen = None
    n_bad_tilenum = 0

    for i, fpath in enumerate(files):
        if i % 500 == 0:
            print(f"    ...{i}/{len(files)}", flush=True)
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
                    exp_key, dt, int(tile_id_vals[k]), float(lon_vals[k]), float(lat_vals[k]),
                    float(obs[k]), float(fcst[k]), float(ana[k]),
                    float(obsvar[k]), float(fcstvar[k]), float(anavar[k]),
                ))

    print(f"[{exp_key}] Assembled {len(rows)} assimilated rows. n_bad_tilenum={n_bad_tilenum}", flush=True)
    df = pd.DataFrame(rows, columns=[
        "experiment", "datetime", "tile_id", "lon", "lat",
        "obs", "fcst", "ana", "obsvar", "fcstvar", "anavar",
    ])
    return df, species_id_seen, exp_dir, tilecoord_file, ana_root


def through_origin_slope(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    denom = np.sum(x * x)
    slope = np.sum(x * y) / denom if denom > 0 else np.nan
    resid = y - slope * x
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return slope, r2


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_dfs = []
    provenance = []
    species_ids = {}

    for exp_key, exp_id in EXPERIMENTS.items():
        df, sp_id, exp_dir, tilecoord_file, ana_root = extract_experiment(exp_key, exp_id)
        all_dfs.append(df)
        species_ids[exp_key] = sp_id
        provenance.append((exp_key, exp_id, exp_dir, tilecoord_file, ana_root))

    df = pd.concat(all_dfs, ignore_index=True)
    df["d"] = df["obs"] - df["fcst"]
    df["K"] = df["fcstvar"] / (df["fcstvar"] + df["obsvar"])
    df["pred_incr"] = df["K"] * df["d"]
    df["actual_incr"] = df["ana"] - df["fcst"]

    df = df.sort_values(["experiment", "datetime", "tile_id"]).reset_index(drop=True)

    summary_lines = []
    summary_lines.append("CYGNSS L1 gain-consistency check (Step 1, OFA-only, no rerun)")
    summary_lines.append(f"species_name = {SPECIES_NAME_WANTED}")
    summary_lines.append("")
    summary_lines.append(f"{'experiment':22s} {'species_id':>10s} {'N':>8s} {'slope(0)':>10s} {'R2(0)':>8s} {'slope(int)':>10s} {'intercept':>10s} {'R2(int)':>8s}")

    results = {}
    for exp_key in EXPERIMENTS:
        sub = df[df["experiment"] == exp_key]
        x = sub["pred_incr"].values
        y = sub["actual_incr"].values
        n = len(x)
        slope0, r2_0 = through_origin_slope(x, y)
        if n > 1 and np.std(x) > 0:
            slope_i, intercept = np.polyfit(x, y, 1)
            yhat = slope_i * x + intercept
            ss_res = np.sum((y - yhat) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2_i = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        else:
            slope_i, intercept, r2_i = np.nan, np.nan, np.nan
        results[exp_key] = dict(n=n, slope0=slope0, r2_0=r2_0, slope_i=slope_i, intercept=intercept, r2_i=r2_i)
        summary_lines.append(
            f"{exp_key:22s} {species_ids[exp_key]:>10d} {n:>8d} {slope0:>10.4f} {r2_0:>8.4f} "
            f"{slope_i:>10.4f} {intercept:>10.4f} {r2_i:>8.4f}"
        )

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    header_lines = ["# CYGNSS L1 gain-consistency check (Step 1 of revised single-obs diagnosis)"]
    header_lines.append(f"# species_name_wanted: {SPECIES_NAME_WANTED}")
    for exp_key, exp_id, exp_dir, tilecoord_file, ana_root in provenance:
        header_lines.append(f"# experiment: {exp_key} = {exp_id}")
        header_lines.append(f"#   path: {exp_dir}")
        header_lines.append(f"#   tilecoord_file: {tilecoord_file}")
        header_lines.append(f"#   ofa_root: {ana_root}")
        header_lines.append(f"#   species_id_resolved: {species_ids[exp_key]}")
    header_lines.append("# formulas: d=obs-fcst; K=fcstvar/(fcstvar+obsvar); pred_incr=K*d; actual_incr=ana-fcst")
    header_lines.append("# obsvar=R (post-scaling); fcstvar=HPH^T through the real nonlinear operator (both saved natively in OFA, no rerun needed)")
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")
    header_lines.append("# columns: experiment,datetime,tile_id,lon,lat,obs,fcst,ana,obsvar,fcstvar,anavar,d,K,pred_incr,actual_incr")

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
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_gain_consistency.pkl"))

    return df, results


if __name__ == "__main__":
    main()
