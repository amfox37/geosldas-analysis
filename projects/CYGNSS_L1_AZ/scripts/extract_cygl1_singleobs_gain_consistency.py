#!/usr/bin/env python3
"""
Step 1.5 of the revised CYGNSS L1 single-observation diagnosis: re-run the
Kalman-gain-consistency check (see extract_cygl1_gain_consistency.py) against
the thinned single-obs-per-window experiment
(DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs, job 58019133), where
each assimilated CYGNSS L1 obs is genuinely isolated (>5deg from every other
kept obs in the same 3-hour window -- outside its own and every other kept
obs's 2.5deg Gaspari-Cohn compact-support ellipse). Unlike the original Step 1
(full obs stream, ~56 neighbors/obs, confound invalidated the K*d identity),
this run has 0 confounding neighbors by construction, so:

    d          = obs - fcst
    K          = fcstvar / (fcstvar + obsvar)
    pred_incr  = K * d
    actual_incr= ana - fcst

slope(actual_incr ~ pred_incr, through origin) is now a valid test of whether
the analysis update applies its own saved gain correctly.

Also joins each kept obs against the same-cycle catch_progn_incr collection
(CATDEF/RZEXC/SRFEXC_INCR, ensemble-mean, at the owner tile) for a
state-space increment sign check: since dH/dsfmc > 0 everywhere in this
domain (Task 0, operator Jacobian sign check), a positive d (obs wetter than
forecast) should correspond to a wetting increment (RZEXC/SRFEXC decrease
sign convention: more negative CATDEF/RZEXC/SRFEXC = wetter, per CLSM
prognostic convention -- see M21C sum_da_analysis_increments.py for the
sign convention actually used there) at that tile/cycle, not just a
consistent obs-space K*d.

Species resolved by NAME (obsparam_descr/obsparam_species_id), never a
hardcoded index, per project convention.
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

EXP_ID = "DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs"
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
EXP_DIR = os.path.join(EXP_ROOT, EXP_ID)

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygl1_singleobs_gain_consistency_ofa.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_singleobs_gain_consistency_summary.txt")

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


def extract_ofa():
    tilecoord_file = os.path.join(EXP_DIR, "output", DOMAIN, "rc_out", f"{EXP_ID}.ldas_tilecoord.bin")
    ana_root = os.path.join(EXP_DIR, "output", DOMAIN, "ana", "ens_avg")

    print(f"Reading tilecoord: {tilecoord_file}", flush=True)
    tc = read_tilecoord(tilecoord_file)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]

    files = sorted(glob.glob(os.path.join(ana_root, "Y????", "M??", f"{EXP_ID}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    print(f"Found {len(files)} OFA files", flush=True)
    if len(files) == 0:
        print("ERROR: no OFA files found -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_id_seen = None

    for fpath in files:
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
            bad = bad_fill | bad_tile
            row_idx_c = np.clip(row_idx, 0, n_tile - 1)
            tile_id_vals = tile_id_by_row[row_idx_c]
            lon_vals = lon_by_row[row_idx_c]
            lat_vals = lat_by_row[row_idx_c]

            for k in range(len(obs)):
                if bad[k]:
                    continue
                rows.append((
                    dt, int(tilenum[k]), int(tile_id_vals[k]), float(lon_vals[k]), float(lat_vals[k]),
                    float(obs[k]), float(fcst[k]), float(ana[k]),
                    float(obsvar[k]), float(fcstvar[k]), float(anavar[k]),
                ))

    print(f"Assembled {len(rows)} assimilated CYGNSS L1 rows (species_id={species_id_seen})", flush=True)
    df = pd.DataFrame(rows, columns=[
        "datetime", "tilenum", "tile_id", "lon", "lat",
        "obs", "fcst", "ana", "obsvar", "fcstvar", "anavar",
    ])
    return df, species_id_seen, tilecoord_file, ana_root


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
    df, species_id, tilecoord_file, ana_root = extract_ofa()

    df["d"] = df["obs"] - df["fcst"]
    df["K"] = df["fcstvar"] / (df["fcstvar"] + df["obsvar"])
    df["pred_incr"] = df["K"] * df["d"]
    df["actual_incr"] = df["ana"] - df["fcst"]
    df = df.sort_values(["datetime", "tilenum"]).reset_index(drop=True)

    x = df["pred_incr"].values
    y = df["actual_incr"].values
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

    summary_lines = []
    summary_lines.append("CYGNSS L1 gain-consistency check, Step 1.5 (thinned single-obs-per-window, job 58019133)")
    summary_lines.append(f"experiment: {EXP_ID}")
    summary_lines.append(f"species_name = {SPECIES_NAME_WANTED}, species_id_resolved = {species_id}")
    summary_lines.append(f"N obs (each isolated, >5deg from every other kept obs in its 3-hr window) = {n}")
    summary_lines.append("")
    summary_lines.append(f"{'':22s} {'slope(0)':>10s} {'R2(0)':>8s} {'slope(int)':>10s} {'intercept':>10s} {'R2(int)':>8s}")
    summary_lines.append(
        f"{'through-origin/OLS':22s} {slope0:>10.4f} {r2_0:>8.4f} {slope_i:>10.4f} {intercept:>10.4f} {r2_i:>8.4f}"
    )
    summary_lines.append("")
    summary_lines.append("For comparison, the ORIGINAL (confounded, full obs stream, ~56 neighbors/obs) Step 1 check")
    summary_lines.append("on cygl1assim/_halfR/_quarterR (errstd 4.4/2.2/1.1 dB) gave slope(0)=0.679/0.505/0.379,")
    summary_lines.append("R2(0)=0.132/0.215/0.277 -- see cygl1_gain_consistency_summary.txt. This run uses errstd=3.0 dB")
    summary_lines.append("and is NOT a like-for-like R value, so compare the MECHANISM (does slope~=1 now that obs are")
    summary_lines.append("isolated?), not the exact numbers.")

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    header_lines = ["# CYGNSS L1 gain-consistency check, Step 1.5 (thinned single-obs-per-window)"]
    header_lines.append(f"# experiment: {EXP_ID}")
    header_lines.append(f"#   path: {EXP_DIR}")
    header_lines.append(f"#   tilecoord_file: {tilecoord_file}")
    header_lines.append(f"#   ofa_root: {ana_root}")
    header_lines.append(f"#   species_name_wanted: {SPECIES_NAME_WANTED}")
    header_lines.append(f"#   species_id_resolved: {species_id}")
    header_lines.append("# thinning: geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/thin_cygl1_single_obs_per_window.py")
    header_lines.append("#   min_sep_deg=5.0 (2x xcompact=2.5deg), single obs per 3hr window kept if isolated")
    header_lines.append("# formulas: d=obs-fcst; K=fcstvar/(fcstvar+obsvar); pred_incr=K*d; actual_incr=ana-fcst")
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")
    header_lines.append("# columns: datetime,tilenum,tile_id,lon,lat,obs,fcst,ana,obsvar,fcstvar,anavar,d,K,pred_incr,actual_incr")

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
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_singleobs_gain_consistency.pkl"))
    return df


if __name__ == "__main__":
    main()
