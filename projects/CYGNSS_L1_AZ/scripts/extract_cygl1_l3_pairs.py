#!/usr/bin/env python3
"""
Extract paired CYGNSS L1 (CYGNSS_L1_DDM3X5_CROP_SCALAR) and CYGNSS L3
(CYGNSS_SM_6hr) observation/forecast/analysis series from the OLv8_M36_all_sensors_AZ_scaled
open-loop / monitoring run's ObsFcstAna files for calendar year 2020.

Experiment path : /discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_all_sensors_AZ_scaled
OFA pattern      : output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y<YYYY>/M<MM>/OLv8_M36_all_sensors_AZ_scaled.ens_avg.ldas_ObsFcstAna.<YYYYMMDD_HHMM>z.nc4
Tilecoord file   : output/SMAP_EASEv2_M36_GLOBAL/rc_out/OLv8_M36_all_sensors_AZ_scaled.ldas_tilecoord.bin

Species are resolved by NAME using each file's own obsparam_descr / obsparam_species_id
variables (never by a hardcoded index), per project convention (species indices are NOT
stable across experiments in this project).

Output: CSV (gzip) with one row per observation of either CYGNSS species, plus a
sanity-check summary printed to stdout (Section 8 of the 2026-08-20 brief).
"""
import glob
import os
import struct
import sys
import hashlib
import subprocess

import numpy as np
import netCDF4 as nc
import pandas as pd

EXP_DIR = "/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_all_sensors_AZ_scaled"
EXP_ID = "OLv8_M36_all_sensors_AZ_scaled"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
TILECOORD_FILE = os.path.join(EXP_DIR, "output", DOMAIN, "rc_out", f"{EXP_ID}.ldas_tilecoord.bin")
ANA_ROOT = os.path.join(EXP_DIR, "output", DOMAIN, "ana", "ens_avg")

SPECIES_NAMES_WANTED = ["CYGNSS_L1_DDM3X5_CROP_SCALAR", "CYGNSS_SM_6hr"]

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygnss_l1_l3_paired_ofa_2020.csv.gz")


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


def main():
    print(f"[1] Reading tilecoord: {TILECOORD_FILE}", flush=True)
    tc = read_tilecoord(TILECOORD_FILE)
    n_tile = tc["N_tile"]
    tile_id_by_row = tc["tile_id"]         # row index (0-based) -> tile_id
    lon_by_row = tc["com_lon"]
    lat_by_row = tc["com_lat"]
    print(f"    N_tile = {n_tile}")

    files = sorted(glob.glob(os.path.join(ANA_ROOT, "Y2020", "M??", f"{EXP_ID}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    print(f"[2] Found {len(files)} OFA files for 2020", flush=True)
    if len(files) == 0:
        print("ERROR: no OFA files found for 2020 -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_name_to_id_seen = {}
    n_bad_tilenum = 0

    for i, fpath in enumerate(files):
        if i % 200 == 0:
            print(f"    ...{i}/{len(files)}", flush=True)
        fname = os.path.basename(fpath)
        # OLv8_M36_all_sensors_AZ_scaled.ens_avg.ldas_ObsFcstAna.20200101_0300z.nc4
        stamp = fname.split(".")[-2]  # 20200101_0300z
        stamp = stamp.rstrip("z")
        yyyymmdd, hhmm = stamp.split("_")
        dt = pd.Timestamp(
            year=int(yyyymmdd[0:4]), month=int(yyyymmdd[4:6]), day=int(yyyymmdd[6:8]),
            hour=int(hhmm[0:2]), minute=int(hhmm[2:4]),
        )

        with nc.Dataset(fpath) as f:
            descr = np.array(f.variables["obsparam_descr"][:])
            spid = np.array(f.variables["obsparam_species_id"][:])
            name_to_id = {str(d): int(s) for d, s in zip(descr, spid)}
            wanted_ids = {}
            for nm in SPECIES_NAMES_WANTED:
                if nm in name_to_id:
                    wanted_ids[name_to_id[nm]] = nm
                    species_name_to_id_seen[nm] = name_to_id[nm]

            species = np.array(f.variables["species"][:])
            mask = np.isin(species, list(wanted_ids.keys()))
            if not np.any(mask):
                continue

            tilenum = np.array(f.variables["tilenum"][:])[mask]
            obs = np.array(f.variables["obs"][:])[mask]
            fcst = np.array(f.variables["fcst"][:])[mask]
            ana = np.array(f.variables["ana"][:])[mask]
            assim_flag = np.array(f.variables["assim_flag"][:])[mask]
            spec_sel = species[mask]

            # fill value handling (1e15)
            for arr in (obs, fcst, ana):
                arr[arr > 1e14] = np.nan

            row_idx = tilenum - 1  # 0-based row index into tilecoord arrays
            bad = (row_idx < 0) | (row_idx >= n_tile)
            n_bad_tilenum += int(bad.sum())
            row_idx_c = np.clip(row_idx, 0, n_tile - 1)

            tile_id_vals = tile_id_by_row[row_idx_c]
            lon_vals = lon_by_row[row_idx_c]
            lat_vals = lat_by_row[row_idx_c]

            for k in range(len(spec_sel)):
                if bad[k]:
                    continue
                sp = int(spec_sel[k])
                rows.append((
                    dt, int(tile_id_vals[k]), float(lon_vals[k]), float(lat_vals[k]),
                    sp, wanted_ids[sp], float(obs[k]), float(fcst[k]), float(ana[k]),
                    int(assim_flag[k]),
                ))

    print(f"[3] Assembled {len(rows)} rows. n_bad_tilenum={n_bad_tilenum}", flush=True)
    if n_bad_tilenum > 0:
        print(f"WARNING: {n_bad_tilenum} observations had out-of-range tilenum and were dropped", file=sys.stderr)

    df = pd.DataFrame(rows, columns=[
        "datetime", "tile_id", "lon", "lat", "species", "species_name",
        "obs", "fcst", "ana", "assim_flag",
    ])
    df = df.sort_values(["species", "datetime", "tile_id"]).reset_index(drop=True)

    print("[4] species_name -> species_id resolved from obsparam (as found in files):")
    for nm, sid in species_name_to_id_seen.items():
        print(f"      {nm} -> {sid}")

    print("[5] Row counts by species:")
    print(df.groupby("species_name").size())

    print(f"[6] Writing CSV: {CSV_PATH}")
    os.makedirs(OUT_DIR, exist_ok=True)

    header_lines = [
        f"# experiment_path: {EXP_DIR}",
        f"# experiment_id: {EXP_ID}",
        f"# ofa_filename_pattern: output/{DOMAIN}/ana/ens_avg/Y<YYYY>/M<MM>/{EXP_ID}.ens_avg.ldas_ObsFcstAna.<YYYYMMDD_HHMM>z.nc4",
        f"# tilecoord_file: {TILECOORD_FILE}",
        "# species_name_to_id_mapping (resolved per-file from obsparam_descr/obsparam_species_id, NOT assumed):",
    ]
    for nm, sid in species_name_to_id_seen.items():
        header_lines.append(f"#   {nm} = {sid}")
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")
    header_lines.append("# columns: datetime,tile_id,lon,lat,species,species_name,obs,fcst,ana,assim_flag")
    header_lines.append("# lon/lat are tilecoord (com_lon/com_lat) tile-center values, NOT the OFA superobs centroid.")

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

    # Save a pickle too for downstream analysis scripts (not part of deliverable, scratch convenience)
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_l3_pairs_2020.pkl"))

    return df


if __name__ == "__main__":
    main()
