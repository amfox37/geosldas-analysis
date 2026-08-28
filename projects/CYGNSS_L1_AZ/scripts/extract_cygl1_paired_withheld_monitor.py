#!/usr/bin/env python3
"""
CYGNSS L1 paired thinning-density experiment: extract the 12 monitor-only
species (everything except CYGNSS_L1_DDM3X5_CROP_SCALAR) from all 4 arms
(OL + 3 DA), May-Oct 2020 only, for the withheld-monitor matched comparison
(cygl1_paired_withheld_monitor_diagnostics.py).

Only obs/fcst are extracted (ana is meaningless for monitor-only obs, since
assim_flag==0 in all four arms for these 12 species -- verified). Species are
resolved by NAME (obsparam_descr/obsparam_species_id) per file, never a
hardcoded index, per project convention (combines the multi-species pattern
of extract_cygl1_l3_pairs.py with the multi-experiment loop pattern of
extract_cygl1_gain_consistency.py).

Experiments (all under /discover/nobackup/projects/land_da/cygl1_operator_test/):
  OL:            OLv8_M36_AZ_paired_monitor             (L1 assim=.false., full L1 stream monitored)
  sparse:        DAv8_M36_AZ_paired_cygl1_sparse
  intermediate:  DAv8_M36_AZ_paired_cygl1_intermediate
  dense:         DAv8_M36_AZ_paired_cygl1_dense

Join key downstream is (species_name, tile_id, datetime) -- tile_id (not
tilenum) is used because it's the stable tilecoord identity shared across
experiments (tilenum is a row index into each experiment's own OFA species
subset ordering and is NOT guaranteed to line up across arms for the same
physical tile/obs). tilenum is included in the output for traceability only.
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

SPECIES_NAMES_MONITOR = [
    "SMOS_fit_Tbh_A", "SMOS_fit_Tbh_D", "SMOS_fit_Tbv_A", "SMOS_fit_Tbv_D",
    "SMAP_L1C_Tbh_A", "SMAP_L1C_Tbh_D", "SMAP_L1C_Tbv_A", "SMAP_L1C_Tbv_D",
    "ASCAT_HSAF_META_SM", "ASCAT_HSAF_METB_SM", "ASCAT_HSAF_METC_SM", "CYGNSS_SM_6hr",
]

ARMS = {
    "OL": "OLv8_M36_AZ_paired_monitor",
    "sparse": "DAv8_M36_AZ_paired_cygl1_sparse",
    "intermediate": "DAv8_M36_AZ_paired_cygl1_intermediate",
    "dense": "DAv8_M36_AZ_paired_cygl1_dense",
}
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
MONTHS_WANTED = ["M05", "M06", "M07", "M08", "M09", "M10"]

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
CSV_PATH = os.path.join(OUT_DIR, "cygl1_paired_withheld_monitor_ofa.csv.gz")

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

    files = []
    for mm in MONTHS_WANTED:
        files.extend(glob.glob(os.path.join(ana_root, "Y2020", mm, f"{exp_id}.ens_avg.ldas_ObsFcstAna.*.nc4")))
    files = sorted(files)
    print(f"[{arm}] Found {len(files)} OFA files (Y2020 M05-M10)", flush=True)
    if len(files) == 0:
        print(f"ERROR: no OFA files found for {arm} -- aborting.", file=sys.stderr)
        sys.exit(1)

    rows = []
    species_name_to_id_seen = {}
    n_bad_tilenum = 0
    n_assim_flag_nonzero = 0

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
            wanted_ids = {}
            for nm in SPECIES_NAMES_MONITOR:
                if nm in name_to_id:
                    wanted_ids[name_to_id[nm]] = nm
                    species_name_to_id_seen[nm] = name_to_id[nm]
            if not wanted_ids:
                continue

            species = np.array(f.variables["species"][:])
            assim_flag = np.array(f.variables["assim_flag"][:])
            mask = np.isin(species, list(wanted_ids.keys()))
            if not np.any(mask):
                continue

            n_assim_flag_nonzero += int((assim_flag[mask] != 0).sum())

            tilenum = np.array(f.variables["tilenum"][:])[mask]
            obs = np.array(f.variables["obs"][:])[mask]
            fcst = np.array(f.variables["fcst"][:])[mask]
            spec_sel = species[mask]

            stacked = np.stack([obs, fcst])
            bad_fill = np.any(stacked > FILL_THRESH, axis=0)

            row_idx = tilenum - 1
            bad_tile = (row_idx < 0) | (row_idx >= n_tile)
            n_bad_tilenum += int(bad_tile.sum())
            bad = bad_fill | bad_tile
            row_idx_c = np.clip(row_idx, 0, n_tile - 1)
            tile_id_vals = tile_id_by_row[row_idx_c]

            for k in range(len(obs)):
                if bad[k]:
                    continue
                sp = int(spec_sel[k])
                rows.append((
                    arm, dt, int(tilenum[k]), int(tile_id_vals[k]),
                    wanted_ids[sp], float(obs[k]), float(fcst[k]),
                ))

    print(f"[{arm}] Assembled {len(rows)} monitor-species rows. n_bad_tilenum={n_bad_tilenum}, "
          f"n_assim_flag_nonzero={n_assim_flag_nonzero} (expect 0 -- all 12 species monitor-only)",
          flush=True)
    df = pd.DataFrame(rows, columns=["arm", "datetime", "tilenum", "tile_id", "species_name", "obs", "fcst"])
    return df, species_name_to_id_seen, n_assim_flag_nonzero


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_dfs = []
    provenance = []
    total_nonzero_assim = 0

    for arm, exp_id in ARMS.items():
        df, sp_map, n_nonzero = extract_arm(arm, exp_id)
        all_dfs.append(df)
        total_nonzero_assim += n_nonzero
        provenance.append((arm, exp_id, sp_map))

    df = pd.concat(all_dfs, ignore_index=True)
    df = df.sort_values(["arm", "species_name", "datetime", "tile_id"]).reset_index(drop=True)

    print("\nRow counts by arm x species:")
    print(df.groupby(["arm", "species_name"]).size().unstack(fill_value=0))

    if total_nonzero_assim > 0:
        print(f"\nWARNING: {total_nonzero_assim} monitor-species rows had assim_flag != 0 "
              f"across all arms -- these species were expected to be monitor-only everywhere.",
              file=sys.stderr)

    header_lines = ["# CYGNSS L1 paired thinning-density experiment: withheld-monitor OFA extraction"]
    header_lines.append("# 12 monitor-only species (all except CYGNSS_L1_DDM3X5_CROP_SCALAR), 4 arms, Y2020 M05-M10")
    for arm, exp_id, sp_map in provenance:
        header_lines.append(f"# arm: {arm} = {exp_id}")
        for nm, sid in sp_map.items():
            header_lines.append(f"#   {nm} -> species_id {sid}")
    header_lines.append(f"# total_assim_flag_nonzero_rows: {total_nonzero_assim} (expect 0)")
    header_lines.append("# join key downstream: (species_name, tile_id, datetime) -- tile_id is the stable")
    header_lines.append("#   tilecoord identity shared across arms; tilenum included for traceability only")
    try:
        git_tag = subprocess.check_output(
            ["git", "-C", "/gpfsm/dnb06/projects/p284/geosldas-analysis", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_tag = "unknown"
    header_lines.append(f"# geosldas-analysis_repo_commit: {git_tag}")
    header_lines.append("# columns: arm,datetime,tilenum,tile_id,species_name,obs,fcst")

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

    print(f"\nWrote {CSV_PATH} and .sha256")
    df.to_pickle(os.path.join(OUT_DIR, "_cygl1_paired_withheld_monitor.pkl"))
    return df


if __name__ == "__main__":
    main()
