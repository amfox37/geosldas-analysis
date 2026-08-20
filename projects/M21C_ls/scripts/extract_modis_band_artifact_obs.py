#!/usr/bin/env python3
"""Extract per-window MODIS SCF observations near 90E from existing ObsFcstAna files.

One-off diagnostic for the longitude-band artifact brief (see
docs/notes/discover_modis_band_artifact.md). Reads already-archived OFA
.bin files (no model rerun), joins tile_id -> lon/lat/i_indg via the
tilecoord file, and writes one row per MODIS observation within a lon/lat
box around the 90E EASEv2 M36 column boundary.

Species IDs (confirmed via obsparam for this experiment, NOT the values
originally assumed in the brief):
    12 = MYD10C1 (Aqua)
    13 = MOD10C1 (Terra)
"""

import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "python" / "io"))
from read_GEOSldas import read_ObsFcstAna, read_tilecoord  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRATCH_DIR = REPO_ROOT / "projects" / "M21C_ls" / "output" / "modis_band_artifact_scratch"
OUTPUT_DIR = REPO_ROOT / "projects" / "M21C_ls" / "output"
OUTPUT_CSV_GZ = OUTPUT_DIR / "modis_obs_by_window_85E_95E.csv.gz"

TILECOORD_PATH = (
    "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/"
    "LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_DAv8_M36.ldas_tilecoord.bin"
)

EXPERIMENT_PATH = (
    "/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/"
    "LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/ (archived: "
    "SMAP_EASEv2_M36_GLOBAL_ana_tars/Y<year>.tar)"
)
OFA_PATTERN = "Y<YYYY>/M<MM>/LS_DAv8_M36.ens_avg.ldas_ObsFcstAna.<YYYYMMDD>_<HHMM>z.bin"
CODE_VERSION_NOTE = (
    "GEOSldas_GridComp v3.2.0 (GEOSldas_release tag v20.2.0, built 2025-12-23) "
    "used for Task 3 code inspection as best-available proxy -- LS_DAv8_M36_v2's "
    "ORIGINAL build (install-Aggressive under the now-deleted "
    "GEOSldas_Landsweeper_v2 tree) is not recoverable; the .bin files extracted "
    "here are LS_DAv8_M36_v2's real output, only the source-code cross-reference "
    "uses the proxy tag."
)

MODIS_SPECIES = {12: "MYD10C1_Aqua", 13: "MOD10C1_Terra"}

LON_MIN, LON_MAX = 85.0, 95.0
LAT_MIN, LAT_MAX = 25.0, 45.0

TARGET_DATES = [
    "20040115", "20040415", "20040715", "20041015",
    "20010715", "20030715", "20050715",
    "20030115", "20050115", "20060115",
]
WINDOWS = ["0000", "0300", "0600", "0900", "1200", "1500", "1800", "2100"]


def find_bin_files():
    files = []
    for date in TARGET_DATES:
        yyyy, mm = date[:4], date[4:6]
        for hhmm in WINDOWS:
            p = list((SCRATCH_DIR / f"Y{yyyy}" / f"M{mm}").glob(
                f"LS_DAv8_M36.ens_avg.ldas_ObsFcstAna.{date}_{hhmm}z.bin"
            ))
            if len(p) != 1:
                raise FileNotFoundError(f"Expected 1 file for {date} {hhmm}z, found {len(p)}")
            files.append(p[0])
    return files


def main():
    print(f"Reading tilecoord: {TILECOORD_PATH}")
    tc = read_tilecoord(TILECOORD_PATH)
    n_tile = tc["N_tile"]
    tile_id = tc["tile_id"]
    com_lon = tc["com_lon"]
    com_lat = tc["com_lat"]
    i_indg = tc["i_indg"]

    assert np.array_equal(tile_id, np.arange(1, n_tile + 1)), (
        "tile_id is not a simple 1..N_tile index; obs_tilenum join assumption is invalid"
    )

    box_mask = (
        (com_lon >= LON_MIN) & (com_lon <= LON_MAX)
        & (com_lat >= LAT_MIN) & (com_lat <= LAT_MAX)
    )
    n_box_tiles = int(box_mask.sum())
    print(f"Tiles in [{LON_MIN},{LON_MAX}]E x [{LAT_MIN},{LAT_MAX}]N box: {n_box_tiles}")
    if n_box_tiles < 500:
        raise RuntimeError(
            f"Only {n_box_tiles} tiles in box; brief says widen if <500 -- stopping, not auto-widening."
        )

    box_tile_ids = set(tile_id[box_mask].tolist())

    bin_files = find_bin_files()
    print(f"Found {len(bin_files)} OFA .bin files to process")

    rows = []
    for fpath in bin_files:
        # filename: ...ldas_ObsFcstAna.<YYYYMMDD>_<HHMM>z.bin
        stem = fpath.name
        date_token = stem.split(".")[-2]  # "<YYYYMMDD>_<HHMM>z"
        date_str, hhmm_z = date_token.split("_")
        window_str = hhmm_z.rstrip("z")

        d = read_ObsFcstAna(str(fpath))
        obs_species = np.asarray(d["obs_species"])
        obs_tilenum = np.asarray(d["obs_tilenum"])
        obs_assim = np.asarray(d["obs_assim"])
        obs_obs = np.asarray(d["obs_obs"])
        obs_fcst = np.asarray(d["obs_fcst"])

        species_mask = np.isin(obs_species, list(MODIS_SPECIES.keys()))
        if not species_mask.any():
            continue

        tilenum_sel = obs_tilenum[species_mask]
        in_box = np.array([t in box_tile_ids for t in tilenum_sel])
        if not in_box.any():
            continue

        idx = np.nonzero(species_mask)[0][in_box]

        for i in idx:
            tnum = int(obs_tilenum[i])
            tc_idx = tnum - 1  # tile_id[tc_idx] == tnum
            rows.append({
                "date": date_str,
                "window": window_str,
                "tile_id": tnum,
                "lon": float(com_lon[tc_idx]),
                "lat": float(com_lat[tc_idx]),
                "i_indg": int(i_indg[tc_idx]),
                "species": int(obs_species[i]),
                "species_name": MODIS_SPECIES[int(obs_species[i])],
                "obs": float(obs_obs[i]),
                "fcst": float(obs_fcst[i]),
                "assim_flag": int(obs_assim[i]),
            })

    print(f"Total MODIS obs rows in box: {len(rows)}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    header_lines = [
        f"# experiment_path: {EXPERIMENT_PATH}",
        f"# ofa_filename_pattern: {OFA_PATTERN}",
        f"# tilecoord_file: {TILECOORD_PATH}",
        f"# code_version_note: {CODE_VERSION_NOTE}",
        f"# species_mapping: 12=MYD10C1(Aqua), 13=MOD10C1(Terra) "
        "(confirmed via obsparam + Section-8 sanity-check reproduction; "
        "differs from the brief's original assumption of 11/12)",
        f"# lon_box: [{LON_MIN},{LON_MAX}]E  lat_box: [{LAT_MIN},{LAT_MAX}]N  n_box_tiles: {n_box_tiles}",
        f"# target_dates: {','.join(TARGET_DATES)}",
        "# lon/lat are tile com_lon/com_lat from the tilecoord file (NOT the OFA file's "
        "own obs_lon/obs_lat, which is a jittering superob center)",
    ]

    cols = ["date", "window", "tile_id", "lon", "lat", "i_indg",
            "species", "species_name", "obs", "fcst", "assim_flag"]

    with gzip.open(OUTPUT_CSV_GZ, "wt") as f:
        for line in header_lines:
            f.write(line + "\n")
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")

    print(f"Wrote {OUTPUT_CSV_GZ}")

    sha256 = hashlib.sha256()
    with open(OUTPUT_CSV_GZ, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256.update(chunk)
    checksum_path = OUTPUT_CSV_GZ.with_suffix(OUTPUT_CSV_GZ.suffix + ".sha256")
    with open(checksum_path, "w") as f:
        f.write(f"{sha256.hexdigest()}  {OUTPUT_CSV_GZ.name}\n")
    print(f"sha256: {sha256.hexdigest()}")
    print(f"Wrote {checksum_path}")


if __name__ == "__main__":
    main()
