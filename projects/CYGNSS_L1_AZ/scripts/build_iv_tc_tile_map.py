#!/usr/bin/env python3
"""Map the AZ domain tile numbering to the M36 cell index used by iv_tc outputs.

The iv_tc step-2/3/4 files key every value by ``idx0``, a 0-based global EASEv2 M36
cell index (``j_indg * 964 + i_indg``).  The CYGNSS_L1_AZ stats and the L1-area
mask use ``tile_axis_index`` 0..908 (one entry per domain tile).  This writes the
bridge between them, plus the L1-domain mask, so TC/IV maps can be drawn on the
same tiles as the O-F and ISMN work.

Output columns: tile_axis_index, tile_id, idx0, i_indg, j_indg, com_lon, com_lat,
n_l1_obs, in_l1_mask (n_l1_obs >= --min-l1-obs).

Run with the repo python, e.g.
  python -I projects/CYGNSS_L1_AZ/scripts/build_iv_tc_tile_map.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
sys.path.insert(0, str(PROJ.parent / "iv_tc"))  # this repo's own iv_tc package
from iv_tc.readers import M36_NX, read_tilecoord  # noqa: E402

LAND_DA = Path("/discover/nobackup/projects/land_da/cygl1_operator_test")
DEFAULT_TILECOORD = (
    LAND_DA / "OLv8_M36_AZ_fixedop/output/SMAP_EASEv2_M36_GLOBAL/rc_out/OLv8_M36_AZ_fixedop.ldas_tilecoord.bin"
)
DEFAULT_TILES = PROJ / "output/operator_skill/da_link/tile_coords.csv"
DEFAULT_COUNTS = PROJ / "output/ismn_fixedop_20200101_20221231/tile_l1_obs_count_2020_2022.csv"
DEFAULT_OUT = PROJ / "output/iv_tc_az/tile_map.csv"


def build(tilecoord_path, tiles_path, counts_path, min_l1_obs):
    tc = read_tilecoord(tilecoord_path)
    tiles = pd.read_csv(tiles_path)
    counts = pd.read_csv(counts_path)

    row_of = {int(t): k for k, t in enumerate(tc["tile_id"])}
    missing = [t for t in tiles.tile_id if int(t) not in row_of]
    if missing:
        raise ValueError(f"{len(missing)} domain tile_ids not in tilecoord file, e.g. {missing[:5]}")
    rows = np.array([row_of[int(t)] for t in tiles.tile_id])

    out = pd.DataFrame(
        {
            "tile_axis_index": tiles.tile_axis_index.to_numpy(),
            "tile_id": tiles.tile_id.to_numpy(),
            "i_indg": tc["i_indg"][rows].astype(int),
            "j_indg": tc["j_indg"][rows].astype(int),
            "com_lon": tc["com_lon"][rows],
            "com_lat": tc["com_lat"][rows],
        }
    )
    out["idx0"] = out.j_indg * M36_NX + out.i_indg

    # Sanity: coordinates in the tilecoord file must agree with the stats bundle.
    dlon = np.abs(out.com_lon - tiles.com_lon.to_numpy()).max()
    dlat = np.abs(out.com_lat - tiles.com_lat.to_numpy()).max()
    if max(dlon, dlat) > 1e-3:
        raise ValueError(f"tile coordinates disagree with stats bundle (max {dlon:.2e} lon, {dlat:.2e} lat)")
    if out.idx0.nunique() != len(out):
        raise ValueError("several domain tiles share an M36 cell; idx0 is not a unique key")

    out = out.merge(counts, left_on="tile_axis_index", right_on="tile_index", how="left").drop(columns="tile_index")
    if out.n_l1_obs.isna().any():
        raise ValueError("L1 obs count missing for some domain tiles")
    out["in_l1_mask"] = out.n_l1_obs >= min_l1_obs
    cols = ["tile_axis_index", "tile_id", "idx0", "i_indg", "j_indg", "com_lon", "com_lat", "n_l1_obs", "in_l1_mask"]
    return out[cols]


def check_against_pairs(tile_map, pair_file):
    idx = np.load(pair_file)["idx0"]
    n_in = int(np.isin(idx, tile_map.idx0).sum())
    print(f"check {pair_file}: {n_in}/{idx.size} pair cells are domain tiles")
    if n_in != idx.size:
        raise SystemExit("some paired cells fall outside the domain tile map")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tilecoord", type=Path, default=DEFAULT_TILECOORD)
    p.add_argument("--tiles", type=Path, default=DEFAULT_TILES)
    p.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    p.add_argument("--min-l1-obs", type=int, default=100)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--check-pairs", type=Path, action="append", default=[], help="step-2 .npz to verify against")
    a = p.parse_args()

    tm = build(a.tilecoord, a.tiles, a.counts, a.min_l1_obs)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    tm.to_csv(a.out, index=False)
    print(f"{len(tm)} tiles, {int(tm.in_l1_mask.sum())} in L1 mask (>= {a.min_l1_obs} obs) -> {a.out}")
    for f in a.check_pairs:
        check_against_pairs(tm, f)


if __name__ == "__main__":
    main()
