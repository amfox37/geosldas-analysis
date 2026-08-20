#!/usr/bin/env python3
"""
Section 8 sanity check: reproduce the OL open-loop temporal-stats numbers for
CYGNSS_L1_DDM3X5_CROP_SCALAR and CYGNSS_SM_6hr from the paired OFA extraction
(cygnss_l1_l3_paired_ofa_2020.csv.gz), and compare against the brief's target
values:

species: O_stdv, F_stdv, OmF_stdv, total N_data
CYGNSS_L1_DDM3X5_CROP_SCALAR: 2.701, 2.738, 3.095, 171483
CYGNSS_SM_6hr: 0.026, 0.029, 0.023, 197211

Method: per-tile O_stdv/F_stdv/OmF_stdv (population, ddof=0) computed over
non-NaN obs & fcst pairs, tiles with N_data >= 20 kept, then N_data-weighted
mean across tiles. Also reports the ddof=1 (sample) variant and a pooled
(variance-combination) variant, since the brief does not specify which
weighting convention was used upstream.
"""
import pandas as pd
import numpy as np
import gzip

CSV_PATH = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygnss_l1_l3_paired_ofa_2020.csv.gz"

TARGETS = {
    "CYGNSS_L1_DDM3X5_CROP_SCALAR": dict(O_stdv=2.701, F_stdv=2.738, OmF_stdv=3.095, N=171483),
    "CYGNSS_SM_6hr": dict(O_stdv=0.026, F_stdv=0.029, OmF_stdv=0.023, N=197211),
}


def read_csv_skip_header(path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    n_hdr = 0
    for line in lines:
        if line.startswith("#"):
            n_hdr += 1
        else:
            break
    from io import StringIO
    return pd.read_csv(StringIO("".join(lines[n_hdr:])))


def per_tile_stats(df, ddof=0):
    out = []
    for tile_id, g in df.groupby("tile_id"):
        valid = g.dropna(subset=["obs", "fcst"])
        n = len(valid)
        if n < 2:
            continue
        o = valid["obs"].values
        fc = valid["fcst"].values
        omf = o - fc
        out.append((tile_id, n, o.std(ddof=ddof), fc.std(ddof=ddof), omf.std(ddof=ddof)))
    return pd.DataFrame(out, columns=["tile_id", "N", "O_stdv", "F_stdv", "OmF_stdv"])


def weighted_mean(vals, weights):
    return float(np.sum(vals * weights) / np.sum(weights))


def pooled_stdv(df, col, nmin=20):
    """Pool all obs (not per-tile stdv averaging) across tiles with N>=20, compute grand stdv."""
    pass


def main():
    df = read_csv_skip_header(CSV_PATH)
    print(f"Loaded {len(df)} rows from {CSV_PATH}")
    print(df['species_name'].value_counts())
    print()

    for sp_name, target in TARGETS.items():
        sub = df[df["species_name"] == sp_name].copy()
        print(f"=== {sp_name} ===  (n_obs total in extraction, incl. tiles<20: {len(sub.dropna(subset=['obs','fcst']))})")

        for ddof in (0, 1):
            ts = per_tile_stats(sub, ddof=ddof)
            keep = ts[ts["N"] >= 20]
            n_total = int(keep["N"].sum())
            o_w = weighted_mean(keep["O_stdv"].values, keep["N"].values)
            f_w = weighted_mean(keep["F_stdv"].values, keep["N"].values)
            omf_w = weighted_mean(keep["OmF_stdv"].values, keep["N"].values)
            print(f"  [N_data-weighted mean of per-tile stdv, ddof={ddof}]")
            print(f"    N_tiles kept (N>=20): {len(keep)} / {len(ts)}")
            print(f"    total N_data        : {n_total}   (target {target['N']})")
            print(f"    O_stdv              : {o_w:.4f}   (target {target['O_stdv']})")
            print(f"    F_stdv              : {f_w:.4f}   (target {target['F_stdv']})")
            print(f"    OmF_stdv            : {omf_w:.4f}   (target {target['OmF_stdv']})")
            r_implied = (o_w**2 + f_w**2 - omf_w**2) / (2 * o_w * f_w)
            print(f"    implied correlation : {r_implied:.4f}")

        # Also: pooled (grand) stdv across all valid obs in kept tiles (alternative convention)
        ts0 = per_tile_stats(sub, ddof=0)
        keep_tiles = set(ts0[ts0["N"] >= 20]["tile_id"])
        pooled = sub[sub["tile_id"].isin(keep_tiles)].dropna(subset=["obs", "fcst"])
        o = pooled["obs"].values
        fc = pooled["fcst"].values
        omf = o - fc
        print(f"  [pooled/grand stdv over all obs in kept tiles, ddof=0]")
        print(f"    N_data   : {len(pooled)}")
        print(f"    O_stdv   : {o.std(ddof=0):.4f}")
        print(f"    F_stdv   : {fc.std(ddof=0):.4f}")
        print(f"    OmF_stdv : {omf.std(ddof=0):.4f}")
        r_implied2 = (o.std()**2 + fc.std()**2 - omf.std()**2) / (2*o.std()*fc.std())
        print(f"    implied correlation : {r_implied2:.4f}")
        print()


if __name__ == "__main__":
    main()
