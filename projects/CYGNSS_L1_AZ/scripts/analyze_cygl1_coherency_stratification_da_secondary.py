#!/usr/bin/env python3
"""
Step 2b/4/5 secondary check: same coherency_state/coherency_ratio stratification and
within-tile control as analyze_cygl1_coherency_stratification.py, but on the DA run
(DAv8_M36_all_sensors_AZ_scaled_cygl1assim) join produced by
extract_cygl1_coherency_join_da_secondary.py.

SECONDARY EVIDENCE ONLY -- this run actually assimilates CYGNSS_L1 (assim_flag==1 for
all matched rows, obsparam_errstd=4.4 dB), so O-F here reflects post-increment feedback
from prior cycles, not pure operator fit. See runs/cygl1_coherency_stratification.md for
how this compares to the OL-run (primary) result.
"""
import gzip
import io

import numpy as np
import pandas as pd
from scipy import stats

CSV_PATH = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygnss_l1_coherency_joined_202001_202002_DAv8_secondary.csv.gz"
OUT_TXT = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygl1_coherency_stratification_DAv8_secondary_summary.txt"
STATE_NAMES = {0: "not_coherent", 1: "coherent", 2: "mixed", 3: "indeterminate"}


def read_csv_skip_header(path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    n_hdr = sum(1 for line in lines if line.startswith("#"))
    return pd.read_csv(io.StringIO("".join(lines[n_hdr:])))


def main():
    df = read_csv_skip_header(CSV_PATH)
    df["omf"] = df["obs_db"] - df["fcst_db"]
    lines = []
    lines.append(f"Loaded {len(df)} joined CYGNSS L1 obs (DA run, secondary), {df['tile_id'].nunique()} distinct tiles.")
    lines.append(f"Pooled O-F: mean={df['omf'].mean():.4f}  median={df['omf'].median():.4f}  sd={df['omf'].std():.4f}")
    lines.append(f"assim_flag value counts: {dict(df['assim_flag'].value_counts())}")

    lines.append("\n=== Bucketed by coherency_state ===")
    for key, sub in df.groupby("coherency_state"):
        lines.append(f"\n-- state={key} ({STATE_NAMES.get(key,key)}) -- N={len(sub)}")
        lines.append(f"   mean O-F={sub['omf'].mean():.4f}  median O-F={sub['omf'].median():.4f}  sd O-F={sub['omf'].std():.4f}")
        lines.append(f"   mean fcst={sub['fcst_db'].mean():.4f}  n_distinct_tiles={sub['tile_id'].nunique()}")

    coh_tiles = set(df.loc[df.coherency_state == 1, "tile_id"].unique())
    noncoh_tiles = set(df.loc[df.coherency_state == 0, "tile_id"].unique())
    both = sorted(coh_tiles & noncoh_tiles)
    lines.append(f"\n=== Within-tile control (DA run, secondary) ===")
    lines.append(f"Tiles with both categories: {len(both)}")
    rows = []
    for tid in both:
        c = df[(df.tile_id == tid) & (df.coherency_state == 1)]
        n = df[(df.tile_id == tid) & (df.coherency_state == 0)]
        rows.append((tid, len(c), len(n), c.omf.mean() - n.omf.mean(),
                     c.omf.abs().mean() - n.omf.abs().mean(), c.omf.std() - n.omf.std()))
    pt = pd.DataFrame(rows, columns=["tile_id", "n_c", "n_n", "signed_diff", "absomf_diff", "sd_diff"])
    if len(pt):
        lines.append(f"signed_diff (coherent-incoherent mean O-F): mean={pt.signed_diff.mean():.4f} median={pt.signed_diff.median():.4f}")
        lines.append(f"absomf_diff: mean={pt.absomf_diff.mean():.4f} median={pt.absomf_diff.median():.4f}  frac_tighter={ (pt.absomf_diff<0).mean():.3f}")
        sdv = pt.sd_diff.dropna()
        lines.append(f"sd_diff: mean={sdv.mean():.4f} median={sdv.median():.4f}  frac_lower_spread={ (sdv<0).mean():.3f}  N={len(sdv)}")
        if len(pt) >= 8:
            w1 = stats.wilcoxon(pt.signed_diff)
            w2 = stats.wilcoxon(pt.absomf_diff)
            lines.append(f"Wilcoxon signed_diff p={w1.pvalue:.4g}   Wilcoxon absomf_diff p={w2.pvalue:.4g}")
        if len(sdv) >= 8:
            w3 = stats.wilcoxon(sdv)
            lines.append(f"Wilcoxon sd_diff p={w3.pvalue:.4g}")
        pt.to_csv("/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/cygl1_coherency_within_tile_pairs_DAv8_secondary.csv", index=False)

    text = "\n".join(lines)
    print(text)
    with open(OUT_TXT, "w") as fo:
        fo.write(text + "\n")


if __name__ == "__main__":
    main()
