#!/usr/bin/env python3
"""
Event-level diagnostics for the gated-dense arm's 22-month extension
(Jan2020-Dec2021), same six-check pattern as cygl1_paired_event_diagnostics_22mo.py.
Full period + monthly-condensed + yearly-split tables, to check whether the
2020 finding ("gate helps substantially but doesn't close the gap to
intermediate's event-level self-consistency: ~70-75% improved/~18-20%
wrong-direction band, no drift within 2020") continues unchanged into 2021.

Input: cygl1_gated_dense_gain_consistency_22mo_ofa.csv.gz
       (extract_cygl1_gated_dense_gain_consistency_22mo.py)
Outputs: cygl1_gated_dense_event_diagnostics_22mo.csv / _summary.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cygl1_paired_event_diagnostics_22mo import verbose_block, condensed_row  # noqa: E402

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
IN_CSV = os.path.join(OUT_DIR, "cygl1_gated_dense_gain_consistency_22mo_ofa.csv.gz")
EVENT_CSV = os.path.join(OUT_DIR, "cygl1_gated_dense_event_diagnostics_22mo.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_gated_dense_event_diagnostics_22mo_summary.txt")

DFLOOR = 0.1  # dB
ARM = "gated_dense"


def main():
    df = pd.read_csv(IN_CSV, comment="#")
    df["datetime"] = pd.to_datetime(df["datetime"])
    print(f"Loaded {len(df)} events, arm={ARM}")

    df["d_f"] = df["obs"] - df["fcst"]
    df["d_a"] = df["obs"] - df["ana"]
    df["delta_h"] = df["ana"] - df["fcst"]
    assert np.allclose(df["d_f"], df["d"]), "d_f mismatch vs precomputed d"
    assert np.allclose(df["delta_h"], df["actual_incr"]), "delta_h mismatch vs precomputed actual_incr"

    well_cond = df["d_f"].abs() >= DFLOOR
    df["gain_proxy"] = np.where(well_cond, df["delta_h"] / df["d_f"], np.nan)
    df["abs_d_f"] = df["d_f"].abs()
    df["abs_d_a"] = df["d_a"].abs()
    df["improved"] = df["abs_d_a"] < df["abs_d_f"]
    df["toward_obs"] = df["d_f"] * df["delta_h"] > 0
    df["toward_obs_or_zero"] = df["d_f"] * df["delta_h"] >= 0
    df["month"] = df["datetime"].dt.strftime("%Y-%m")

    df = df.sort_values(["datetime", "tilenum"]).reset_index(drop=True)

    lines = []
    lines.append("CYGNSS L1 gated-dense arm: event-level diagnostics, 22-month extension")
    lines.append("arm: gated_dense (full/unthinned stream, hard-gated binary 406206a)")
    lines.append("period: Jan2020-Dec2021 (full calendar years, no spin-up exclusion)")
    lines.append(f"gain_proxy computed only for |d_f| >= {DFLOOR} dB (well-conditioned)")
    lines.append("")

    lines.append("=" * 100)
    lines.append("FULL-PERIOD (Jan2020-Dec2021 pooled) -- verbose six-check block")
    lines.append("=" * 100)
    lines.extend(verbose_block(df, f"arm={ARM}"))

    lines.append("=" * 100)
    lines.append("MONTHLY -- condensed (watch for any drift/degradation as cycles accumulate)")
    lines.append("=" * 100)
    for month, g in df.groupby("month"):
        lines.append(condensed_row(g, month))
    lines.append("")

    lines.append("=" * 100)
    lines.append("YEARLY SPLIT -- condensed")
    lines.append("=" * 100)
    df["year"] = df["datetime"].dt.year
    for year, g in df.groupby("year"):
        lines.append(condensed_row(g, str(year)))
    lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    event_cols = ["arm", "datetime", "tilenum", "tile_id", "lon", "lat", "obs", "fcst", "ana",
                  "obsvar", "fcstvar", "K", "d_f", "d_a", "delta_h", "gain_proxy",
                  "improved", "toward_obs", "month"]
    df[event_cols].to_csv(EVENT_CSV, index=False)
    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {EVENT_CSV} and {SUMMARY_PATH}")
    return df


if __name__ == "__main__":
    main()
