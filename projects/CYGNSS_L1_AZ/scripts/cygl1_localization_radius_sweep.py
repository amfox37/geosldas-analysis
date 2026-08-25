#!/usr/bin/env python3
"""
Follow-up to the paired thinning-density experiment: instead of thinning obs
(as DA-sparse/intermediate/dense did), test whether tightening localization
(xcompact/ycompact) on the FULL/unthinned CYGNSS L1 obs stream can reproduce
DA-intermediate's local-interaction-count profile (median=1, p90=2, p99=3,
achieved via min_sep=2.40deg thinning at the actual xcompact=1.25deg).

Constraint (verified against clsm_ensupd_upd_routines.F90 this session):
update_type=12 in this project's LDASsa_SPECIAL_inputs_ensupd.nml puts
xcompact/ycompact under check_compact()'s hard floor of 2x the largest
relevant xcorr/ycorr. The binding floor here is the forcing-perturbation
correlation scale (pcp/sw/lw, xcorr=ycorr=0.5deg) -> floor = 1.0deg. Below
that the GEOSldas binary aborts (ldas_abort) before running at all. So only
xcompact in [1.0, 1.25] is actually runnable; this sweep asks whether the
lower end of that narrow, code-enforced range gets anywhere close to
intermediate's target on the full obs stream.

Reuses (imports, does not duplicate) load_all_obs/local_interaction_counts
from thin_cygl1_nested_density_6mo.py -- monkeypatches module-level
XCOMPACT_DEG/YCOMPACT_DEG per sweep point (the function reads those globals,
so no signature change needed).
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_cygl1_nested_density_6mo as thin  # noqa: E402

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
OUT_CSV = os.path.join(OUT_DIR, "cygl1_localization_radius_sweep.csv")

# xcompact/ycompact values to test on the FULL (dense) obs stream.
# 1.25 = current/actual value (reference point, matches the paired-experiment dense arm)
# 1.00 = hard floor from check_compact() against forcing-perturbation xcorr=0.5deg
SWEEP_DEG = [1.25, 1.20, 1.10, 1.05, 1.00]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading full Jan-Oct 2020 CYGNSS L1 obs stream...", flush=True)
    t0 = time.time()
    df, src_paths, dates = thin.load_all_obs()
    print(f"  loaded {len(df)} obs in {time.time()-t0:.0f}s", flush=True)

    dense_mask = np.ones(len(df), dtype=bool)

    rows = []
    for deg in SWEEP_DEG:
        thin.XCOMPACT_DEG = deg
        thin.YCOMPACT_DEG = deg
        print(f"\nxcompact=ycompact={deg:.2f}deg (interaction threshold = 2x = {2*deg:.2f}deg)...",
              flush=True)
        t0 = time.time()
        counts = thin.local_interaction_counts(df, dense_mask)
        dt = time.time() - t0
        c = counts[dense_mask]
        median_c = float(np.median(c))
        p90_c = float(np.percentile(c, 90))
        p99_c = float(np.percentile(c, 99))
        frac_zero = float((c == 0).mean())
        print(f"  [{dt:.0f}s] median={median_c:.2f}, p90={p90_c:.2f}, p99={p99_c:.2f}, "
              f"max={c.max()}, frac_zero={frac_zero:.3f}", flush=True)
        rows.append((deg, 2 * deg, len(c), median_c, p90_c, p99_c, int(c.max()), frac_zero))

    import pandas as pd
    result = pd.DataFrame(rows, columns=["xcompact_ycompact_deg", "interaction_threshold_deg",
                                          "n_obs", "median_local_count", "p90_local_count",
                                          "p99_local_count", "max_local_count", "frac_zero"])
    result.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")
    print(result.to_string(index=False))
    print("\nFor reference, DA-intermediate's target (via thinning, at the ACTUAL xcompact=1.25deg): "
          "median=1, p90=2, p99=3 (locked-in calibration).")
    print("check_compact() hard floor in this project's config: xcompact/ycompact >= 1.00deg "
          "(2x forcing-perturbation xcorr=0.5deg) -- below that GEOSldas aborts before running.")


if __name__ == "__main__":
    main()
