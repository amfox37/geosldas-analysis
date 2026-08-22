#!/usr/bin/env python3
"""
Follow-up to extract_cygl1_singleobs_v2_gain_consistency.py: for the
genuinely-isolated CYGNSS L1 obs from the RE-THINNED, window-bug-fixed rerun
(Step 1.5 v2, DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs_v2, job
58019976), join the same-cycle catch_progn_incr collection at the owner tile
to see the actual state-space increment split.

Purpose: the obs-space gain-consistency check found actual_incr consistently
~50% of the naive single-state-variable K*d prediction (slope=0.50, R2=0.96),
with correct sign in 46/47 obs. Leading hypothesis: N_state=2 (srfexc+rzexc,
confirmed in obs_error_variance_diagnosis.md 2026-08-20 entry via get_select_
species/found_Tb_obs code read) means the true update spreads correction
across two state variables in a way the naive scalar K*d doesn't capture,
which could show up as an obs-space undercount even with a "correct" gain.
This script pulls the actual per-state-variable increments to check that
directly -- and, as a side effect, empirically confirms/falsifies N_state=2
by checking whether CATDEF/snow/temperature increments are exactly zero
(non-peat tiles) at every kept-obs tile/cycle.

Sign convention (CLSM prognostics), confirmed 2026-08-21 against
M21C_ls/scripts/sum_da_analysis_increments.py's water-budget formula
(dWTOT_from_increments = -CATDEF_INCR + RZEXC_INCR + SRFEXC_INCR): CATDEF is
a deficit (more positive = drier), but RZEXC/SRFEXC are "excess" terms with
the OPPOSITE sign convention (more POSITIVE = wetter) -- confirmed empirically
below too (0/47 "correct" under the initially-assumed opposite convention is
itself the tell: exactly 0, not noisy, means the assumed sign was simply
inverted). Since dH/dsfmc > 0 everywhere in this domain (Task 0, operator
Jacobian sign check: wetter soil -> higher dB), a positive d=obs-fcst (obs
wetter/higher-dB than fcst) should correspond to a wetting increment, i.e.
POSITIVE RZEXC_INCR/SRFEXC_INCR, at that tile/cycle.
"""
import glob
import os
import gzip

import numpy as np
import netCDF4 as nc
import pandas as pd

EXP_ID = "DAv8_M36_all_sensors_AZ_scaled_cygl1assim_singleobs_v2"
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
EXP_DIR = os.path.join(EXP_ROOT, EXP_ID)
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"

OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"
GAIN_CSV = os.path.join(OUT_DIR, "cygl1_singleobs_v2_gain_consistency_ofa.csv.gz")
OUT_CSV = os.path.join(OUT_DIR, "cygl1_singleobs_v2_state_incr.csv.gz")
SUMMARY_PATH = os.path.join(OUT_DIR, "cygl1_singleobs_v2_state_incr_summary.txt")

# All prognostic increment fields in the collection, so we can empirically
# confirm N_state=2 (everything except SRFEXC/RZEXC should be exactly zero).
INCR_FIELDS = [
    "CAPAC_INCR", "CATDEF_INCR", "RZEXC_INCR", "SRFEXC_INCR",
    "GHTCNT1_INCR", "GHTCNT2_INCR", "GHTCNT3_INCR", "GHTCNT4_INCR", "GHTCNT5_INCR", "GHTCNT6_INCR",
    "TCFSAT_INCR", "TCFTRN_INCR", "TCFWLT_INCR",
    "QCFSAT_INCR", "QCFTRN_INCR", "QCFWLT_INCR",
    "WESNN1_INCR", "WESNN2_INCR", "WESNN3_INCR",
    "SNDZN1_INCR", "SNDZN2_INCR", "SNDZN3_INCR",
    "HTSNNN1_INCR", "HTSNNN2_INCR", "HTSNNN3_INCR",
]


def incr_fname(dt):
    return f"{EXP_ID}.catch_progn_incr.{dt.strftime('%Y%m%d')}_{dt.strftime('%H%M')}z.nc4"


def main():
    df = pd.read_csv(GAIN_CSV, comment="#")
    df["datetime"] = pd.to_datetime(df["datetime"])
    print(f"Loaded {len(df)} isolated obs from {GAIN_CSV}", flush=True)

    incr_root = os.path.join(EXP_DIR, "output", DOMAIN, "cat", "ens_avg")

    rows = []
    missing = []
    for _, r in df.iterrows():
        dt = r["datetime"]
        yyyy, mm = dt.strftime("%Y"), dt.strftime("%m")
        fpath = os.path.join(incr_root, f"Y{yyyy}", f"M{mm}", incr_fname(dt))
        if not os.path.exists(fpath):
            missing.append((dt, r["tilenum"], fpath))
            continue
        row_idx = int(r["tilenum"]) - 1
        with nc.Dataset(fpath) as f:
            vals = {}
            for field in INCR_FIELDS:
                v = float(np.array(f.variables[field][0, row_idx]))
                if v > 1e14:
                    v = np.nan
                vals[field] = v
        out = {
            "datetime": dt, "tilenum": int(r["tilenum"]), "tile_id": int(r["tile_id"]),
            "d": float(r["d"]), "K": float(r["K"]), "pred_incr": float(r["pred_incr"]),
            "actual_incr": float(r["actual_incr"]),
        }
        out.update(vals)
        rows.append(out)

    if missing:
        print(f"WARNING: {len(missing)} obs had no matching catch_progn_incr file:")
        for dt, tn, fp in missing:
            print(f"  {dt} tile {tn}: {fp}")

    out_df = pd.DataFrame(rows)
    print(f"Joined {len(out_df)} / {len(df)} obs to catch_progn_incr", flush=True)

    # Empirical N_state check: are all non-srfexc/rzexc increments exactly zero?
    other_fields = [f for f in INCR_FIELDS if f not in ("SRFEXC_INCR", "RZEXC_INCR")]
    other_nonzero = (out_df[other_fields].abs() > 0).any(axis=1)
    n_other_nonzero = int(other_nonzero.sum())

    srfexc = out_df["SRFEXC_INCR"]
    rzexc = out_df["RZEXC_INCR"]

    # Sign check: wetting increment (POSITIVE SRFEXC/RZEXC, per sum_da_analysis_increments.py
    # convention) should coincide with positive d (obs wetter than fcst).
    expected_sign = np.sign(out_df["d"])  # expected sign of SRFEXC/RZEXC_INCR
    srfexc_sign_ok = np.sign(srfexc) == expected_sign
    rzexc_sign_ok = np.sign(rzexc) == expected_sign

    summary_lines = []
    summary_lines.append("CYGNSS L1 state-space increment check, Step 1.5 v2 follow-up (window-bug-fixed)")
    summary_lines.append(f"experiment: {EXP_ID}")
    summary_lines.append(f"N obs joined = {len(out_df)} (of {len(df)} candidates, {len(missing)} missing)")
    summary_lines.append("")
    summary_lines.append(f"N_state=2 empirical check: obs with any nonzero increment OUTSIDE srfexc/rzexc = {n_other_nonzero} / {len(out_df)}")
    summary_lines.append("  (CONFIRMS N_state=2 if this is 0 -- catdef/snow/temperature/canopy never touched)")
    summary_lines.append("")
    summary_lines.append(f"SRFEXC_INCR: mean={srfexc.mean():.5f}, std={srfexc.std():.5f}, |mean|={srfexc.abs().mean():.5f} kg/m2")
    summary_lines.append(f"RZEXC_INCR:  mean={rzexc.mean():.5f}, std={rzexc.std():.5f}, |mean|={rzexc.abs().mean():.5f} kg/m2")
    summary_lines.append(f"ratio |SRFEXC_INCR| / |RZEXC_INCR| (mean of ratios, obs with nonzero rzexc) = "
                          f"{(srfexc.abs() / rzexc.abs().replace(0, np.nan)).mean():.3f}")
    summary_lines.append("")
    summary_lines.append(f"Sign check (expect SRFEXC/RZEXC_INCR sign = +sign(d), since wetter=more positive [RZEXC/SRFEXC] and dH/dsfmc>0):")
    summary_lines.append(f"  SRFEXC_INCR correct sign: {int(srfexc_sign_ok.sum())} / {len(out_df)} ({srfexc_sign_ok.mean():.3f})")
    summary_lines.append(f"  RZEXC_INCR  correct sign: {int(rzexc_sign_ok.sum())} / {len(out_df)} ({rzexc_sign_ok.mean():.3f})")
    summary_lines.append("")
    # zero-d edge cases
    nonzero_d = out_df["d"].abs() > 1e-6
    summary_lines.append(f"(sign checks restricted to nonzero d: {int(nonzero_d.sum())}/{len(out_df)} obs have |d|>1e-6)")
    srfexc_sign_ok_nz = srfexc_sign_ok[nonzero_d]
    rzexc_sign_ok_nz = rzexc_sign_ok[nonzero_d]
    summary_lines.append(f"  SRFEXC_INCR correct sign (|d|>0 only): {int(srfexc_sign_ok_nz.sum())} / {len(srfexc_sign_ok_nz)} ({srfexc_sign_ok_nz.mean():.3f})")
    summary_lines.append(f"  RZEXC_INCR  correct sign (|d|>0 only): {int(rzexc_sign_ok_nz.sum())} / {len(rzexc_sign_ok_nz)} ({rzexc_sign_ok_nz.mean():.3f})")

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    with gzip.open(OUT_CSV, "wt") as fo:
        fo.write("# CYGNSS L1 state-space increment check, Step 1.5 v2 follow-up (window-bug-fixed)\n")
        fo.write(f"# experiment: {EXP_ID}\n")
        out_df.to_csv(fo, index=False)

    with open(SUMMARY_PATH, "w") as fo:
        fo.write(summary_text + "\n")

    print(f"\nWrote {OUT_CSV} and {SUMMARY_PATH}")
    return out_df


if __name__ == "__main__":
    main()
