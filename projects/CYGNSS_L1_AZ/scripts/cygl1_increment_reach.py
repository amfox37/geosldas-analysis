#!/usr/bin/env python3
"""
Plan B diagnostic: does the CygL1 analysis increment reach state variables it
has no direct physical information about (TSURF, TSOIL1, RZMC, PRMC), not just
the near-surface moisture (SFMC) it's plausibly informative for?

Since CygL1 is the ONLY assim=.true. species in this arm, any non-zero
ANA-FCST anywhere in the domain at a given cycle is necessarily caused by a
nearby CygL1 update (no other species can move the ensemble mean). So for
every CygL1 assim event we pull the same tile/cycle's inst3_1d_lndfcstana_Nt
increments (SFMC/RZMC/PRMC/TSURF/TSOIL1, ANA-FCST) and correlate them against
CygL1's own obs-fcst residual (d_f). A real, non-zero correlation for
TSURF/TSOIL1 (which the CygL1 radar-backscatter operator has no business
informing) would be direct evidence of a spurious/leaky cross-covariance
channel -- a candidate mechanism for the persistent local Tb harm found in
cygl1_local_skill_decay.py's radius-based event study.

Usage:
  cygl1_increment_reach.py --da-expid DAv8_M36_AZ_paired_cygl1_dense075_coh05 \\
                            --arm-tag dense075_coh05 --start 20200101 --end 20220101
"""
import argparse
import glob
import gzip
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_cygl1_paired_gain_consistency import read_tilecoord  # noqa: E402

DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
EXP_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test"
OUT_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output"

CYGL1_DESCR = "CYGNSS_L1_DDM3X5_CROP_SCALAR"
FILL_THRESH = 1e14
DEFAULT_SPECIES_DESCR = CYGL1_DESCR
DFLOOR = 0.1  # dB, matches the project's existing gain_proxy well-conditioned floor

STATE_VARS = ["SFMC", "RZMC", "PRMC", "TSURF", "TSOIL1"]


def month_range(start, end):
    ym = []
    y, m = start.year, start.month
    while (y, m) < (end.year, end.month):
        ym.append((y, m))
        m += 1
        if m == 13:
            m = 1
            y += 1
    return ym


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--da-expid", default="DAv8_M36_AZ_paired_cygl1_dense075_coh05")
    ap.add_argument("--arm-tag", default="dense075_coh05")
    ap.add_argument("--start", default="20200101")
    ap.add_argument("--end", default="20220101")
    ap.add_argument("--species-descr", default=DEFAULT_SPECIES_DESCR,
                     help="comma-separated obsparam descr string(s) to pool as the "
                          "'own species' assim event, e.g. multiple SMAP orbit/pol species")
    ap.add_argument("--dfloor", type=float, default=DFLOOR,
                     help="well-conditioned |d_f| floor, in the species' own obs units "
                          "(default 0.1 dB, calibrated for CygL1 -- override for other units, "
                          "e.g. m3/m3 for SM species)")
    args = ap.parse_args()
    wanted_descrs = set(args.species_descr.split(","))

    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)

    exp_dir = os.path.join(EXP_ROOT, args.da_expid)
    tilecoord_file = os.path.join(exp_dir, "output", DOMAIN, "rc_out", f"{args.da_expid}.ldas_tilecoord.bin")
    ofa_root = os.path.join(exp_dir, "output", DOMAIN, "ana", "ens_avg")
    incr_root = os.path.join(exp_dir, "output", DOMAIN, "cat", "ens_avg")

    tc = read_tilecoord(tilecoord_file)
    n_tile = tc["N_tile"]

    rows = []
    n_cycles = 0
    n_missing_incr = 0

    for yy, mm in month_range(start, end):
        ofa_files = sorted(glob.glob(os.path.join(ofa_root, f"Y{yy}", f"M{mm:02d}", f"{args.da_expid}.ens_avg.ldas_ObsFcstAna.*.nc4")))
        print(f"[{yy}-{mm:02d}] {len(ofa_files)} OFA files", flush=True)
        for fpath in ofa_files:
            fname = os.path.basename(fpath)
            stamp = fname.split(".")[-2]  # e.g. 20200101_0300z
            with nc.Dataset(fpath) as f:
                descr = np.array(f.variables["obsparam_descr"][:])
                spid = np.array(f.variables["obsparam_species_id"][:])
                name_to_id = {str(d): int(s) for d, s in zip(descr, spid)}
                wanted_ids = {name_to_id[d] for d in wanted_descrs if d in name_to_id}
                if not wanted_ids:
                    continue
                species = np.array(f.variables["species"][:])
                assim_flag = np.array(f.variables["assim_flag"][:])
                mask = np.isin(species, list(wanted_ids)) & (assim_flag == 1)
                if not np.any(mask):
                    continue
                tilenum = np.array(f.variables["tilenum"][:])[mask]
                obs = np.array(f.variables["obs"][:])[mask]
                fcst = np.array(f.variables["fcst"][:])[mask]
                bad = (np.abs(obs) > FILL_THRESH) | (np.abs(fcst) > FILL_THRESH)
                row_idx = np.clip(tilenum - 1, 0, n_tile - 1)
                bad = bad | (tilenum - 1 < 0) | (tilenum - 1 >= n_tile)

            incr_path = os.path.join(incr_root, f"Y{yy}", f"M{mm:02d}", f"{args.da_expid}.inst3_1d_lndfcstana_Nt.{stamp}.nc4")
            if not os.path.exists(incr_path):
                n_missing_incr += 1
                continue
            n_cycles += 1
            with nc.Dataset(incr_path) as fi:
                incr = {}
                for v in STATE_VARS:
                    fcst_arr = np.array(fi.variables[f"{v}_FCST"][0, :])
                    ana_arr = np.array(fi.variables[f"{v}_ANA"][0, :])
                    incr[v] = ana_arr - fcst_arr

            for k in range(len(tilenum)):
                if bad[k]:
                    continue
                ri = row_idx[k]
                rows.append((float(obs[k] - fcst[k]), *[float(incr[v][ri]) for v in STATE_VARS]))

    print(f"n_cycles_with_cygl1_events={n_cycles}, n_missing_incr_files={n_missing_incr}", flush=True)
    df = pd.DataFrame(rows, columns=["d_f"] + STATE_VARS)
    print(f"Total CygL1 assim events joined to increments: {len(df)}")

    well = df[df["d_f"].abs() >= args.dfloor]
    print(f"Well-conditioned (|d_f|>={args.dfloor}): {len(well)}")

    lines = []
    lines.append(f"CYGNSS L1 increment-reach diagnostic: {args.arm_tag}")
    lines.append(f"period: {args.start}-{args.end}, N events = {len(df)}, well-conditioned N = {len(well)}")
    lines.append("")
    lines.append("Does the CygL1 increment reach state variables beyond near-surface moisture?")
    lines.append(f"{'var':10s}{'mean_incr':>14s}{'std_incr':>14s}{'corr(d_f)':>14s}{'corr(|d_f|)':>14s}")
    for v in STATE_VARS:
        mean_i = well[v].mean()
        std_i = well[v].std()
        corr = well[["d_f", v]].corr().iloc[0, 1]
        corr_abs = np.corrcoef(well["d_f"].abs(), well[v].abs())[0, 1]
        lines.append(f"{v:10s}{mean_i:>14.6f}{std_i:>14.6f}{corr:>14.4f}{corr_abs:>14.4f}")

    summary_text = "\n".join(lines)
    print()
    print(summary_text)

    os.makedirs(OUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUT_DIR, f"cygl1_increment_reach_{args.arm_tag}.csv.gz")
    summary_path = os.path.join(OUT_DIR, f"cygl1_increment_reach_{args.arm_tag}_summary.txt")
    with gzip.open(csv_path, "wt") as fo:
        fo.write(f"# {args.da_expid}, {args.start}-{args.end}\n")
        df.to_csv(fo, index=False)
    with open(summary_path, "w") as fo:
        fo.write(summary_text + "\n")
    print(f"\nWrote {csv_path} and {summary_path}")


if __name__ == "__main__":
    main()
