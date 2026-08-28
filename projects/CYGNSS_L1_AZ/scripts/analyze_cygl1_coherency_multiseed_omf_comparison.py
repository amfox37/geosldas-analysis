#!/usr/bin/env python3
"""
CYGNSS L1 coherency-screening experiment, multi-seed re-draw: same aggregation
as analyze_cygl1_coherency_omf_comparison.py (N-weighted pooled bias/stdv per
species group: Tb 1-8, SM 9-12, CygL1 13; skill delta = OmF(OL-xmask) -
OmF(arm)), but across Arm A + THREE independent Arm-B random-matched draws
(seed1=20260827, seed2=20260828, seed3=20260829), to check whether Arm A's
CygL1-skill advantage over Arm B is a real effect or a single-seed artifact.

Does not touch analyze_cygl1_coherency_omf_comparison.py or its outputs.
"""
import pickle
import numpy as np

STATS_DIR = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/postproc_paired_density/stats_output/"

SPECIES_NAMES = [
    "SMOS_Tbh_A", "SMOS_Tbh_D", "SMOS_Tbv_A", "SMOS_Tbv_D",
    "SMAP_Tbh_A", "SMAP_Tbh_D", "SMAP_Tbv_A", "SMAP_Tbv_D",
    "ASCAT_A", "ASCAT_B", "ASCAT_C", "CYGNSS_SM_6hr",
    "CYGNSS_L1",
]
GROUPS = {"Tb (K)": list(range(0, 8)), "SM (m3/m3)": list(range(8, 12)), "CygL1 (dB)": [12]}

ARM_FILES = {
    "A":     ("spatial_stats_DA_coherency_screened_202001_202012.pkl",
              "spatial_stats_OL_paired_monitor_xmask_coherency_A_202001_202012.pkl"),
    "B_s1":  ("spatial_stats_DA_coherency_randmatch_202001_202012.pkl",
              "spatial_stats_OL_paired_monitor_xmask_coherency_B_202001_202012.pkl"),
    "B_s2":  ("spatial_stats_DA_coherency_randmatch_seed2_202001_202012.pkl",
              "spatial_stats_OL_paired_monitor_xmask_coherency_B_seed2_202001_202012.pkl"),
    "B_s3":  ("spatial_stats_DA_coherency_randmatch_seed3_202001_202012.pkl",
              "spatial_stats_OL_paired_monitor_xmask_coherency_B_seed3_202001_202012.pkl"),
}
B_ARMS = ["B_s1", "B_s2", "B_s3"]


def load(filename):
    with open(STATS_DIR + filename, "rb") as f:
        return pickle.load(f)


def group_bias_stdv(d, indices, n_months=12):
    N = d["N_data"][:n_months, indices]
    mean = d["OmF_mean"][:n_months, indices]
    stdv = d["OmF_stdv"][:n_months, indices]
    valid = (N > 0) & ~np.isnan(mean)
    if not valid.any():
        return np.nan, np.nan, 0
    Nf, meanf, stdvf = N[valid], mean[valid], stdv[valid]
    Ntot = Nf.sum()
    pooled_mean = (Nf * meanf).sum() / Ntot
    pooled_var = (Nf * (stdvf ** 2 + meanf ** 2)).sum() / Ntot - pooled_mean ** 2
    pooled_stdv = np.sqrt(max(pooled_var, 0))
    return pooled_mean, pooled_stdv, int(Ntot)


def main():
    data = {arm: (load(own_fn), load(xmask_fn)) for arm, (own_fn, xmask_fn) in ARM_FILES.items()}

    # per-arm, per-group: bias/stdv of own OmF, OL-xmask OmF, and the skill deltas
    results = {}  # results[arm][group] = dict(...)
    for arm, (own, xmask) in data.items():
        results[arm] = {}
        for gname, idxs in GROUPS.items():
            m, s, n = group_bias_stdv(own, idxs)
            mX, sX, nX = group_bias_stdv(xmask, idxs)
            d_bias = abs(mX) - abs(m)
            d_stdv = sX - s
            results[arm][gname] = dict(mean=m, stdv=s, n=n, meanX=mX, stdvX=sX, nX=nX,
                                        d_bias=d_bias, d_stdv=d_stdv)

    print("=" * 115)
    print("PER-ARM, PER-GROUP: own OmF, OL-xmask OmF, and skill deltas (Delta_stdv = stdv(OL-xmask) - stdv(arm))")
    print("=" * 115)
    for gname in GROUPS:
        print(f"\n--- {gname} ---")
        hdr = f"{'arm':10s} {'OmF(arm) mean/stdv':>24s} {'OmF(OL-xmask) mean/stdv':>26s} {'Delta_bias':>12s} {'Delta_stdv':>12s} {'N(arm)':>10s}"
        print(hdr)
        for arm in ["A"] + B_ARMS:
            r = results[arm][gname]
            own_str = f"{r['mean']:+.4f}/{r['stdv']:.4f}"
            xmask_str = f"{r['meanX']:+.4f}/{r['stdvX']:.4f}"
            print(f"{arm:10s} {own_str:>24s} {xmask_str:>26s} "
                  f"{r['d_bias']:+12.4f} {r['d_stdv']:+12.4f} {r['n']:>10d}")

    print()
    print("=" * 115)
    print("SEED-SPREAD SUMMARY: Arm A vs. mean/std/range of Delta_stdv across the 3 Arm-B seeds")
    print("=" * 115)
    verdicts = {}
    for gname in GROUPS:
        a_val = results["A"][gname]["d_stdv"]
        b_vals = np.array([results[arm][gname]["d_stdv"] for arm in B_ARMS])
        b_mean, b_std = b_vals.mean(), b_vals.std(ddof=1)
        b_min, b_max = b_vals.min(), b_vals.max()
        outside = (a_val < b_min) or (a_val > b_max)
        verdicts[gname] = outside
        print(f"\n{gname}:")
        print(f"  Delta_stdv: Arm A = {a_val:+.4f}")
        print(f"  Delta_stdv: Arm B seed1/2/3 = {b_vals[0]:+.4f} / {b_vals[1]:+.4f} / {b_vals[2]:+.4f}")
        print(f"  Delta_stdv: Arm B mean±std = {b_mean:+.4f} ± {b_std:.4f}  (range {b_min:+.4f} to {b_max:+.4f})")
        print(f"  Verdict: Arm A is {'OUTSIDE' if outside else 'WITHIN'} the 3-seed Arm-B range "
              f"({'real, seed-independent advantage' if outside else 'NOT distinguishable from Arm-B seed-to-seed noise'})")

    print()
    print("=" * 115)
    print("BOTTOM LINE")
    print("=" * 115)
    for gname in GROUPS:
        tag = "REAL (outside B-seed spread)" if verdicts[gname] else "WITHIN NOISE (overlaps B-seed spread)"
        print(f"  {gname:14s}: Arm A's advantage over Arm B is {tag}")

    return results


if __name__ == "__main__":
    main()
