#!/usr/bin/env python3
"""
CYGNSS L1 coherency-screening experiment: aggregate the postproc monthly
spatial_stats pkls (OmF_mean/OmF_stdv/N_data per species per month, Jan-Dec
2020) into a single N-weighted bias/stdv number per species (and per species
group: Tb 1-8, SM 9-12, CygL1 13), for:

  DA_coherency_screened / DA_coherency_randmatch          (own OmF, own obs)
  OL_paired_monitor_xmask_coherency_A / _B                (OL xmasked to each arm's obs)
  DA_gated_dense (Jan-Dec 2020 slice)                      (benchmark, not obs-count-matched)

Reports bias+stdv per feedback_da_evaluation_standard (mean, not just spread --
coherency_ratio is known to track a systematic bias, not just noise).
Prints skill deltas Delta_A = OmF(xmask_A) - OmF(A), Delta_B = OmF(xmask_B) - OmF(B).
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


def load(tag, filename):
    with open(STATS_DIR + filename, "rb") as f:
        d = pickle.load(f)
    return d


def pooled_bias_stdv(d, species_idx, n_months=12):
    """N-weighted pooled mean/std across the first n_months months for one species index."""
    N = d["N_data"][:n_months, species_idx]
    mean = d["OmF_mean"][:n_months, species_idx]
    stdv = d["OmF_stdv"][:n_months, species_idx]
    valid = (N > 0) & ~np.isnan(mean)
    if not valid.any():
        return np.nan, np.nan, 0
    N, mean, stdv = N[valid], mean[valid], stdv[valid]
    Ntot = N.sum()
    pooled_mean = (N * mean).sum() / Ntot
    # pooled variance from sub-sample means/vars (treats each month as a sub-sample)
    pooled_var = (N * (stdv ** 2 + mean ** 2)).sum() / Ntot - pooled_mean ** 2
    pooled_stdv = np.sqrt(max(pooled_var, 0))
    return pooled_mean, pooled_stdv, int(Ntot)


def group_bias_stdv(d, indices, n_months=12):
    """Nobs-weighted pooled bias/stdv across a group of species (valid within
    each unit group only -- Tb[K]/SM[m3m3]/CygL1[dB] are never mixed)."""
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
    files = {
        "A (screened, own OF)":     "spatial_stats_DA_coherency_screened_202001_202012.pkl",
        "B (randmatch, own OF)":    "spatial_stats_DA_coherency_randmatch_202001_202012.pkl",
        "OL-xmask-A":               "spatial_stats_OL_paired_monitor_xmask_coherency_A_202001_202012.pkl",
        "OL-xmask-B":               "spatial_stats_OL_paired_monitor_xmask_coherency_B_202001_202012.pkl",
        "gated-dense (own OF)":     "spatial_stats_DA_gated_dense_202001_202012.pkl",
    }
    data = {tag: load(tag, fn) for tag, fn in files.items()}

    print("=" * 100)
    print("PER-SPECIES bias (OmF_mean) and stdv (OmF_stdv), Jan-Dec 2020, N-weighted pooled over 12 months")
    print("=" * 100)
    hdr = f"{'species':26s}" + "".join(f"{tag:>26s}" for tag in files)
    print(hdr)
    for i, name in enumerate(SPECIES_NAMES):
        row = f"{name:26s}"
        for tag in files:
            mean, stdv, n = pooled_bias_stdv(data[tag], i)
            row += f"{f'{mean:+.4f}/{stdv:.4f}(N={n})':>26s}"
        print(row)

    print()
    print("=" * 100)
    print("GROUPED bias/stdv (Nobs-weighted within unit group), Jan-Dec 2020")
    print("=" * 100)
    print(f"{'group':16s}" + "".join(f"{tag:>26s}" for tag in files))
    group_results = {}
    for gname, idxs in GROUPS.items():
        row = f"{gname:16s}"
        group_results[gname] = {}
        for tag in files:
            mean, stdv, n = group_bias_stdv(data[tag], idxs)
            group_results[gname][tag] = (mean, stdv, n)
            row += f"{f'{mean:+.4f}/{stdv:.4f}(N={n})':>26s}"
        print(row)

    print()
    print("=" * 100)
    print("SKILL DELTAS: Delta = OmF(OL-xmask) - OmF(arm)  [positive = assim reduces OmF vs its own OL baseline]")
    print("(bias delta = |bias(OL-xmask)| - |bias(arm)| since bias sign is not meaningful to subtract directly;")
    print(" stdv delta = stdv(OL-xmask) - stdv(arm), positive = assim reduces spread)")
    print("=" * 100)
    for gname in GROUPS:
        mA, sA, nA = group_results[gname]["A (screened, own OF)"]
        mB, sB, nB = group_results[gname]["B (randmatch, own OF)"]
        mXA, sXA, nXA = group_results[gname]["OL-xmask-A"]
        mXB, sXB, nXB = group_results[gname]["OL-xmask-B"]
        d_bias_A = abs(mXA) - abs(mA)
        d_bias_B = abs(mXB) - abs(mB)
        d_stdv_A = sXA - sA
        d_stdv_B = sXB - sB
        print(f"{gname}:")
        print(f"  Arm A: OmF(OL-xmask-A)={mXA:+.4f}/{sXA:.4f}  OmF(A)={mA:+.4f}/{sA:.4f}  "
              f"Delta_bias_A={d_bias_A:+.4f}  Delta_stdv_A={d_stdv_A:+.4f}  "
              f"({'A REDUCES stdv' if d_stdv_A>0 else 'A INCREASES stdv'})")
        print(f"  Arm B: OmF(OL-xmask-B)={mXB:+.4f}/{sXB:.4f}  OmF(B)={mB:+.4f}/{sB:.4f}  "
              f"Delta_bias_B={d_bias_B:+.4f}  Delta_stdv_B={d_stdv_B:+.4f}  "
              f"({'B REDUCES stdv' if d_stdv_B>0 else 'B INCREASES stdv'})")
        print(f"  A vs B: Delta_stdv_A - Delta_stdv_B = {d_stdv_A - d_stdv_B:+.4f}  "
              f"({'A beats B (coherency screen adds skill beyond count)' if d_stdv_A > d_stdv_B + 1e-6 else 'A does NOT beat B'})")
        print()

    print("=" * 100)
    print("A/B vs gated-dense (own OmF, NOT obs-count-matched -- caveat per design spec)")
    print("=" * 100)
    for gname in GROUPS:
        mA, sA, nA = group_results[gname]["A (screened, own OF)"]
        mB, sB, nB = group_results[gname]["B (randmatch, own OF)"]
        mGD, sGD, nGD = group_results[gname]["gated-dense (own OF)"]
        print(f"{gname}: gated-dense={mGD:+.4f}/{sGD:.4f}(N={nGD})  A={mA:+.4f}/{sA:.4f}(N={nA})  "
              f"B={mB:+.4f}/{sB:.4f}(N={nB})")


if __name__ == "__main__":
    main()
