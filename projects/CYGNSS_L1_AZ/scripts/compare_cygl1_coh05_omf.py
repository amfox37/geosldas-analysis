#!/usr/bin/env python3
"""
CYGNSS L1 paired thinning-density + coherency-screening (coh05) family:
compares each DA arm's own species-13 (CYGNSS L1) O-F std, plus the pooled
monitor Tb/SM O-F std, against its own per-arm OL cross-mask baseline (same
obs population, via use_obs=True on the DA arm -- see
run_cygl1_paired_OL_xmask_coh05*.py), which is the correct "does DA beat OL"
test per feedback_da_evaluation_standard for species 13 (CYGNSS L1), and
equivalent to the plain OL baseline for species 1-12 (monitor obs
populations coincide for OL and every DA arm in this family).

All % changes reported as (DA-OL)/OL for every species including CygL1, per
feedback_omf_stdv_sign_convention -- no sign flip for CygL1.

Each period below reflects how far that arm's run was actually extended:
dense_coh05 was never extended past 6mo, intermediate_coh05 stopped at 12mo,
dense075_coh05 reached the full 24mo record (Jan 2020-Dec 2021).

Consolidates three prior near-duplicate scripts that differed only in which
arms/period they covered (compare_cygl1_coh05_omf.py [orig 6mo, all 3 arms,
pre-cross-mask methodology], compare_cygl1_coh05_omf_12mo.py [12mo,
intermediate+dense075], compare_cygl1_coh05_omf_24mo.py [24mo, dense075
only]) -- merged 2026-09-08, now driven by CLI args instead of separate
files, and all periods use the same OL-cross-mask methodology (the 6mo-era
"pre-coh05 baseline" comparison, made obsolete once the OL-xmask companion
stats existed for all 3 arms, is dropped here).

Inputs: spatial_stats_{DA,OL_paired_monitor_xmask}_<arm>_<suffix>.pkl, built
by run_cygl1_paired_density_coh05*.py / run_cygl1_paired_OL_xmask_coh05*.py
(postproc_ObsFcstAna toolkit).

Usage: compare_cygl1_coh05_omf.py [--period {6mo,12mo,24mo,all}]
"""
import argparse
import pickle
import numpy as np

STATS = "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/output/postproc_paired_density/stats_output/"

SPECIES_NAMES = [
    "SMOS_Tbh_A", "SMOS_Tbh_D", "SMOS_Tbv_A", "SMOS_Tbv_D",
    "SMAP_Tbh_A", "SMAP_Tbh_D", "SMAP_Tbv_A", "SMAP_Tbv_D",
    "ASCAT_A", "ASCAT_B", "ASCAT_C", "CYGNSS_SM_6hr",
    "CYGNSS_L1",
]
TB_IDX = list(range(0, 8))
SM_IDX = list(range(8, 12))
CYGL1_IDX = 12

PERIODS = {
    "6mo":  dict(suffix="202001_202006", label="Jan-Jun 2020 (6-month)",
                 arms=["dense_coh05", "dense075_coh05", "intermediate_coh05"]),
    "12mo": dict(suffix="202001_202012", label="Jan-Dec 2020 (12-month)",
                 arms=["intermediate_coh05", "dense075_coh05"]),
    "24mo": dict(suffix="202001_202112", label="Jan 2020-Dec 2021 (24-month)",
                 arms=["dense075_coh05"]),
}


def pooled(N, mean, stdv):
    valid = (N > 0) & ~np.isnan(mean)
    if not valid.any():
        return np.nan, np.nan, 0
    N, mean, stdv = N[valid], mean[valid], stdv[valid]
    Ntot = N.sum()
    pooled_mean = (N * mean).sum() / Ntot
    pooled_var = (N * (stdv ** 2 + mean ** 2)).sum() / Ntot - pooled_mean ** 2
    return pooled_mean, np.sqrt(max(pooled_var, 0)), int(Ntot)


def group_pool(N, OmF_mean, OmF_stdv, idxs):
    Ng = N[:, idxs]; Mg = OmF_mean[:, idxs]; Sg = OmF_stdv[:, idxs]
    valid = (Ng > 0) & ~np.isnan(Mg)
    if not valid.any():
        return np.nan, np.nan, 0
    Ngf, Mgf, Sgf = Ng[valid], Mg[valid], Sg[valid]
    Ntot = Ngf.sum()
    pm = (Ngf * Mgf).sum() / Ntot
    ps = np.sqrt(max((Ngf * (Sgf ** 2 + Mgf ** 2)).sum() / Ntot - pm ** 2, 0))
    return pm, ps, int(Ntot)


def load(fname):
    with open(STATS + fname, "rb") as f:
        d = pickle.load(f)
    N = np.array(d["N_data"])
    OmF_mean = np.array(d["OmF_mean"])
    OmF_stdv = np.array(d["OmF_stdv"])

    per_species = {}
    for i, name in enumerate(SPECIES_NAMES):
        m, s, n = pooled(N[:, i], OmF_mean[:, i], OmF_stdv[:, i])
        per_species[name] = {"OmF_mean": m, "OmF_stdv": s, "N": n}

    per_group = {}
    for gname, idxs in [("Tb (K)", TB_IDX), ("SM (m3/m3)", SM_IDX)]:
        m, s, n = group_pool(N, OmF_mean, OmF_stdv, idxs)
        per_group[gname] = {"OmF_stdv": s}

    return per_species, per_group


def run_period(period_key):
    cfg = PERIODS[period_key]
    suffix, arms = cfg["suffix"], cfg["arms"]

    print("=" * 100)
    print(f"{cfg['label']}, Arizona-box paired experiment:")
    print("DA vs its own per-arm OL cross-mask baseline (same obs population)")
    print("=" * 100)
    header = f"{'':22s}{'Tb OmF_stdv (K)':>18s}{'% vs OL':>10s}{'SM OmF_stdv':>14s}{'% vs OL':>10s}" \
             f"{'CygL1 OmF_stdv (dB)':>22s}{'% vs OL':>10s}{'N CygL1':>10s}"
    print(header)

    for arm in arms:
        da_sp, da_grp = load(f"spatial_stats_DA_paired_{arm}_{suffix}.pkl")
        ol_sp, ol_grp = load(f"spatial_stats_OL_paired_monitor_xmask_{arm}_{suffix}.pkl")

        tb_da = da_grp["Tb (K)"]["OmF_stdv"]; tb_ol = ol_grp["Tb (K)"]["OmF_stdv"]
        sm_da = da_grp["SM (m3/m3)"]["OmF_stdv"]; sm_ol = ol_grp["SM (m3/m3)"]["OmF_stdv"]
        cyg_da = da_sp["CYGNSS_L1"]["OmF_stdv"]; cyg_ol = ol_sp["CYGNSS_L1"]["OmF_stdv"]
        n_cyg = da_sp["CYGNSS_L1"]["N"]

        tb_pct = 100 * (tb_da - tb_ol) / tb_ol
        sm_pct = 100 * (sm_da - sm_ol) / sm_ol
        cyg_pct = 100 * (cyg_da - cyg_ol) / cyg_ol

        print(f"{arm:22s}{tb_da:18.4f}{tb_pct:+10.2f}{sm_da:14.5f}{sm_pct:+10.2f}"
              f"{cyg_da:22.4f}{cyg_pct:+10.2f}{n_cyg:10d}")
        print(f"{'  (OL baseline)':22s}{tb_ol:18.4f}{'--':>10s}{sm_ol:14.5f}{'--':>10s}"
              f"{cyg_ol:22.4f}{'--':>10s}{ol_sp['CYGNSS_L1']['N']:10d}")

    print()
    print("Per-species monitor Tb/SM detail (% vs each arm's own OL cross-mask baseline):")
    for arm in arms:
        da_sp, _ = load(f"spatial_stats_DA_paired_{arm}_{suffix}.pkl")
        ol_sp, _ = load(f"spatial_stats_OL_paired_monitor_xmask_{arm}_{suffix}.pkl")
        print(f"\n-- {arm} --")
        for name in SPECIES_NAMES[:12]:
            m = da_sp[name]["OmF_stdv"]; o = ol_sp[name]["OmF_stdv"]
            pct = 100 * (m - o) / o if o and not np.isnan(o) else np.nan
            print(f"  {name:14s} OmF_stdv={m:8.4f}  OL={o:8.4f}  {pct:+7.2f}%")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--period", choices=list(PERIODS) + ["all"], default="all")
    args = ap.parse_args()

    periods = list(PERIODS) if args.period == "all" else [args.period]
    for p in periods:
        run_period(p)


if __name__ == "__main__":
    main()
