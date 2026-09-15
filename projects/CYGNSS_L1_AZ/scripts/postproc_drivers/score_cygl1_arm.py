#!/usr/bin/env python3
"""
Generalized single-arm scorer for the CYGNSS L1 paired thinning-density/
coherency-screening family: computes own-arm O-F stats AND the OL
cross-mask baseline (per feedback_da_evaluation_standard) for ONE DA arm,
through whatever months are actually complete so far -- then prints the
DA-vs-OL comparison table.

Replaces writing a new hardcoded run_cygl1_paired_density_<arm>.py /
run_cygl1_paired_OL_xmask_<arm>.py / compare_cygl1_<arm>_omf.py trio per
arm (see feedback_avoid_per_run_script_copies) -- one script, arm passed
via CLI args, safe to rerun as a run progresses (save_monthly_sums() is
idempotent, only computes new months).

Usage:
  score_cygl1_arm.py --da-expid DAv8_M36_AZ_paired_cygl1_dense050_coh05 \\
                      --arm-tag dense050_coh05
"""
import sys
import argparse
import warnings
warnings.filterwarnings("ignore")
import os
import glob
import pickle
import numpy as np

from datetime               import datetime, timedelta
from dateutil.relativedelta import relativedelta

TOOLKIT_DIR = ('/gpfsm/dnb06/projects/p284/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/'
               'GEOSldas_App/util/postproc/ObsFcstAna_stats')
SHARED_PYTHON_DIR = ('/gpfsm/dnb06/projects/p284/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/'
                      'GEOSldas_App/util/shared/python')
sys.path.append(SHARED_PYTHON_DIR)
sys.path.append(TOOLKIT_DIR)

from run_cygl1_paired_density import load_exp, EXPDIR, DOMAIN, SUM_ROOT, OUT_PATH
from postproc_ObsFcstAna    import postproc_ObsFcstAna

OL_EXPID = 'OLv8_M36_AZ_paired_monitor'

SPECIES_NAMES = [
    "SMOS_Tbh_A", "SMOS_Tbh_D", "SMOS_Tbv_A", "SMOS_Tbv_D",
    "SMAP_Tbh_A", "SMAP_Tbh_D", "SMAP_Tbv_A", "SMAP_Tbv_D",
    "ASCAT_A", "ASCAT_B", "ASCAT_C", "CYGNSS_SM_6hr",
    "CYGNSS_L1",
]
TB_IDX = list(range(0, 8))
SM_IDX = list(range(8, 12))


def last_complete_month(expid, start_time, end_time):
    base = EXPDIR + expid + '/output/' + DOMAIN + '/ana/ens_avg/'
    month = start_time
    last_good = None
    while True:
        mo_dir = base + month.strftime('Y%Y/M%m') + '/'
        n_expected = (month + relativedelta(months=1) - month).days * 8
        n_found = len(glob.glob(mo_dir + '*ldas_ObsFcstAna*'))
        if not os.path.isdir(mo_dir) or n_found < n_expected - 8:
            break
        last_good = month
        month = month + relativedelta(months=1)
        if month >= end_time:
            break
    if last_good is None:
        return None
    return last_good + relativedelta(months=1)


def run_postproc(exp_list, exptag, start_time, exp_end_time):
    sum_path = SUM_ROOT + exptag + '/'
    os.makedirs(sum_path, exist_ok=True)
    postproc = postproc_ObsFcstAna(exp_list, start_time, exp_end_time, sum_path=sum_path)
    postproc.save_monthly_sums()

    suffix = start_time.strftime('%Y%m') + '_' + (exp_end_time + relativedelta(months=-1)).strftime('%Y%m')
    pkl_file = OUT_PATH + 'spatial_stats_' + exptag + '_' + suffix + '.pkl'
    stats_spatial = postproc.calc_spatial_stats_from_sums()
    with open(pkl_file, 'wb') as f:
        pickle.dump(stats_spatial, f)
    print(f'wrote {pkl_file}')

    nc4_suffix = start_time.strftime('%Y%m%d') + '_' + (exp_end_time + timedelta(days=-1)).strftime('%Y%m%d')
    nc4_file = OUT_PATH + 'temporal_stats_' + exptag + '_' + nc4_suffix + '.nc4'
    postproc.calc_temporal_stats_from_sums(write_to_nc=True, fout_stats=nc4_file)
    print(f'wrote {nc4_file}')
    return pkl_file


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


def load_stats(pkl_file):
    with open(pkl_file, 'rb') as f:
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


def print_comparison(arm_tag, da_pkl, ol_pkl):
    da_sp, da_grp = load_stats(da_pkl)
    ol_sp, ol_grp = load_stats(ol_pkl)

    print()
    print("=" * 100)
    print(f"{arm_tag}: DA vs its own OL cross-mask baseline (same obs population)")
    print("=" * 100)
    header = f"{'':22s}{'Tb OmF_stdv (K)':>18s}{'% vs OL':>10s}{'SM OmF_stdv':>14s}{'% vs OL':>10s}" \
             f"{'CygL1 OmF_stdv (dB)':>22s}{'% vs OL':>10s}{'N CygL1':>10s}"
    print(header)

    tb_da = da_grp["Tb (K)"]["OmF_stdv"]; tb_ol = ol_grp["Tb (K)"]["OmF_stdv"]
    sm_da = da_grp["SM (m3/m3)"]["OmF_stdv"]; sm_ol = ol_grp["SM (m3/m3)"]["OmF_stdv"]
    cyg_da = da_sp["CYGNSS_L1"]["OmF_stdv"]; cyg_ol = ol_sp["CYGNSS_L1"]["OmF_stdv"]
    n_cyg = da_sp["CYGNSS_L1"]["N"]

    tb_pct = 100 * (tb_da - tb_ol) / tb_ol
    sm_pct = 100 * (sm_da - sm_ol) / sm_ol
    cyg_pct = 100 * (cyg_da - cyg_ol) / cyg_ol

    print(f"{arm_tag:22s}{tb_da:18.4f}{tb_pct:+10.2f}{sm_da:14.5f}{sm_pct:+10.2f}"
          f"{cyg_da:22.4f}{cyg_pct:+10.2f}{n_cyg:10d}")
    print(f"{'  (OL baseline)':22s}{tb_ol:18.4f}{'--':>10s}{sm_ol:14.5f}{'--':>10s}"
          f"{cyg_ol:22.4f}{'--':>10s}{ol_sp['CYGNSS_L1']['N']:10d}")

    print()
    print("Per-species monitor Tb/SM detail (% vs OL cross-mask baseline):")
    for name in SPECIES_NAMES[:12]:
        m = da_sp[name]["OmF_stdv"]; o = ol_sp[name]["OmF_stdv"]
        pct = 100 * (m - o) / o if o and not np.isnan(o) else np.nan
        print(f"  {name:14s} OmF_stdv={m:8.4f}  OL={o:8.4f}  {pct:+7.2f}%")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--da-expid', required=True)
    ap.add_argument('--arm-tag', required=True, help='e.g. dense050_coh05')
    ap.add_argument('--ol-expid', default=OL_EXPID)
    ap.add_argument('--start', default='20200101', help='yyyymmdd, default 20200101')
    ap.add_argument('--cap-end', default=None, help='yyyymmdd, optional cap on how far to score')
    args = ap.parse_args()

    start_time = datetime.strptime(args.start, '%Y%m%d')
    cap_end = datetime.strptime(args.cap_end, '%Y%m%d') if args.cap_end else datetime(2099, 1, 1)

    own_end = last_complete_month(args.da_expid, start_time, cap_end)
    if own_end is None or own_end <= start_time:
        print(f'{args.arm_tag}: no complete months yet, nothing to score')
        return
    own_exp = load_exp(args.da_expid, exptag='DA_paired_' + args.arm_tag)
    da_pkl = run_postproc([own_exp], 'DA_paired_' + args.arm_tag, start_time, own_end)

    xmask_end = min(
        last_complete_month(args.ol_expid, start_time, cap_end),
        own_end,
    )
    xmask_tag = 'OL_paired_monitor_xmask_' + args.arm_tag
    ol_main = load_exp(args.ol_expid, exptag=xmask_tag)
    da_sup = load_exp(args.da_expid)
    da_sup['use_obs'] = True
    ol_pkl = run_postproc([ol_main, da_sup], xmask_tag, start_time, xmask_end)

    print_comparison(args.arm_tag, da_pkl, ol_pkl)


if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
