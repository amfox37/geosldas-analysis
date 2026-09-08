#!/usr/bin/env python3

# 24-month extension of run_cygl1_paired_density_coh05_12mo.py: dense075_coh05
# was restarted from BEG_DATE 2026-09-07 (job 58274545, superseding
# 58274014/58274452) with a corrected HISTORY.rc collection set, and run all
# the way through END_DATE 20220101 (cap_restart confirmed 20220101 000000,
# full resubmission chain COMPLETED, zero LDAS ERROR/forrtl in
# GEOSldas_log_txt). intermediate_coh05 was NOT extended past its existing
# 12-month (Jan-Dec 2020) record, so it is intentionally omitted here --
# only dense075_coh05 is processed.
#
# save_monthly_sums() skips months whose sums nc4 already exists, so this
# only computes the new Jan-Dec 2021 months and appends into the same
# sum_path directory already holding Jan-Dec 2020.
#
# Relocated 2026-09-08 from the hsaf_cdr_test obsfcstana-nc4-postproc
# worktree (its .git pointer is stale/broken -- points at a moved
# GEOSldas_develop path -- and it's the wrong repo for our own driver
# scripts anyway) into geosldas-analysis, this project's actual analysis
# repo. The toolkit itself (postproc_ObsFcstAna.py, shared/python/) stays
# put in that same worktree, referenced via absolute sys.path -- tried
# repointing at the "official" dnb34 GEOSldas_cygnss_operator checkout
# instead, but that copy is an older/divergent version (no NC4
# ObsFcstAna support, different obs_param format) and crashes; see the
# NOTE below. Only sys.path changed here, no driver logic changed. Still
# depends on run_cygl1_paired_density.py, which has NOT been moved yet
# (still lives in the old worktree location) -- see PARENT_DIR below.

import sys

# NOTE (2026-09-08): the dnb34 GEOSldas_cygnss_operator checkout's own copy of
# this toolkit is NOT usable here -- it's an older/divergent version missing
# NC4 ObsFcstAna support (read_ObsFcstAna_nc4) and using a different
# obs_param text format, confirmed via a real crash when tried. The toolkit
# copy that actually works for this project is the one physically sitting in
# the (git-unregistered) hsaf_cdr_test worktree below -- it's ahead of the
# "official" checkout, apparently never merged back upstream. Pointing here
# is a read-only import dependency, not something we're committing there.
TOOLKIT_DIR = ('/gpfsm/dnb06/projects/p284/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/'
               'GEOSldas_App/util/postproc/ObsFcstAna_stats')
SHARED_PYTHON_DIR = ('/gpfsm/dnb06/projects/p284/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/'
                      'GEOSldas_App/util/shared/python')
# run_cygl1_paired_density.py has not been relocated yet (still in the old
# worktree) -- keep this on the path until it is.
PARENT_DIR = TOOLKIT_DIR

sys.path.append(SHARED_PYTHON_DIR)
sys.path.append(TOOLKIT_DIR)
sys.path.append(PARENT_DIR)

import warnings;  warnings.filterwarnings("ignore")
import os
import pickle

from datetime               import timedelta
from dateutil.relativedelta import relativedelta

from run_cygl1_paired_density import load_exp, EXPDIR, DOMAIN, SUM_ROOT, OUT_PATH, start_time
from postproc_ObsFcstAna    import postproc_ObsFcstAna

end_time = start_time + relativedelta(months=24)   # 2020-01-01 -> 2022-01-01


def last_complete_month_local(expid):
    import glob
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


def process(expid, exptag):
    exp_end_time = last_complete_month_local(expid)
    if exp_end_time is None or exp_end_time <= start_time:
        print(f'\n=== {exptag}: SKIPPED, no complete months yet ===')
        return

    exp = load_exp(expid, exptag=exptag)
    sum_path = SUM_ROOT + exptag + '/'
    os.makedirs(sum_path, exist_ok=True)

    print(f'\n=== {exptag} -> {sum_path} (through '
          f'{(exp_end_time - relativedelta(months=1)).strftime("%Y%m")}) ===')

    postproc = postproc_ObsFcstAna([exp], start_time, exp_end_time, sum_path=sum_path)
    postproc.save_monthly_sums()

    pkl_file = OUT_PATH + 'spatial_stats_' + exptag + '_' + start_time.strftime('%Y%m') + \
        '_' + (exp_end_time + timedelta(days=-1)).strftime('%Y%m') + '.pkl'
    stats_spatial = postproc.calc_spatial_stats_from_sums()
    with open(pkl_file, 'wb') as f:
        pickle.dump(stats_spatial, f)
    print(f'wrote {pkl_file}')

    nc4_file = OUT_PATH + 'temporal_stats_' + exptag + '_' + start_time.strftime('%Y%m%d') + \
        '_' + (exp_end_time + timedelta(days=-1)).strftime('%Y%m%d') + '.nc4'
    postproc.calc_temporal_stats_from_sums(write_to_nc=True, fout_stats=nc4_file)
    print(f'wrote {nc4_file}')


def main():
    process('DAv8_M36_AZ_paired_cygl1_dense075_coh05', 'DA_paired_dense075_coh05')


if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
