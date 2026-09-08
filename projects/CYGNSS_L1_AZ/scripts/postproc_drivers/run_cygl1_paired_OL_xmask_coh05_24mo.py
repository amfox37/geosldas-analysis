#!/usr/bin/env python3

# 24-month extension of run_cygl1_paired_OL_xmask_coh05_12mo.py, matching
# run_cygl1_paired_density_coh05_24mo.py's period (Jan 2020-Dec 2021). Only
# dense075_coh05 was extended past 12mo (intermediate_coh05 stays at Jan-Dec
# 2020), so it is the only arm processed here.
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
import glob
import pickle

from datetime               import timedelta
from dateutil.relativedelta import relativedelta

from run_cygl1_paired_density import load_exp, EXPDIR, DOMAIN, SUM_ROOT, OUT_PATH, start_time
from postproc_ObsFcstAna    import postproc_ObsFcstAna

end_time = start_time + relativedelta(months=24)   # 2022-01-01


def last_complete_month_local(expid):
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


ARMS = {
    'DAv8_M36_AZ_paired_cygl1_dense075_coh05': 'OL_paired_monitor_xmask_dense075_coh05',
}


def process(da_expid, exptag):
    ol_main = load_exp('OLv8_M36_AZ_paired_monitor', exptag=exptag)
    da_sup  = load_exp(da_expid)
    da_sup['use_obs'] = True

    exp_end_time = min(last_complete_month_local('OLv8_M36_AZ_paired_monitor'),
                        last_complete_month_local(da_expid))
    if exp_end_time is None or exp_end_time <= start_time:
        print(f'\n=== {exptag}: SKIPPED, no complete months yet ===')
        return

    sum_path = SUM_ROOT + exptag + '/'
    os.makedirs(sum_path, exist_ok=True)

    print(f'\n=== {exptag} -> {sum_path} (through '
          f'{(exp_end_time - relativedelta(months=1)).strftime("%Y%m")}) ===')

    postproc = postproc_ObsFcstAna([ol_main, da_sup], start_time, exp_end_time, sum_path=sum_path)
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
    for da_expid, exptag in ARMS.items():
        process(da_expid, exptag)


if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
