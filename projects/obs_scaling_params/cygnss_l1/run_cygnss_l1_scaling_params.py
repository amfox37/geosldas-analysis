#!/usr/bin/env python3

# Compute CYGNSS L1 scaling parameters (climatological mean and stdv of obs and model
# forecasts) on the model tiles of an EASEv2 M36 GEOSldas experiment, from the
# ObsFcstAna files of an unscaled (monitoring) run.  See get_obs_scaling_owner_tile.py.
#
# Output goes to <exp-path>/<exp-run>/output/<domain>/stats/<out-dir>/.
#
# Example:
#
#   python cygnss_l1/run_cygnss_l1_scaling_params.py --exp-path /path/to/exps \
#       --exp-run OLv8_M36_AZ_fixedop --start 2020-01 --end 2022-12
#
# This method is not part of the GEOSldas_GridComp obs_scaling_params scripts; it
# reuses their readers in ../obs_scaling_utils.py.
#
# amfox - Oct 2026

import argparse
import os
import sys

import numpy as np

from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from obs_scaling_utils          import read_obs_param, get_year_bounds
from get_obs_scaling_owner_tile import get_obs_scaling_owner_tile, cygnss_l1_descr

def parse_args(argv=None):

    p = argparse.ArgumentParser(description='Compute CYGNSS L1 scaling parameters on the model tiles '
                                            '(EASEv2 M36) from GEOSldas ObsFcstAna files.')

    p.add_argument('--exp-path',    required=True, help='directory that contains the experiment')
    p.add_argument('--exp-run',     required=True, help='experiment ID')
    p.add_argument('--domain',      default='SMAP_EASEv2_M36_GLOBAL')
    p.add_argument('--start',       required=True, help='first month of the period, YYYY-MM')
    p.add_argument('--end',         required=True, help='last month of the period, YYYY-MM')
    p.add_argument('--prefix',      default='CYGNSS_L1_owner_zscore_', help='output file name prefix')
    p.add_argument('--out-dir',     default='cygnss_l1_z_score_clim', help='output subdirectory of output/<domain>/stats/')
    p.add_argument('--descr',       default=cygnss_l1_descr, help='CYGNSS L1 obs descr in ldas_obsparam')
    p.add_argument('--run-months',  default='',
                   help='comma-separated months to process (default: 1-12 plus the months needed to wrap the window)')
    p.add_argument('--window-days', type=int,   default=75,    help='length of moving window [days]')
    p.add_argument('--ndata-min',   type=int,   default=20,    help='minimum number of obs in window')
    p.add_argument('--std-epsilon', type=float, default=1e-6,  help='minimum obs and model stdv')
    p.add_argument('--dt-assim',    type=int,   default=10800, help='ObsFcstAna file interval [s]')
    p.add_argument('--t0-assim',    type=int,   default=0,     help='time of first ObsFcstAna file of the day [s]')
    p.add_argument('--obsparam',    default='',
                   help='ldas_obsparam file (default: rc_out/Y<start>/M<start>/<exp-run>.ldas_obsparam.<start>01_0000z.txt)')
    p.add_argument('--tilecoord',   default='', help='ldas_tilecoord.bin file (default: from rc_out/)')
    p.add_argument('--tilegrids',   default='', help='ldas_tilegrids.bin file (default: from rc_out/)')
    p.add_argument('--no-overwrite', dest='overwrite', action='store_false')

    g = p.add_argument_group('obs keep-list')
    g.add_argument('--keep-list', default='',
                   help='parquet file of obs (T = window center [s since 2020-01-01], float32 lon/lat, '
                        'coherency_ratio) from CYGNSS_L1_AZ/scripts/cygl1_coherency_filter.py annotate; '
                        'only obs whose (T, lon, lat) has keep-lo <= coherency_ratio <= keep-hi are used')
    g.add_argument('--keep-lo', type=float, default=0.40)
    g.add_argument('--keep-hi', type=float, default=2.16)

    args = p.parse_args(argv)

    for name in ('start', 'end'):
        try:
            setattr(args, name, datetime.strptime(getattr(args, name), '%Y-%m'))
        except ValueError:
            p.error(f'--{name} must be YYYY-MM, not {getattr(args, name)}')
    if args.start > args.end:
        p.error('--start must not be after --end')

    try:
        args.run_months = [int(x) for x in args.run_months.split(',') if x.strip()]
    except ValueError:
        p.error(f'--run-months must be comma-separated integers, not {args.run_months}')
    if any(m < 1 or m > 12 for m in args.run_months):
        p.error(f'--run-months must be between 1 and 12: {args.run_months}')

    return args

def main(argv=None):

    args = parse_args(argv)

    start_year, start_month = args.start.year, args.start.month
    end_year,   end_month   = args.end.year,   args.end.month

    if args.run_months:
        run_months = args.run_months
    else:
        run_months = list(range(1, 13)) + list(range(1, int(np.ceil(args.window_days / 30)) + 1))

    first_year, last_year = get_year_bounds(run_months, start_year, start_month, end_year, end_month)

    obsparam_file = args.obsparam or (
        f'{args.exp_path}/{args.exp_run}/output/{args.domain}/rc_out/Y{start_year:04d}/M{start_month:02d}/'
        f'{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt')

    obs_param = read_obs_param(obsparam_file)

    print('CYGNSS L1 scaling configuration:')
    for key, value in vars(args).items():
        if isinstance(value, datetime):
            value = value.strftime('%Y-%m')
        print(f'  {key:12s}: {value}')

    obs_keep = make_keep_filter(args.keep_list, args.keep_lo, args.keep_hi) if args.keep_list else None

    fname = get_obs_scaling_owner_tile(args.exp_path, args.exp_run, args.domain, obs_param, args.descr,
                                       run_months, first_year, last_year, args.dt_assim, args.t0_assim,
                                       args.window_days, args.ndata_min, args.prefix, args.out_dir,
                                       std_epsilon=args.std_epsilon, tilecoord_file=args.tilecoord,
                                       tilegrids_file=args.tilegrids, overwrite=args.overwrite,
                                       obs_keep=obs_keep)

    print(f'Wrote {fname}')

# ----------------------------------------------------------------------------
#
# obs_keep(time, lon, lat): True where the obs (same ObsFcstAna time, exact float32
# lon/lat) is in the keep-list with keep_lo <= coherency_ratio <= keep_hi

def make_keep_filter(fname, keep_lo, keep_hi):

    import pandas as pd

    def lonlat_key(lon, lat):
        return ((np.asarray(lon, np.float32).view(np.uint32).astype(np.int64) << 32) |
                 np.asarray(lat, np.float32).view(np.uint32).astype(np.int64))

    d = pd.read_parquet(fname, columns=['T', 'lon', 'lat', 'coherency_ratio'])
    d = d[(d.coherency_ratio >= keep_lo) & (d.coherency_ratio <= keep_hi)]

    key  = lonlat_key(d.lon.to_numpy(), d.lat.to_numpy())
    keep = {int(t): set(k) for t, k in pd.Series(key).groupby(d['T'].to_numpy()).agg(set).items()}

    print(f'  keep-list: {len(d)} obs with {keep_lo} <= coherency_ratio <= {keep_hi} from {fname}')

    epoch = datetime(2020, 1, 1)

    def obs_keep(time, lon, lat):
        s = keep.get(int((time - epoch).total_seconds()), set())
        return np.fromiter((int(k) in s for k in lonlat_key(lon, lat)), dtype=bool, count=len(lon))

    return obs_keep

if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
