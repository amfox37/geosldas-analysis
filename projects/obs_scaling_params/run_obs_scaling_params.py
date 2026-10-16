#!/usr/bin/env python3

# Compute observation scaling parameters (climatological mean and stdv of obs and
# model forecasts) from the ObsFcstAna files of a GEOSldas experiment.
#
#   --method latlon : soil moisture obs (e.g., ASCAT, CYGNSS) on a global lat/lon grid,
#                     nc4 output for GEOSldas scale_obs_sfmc_zscore
#   --method tb     : Tb obs (e.g., SMAP, SMOS) in tile space, one orbit per run,
#                     binary output for GEOSldas scale_obs_Tb_zscore
#
# Output goes to <exp-path>/<exp-run>/output/<domain>/stats/<out-dir>/.
#
# Examples:
#
#   python run_obs_scaling_params.py --method latlon --exp-path /path/to/exps --exp-run OL1 \
#       --domain CF0360x6C_GLOBAL --start 2019-02 --end 2020-08 \
#       --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM \
#       --prefix CF0360_OL1_HSAF_zscore_stats_
#
#   python run_obs_scaling_params.py --method tb --exp-path /path/to/exps --exp-run OL1 \
#       --domain SMAP_EASEv2_M36_GLOBAL --start 2015-04 --end 2021-03 \
#       --descr-contains SMAP_L1C --orbit D --prefix OL1_SMAP_zscore_stats_
#
# amfox - Oct 2026

import argparse
import csv

import numpy as np

from obs_scaling_utils      import read_obs_param, get_year_bounds
from get_obs_scaling_latlon import get_obs_scaling_latlon
from get_obs_scaling_tb     import get_obs_scaling_tb

def parse_args(argv=None):

    p = argparse.ArgumentParser(description='Compute observation scaling parameters from GEOSldas ObsFcstAna files.')

    p.add_argument('--method',      required=True, choices=['latlon', 'tb'])
    p.add_argument('--exp-path',    required=True, help='directory that contains the experiment')
    p.add_argument('--exp-run',     required=True, help='experiment ID')
    p.add_argument('--domain',      required=True, help='e.g. SMAP_EASEv2_M36_GLOBAL')
    p.add_argument('--start',       required=True, help='first month of the period, YYYY-MM')
    p.add_argument('--end',         required=True, help='last month of the period, YYYY-MM')
    p.add_argument('--prefix',      required=True, help='output file name prefix')
    p.add_argument('--out-dir',     default='z_score_clim', help='output subdirectory of output/<domain>/stats/')
    p.add_argument('--run-months',  default='',
                   help='comma-separated months to process (default: 1-12 plus the months needed to wrap the window)')
    p.add_argument('--window-days', type=int, default=75, help='length of moving window [days]')
    p.add_argument('--ndata-min',   type=int, default=20, help='minimum number of obs in window')
    p.add_argument('--dt-assim',    type=int, default=10800, help='ObsFcstAna file interval [s]')
    p.add_argument('--t0-assim',    type=int, default=0,     help='time of first ObsFcstAna file of the day [s]')
    p.add_argument('--obsparam',    default='',
                   help='ldas_obsparam file (default: rc_out/Y<start>/M<start>/<exp-run>.ldas_obsparam.<start>01_0000z.txt)')
    p.add_argument('--each-doy',    action='store_true', help='also write one file per DOY')

    g = p.add_argument_group('latlon method')
    g.add_argument('--species',              default='', help='comma-separated obs descr from ldas_obsparam')
    g.add_argument('--grid-resolution',      type=float, default=0.25, help='lat/lon grid spacing [deg]')
    g.add_argument('--no-combine-species',   dest='combine_species', action='store_false',
                   help='one set of statistics per species (default: all species combined)')
    g.add_argument('--each-pentad',          action='store_true', help='also write one file per pentad')
    g.add_argument('--no-all-pentads',       dest='all_pentads', action='store_false',
                   help='do not write the file with all pentads')
    g.add_argument('--exclude-tiles-csv',    default='',
                   help="CSV file with a 'tilenum' column and a column (see --exclude-tiles-column) that is "
                        "true for tiles whose obs are not used, e.g. to apply a QC that postdates the ObsFcstAna files")
    g.add_argument('--exclude-tiles-column', default='reject')

    g = p.add_argument_group('tb method')
    g.add_argument('--descr-contains',  default='SMAP_L1C', help='select Tb species whose descr contains this string')
    g.add_argument('--orbit',           choices=['A', 'D'], default='D')
    g.add_argument('--angles',          default='40', help='comma-separated incidence angles [deg]')
    g.add_argument('--no-convert-grid', dest='convert_grid', action='store_false',
                   help='use all model tiles (default: EASEv2 M36 administering tiles)')

    args = p.parse_args(argv)

    if args.method == 'latlon' and not args.species:
        p.error('--species is required for --method latlon')

    return args

def main(argv=None):

    args = parse_args(argv)

    start_year, start_month = [int(x) for x in args.start.split('-')]
    end_year,   end_month   = [int(x) for x in args.end.split('-')]

    if args.run_months:
        run_months = [int(x) for x in args.run_months.split(',')]
    else:
        run_months = list(range(1, 13)) + list(range(1, int(np.ceil(args.window_days / 30)) + 1))

    first_year, last_year = get_year_bounds(run_months, start_year, start_month, end_year, end_month)

    obsparam_file = args.obsparam or (
        f'{args.exp_path}/{args.exp_run}/output/{args.domain}/rc_out/Y{start_year:04d}/M{start_month:02d}/'
        f'{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt')

    obs_param = read_obs_param(obsparam_file)

    print('Observation scaling configuration:')
    for key, value in vars(args).items():
        print(f'  {key:22s}: {value}')

    if args.method == 'latlon':

        species_names = args.species.split(',')
        species       = []
        for name in species_names:
            ids = [int(p['species']) for p in obs_param if p['descr'] == name]
            if not ids:
                raise ValueError(f'species {name} not found in {obsparam_file}')
            species += ids
        species = sorted(set(species))
        print(f'  species ids             : {species}')

        exclude_tile_mask = None
        if args.exclude_tiles_csv:
            exclude_tile_mask = read_exclude_tile_mask(args.exclude_tiles_csv, args.exclude_tiles_column)

        get_obs_scaling_latlon(args.exp_path, args.exp_run, args.domain, species, species_names,
                               args.combine_species, run_months, first_year, last_year,
                               args.dt_assim, args.t0_assim, args.grid_resolution,
                               args.window_days, args.ndata_min, args.prefix, args.out_dir,
                               args.each_doy, args.each_pentad, args.all_pentads,
                               exclude_tile_mask=exclude_tile_mask)

    else:

        angles = [float(x) for x in args.angles.split(',')]
        orbit  = 1 if args.orbit == 'A' else 2

        get_obs_scaling_tb(args.exp_path, args.exp_run, args.domain, obs_param, args.descr_contains,
                           orbit, angles, run_months, first_year, last_year,
                           args.dt_assim, args.t0_assim, args.window_days, args.ndata_min,
                           args.prefix, args.out_dir, convert_grid=args.convert_grid,
                           print_each_doy=args.each_doy)

# ----------------------------------------------------------------------------
#
# read CSV file with tilenum and exclude flag; returns boolean array indexed by
# tilenum (1-based; index 0 unused), True = use the obs, False = drop the obs

def read_exclude_tile_mask(fname, column):

    tilenum = []
    reject  = []
    with open(fname, newline='') as fid:
        reader = csv.DictReader(fid)
        if column not in (reader.fieldnames or []):
            raise ValueError(f'column {column} not found in {fname} (columns: {reader.fieldnames})')
        for row in reader:
            tilenum.append(int(row['tilenum']))
            reject.append(str(row[column]).strip().lower() in ('true', '1', 'yes'))

    tilenum = np.array(tilenum)
    reject  = np.array(reject, dtype=bool)

    keep = np.ones(tilenum.max() + 1, dtype=bool)
    keep[tilenum[reject]] = False

    print(f'  excluding obs at {reject.sum()} tiles listed in {fname}')

    return keep

if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
