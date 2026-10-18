
# Utilities for the GEOSldas observation scaling scripts:
#  - readers for ldas_obsparam, ldas_tilecoord, ldas_tilegrids and nc4 ObsFcstAna files
#  - helpers for the moving-window loop over the climatological year
#
# The readers return dicts with the same keys as GEOSldas_App/util/shared/python/read_GEOSldas.py.
#
# amfox - Oct 2026

import os
import struct

import numpy as np

from datetime import datetime
from netCDF4  import Dataset

nodata = -9999.

# ----------------------------------------------------------------------------
#
# reader for GEOSldas obsparam file (ascii)
#
# One field per line; string fields are quoted with ' or ".

def read_obs_param(fname):

    print(f"Reading {fname}")

    def get_str(fid):
        return fid.readline().strip().strip("'\"")

    def get_float(fid):
        return float(fid.readline().strip())

    with open(fname, 'r') as fid:
        N_obs_param = int(fid.readline().strip())

        obs_param = []
        for _ in range(N_obs_param):
            param = {}
            param['descr']          = get_str(  fid)
            param['species']        = get_float(fid)
            param['orbit']          = get_float(fid)
            param['pol']            = get_float(fid)
            param['N_ang']          = int(get_float(fid))
            param['ang']            = np.array([float(x) for x in fid.readline().split()])
            param['freq']           = get_float(fid)
            param['FOV']            = get_float(fid)
            param['FOV_units']      = get_str(  fid)
            param['assim']          = get_str(  fid)
            param['scale']          = get_str(  fid)
            param['getinnov']       = get_str(  fid)
            param['RTM_ID']         = get_float(fid)
            param['bias_Npar']      = get_float(fid)
            param['bias_trel']      = get_float(fid)
            param['bias_tcut']      = get_float(fid)
            param['nodata']         = get_float(fid)
            param['varname']        = get_str(  fid)
            param['units']          = get_str(  fid)
            param['fcstvarname']    = get_str(  fid)
            param['fcstunits']      = get_str(  fid)
            param['path']           = get_str(  fid)
            param['name']           = get_str(  fid)
            param['maskpath']       = get_str(  fid)
            param['maskname']       = get_str(  fid)
            param['scalepath']      = get_str(  fid)
            param['scalename']      = get_str(  fid)
            param['flistpath']      = get_str(  fid)
            param['flistname']      = get_str(  fid)
            param['errstd']         = get_float(fid)
            param['std_normal_max'] = get_float(fid)
            param['zeromean']       = get_str(  fid)
            param['coarsen_pert']   = get_str(  fid)
            param['xcorr']          = get_float(fid)
            param['ycorr']          = get_float(fid)
            param['adapt']          = get_float(fid)

            obs_param.append(param)

    return obs_param

# ----------------------------------------------------------------------------
#
# reader for GEOSldas tilecoord file (little-endian Fortran sequential binary)

def read_tilecoord(fname):

    print(f"Reading {fname}")

    int_fields = ['tile_id', 'typ', 'pfaf', 'i_indg', 'j_indg']
    fields     = ['tile_id', 'typ', 'pfaf', 'com_lon', 'com_lat', 'min_lon', 'max_lon',
                  'min_lat', 'max_lat', 'i_indg', 'j_indg', 'frac_cell', 'frac_pfaf',
                  'area', 'elev']

    tile_coord = {}

    with open(fname, 'rb') as ifp:
        ifp.read(4)
        tile_coord['N_tile'] = struct.unpack('<i', ifp.read(4))[0]
        ifp.read(4)

        Nt = tile_coord['N_tile']

        for field in fields:
            dtype = '<i4' if field in int_fields else '<f4'
            ifp.read(4)
            tile_coord[field] = np.frombuffer(ifp.read(Nt*4), dtype=dtype)
            ifp.read(4)

    return tile_coord

# ----------------------------------------------------------------------------
#
# reader for GEOSldas tilegrids file (little-endian Fortran sequential binary);
# returns the global and the domain tile grid

def read_tilegrids(fname):

    print(f"Reading {fname}")

    tile_grids = []

    with open(fname, 'rb') as ifp:
        for _ in range(2):
            ifp.read(4)
            tile_grid = {}
            tile_grid['gridtype'] = ifp.read(40).decode('ascii').strip('\x00 ')
            ints   = struct.unpack('<7i', ifp.read(28))
            floats = struct.unpack('<6f', ifp.read(24))
            ifp.read(4)
            for key, value in zip(['ind_base', 'i_dir', 'j_dir', 'N_lon', 'N_lat', 'i_offg', 'j_offg'], ints):
                tile_grid[key] = value
            for key, value in zip(['ll_lon', 'll_lat', 'ur_lon', 'ur_lat', 'dlon', 'dlat'], floats):
                tile_grid[key] = value
            tile_grids.append(tile_grid)

    return tile_grids[0], tile_grids[1]

# ----------------------------------------------------------------------------
#
# reader for GEOSldas ObsFcstAna file (nc4)
#
# Returns None if the file does not exist.  No-data values (-9999 or _FillValue)
# are set to NaN for obs, obsvar, fcst, fcstvar, ana and anavar.

def read_ObsFcstAna_nc4(fname):

    if not os.path.exists(fname):
        return None

    with Dataset(fname, 'r') as ncid:

        def get_var(name, dtype):
            return np.asarray(np.ma.filled(ncid.variables[name][:], np.nan)).astype(dtype)

        def get_data(name):
            data = get_var(name, np.float32)
            data[data == nodata] = np.nan
            return data

        ofa = {
            'obs_assim'  : get_var( 'assim_flag', bool    ),
            'obs_species': get_var( 'species',    np.int32),
            'obs_tilenum': get_var( 'tilenum',    np.int32),
            'obs_lon'    : get_var( 'lon',        np.float32),
            'obs_lat'    : get_var( 'lat',        np.float32),
            'obs_obs'    : get_data('obs'    ),
            'obs_obsvar' : get_data('obsvar' ),
            'obs_fcst'   : get_data('fcst'   ),
            'obs_fcstvar': get_data('fcstvar'),
            'obs_ana'    : get_data('ana'    ),
            'obs_anavar' : get_data('anavar' ),
        }

    return ofa

# ----------------------------------------------------------------------------
#
# helpers for the loop over the climatological year

# Days per month.  The climatological year is always the non-leap year 2014.

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def get_ObsFcstAna_fname(inpath, exp_run, year, month, day, hour, minute):

    return (inpath + '/ana/ens_avg/Y' + f'{year:04d}' + '/M' + f'{month:02d}' + '/' + exp_run +
            '.ens_avg.ldas_ObsFcstAna.' + f'{year:04d}{month:02d}{day:02d}_{hour:02d}{minute:02d}z.nc4')

def count_ObsFcstAna_files(inpath, exp_run, run_months, first_year, last_year, cycles):

    # number of ObsFcstAna files that exist for run_months, years and cycles (seconds of day)

    n = 0
    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month[month-1] + 1):
            for seconds in cycles:
                for year in range(first_year[imonth], last_year[imonth] + 1):
                    fname = get_ObsFcstAna_fname(inpath, exp_run, year, month, day,
                                                 seconds // 3600, (seconds % 3600) // 60)
                    n += os.path.exists(fname)
    return n

def get_year_bounds(run_months, start_year, start_month, end_year, end_month):

    # first and last year to read for each entry of run_months, so that only
    # months between (start_year, start_month) and (end_year, end_month) are used

    first_year = []
    last_year  = []
    for month in run_months:
        first_year.append(start_year + 1 if month < start_month else start_year)
        last_year.append( end_year   - 1 if month > end_month   else end_year  )

    return first_year, last_year

def get_doy_pentad_range(run_months, first_year, last_year):

    # first/last DOY and pentad covered by run_months, for the output file names

    run_months = np.asarray(run_months)
    first_year = np.asarray(first_year)
    last_year  = np.asarray(last_year)

    first_month = int(run_months[first_year == first_year.min()].min())
    last_month  = int(run_months[last_year  == last_year.max() ].max())

    cumulative_days = np.concatenate(([0], np.cumsum(days_in_month)))

    doy_start    = int(cumulative_days[first_month-1] + 1)
    doy_end      = int(cumulative_days[last_month])
    pentad_start = int(np.ceil( doy_start / 5))
    pentad_end   = int(np.floor(doy_end   / 5))

    return doy_start, doy_end, pentad_start, pentad_end

def days_since_1950(times):

    return np.array([(t - datetime(1950, 1, 1)).total_seconds() / 86400. for t in times])

# ====================== EOF =========================================================
