
# Compute climatological mean and stdv of CYGNSS L1 observations and model forecasts
# from GEOSldas ObsFcstAna files on the model tiles ("owner tiles"), for use in the
# z-score scaling of CYGNSS L1 observations before assimilation.
#
# For each day of the climatological year, the statistics are computed from all
# obs/fcst pairs within a moving window of w_days days (across all years), and
# stored at the center of the window.  One nc4 file with all 73 pentads is written,
# with the tile statistics placed on the dense global EASEv2 M36 grid (x, y = 0-based
# global column and row of each model tile).
#
# Requires an EASEv2 M36 tile space (one tile per grid cell), an unscaled
# (monitoring) run, and at most one CYGNSS L1 obs per tile and ObsFcstAna file.
#
# This method is not part of the GEOSldas_GridComp obs_scaling_params scripts; it
# reuses their readers in ../obs_scaling_utils.py.
#
# amfox - Oct 2026

import os

import numpy as np

from datetime          import datetime, timedelta
from netCDF4           import Dataset

from obs_scaling_utils import read_ObsFcstAna_nc4, read_tilecoord, read_tilegrids
from obs_scaling_utils import get_ObsFcstAna_fname, get_doy_pentad_range, days_in_month
from obs_scaling_utils import days_since_1950, nodata, count_ObsFcstAna_files

N_pentads = 73

cygnss_l1_descr   = 'CYGNSS_L1_DDM3X5_CROP_SCALAR'
cygnss_l1_varname = 'cygl1scal'
cygnss_l1_units   = 'dB'

# sums in the moving window, and fields in the output data array

field_names = ('o_mean', 'o_std', 'm_mean', 'm_std', 'n_data', 'm_min', 'm_max')

def get_obs_scaling_owner_tile(exp_path, exp_run, domain, obs_param, descr,
                               run_months, first_year, last_year, dt_assim, t0_assim,
                               w_days, ndata_min, prefix, out_dir, std_epsilon=1e-6,
                               tilecoord_file='', tilegrids_file='', overwrite=True, obs_keep=None):

    # obs_keep, if given, is called as obs_keep(time, lon, lat) with the CYGNSS L1 obs of
    # each ObsFcstAna file and returns a boolean mask of the obs to use (e.g., to apply
    # an obs QC that the assimilation run applies but the OL run did not).
    #
    # Returns the name of the output file.

    if w_days % 5 or w_days % 10 == 0:
        raise ValueError('w_days must be an odd number of pentads (for example, 75)')
    if ndata_min < 1:
        raise ValueError('ndata_min must be at least 1')
    if std_epsilon < 0:
        raise ValueError('std_epsilon cannot be negative')
    if dt_assim <= 0:
        raise ValueError('dt_assim must be positive')
    if len(run_months) != len(first_year) or len(run_months) != len(last_year):
        raise ValueError('run_months, first_year and last_year must have equal length')

    species = select_cygnss_l1_species(obs_param, descr)

    inpath = exp_path + '/' + exp_run + '/output/' + domain

    t0_assim = t0_assim % dt_assim
    cycles   = range(t0_assim, 86400, dt_assim)

    if count_ObsFcstAna_files(inpath, exp_run, run_months, first_year, last_year, cycles) == 0:
        raise ValueError(f'no ObsFcstAna files for the requested period in {inpath}/ana/ens_avg')

    tile_coord     = read_tilecoord(tilecoord_file or inpath + '/rc_out/' + exp_run + '.ldas_tilecoord.bin')
    global_grid, _ = read_tilegrids(tilegrids_file or inpath + '/rc_out/' + exp_run + '.ldas_tilegrids.bin')

    tile_i, tile_j = get_owner_grid_ind(tile_coord, global_grid)

    N_tile = tile_coord['N_tile']

    # sums for each tile and day of the moving window:
    #  sum(obs), sum(obs^2), sum(fcst), sum(fcst^2), N, min(fcst), max(fcst)

    ring    = np.zeros((len(field_names), N_tile, w_days))
    ring[5] =  np.inf
    ring[6] = -np.inf

    data_out    = np.full((len(field_names), N_tile, N_pentads), np.nan)
    start_times = [datetime(2014, 1, 1) for _ in range(N_pentads)]
    end_times   = [datetime(2014, 1, 1) for _ in range(N_pentads)]

    n_read = 0                          # number of ObsFcstAna files read
    n_obs  = 0                          # number of CYGNSS L1 obs used

    count = 0

    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month[month-1] + 1):

            if count < w_days:
                count += 1
            idx = count - 1             # index of this day in the moving window

            for seconds in cycles:
                hour   = seconds // 3600
                minute = (seconds % 3600) // 60

                for year in range(first_year[imonth], last_year[imonth] + 1):

                    fname = get_ObsFcstAna_fname(inpath, exp_run, year, month, day, hour, minute)
                    ofa   = read_ObsFcstAna_nc4(fname)
                    if ofa is None:
                        continue
                    n_read += 1

                    check_ObsFcstAna_species_metadata(fname, species)

                    mask = ofa['obs_species'] == species['species']
                    if not np.any(mask):
                        continue

                    tiles = ofa['obs_tilenum'][mask].astype(np.int64) - 1

                    if np.any(tiles < 0) or np.any(tiles >= N_tile):
                        raise ValueError(f"species {species['species']}: tile number out of range in {fname}")
                    uniq, mult = np.unique(tiles, return_counts=True)
                    if np.any(mult > 1):
                        details = ', '.join(f'tile {t + 1} x{m}' for t, m in zip(uniq[mult > 1], mult[mult > 1]))
                        raise ValueError(f'more than one CYGNSS L1 obs per tile in {fname}: {details}')

                    obs  = ofa['obs_obs' ][mask].astype(np.float64)
                    fcst = ofa['obs_fcst'][mask].astype(np.float64)
                    ok   = np.isfinite(obs) & np.isfinite(fcst)
                    if obs_keep is not None:
                        ok &= obs_keep(datetime(year, month, day, hour, minute),
                                       ofa['obs_lon'][mask], ofa['obs_lat'][mask])

                    n_obs += np.count_nonzero(ok)

                    tiles = tiles[ok]
                    obs   = obs[ok]
                    fcst  = fcst[ok]

                    np.add.at(     ring[0, :, idx], tiles, obs    )
                    np.add.at(     ring[1, :, idx], tiles, obs**2 )
                    np.add.at(     ring[2, :, idx], tiles, fcst   )
                    np.add.at(     ring[3, :, idx], tiles, fcst**2)
                    np.add.at(     ring[4, :, idx], tiles, 1.     )
                    np.minimum.at( ring[5, :, idx], tiles, fcst   )
                    np.maximum.at( ring[6, :, idx], tiles, fcst   )

            if count < w_days:
                continue

            # The window is full: compute the statistics, store them for the center of
            # the window if it is the center of a pentad, and shift the window by one day.

            end_time   = datetime(2014, month, day, hour, minute)
            start_time = end_time - timedelta(seconds=w_days*86400)
            doy        = (end_time - timedelta(seconds=w_days*86400//2)).timetuple().tm_yday

            if (doy + 2) % 5 == 0:
                pentad = (doy + 2) // 5
                data_out[:, :, pentad-1] = get_window_stats_owner_tile(ring, ndata_min, std_epsilon)
                start_times[pentad-1]    = start_time
                end_times[  pentad-1]    = end_time

            ring[...,  :-1] = ring[..., 1:]
            ring[:5, :, -1] = 0.
            ring[5,  :, -1] =  np.inf
            ring[6,  :, -1] = -np.inf

    print(f'Read {n_read} ObsFcstAna files, {n_obs} CYGNSS L1 obs used')

    if n_obs == 0:
        raise ValueError('no CYGNSS L1 obs used (output file not written)')

    d1, d2, _, _ = get_doy_pentad_range(run_months, first_year, last_year)

    fname_out = (inpath + '/stats/' + out_dir + '/' + prefix +
                 f'{min(first_year)}_doy{d1}_{max(last_year)}_doy{d2}' +
                 f'_W_{w_days}d_Nmin_{ndata_min}_sp_{descr}_all_pentads.nc4')

    write_owner_grid_nc4(fname_out, data_out, tile_i, tile_j, global_grid, start_times, end_times,
                         w_days, ndata_min, std_epsilon, species['species'], descr, overwrite=overwrite)

    return fname_out

# ----------------------------------------------------------------------------
#
# CYGNSS L1 species: exactly one obsparam entry with descr, unscaled, in dB

def select_cygnss_l1_species(obs_param, descr=cygnss_l1_descr):

    matches = [p for p in obs_param if p['descr'] == descr]
    if len(matches) != 1:
        raise ValueError(f'expected one ldas_obsparam entry with descr {descr}; '
                         f'found species {[p["species"] for p in matches]}')
    p = matches[0]

    actual   = (p['varname'], p['units'], p['fcstvarname'], p['fcstunits'])
    expected = (cygnss_l1_varname, cygnss_l1_units, cygnss_l1_varname, cygnss_l1_units)
    if actual != expected:
        raise ValueError(f'unexpected CYGNSS L1 variable metadata {actual}; expected {expected}')
    if is_true(p['scale']):
        raise ValueError('CYGNSS L1 scaling parameters require an unscaled (monitoring) run')

    print(f"CYGNSS L1 species {int(p['species'])}: {p['descr']}")

    return {'species': int(p['species']), 'descr': p['descr']}

def is_true(value):

    value = value.strip().lower()
    if value in ('t', 'true', '.true.', '1'):
        return True
    if value in ('f', 'false', '.false.', '0'):
        return False
    raise ValueError(f'unrecognized logical value: {value}')

# ----------------------------------------------------------------------------
#
# check the obsparam metadata stored in the ObsFcstAna file (obsparam_* variables)
# for the CYGNSS L1 species: same descr, dB, and obs not scaled

def check_ObsFcstAna_species_metadata(fname, species):

    names = ['obsparam_species_id', 'obsparam_scale', 'obsparam_varname', 'obsparam_units',
             'obsparam_fcstvarname', 'obsparam_fcstunits', 'obsparam_descr']

    with Dataset(fname, 'r') as ncid:
        if not all(name in ncid.variables for name in names):
            raise ValueError(f'ObsFcstAna file has no obsparam metadata: {fname}')
        meta = {name: ncid.variables[name][:] for name in names}

    ind = np.flatnonzero(np.asarray(meta['obsparam_species_id']) == species['species'])
    if len(ind) != 1:
        raise ValueError(f"expected one obsparam entry for species {species['species']} in {fname}; "
                         f"found {len(ind)}")
    k = ind[0]

    actual   = tuple(str(meta[name][k]) for name in names[2:])
    expected = (cygnss_l1_varname, cygnss_l1_units, cygnss_l1_varname, cygnss_l1_units, species['descr'])
    if actual != expected:
        raise ValueError(f'CYGNSS L1 metadata mismatch in {fname}: {actual}; expected {expected}')
    if bool(meta['obsparam_scale'][k]):
        raise ValueError(f'CYGNSS L1 obs are already scaled in {fname}')

# ----------------------------------------------------------------------------
#
# 0-based column (i) and row (j) of each model tile in the global EASEv2 M36 grid

def get_owner_grid_ind(tile_coord, global_grid):

    if global_grid['gridtype'] != 'EASEv2_M36':
        raise ValueError(f"CYGNSS L1 owner-tile scaling requires EASEv2_M36, not {global_grid['gridtype']}")

    i = tile_coord['i_indg'].astype(np.int64) - global_grid['i_offg'] - global_grid['ind_base']
    j = tile_coord['j_indg'].astype(np.int64) - global_grid['j_offg'] - global_grid['ind_base']

    if np.any(i < 0) or np.any(i >= global_grid['N_lon']):
        raise ValueError('tilecoord i_indg outside the global tile grid')
    if np.any(j < 0) or np.any(j >= global_grid['N_lat']):
        raise ValueError('tilecoord j_indg outside the global tile grid')
    if np.unique(j * global_grid['N_lon'] + i).size != i.size:
        raise ValueError('more than one model tile in a global grid cell')

    return i, j

# ----------------------------------------------------------------------------
#
# statistics over the moving window, (7, N_tile):
#   o_mean, o_std, m_mean, m_std, N, m_min, m_max
# (all NaN if N < ndata_min or if the obs or model stdv is <= std_epsilon)

def get_window_stats_owner_tile(ring, ndata_min, std_epsilon):

    N = ring[4].sum(axis=-1)

    def mean(x):
        return np.divide(x, N, out=np.full(N.shape, np.nan), where=N > 0)

    o_mean = mean(ring[0].sum(axis=-1))
    m_mean = mean(ring[2].sum(axis=-1))
    o_std  = np.sqrt(np.maximum(mean(ring[1].sum(axis=-1)) - o_mean**2, 0.))
    m_std  = np.sqrt(np.maximum(mean(ring[3].sum(axis=-1)) - m_mean**2, 0.))

    m_min = ring[5].min(axis=-1)
    m_max = ring[6].max(axis=-1)
    m_min[~np.isfinite(m_min)] = np.nan
    m_max[~np.isfinite(m_max)] = np.nan

    data  = np.stack((o_mean, o_std, m_mean, m_std, N, m_min, m_max))
    valid = (N >= ndata_min) & (o_std > std_epsilon) & (m_std > std_epsilon)
    data[:, ~valid] = np.nan

    return data

# ----------------------------------------------------------------------------
#
# write the tile statistics (7, N_tile, N_pentads) on the dense global grid (nc4);
# m_min and m_max are written as the min/max over all pentads

def write_owner_grid_nc4(fname, data, tile_i, tile_j, grid, start_times, end_times,
                         w_days, ndata_min, std_epsilon, species_id, descr, overwrite=True):

    if os.path.exists(fname) and not overwrite:
        raise FileExistsError(f'{fname} exists and overwrite is False')
    os.makedirs(os.path.dirname(fname), exist_ok=True)

    N_lon     = grid['N_lon']
    N_lat     = grid['N_lat']
    N_pentads = data.shape[2]

    long_names = {
        'o_mean': 'Mean unscaled CYGNSS L1 observation',
        'o_std' : 'Population standard deviation of unscaled CYGNSS L1 observations',
        'm_mean': 'Mean CYGNSS L1 model forecast equivalent',
        'm_std' : 'Population standard deviation of CYGNSS L1 model forecast equivalents',
        'n_data': 'Number of paired observation and forecast samples',
        'm_min' : 'Minimum CYGNSS L1 model forecast equivalent over all pentads',
        'm_max' : 'Maximum CYGNSS L1 model forecast equivalent over all pentads',
    }

    def to_dense(values):
        dense = np.full((N_lat, N_lon), nodata)
        valid = np.isfinite(values)
        dense[tile_j[valid], tile_i[valid]] = values[valid]
        return dense

    with Dataset(fname, 'w', format='NETCDF4') as ncid:

        ncid.setncatts({
            'title'                         : 'CYGNSS L1 owner-tile z-score scaling climatology',
            'schema_version'                : '1.0',
            'gridtype'                      : grid['gridtype'],
            'source_ind_base'               : grid['ind_base'],
            'source_i_offg'                 : grid['i_offg'],
            'source_j_offg'                 : grid['j_offg'],
            'window_days'                   : w_days,
            'ndata_min'                     : ndata_min,
            'std_epsilon'                   : std_epsilon,
            'standard_deviation_convention' : 'population (ddof=0)',
            'species_id'                    : species_id,
            'species_description'           : descr,
            'spatial_key'                   : 'global owner tile i_indg/j_indg',
            'sampling_contract'             : 'at most one selected observation per owner tile and cycle',
        })

        ncid.createDimension('pentad', N_pentads)
        ncid.createDimension('y',      N_lat)
        ncid.createDimension('x',      N_lon)

        var = ncid.createVariable('pentad', 'i4', ('pentad',))
        var[:] = np.arange(1, N_pentads + 1, dtype=np.int32)
        var.setncatts({'long_name': 'pentad of year', 'units': '1', 'axis': 'T'})

        var = ncid.createVariable('x', 'i4', ('x',))
        var[:] = np.arange(N_lon, dtype=np.int32)
        var.setncatts({'long_name': 'zero-based global EASE-grid column index', 'units': '1', 'axis': 'X'})

        var = ncid.createVariable('y', 'i4', ('y',))
        var[:] = np.arange(N_lat, dtype=np.int32)
        var.setncatts({'long_name': 'zero-based global EASE-grid row index', 'units': '1', 'axis': 'Y'})

        for name, times in (('start_time', start_times), ('end_time', end_times)):
            var = ncid.createVariable(name, 'f8', ('pentad',))
            var[:] = days_since_1950(times)
            var.setncatts({'units': 'days since 1950-01-01 00:00:00', 'calendar': 'standard'})

        chunks = (1, min(64, N_lat), min(128, N_lon))

        for k, name in enumerate(field_names[:5]):
            var = ncid.createVariable(name, 'f8', ('pentad', 'y', 'x'), fill_value=nodata,
                                      zlib=True, complevel=5, chunksizes=chunks)
            var.setncatts({'long_name': long_names[name], 'units': '1' if name == 'n_data' else 'dB'})
            for ip in range(N_pentads):
                var[ip, :, :] = to_dense(data[k, :, ip])

        for k, name in ((5, 'm_min'), (6, 'm_max')):
            var = ncid.createVariable(name, 'f8', ('y', 'x'), fill_value=nodata,
                                      zlib=True, complevel=5, chunksizes=chunks[1:])
            var.setncatts({'long_name': long_names[name], 'units': 'dB'})
            finite = np.isfinite(data[k])
            if name == 'm_min':
                values = np.where(finite, data[k],  np.inf).min(axis=1)
            else:
                values = np.where(finite, data[k], -np.inf).max(axis=1)
            values[~finite.any(axis=1)] = np.nan
            var[:, :] = to_dense(values)

# ====================== EOF =========================================================
