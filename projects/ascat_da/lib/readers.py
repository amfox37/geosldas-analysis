"""
Readers for Legacy BUFR ASCAT, H121 CDR NetCDF, and GEOSldas OFA files.

Each reader returns a dict of 1-D numpy arrays keyed by field name, plus
a 'window' field (int 0-7, GEOSldas 3-hour analysis cycle index) and a
'cycle' field when timestamps are available.

GEOSldas ASCAT window convention:
    Cycles are centered on analysis times 0000, 0300, ..., 2100 UTC.
    For dtstep_assim=3h, cycle 1 (0300z) covers approximately 0130-0430.
    Legacy BUFR keeps [low, up); HSAF/H121 keeps (low, up].

Index convention:
    window=0 is the 0000z cycle
    window=1 is the 0300z cycle
    ...
    window=7 is the 2100z cycle

    cycle is an integer preserving day rollover relative to the requested date:
      0..7 are cycles on the requested date; 8 is next-day 0000z.
"""

import numpy as np
import glob
import os
from datetime import datetime, timedelta, timezone

import eccodes as ecc
import netCDF4 as nc

from .qc import apply_bufr_qc, apply_h121_qc, QC_DEFAULT_BUFR, QC_DEFAULT_H121


# ── Helpers ───────────────────────────────────────────────────────────────────

def _hour_to_geos_cycle(hour_utc, interval='left_closed'):
    """Fractional UTC hour -> GEOSldas centered 3-hour cycle.

    Parameters
    ----------
    hour_utc : array-like
        Fractional hour in UTC on the requested data day.
    interval : {'left_closed', 'right_closed'}
        GEOSldas legacy BUFR uses [low, up); HSAF/H121 uses (low, up].

    Returns
    -------
    cycle : ndarray of int
        Cycle index relative to the requested date. Values 0..7 are same-day
        cycles; 8 is the next day's 0000z cycle.
    window : ndarray of int
        Cycle index modulo 8.
    """
    hour = np.asarray(hour_utc, float)
    eps = 1e-9
    if interval == 'right_closed':
        shifted = hour + 1.5 - eps
    elif interval == 'left_closed':
        shifted = hour + 1.5
    else:
        raise ValueError("interval must be 'left_closed' or 'right_closed'")
    cycle = np.floor(shifted / 3.0).astype(int)
    return cycle, cycle % 8


def _inbox(lat, lon, domain):
    """Boolean mask: obs inside (lat0, lon0, lat1, lon1) bounding box."""
    if domain is None:
        return np.ones(len(lat), dtype=bool)
    lat0, lon0, lat1, lon1 = domain
    return ((lat >= lat0) & (lat <= lat1) &
            (lon >= lon0) & (lon <= lon1))


# ── Legacy BUFR reader ────────────────────────────────────────────────────────

_BUFR_KEYS = {
    'ssm':  'surfaceSoilMoisture',
    'smpf': 'soilMoistureProcessingFlag',
    'smcf': 'soilMoistureCorrectionFlag',
    'tpcx': 'topographicComplexity',
    'iwfr': 'inundationAndWetlandFraction',
    'lat':  'latitude',
    'lon':  'longitude',
    'hour': 'hour',
    'minute': 'minute',
}


def read_bufr(bufr_dir, date, prefix, domain=None, qc=None):
    """Read one day of Legacy BUFR ASCAT files after QC.

    Parameters
    ----------
    bufr_dir : str
        Directory containing .bfr files for the given date.
    date : datetime
        Date to read (matches filenames by YYYYMMDD).
    prefix : str
        Filename prefix, e.g. 'M02-ASCA-ASCSMO02-NA-5.0-'.
    domain : (lat0, lon0, lat1, lon1) or None
        If given, restrict to bounding box.
    qc : dict or None
        QC config; defaults to QC_DEFAULT_BUFR.

    Returns
    -------
    dict with arrays: lat, lon, ssm, window
    """
    if qc is None:
        qc = QC_DEFAULT_BUFR

    out = {k: [] for k in ('lat', 'lon', 'ssm', 'cycle', 'window')}

    pattern = os.path.join(bufr_dir, f'{prefix}{date.strftime("%Y%m%d")}*.bfr')
    for fpath in sorted(glob.glob(pattern)):
        with open(fpath, 'rb') as fh:
            bufr = ecc.codes_bufr_new_from_file(fh)
            while bufr is not None:
                ecc.codes_set(bufr, 'unpack', 1)
                try:
                    arrays = {}
                    for k, eckey in _BUFR_KEYS.items():
                        try:
                            arrays[k] = ecc.codes_get_array(bufr, eckey).astype(float)
                        except Exception:
                            arrays[k] = np.array([np.nan])
                    n = max(len(arrays[k]) for k in ('lat', 'lon', 'ssm'))
                    for k in arrays:
                        if len(arrays[k]) == 1:
                            arrays[k] = np.repeat(arrays[k], n)

                    mask = apply_bufr_qc(arrays, qc)
                    mask &= _inbox(arrays['lat'], arrays['lon'], domain)

                    hour_frac = arrays['hour'] + arrays['minute'] / 60.
                    out['lat'].extend(arrays['lat'][mask].tolist())
                    out['lon'].extend(arrays['lon'][mask].tolist())
                    out['ssm'].extend(arrays['ssm'][mask].tolist())
                    cyc, win = _hour_to_geos_cycle(hour_frac[mask], interval='left_closed')
                    out['cycle'].extend(cyc.tolist())
                    out['window'].extend(win.tolist())
                except Exception:
                    pass
                ecc.codes_release(bufr)
                bufr = ecc.codes_bufr_new_from_file(fh)

    return {k: np.array(v) for k, v in out.items()}


# ── H121 CDR NetCDF reader ────────────────────────────────────────────────────

def read_h121(h121_dir, date, domain=None, qc=None):
    """Read one day of H121 CDR NetCDF files after QC.

    Parameters
    ----------
    h121_dir : str
        Directory containing H121 .nc files for the given date.
    date : datetime
        Date to read (matched by YYYYMMDD in filename).
    domain : (lat0, lon0, lat1, lon1) or None
    qc : dict or None
        QC config; defaults to QC_DEFAULT_H121.

    Returns
    -------
    dict with arrays: lat, lon, ssm, subsfc, window
    """
    if qc is None:
        qc = QC_DEFAULT_H121

    out = {k: [] for k in ('lat', 'lon', 'ssm', 'subsfc', 'cycle', 'window')}

    pattern = os.path.join(h121_dir, f'*_{date.strftime("%Y%m%d")}*.nc')
    for fpath in sorted(glob.glob(pattern)):
        ds = nc.Dataset(fpath)
        try:
            lat    = ds.variables['latitude'][:].filled(np.nan)
            lon    = ds.variables['longitude'][:].filled(np.nan)
            ssm    = ds.variables['surface_soil_moisture'][:].filled(np.nan)
            def _load(name, fill):
                v = ds.variables[name][:]
                return v.data.astype(float) if hasattr(v, 'data') else np.asarray(v, float)

            sf     = _load('surface_flag', 255)
            pf     = _load('processing_flag', 255)
            cf     = _load('correction_flag', 255)
            wf     = _load('wetland_fraction', -128)
            tc     = _load('topographic_complexity', -128)
            subsfc = _load('subsurface_scattering_probability', -128)
            sens   = _load('surface_soil_moisture_sensitivity', np.inf)
            bsflag = _load('backscatter40_flag', 255)
            # fill value (-128) means the probability was not computed for that obs,
            # which is treated as 0% chance of subsurface contamination.
            fill_s = float(getattr(
                ds.variables['subsurface_scattering_probability'], '_FillValue', -128))
            subsfc[subsfc == fill_s] = 0.0

            # timestamp from file — H121 filenames encode hour as HH in ...YYYYMMDDTHHMMSS...
            # fall back to netCDF time variable if available
            try:
                t_var = ds.variables['time']
                t_units = getattr(t_var, 'units', '')
                times = nc.num2date(t_var[:], t_units)
                hour_frac = np.array([t.hour + t.minute / 60. for t in times],
                                     dtype=float)
            except Exception:
                # parse from filename: ..._YYYYMMDDTHHMMSS_...
                fname = os.path.basename(fpath)
                try:
                    tok = [p for p in fname.split('_') if 'T' in p and len(p) >= 13][0]
                    hh = int(tok[9:11])
                    mm = int(tok[11:13])
                    hour_frac = np.full(len(lat), hh + mm / 60., dtype=float)
                except Exception:
                    hour_frac = np.full(len(lat), np.nan)

            arrays = dict(ssm=ssm, sf=sf, pf=pf, cf=cf, tc=tc, wf=wf, subsfc=subsfc, sens=sens, bsflag=bsflag)
            mask = apply_h121_qc(arrays, qc)
            mask &= _inbox(lat, lon, domain)

            out['lat'].extend(lat[mask].tolist())
            out['lon'].extend(lon[mask].tolist())
            out['ssm'].extend(ssm[mask].tolist())
            out['subsfc'].extend(subsfc[mask].tolist())
            cyc, win = _hour_to_geos_cycle(hour_frac[mask], interval='right_closed')
            out['cycle'].extend(cyc.tolist())
            out['window'].extend(win.tolist())
        finally:
            ds.close()

    return {k: np.array(v) for k, v in out.items()}


# ── GEOSldas OFA reader ───────────────────────────────────────────────────────

def read_ofa(ofa_dir, date, species=None, domain=None):
    """Read GEOSldas ObsFcstAna files for one day.

    Parameters
    ----------
    ofa_dir : str
        Directory containing .nc4 OFA files (one per 3-hour window).
    date : datetime
        Date to read (matched by YYYYMMDD in filename).
    species : int, list of int, or None
        Filter to these species IDs. None = all species.
    domain : (lat0, lon0, lat1, lon1) or None

    Returns
    -------
    dict with arrays: lat, lon, obs, fcst, ana, innov, incr,
                      assim_flag, species, window (0–7)
    """
    if isinstance(species, int):
        species = [species]

    fields = ('lat', 'lon', 'obs', 'fcst', 'ana', 'assim_flag', 'species')
    out = {f: [] for f in fields}
    out['innov']  = []
    out['incr']   = []
    out['window'] = []

    pattern = os.path.join(ofa_dir, f'*{date.strftime("%Y%m%d")}*.nc4')
    for fpath in sorted(glob.glob(pattern)):
        # window index from filename timestamp e.g. ...20200101_0300z...
        fname = os.path.basename(fpath)
        try:
            tok = [p for p in fname.replace('.', '_').split('_') if p.endswith('z')][0]
            hh = int(tok[:2])
            win = int(hh // 3) % 8
        except Exception:
            win = -1

        ds = nc.Dataset(fpath)
        try:
            sp_arr   = ds.variables['species'][:]
            lat_arr  = ds.variables['lat'][:]
            lon_arr  = ds.variables['lon'][:]
            obs_arr  = ds.variables['obs'][:]
            fcst_arr = ds.variables['fcst'][:]
            ana_arr  = ds.variables['ana'][:]
            af_arr   = ds.variables['assim_flag'][:]

            mask = np.ones(len(sp_arr), dtype=bool)
            if species is not None:
                mask &= np.isin(sp_arr, species)
            mask &= _inbox(lat_arr, lon_arr, domain)

            n = mask.sum()
            out['lat'].extend(lat_arr[mask].tolist())
            out['lon'].extend(lon_arr[mask].tolist())
            out['obs'].extend(obs_arr[mask].tolist())
            out['fcst'].extend(fcst_arr[mask].tolist())
            out['ana'].extend(ana_arr[mask].tolist())
            out['innov'].extend((obs_arr[mask] - fcst_arr[mask]).tolist())
            out['incr'].extend((ana_arr[mask] - fcst_arr[mask]).tolist())
            out['assim_flag'].extend(af_arr[mask].tolist())
            out['species'].extend(sp_arr[mask].tolist())
            out['window'].extend([win] * n)
        finally:
            ds.close()

    return {k: np.array(v) for k, v in out.items()}
