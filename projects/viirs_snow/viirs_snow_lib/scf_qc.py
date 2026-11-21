"""
Python mirror of the snow cover fraction (SCF) QC in the GEOSldas readers.

Mirrors (GEOSldas_GridComp, GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90):

- ``read_JPSS_SCF_hdf``  -- VIIRS VNP10C1 / VJ110C1 / VJ210C1, version 2, hdf5
  (branch ``feature/amfox/jpss_hdf5_reader``, checked at commit 6d0ed32)
- ``read_MODIS_SCF_hdf`` -- MODIS MOD10C1 / MYD10C1, version 6.1, hdf4
  (develop, unchanged since reichle 2023)

Both products are daily on the 0.05 deg climate modeling grid (CMG): 7200 x 3600
cells, row 0 at the north pole, column 0 at -180 deg.

When a threshold changes in the Fortran reader, change it here too, and rerun
``scripts/viirs_obs_study/check_qc_vs_fortran.py`` (see the README there).
Last check: GridComp 6d0ed32, 9 files (VIIRS J1/J2/S-NPP 2020 and 2024, MODIS Aqua
2020): same cells in the same order, bit-identical SCF.

QC dicts are plain dicts; override single thresholds with e.g.
``{**QC_VIIRS, 'basic_qa_keep': (0, 3)}``.

VIIRS ``Basic_QA`` key (CMG file attribute ``key_to_QA_values`` and VNP10C1 v2
user guide): 0 = good, 1 = poor, 2 = bad, 3 = other (low illumination,
70 <= SZA < 85 deg).  The reader keeps 0, 1 and 3.
"""

import numpy as np

# ── CMG grid ──────────────────────────────────────────────────────────────────

CMG_N_LON = 7200
CMG_N_LAT = 3600
CMG_DLON = 0.05
CMG_DLAT = 0.05


def cmg_lon_lat():
    """Cell-center lon (7200,) and lat (3600,) of the CMG grid, as in the readers."""
    lon = -180. + 0.5 * CMG_DLON + CMG_DLON * np.arange(CMG_N_LON)
    lat = 90. - 0.5 * CMG_DLAT - CMG_DLAT * np.arange(CMG_N_LAT)
    return lon, lat


# ── QC configs (values as in the Fortran readers) ─────────────────────────────

QC_VIIRS = dict(
    snow_cover_max=100,          # Snow_Cover <= 100: excludes night, lake, ocean, fill, no decision
    clear_index_min=20,          # Clear_Index > 20: sufficiently clear
    clear_index_max=100,         # Clear_Index <= 100: excludes flag values (Antarctica = 243)
    basic_qa_keep=(0, 1, 3),     # Basic_QA classes kept (2 = "bad" dropped)
    scf_cap=1.0,                 # SCF = Snow_Cover/Clear_Index capped at 1 (VIIRS can reach ~1.05)
)

QC_MODIS = dict(
    snow_cover_max=100,          # Day_CMG_Snow_Cover <= 100
    clear_index_min=20,          # Day_CMG_Clear_Index > 20
    snow_spatial_qa_max=2,       # Snow_Spatial_QA <= 2 (0=best, 1=good, 2=OK; 3=poor, 4=other dropped)
    scf_cap=None,                # MODIS reader does not cap (MODIS SCF never exceeds 1)
)

VIIRS_FIELD_PATH = 'HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields'
VIIRS_FIELDS = ('Snow_Cover', 'Clear_Index', 'Basic_QA')
MODIS_FIELDS = ('Day_CMG_Snow_Cover', 'Day_CMG_Clear_Index', 'Snow_Spatial_QA')


# ── Readers (raw integer fields, full grid) ───────────────────────────────────

def read_viirs_cmg(path):
    """Raw VIIRS CMG fields as int16 arrays (3600, 7200): Snow_Cover, Clear_Index, Basic_QA."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        grp = ds[VIIRS_FIELD_PATH]
        out = {}
        for name in VIIRS_FIELDS:
            var = grp[name]
            var.set_auto_maskandscale(False)
            out[name] = var[:].astype(np.int16)
        out['LocalGranuleID'] = getattr(ds, 'LocalGranuleID', '')
    return out


def read_modis_cmg(path):
    """Raw MODIS CMG fields as int16 arrays (3600, 7200): Day_CMG_Snow_Cover, Day_CMG_Clear_Index, Snow_Spatial_QA."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        out = {}
        for name in MODIS_FIELDS:
            var = ds[name]
            var.set_auto_maskandscale(False)
            out[name] = var[:].astype(np.int16)
    return out


# ── QC ────────────────────────────────────────────────────────────────────────

def viirs_scf(raw, qc=QC_VIIRS):
    """
    SCF (float32, NaN where rejected) and the keep mask, as in read_JPSS_SCF_hdf.

    SCF = Snow_Cover / Clear_Index in single precision (as in the Fortran reader),
    then capped at qc['scf_cap'].
    """
    sc, ci, qa = raw['Snow_Cover'], raw['Clear_Index'], raw['Basic_QA']
    keep = ((sc <= qc['snow_cover_max']) &
            (ci > qc['clear_index_min']) &
            (ci <= qc['clear_index_max']) &
            np.isin(qa, qc['basic_qa_keep']))
    scf = np.full(sc.shape, np.nan, dtype=np.float32)
    scf[keep] = sc[keep].astype(np.float32) / ci[keep].astype(np.float32)
    if qc['scf_cap'] is not None:
        np.minimum(scf, np.float32(qc['scf_cap']), out=scf, where=keep)
    return scf, keep


def modis_scf(raw, qc=QC_MODIS):
    """SCF (float32, NaN where rejected) and the keep mask, as in read_MODIS_SCF_hdf."""
    sc, ci, qa = raw['Day_CMG_Snow_Cover'], raw['Day_CMG_Clear_Index'], raw['Snow_Spatial_QA']
    keep = ((sc <= qc['snow_cover_max']) &
            (ci > qc['clear_index_min']) &
            (qa <= qc['snow_spatial_qa_max']))
    scf = np.full(sc.shape, np.nan, dtype=np.float32)
    scf[keep] = sc[keep].astype(np.float32) / ci[keep].astype(np.float32)
    if qc['scf_cap'] is not None:
        np.minimum(scf, np.float32(qc['scf_cap']), out=scf, where=keep)
    return scf, keep


def viirs_qc_without_qa(raw, qc=QC_VIIRS):
    """
    VIIRS keep mask with every reader test except Basic_QA, plus the QA grid.

    Used to study the QA classes separately: ``keep & (qa == k)`` gives the cells
    that the reader would keep if class k were accepted.
    """
    sc, ci = raw['Snow_Cover'], raw['Clear_Index']
    keep = ((sc <= qc['snow_cover_max']) &
            (ci > qc['clear_index_min']) &
            (ci <= qc['clear_index_max']))
    return keep, raw['Basic_QA']
