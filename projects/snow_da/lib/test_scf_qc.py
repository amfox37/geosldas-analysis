"""
Unit tests for lib/scf_qc.py (synthetic cells; no data files needed).

run:  /Users/amfox/mamba/envs/regrid/bin/python -m pytest projects/snow_da/lib/test_scf_qc.py
The check against the Fortran readers on real files is scripts/viirs_obs_study/check_qc_vs_fortran.py.
"""

import numpy as np

from .scf_qc import (QC_VIIRS, cmg_lon_lat, modis_scf, viirs_qc_without_qa, viirs_scf)


def _viirs(sc, ci, qa):
    a = lambda x: np.array([x], dtype=np.int16)
    return {'Snow_Cover': a(sc), 'Clear_Index': a(ci), 'Basic_QA': a(qa)}


def _modis(sc, ci, qa):
    a = lambda x: np.array([x], dtype=np.int16)
    return {'Day_CMG_Snow_Cover': a(sc), 'Day_CMG_Clear_Index': a(ci), 'Snow_Spatial_QA': a(qa)}


def test_viirs_keeps_qa_0_1_3_and_drops_2():
    sc = [50, 50, 50, 50]; ci = [80, 80, 80, 80]; qa = [0, 1, 2, 3]
    scf, keep = viirs_scf(_viirs(sc, ci, qa))
    assert keep.tolist() == [[True, True, False, True]]
    assert np.isnan(scf[0, 2])


def test_viirs_thresholds():
    # Snow_Cover > 100 (flag), Clear_Index == 20 (not > 20), Clear_Index == 243 (Antarctica flag)
    scf, keep = viirs_scf(_viirs([101, 50, 100], [80, 20, 243], [0, 0, 0]))
    assert not keep.any()


def test_viirs_cap_at_one_and_single_precision():
    scf, keep = viirs_scf(_viirs([100, 33], [99, 97], [0, 0]))
    assert scf[0, 0] == np.float32(1.0)                       # 100/99 > 1 -> capped
    assert scf[0, 1] == np.float32(33) / np.float32(97)       # real*4 divide as in Fortran


def test_viirs_override_drops_class_1():
    scf, keep = viirs_scf(_viirs([50, 50], [80, 80], [0, 1]), {**QC_VIIRS, 'basic_qa_keep': (0, 3)})
    assert keep.tolist() == [[True, False]]


def test_viirs_qc_without_qa_ignores_qa_only():
    keep, qa = viirs_qc_without_qa(_viirs([50, 50, 101], [80, 80, 80], [2, 1, 0]))
    assert keep.tolist() == [[True, True, False]]
    assert qa.tolist() == [[2, 1, 0]]


def test_modis_qc_and_no_cap():
    scf, keep = modis_scf(_modis([50, 50, 50, 100], [80, 80, 20, 80], [2, 3, 0, 0]))
    assert keep.tolist() == [[True, False, False, True]]
    assert scf[0, 3] == np.float32(100) / np.float32(80)     # MODIS reader does not cap


def test_cmg_grid():
    lon, lat = cmg_lon_lat()
    assert lon.size == 7200 and lat.size == 3600
    assert np.isclose(lon[0], -179.975) and np.isclose(lon[-1], 179.975)
    assert np.isclose(lat[0], 89.975) and np.isclose(lat[-1], -89.975)
