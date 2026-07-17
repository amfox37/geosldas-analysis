from datetime import date
from pathlib import Path
import struct
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

netCDF4 = pytest.importorskip("netCDF4")

from iv_tc.readers import (  # noqa: E402
    _ease2_m36_lon_lat_for_ij,
    read_ascat_h119_h120_model_pair,
    read_ascat_h119_h120_sparse,
    read_ascat_h121_model_pair,
    read_cygnss_l3_model_pair,
    read_cygnss_l3_sparse,
    read_smap_l3_model_pair,
    read_smap_l3_sparse,
    read_smosic_model_pair,
    read_smosic_sparse,
    read_tilecoord,
    representative_tiles_from_tilecoord,
)


DISCOVER_SAMPLE = PROJECT_ROOT / "test_data" / "inputs" / "discover_sample"


def test_read_smosic_sparse_preserves_indices_values_and_qc(tmp_path):
    path = tmp_path / "smos_ic_sm_m36_20200102.nc"
    _write_smosic_sparse(
        path,
        idx=[12, 13, 14],
        values=[0.21, 0.32, 0.43],
        coverage=[0.4, 0.6, 0.8],
    )

    obs = read_smosic_sparse(path)

    assert obs.date == date(2020, 1, 2)
    assert obs.sensor == "SMOS-IC"
    assert obs.units == "m3 m-3"
    np.testing.assert_array_equal(obs.idx, np.array([12, 13, 14]))
    np.testing.assert_allclose(obs.values, np.array([0.21, 0.32, 0.43]), rtol=1e-6)
    assert obs.qc_summary["n_points"] == 3
    assert obs.qc_summary["coverage_min"] == pytest.approx(0.4)
    assert obs.qc_summary["coverage_mean"] == pytest.approx(0.6)
    assert obs.qc_summary["coverage_max"] == pytest.approx(0.8)
    assert obs.qc_summary["qc_tb_rmse_max"] == pytest.approx(8.0)


def test_read_tilecoord_and_representative_tiles_use_fraction_then_tile_id(tmp_path):
    path = tmp_path / "test.ldas_tilecoord.bin"
    _write_tilecoord(path)

    tilecoord = read_tilecoord(path)
    reps = representative_tiles_from_tilecoord(tilecoord, nx=10)

    np.testing.assert_array_equal(tilecoord["N_tile"], np.array(5, dtype=np.int32))
    np.testing.assert_array_equal(reps.m36_linear, np.array([12, 13, 14]))
    np.testing.assert_array_equal(reps.tile_index, np.array([3, 2, 4]))


def test_read_smosic_model_pair_matches_representative_model_tiles(tmp_path):
    smosic_path = tmp_path / "smos_ic_sm_m36_20200102.nc"
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    model_path = tmp_path / "model.nc4"

    _write_smosic_sparse(
        smosic_path,
        idx=[12, 13, 99],
        values=[0.21, 0.32, 0.99],
        coverage=[0.4, 0.6, 1.0],
    )
    _write_tilecoord(tilecoord_path)
    _write_model(model_path, values=[0.01, 0.02, 0.33, 0.14, 0.55])

    pair = read_smosic_model_pair(smosic_path, model_path, tilecoord_path, run="OL", nx=10)

    assert pair.date == date(2020, 1, 2)
    assert pair.sensor == "SMOS-IC"
    assert pair.run == "OL"
    assert pair.obs_units == "m3 m-3"
    assert pair.model_units == "m3 m-3"
    np.testing.assert_array_equal(pair.idx, np.array([12, 13]))
    np.testing.assert_allclose(pair.obs, np.array([0.21, 0.32]), rtol=1e-6)
    np.testing.assert_allclose(pair.model, np.array([0.14, 0.33]), rtol=1e-6)


def test_read_smap_l3_sparse_matches_matlab_qc_and_am_pm_mean(tmp_path):
    path = tmp_path / "SMAP_L3_SM_P_20200102_R19240_001.h5"
    _write_smap_l3(path)

    obs = read_smap_l3_sparse(path, nx=5)

    assert obs.date == date(2020, 1, 2)
    assert obs.sensor == "SMAP L3"
    assert obs.units == "m3 m-3"
    np.testing.assert_array_equal(obs.idx, np.array([7, 8]))
    np.testing.assert_allclose(obs.values, np.array([0.30, 0.50]), rtol=1e-6)
    assert obs.qc_summary["am_raw_finite"] == 3
    assert obs.qc_summary["am_kept"] == 1
    assert obs.qc_summary["pm_raw_finite"] == 3
    assert obs.qc_summary["pm_kept"] == 2
    assert obs.qc_summary["n_points"] == 2


def test_read_smap_l3_model_pair_matches_representative_model_tiles(tmp_path):
    smap_path = tmp_path / "SMAP_L3_SM_P_20200102_R19240_001.h5"
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    model_path = tmp_path / "model.nc4"

    _write_smap_l3(smap_path)
    _write_tilecoord(tilecoord_path)
    _write_model(model_path, values=[0.01, 0.02, 0.33, 0.14, 0.55])

    pair = read_smap_l3_model_pair(smap_path, model_path, tilecoord_path, run="OL", nx=10)

    assert pair.date == date(2020, 1, 2)
    assert pair.sensor == "SMAP L3"
    assert pair.run == "OL"
    assert pair.obs_units == "m3 m-3"
    assert pair.model_units == "m3 m-3"
    np.testing.assert_array_equal(pair.idx, np.array([12, 13]))
    np.testing.assert_allclose(pair.obs, np.array([0.30, 0.50]), rtol=1e-6)
    np.testing.assert_allclose(pair.model, np.array([0.14, 0.33]), rtol=1e-6)


def test_read_cygnss_l3_sparse_uses_daily_field_qc_and_spatial_orientation(tmp_path):
    cygnss_path = (
        tmp_path
        / "cyg.ddmi.s20200102-030000-e20200102-210000."
        "l3.grid-soil-moisture-36km.a32.d33.nc"
    )
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    target_i = np.array([480, 481, 480, 481], dtype=np.int64)
    target_j = np.array([200, 200, 201, 201], dtype=np.int64)
    lon, lat = _ease2_m36_lon_lat_for_ij(target_i, target_j)

    _write_tilecoord_cells(tilecoord_path, target_i, target_j)
    _write_cygnss_l3(cygnss_path, lon, lat)

    obs = read_cygnss_l3_sparse(cygnss_path, tilecoord_path)

    assert obs.date == date(2020, 1, 2)
    assert obs.sensor == "CYGNSS L3"
    assert obs.units == "m3 m-3"
    np.testing.assert_array_equal(obs.idx, np.array([193280, 194244, 193281]))
    # SM_subdaily is intentionally 0.9 everywhere in this fixture. These are
    # SM_daily values after the daily sigma mask and x/y reorientation.
    np.testing.assert_allclose(obs.values, np.array([0.10, 0.30, 0.20]), rtol=1e-6)
    assert obs.qc_summary["raw_finite"] == 4
    assert obs.qc_summary["kept"] == 3
    assert obs.qc_summary["source_variable"] == "SM_daily"
    assert obs.qc_summary["sigma_variable"] == "SIGMA_daily"
    assert obs.qc_summary["fill_nearest"] == 0


def test_read_cygnss_l3_model_pair_matches_representative_model_tiles(tmp_path):
    cygnss_path = (
        tmp_path
        / "cyg.ddmi.s20200102-030000-e20200102-210000."
        "l3.grid-soil-moisture-36km.a32.d33.nc"
    )
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    model_path = tmp_path / "model.nc4"
    target_i = np.array([480, 481, 480, 481], dtype=np.int64)
    target_j = np.array([200, 200, 201, 201], dtype=np.int64)
    lon, lat = _ease2_m36_lon_lat_for_ij(target_i, target_j)

    _write_tilecoord_cells(tilecoord_path, target_i, target_j)
    _write_cygnss_l3(cygnss_path, lon, lat)
    _write_model(model_path, values=[0.11, 0.22, 0.33, 0.44])

    pair = read_cygnss_l3_model_pair(cygnss_path, model_path, tilecoord_path, run="OL")

    assert pair.date == date(2020, 1, 2)
    assert pair.sensor == "CYGNSS L3"
    assert pair.run == "OL"
    np.testing.assert_array_equal(pair.idx, np.array([193280, 194244, 193281]))
    np.testing.assert_allclose(pair.obs, np.array([0.10, 0.30, 0.20]), rtol=1e-6)
    np.testing.assert_allclose(pair.model, np.array([0.11, 0.33, 0.22]), rtol=1e-6)


def test_copied_discover_sample_pairs_cygnss_l3_with_ol_model():
    cygnss_path = (
        DISCOVER_SAMPLE
        / "CYGNSS"
        / "Y2018"
        / "M10"
        / "cyg.ddmi.s20181015-030000-e20181015-210000."
        "l3.grid-soil-moisture-36km.a32.d33.nc"
    )
    model_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "cat"
        / "ens_avg"
        / "Y2018"
        / "M10"
        / "OLv7_M36_MULTI_type_13_H121.tavg24_1d_lnd_Nt.20181015_1200z.nc4"
    )
    tilecoord_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "rc_out"
        / "OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin"
    )
    missing = [path for path in (cygnss_path, model_path, tilecoord_path) if not path.exists()]
    if missing:
        pytest.skip(f"Discover sample fixture is not present: {missing[0]}")

    pair = read_cygnss_l3_model_pair(cygnss_path, model_path, tilecoord_path, run="OL")

    assert pair.date == date(2018, 10, 15)
    assert pair.sensor == "CYGNSS L3"
    assert pair.run == "OL"
    assert pair.idx.size > 1000
    assert pair.idx.size == pair.obs.size == pair.model.size
    assert np.isfinite(pair.obs).all()
    assert np.isfinite(pair.model).all()
    assert pair.obs_units == "m3 m-3"
    assert pair.model_units == "m3 m-3"


def test_read_ascat_h119_h120_sparse_uses_matlab_qc_and_linear_interpolation(tmp_path):
    mat_path = tmp_path / "ASCAT_HSAF_H119_SM_20200102_AD.mat"
    aux_path = tmp_path / "TUW_WARP5_grid_info_2_2.nc"
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    target_i = np.array([480, 481, 480, 481], dtype=np.int64)
    target_j = np.array([200, 200, 201, 201], dtype=np.int64)
    lon, lat = _ease2_m36_lon_lat_for_ij(target_i, target_j)

    _write_tilecoord_cells(tilecoord_path, target_i, target_j)
    _write_ascat_aux(aux_path, lon, lat)
    _write_ascat_h119_h120_mat(mat_path, sm=[10.0, 20.0, 30.0, 40.0], conf=[0, 0, 0, 1])

    obs = read_ascat_h119_h120_sparse(mat_path, aux_path, tilecoord_path, fill_nearest=False)

    assert obs.date == date(2020, 1, 2)
    assert obs.sensor == "ASCAT H119/H120"
    assert obs.units == "percent saturation"
    np.testing.assert_array_equal(obs.idx, np.array([193280, 194244, 193281]))
    np.testing.assert_allclose(obs.values, np.array([10.0, 30.0, 20.0]), rtol=1e-6)
    assert obs.qc_summary["raw_finite"] == 4
    assert obs.qc_summary["kept"] == 3
    assert obs.qc_summary["n_targets"] == 4
    assert obs.qc_summary["n_points"] == 3


def test_read_ascat_h119_h120_sparse_default_leaves_gaps_outside_hull_as_nan(tmp_path):
    mat_path = tmp_path / "ASCAT_HSAF_H119_SM_20200102_AD.mat"
    aux_path = tmp_path / "TUW_WARP5_grid_info_2_2.nc"
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"

    gpi_i = np.array([480, 481, 480, 481], dtype=np.int64)
    gpi_j = np.array([200, 200, 201, 201], dtype=np.int64)
    gpi_lon, gpi_lat = _ease2_m36_lon_lat_for_ij(gpi_i, gpi_j)

    # One target tile sits far outside the GPI convex hull.
    target_i = np.array([480, 481, 480, 481, 700], dtype=np.int64)
    target_j = np.array([200, 200, 201, 201, 350], dtype=np.int64)

    _write_tilecoord_cells(tilecoord_path, target_i, target_j)
    _write_ascat_aux(aux_path, gpi_lon, gpi_lat)
    _write_ascat_h119_h120_mat(mat_path, sm=[10.0, 20.0, 30.0, 40.0], conf=[0, 0, 0, 0])

    obs = read_ascat_h119_h120_sparse(mat_path, aux_path, tilecoord_path)

    assert obs.qc_summary["n_targets"] == 5
    assert obs.qc_summary["fill_nearest"] == 0
    far_idx = 700 + 350 * 964
    assert far_idx not in obs.idx.tolist()
    assert obs.qc_summary["n_points"] < obs.qc_summary["n_targets"]


def test_read_ascat_h119_h120_model_pair_matches_representative_model_tiles(tmp_path):
    mat_path = tmp_path / "ASCAT_HSAF_H119_SM_20200102_AD.mat"
    aux_path = tmp_path / "TUW_WARP5_grid_info_2_2.nc"
    tilecoord_path = tmp_path / "test.ldas_tilecoord.bin"
    model_path = tmp_path / "model.nc4"
    target_i = np.array([480, 481, 480, 481], dtype=np.int64)
    target_j = np.array([200, 200, 201, 201], dtype=np.int64)
    lon, lat = _ease2_m36_lon_lat_for_ij(target_i, target_j)

    _write_tilecoord_cells(tilecoord_path, target_i, target_j)
    _write_ascat_aux(aux_path, lon, lat)
    _write_ascat_h119_h120_mat(mat_path, sm=[10.0, 20.0, 30.0, 40.0], conf=[0, 0, 0, 1])
    _write_model(model_path, values=[0.11, 0.22, 0.33, 0.44])

    pair = read_ascat_h119_h120_model_pair(
        mat_path,
        aux_path,
        model_path,
        tilecoord_path,
        run="OL",
        fill_nearest=False,
    )

    assert pair.date == date(2020, 1, 2)
    assert pair.sensor == "ASCAT H119/H120"
    assert pair.run == "OL"
    assert pair.obs_units == "percent saturation"
    assert pair.model_units == "m3 m-3"
    np.testing.assert_array_equal(pair.idx, np.array([193280, 194244, 193281]))
    np.testing.assert_allclose(pair.obs, np.array([10.0, 30.0, 20.0]), rtol=1e-6)
    np.testing.assert_allclose(pair.model, np.array([0.11, 0.33, 0.22]), rtol=1e-6)


def test_copied_discover_sample_pairs_smosic_with_ol_model():
    smosic_path = DISCOVER_SAMPLE / "SMOS_IC" / "preprocessed_m36_daily" / "smos_ic_sm_m36_20181015.nc"
    model_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "cat"
        / "ens_avg"
        / "Y2018"
        / "M10"
        / "OLv7_M36_MULTI_type_13_H121.tavg24_1d_lnd_Nt.20181015_1200z.nc4"
    )
    tilecoord_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "rc_out"
        / "OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin"
    )
    missing = [path for path in (smosic_path, model_path, tilecoord_path) if not path.exists()]
    if missing:
        pytest.skip(f"Discover sample fixture is not present: {missing[0]}")

    pair = read_smosic_model_pair(smosic_path, model_path, tilecoord_path, run="OL")

    assert pair.date == date(2018, 10, 15)
    assert pair.sensor == "SMOS-IC"
    assert pair.run == "OL"
    assert pair.idx.size > 1000
    assert pair.idx.size == pair.obs.size == pair.model.size
    assert np.isfinite(pair.obs).all()
    assert np.isfinite(pair.model).all()
    assert pair.obs_units == "m3 m-3"
    assert pair.model_units == "m3 m-3"


def test_copied_discover_sample_pairs_smap_l3_with_ol_model():
    smap_path = DISCOVER_SAMPLE / "SPL3SMP_v009" / "Y2018" / "SMAP_L3_SM_P_20181015_R19240_001.h5"
    model_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "cat"
        / "ens_avg"
        / "Y2018"
        / "M10"
        / "OLv7_M36_MULTI_type_13_H121.tavg24_1d_lnd_Nt.20181015_1200z.nc4"
    )
    tilecoord_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "rc_out"
        / "OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin"
    )
    missing = [path for path in (smap_path, model_path, tilecoord_path) if not path.exists()]
    if missing:
        pytest.skip(f"Discover sample fixture is not present: {missing[0]}")

    pair = read_smap_l3_model_pair(smap_path, model_path, tilecoord_path, run="OL")

    assert pair.date == date(2018, 10, 15)
    assert pair.sensor == "SMAP L3"
    assert pair.run == "OL"
    assert pair.idx.size > 1000
    assert pair.idx.size == pair.obs.size == pair.model.size
    assert np.isfinite(pair.obs).all()
    assert np.isfinite(pair.model).all()
    assert pair.obs_units == "m3 m-3"
    assert pair.model_units == "m3 m-3"


def test_copied_discover_sample_pairs_h119_h120_with_ol_model():
    mat_path = (
        DISCOVER_SAMPLE
        / "ASCAT_HSAF"
        / "H119_H120_processed"
        / "Y2018"
        / "M10"
        / "ASCAT_HSAF_H119_SM_20181015_AD.mat"
    )
    aux_path = DISCOVER_SAMPLE / "ASCAT_HSAF" / "Auxiliary" / "TUW_WARP5_grid_info_2_2.nc"
    model_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "cat"
        / "ens_avg"
        / "Y2018"
        / "M10"
        / "OLv7_M36_MULTI_type_13_H121.tavg24_1d_lnd_Nt.20181015_1200z.nc4"
    )
    tilecoord_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "rc_out"
        / "OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin"
    )
    missing = [path for path in (mat_path, aux_path, model_path, tilecoord_path) if not path.exists()]
    if missing:
        pytest.skip(f"Discover sample fixture is not present: {missing[0]}")

    pair = read_ascat_h119_h120_model_pair(
        mat_path,
        aux_path,
        model_path,
        tilecoord_path,
        run="OL",
        method="nearest",
    )

    assert pair.date == date(2018, 10, 15)
    assert pair.sensor == "ASCAT H119/H120"
    assert pair.run == "OL"
    assert pair.idx.size > 1000
    assert pair.idx.size == pair.obs.size == pair.model.size
    assert np.isfinite(pair.obs).all()
    assert np.isfinite(pair.model).all()
    assert np.nanmin(pair.obs) >= 0.0
    assert np.nanmax(pair.obs) <= 100.0
    assert pair.obs_units == "percent saturation"
    assert pair.model_units == "m3 m-3"


def test_copied_discover_sample_pairs_h121_with_ol_model(tmp_path):
    h121_source_dir = DISCOVER_SAMPLE / "ASCAT_SSM_CDR" / "H121" / "metop_a" / "Y2018" / "M10"
    model_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "cat"
        / "ens_avg"
        / "Y2018"
        / "M10"
        / "OLv7_M36_MULTI_type_13_H121.tavg24_1d_lnd_Nt.20181015_1200z.nc4"
    )
    tilecoord_path = (
        DISCOVER_SAMPLE
        / "runs"
        / "OL"
        / "output"
        / "SMAP_EASEv2_M36_GLOBAL"
        / "rc_out"
        / "OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin"
    )
    h121_files = sorted(h121_source_dir.glob("*.nc"))
    missing = [path for path in (model_path, tilecoord_path) if not path.exists()]
    if missing or not h121_files:
        skip_path = missing[0] if missing else h121_source_dir
        pytest.skip(f"Discover sample fixture is not present: {skip_path}")

    h121_test_dir = tmp_path / "h121"
    h121_test_dir.mkdir()
    (h121_test_dir / h121_files[0].name).symlink_to(h121_files[0].resolve())

    pair = read_ascat_h121_model_pair(h121_test_dir, date(2018, 10, 15), model_path, tilecoord_path, run="OL")

    assert pair.date == date(2018, 10, 15)
    assert pair.sensor == "ASCAT H121"
    assert pair.run == "OL"
    assert pair.idx.size > 100
    assert pair.idx.size == pair.obs.size == pair.model.size
    assert np.isfinite(pair.obs).all()
    assert np.isfinite(pair.model).all()
    assert pair.obs_units == "percent saturation"
    assert pair.model_units == "m3 m-3"


def _write_smosic_sparse(path, idx, values, coverage):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("n_points", len(idx))
        idx_var = ds.createVariable("idx_EASEv2_lonxlat", "i4", ("n_points",))
        sm_var = ds.createVariable("sm_obs", "f4", ("n_points",))
        cov_var = ds.createVariable("coverage_frac", "f4", ("n_points",))

        idx_var[:] = np.asarray(idx, dtype=np.int32)
        sm_var[:] = np.asarray(values, dtype=np.float32)
        cov_var[:] = np.asarray(coverage, dtype=np.float32)

        sm_var.units = "m3 m-3"
        ds.date = "20200102"
        ds.qc_scene_flag_max = 1
        ds.qc_tb_rmse_max = 8.0
        ds.qc_sm_min = 0.0
        ds.qc_sm_max = 1.0
        ds.min_coverage_frac = 0.05


def _write_smap_l3(path):
    fill = -9999.0
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("y", 2)
        ds.createDimension("x", 5)

        am = ds.createGroup("Soil_Moisture_Retrieval_Data_AM")
        _write_smap_group(
            am,
            sm_name="soil_moisture",
            qf_name="retrieval_qual_flag",
            suffix="",
            sm=[
                [fill, fill, fill, fill, fill],
                [fill, fill, 0.20, 0.30, 0.40],
            ],
            qf=[
                [0, 0, 0, 0, 0],
                [0, 0, 0, 4, 0],
            ],
            surface_status=np.zeros((2, 5), dtype=np.uint16),
            surface_flag=[
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 32],
            ],
        )

        pm = ds.createGroup("Soil_Moisture_Retrieval_Data_PM")
        surface_status_pm = np.zeros((2, 5), dtype=np.uint16)
        surface_status_pm[0, 1] = 1
        _write_smap_group(
            pm,
            sm_name="soil_moisture_dca_pm",
            qf_name="retrieval_qual_flag_dca_pm",
            suffix="_pm",
            sm=[
                [fill, 0.10, fill, fill, fill],
                [fill, fill, 0.40, 0.50, fill],
            ],
            qf=np.zeros((2, 5), dtype=np.uint16),
            surface_status=surface_status_pm,
            surface_flag=np.zeros((2, 5), dtype=np.uint16),
        )


def _write_smap_group(group, sm_name, qf_name, suffix, sm, qf, surface_status, surface_flag):
    sm_var = group.createVariable(sm_name, "f4", ("y", "x"), fill_value=-9999.0)
    sm_var.units = "cm**3/cm**3"
    sm_var[:] = np.asarray(sm, dtype=np.float32)

    qf_var = group.createVariable(qf_name, "u2", ("y", "x"), fill_value=np.uint16(65534))
    qf_var[:] = np.asarray(qf, dtype=np.uint16)

    status_var = group.createVariable(f"grid_surface_status{suffix}", "u2", ("y", "x"), fill_value=np.uint16(65534))
    status_var[:] = np.asarray(surface_status, dtype=np.uint16)

    flag_var = group.createVariable(f"surface_flag{suffix}", "u2", ("y", "x"), fill_value=np.uint16(65534))
    flag_var[:] = np.asarray(surface_flag, dtype=np.uint16)

    rfi_h_var = group.createVariable(f"tb_qual_flag_h{suffix}", "u2", ("y", "x"), fill_value=np.uint16(65534))
    rfi_v_var = group.createVariable(f"tb_qual_flag_v{suffix}", "u2", ("y", "x"), fill_value=np.uint16(65534))
    rfi_h_var[:] = np.zeros((2, 5), dtype=np.uint16)
    rfi_v_var[:] = np.zeros((2, 5), dtype=np.uint16)


def _write_cygnss_l3(path, lon, lat):
    # lon/lat are supplied in x/y order. SM_daily is intentionally stored in
    # the opposite y/x order behind a singleton time dimension.
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("x", 2)
        ds.createDimension("y", 2)
        ds.createDimension("time", 1)
        lon_var = ds.createVariable("longitude", "f8", ("x", "y"))
        lat_var = ds.createVariable("latitude", "f8", ("x", "y"))
        sm_var = ds.createVariable("SM_daily", "f4", ("time", "y", "x"), fill_value=-9999.0)
        sigma_var = ds.createVariable("SIGMA_daily", "f4", ("time", "y", "x"), fill_value=-9999.0)
        subdaily = ds.createVariable("SM_subdaily", "f4", ("time", "y", "x"), fill_value=-9999.0)

        lon_var[:] = np.asarray(lon, dtype=np.float64).reshape(2, 2)
        lat_var[:] = np.asarray(lat, dtype=np.float64).reshape(2, 2)
        daily_xy = np.array([[0.10, 0.30], [0.20, 0.40]], dtype=np.float32)
        sigma_xy = np.array([[0.01, 0.01], [0.01, np.nan]], dtype=np.float32)
        sm_var[0, :, :] = daily_xy.T
        sigma_var[0, :, :] = sigma_xy.T
        subdaily[0, :, :] = np.full((2, 2), 0.9, dtype=np.float32)
        sm_var.units = "m3/m3"
        sm_var.valid_min = 0.0
        sm_var.valid_max = 1.0


def _write_ascat_aux(path, lon, lat):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("gpi", len(lon) + 1)
        land_flag = ds.createVariable("land_flag", "i4", ("gpi",))
        lon_var = ds.createVariable("lon", "f8", ("gpi",))
        lat_var = ds.createVariable("lat", "f8", ("gpi",))

        land_flag[:] = np.asarray([*([1] * len(lon)), 0], dtype=np.int32)
        lon_var[:] = np.asarray([*lon, 0.0], dtype=np.float64)
        lat_var[:] = np.asarray([*lat, 0.0], dtype=np.float64)


def _write_ascat_h119_h120_mat(path, sm, conf):
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.savemat(
        path,
        {
            "sm_tile": np.asarray(sm, dtype=np.float64),
            "conf_flag_tile": np.asarray(conf, dtype=np.int16),
        },
    )


def _write_model(path, values):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("tile", len(values))
        var = ds.createVariable("SFMC", "f4", ("time", "tile"))
        var.units = "m3 m-3"
        var[0, :] = np.asarray(values, dtype=np.float32)


def _write_tilecoord_cells(path, i_indg, j_indg):
    i_indg = list(i_indg)
    j_indg = list(j_indg)
    n_tile = len(i_indg)
    fields = {
        "tile_id": list(range(1, n_tile + 1)),
        "typ": [1] * n_tile,
        "pfaf": [1] * n_tile,
        "com_lon": [0.0] * n_tile,
        "com_lat": [0.0] * n_tile,
        "min_lon": [0.0] * n_tile,
        "max_lon": [1.0] * n_tile,
        "min_lat": [0.0] * n_tile,
        "max_lat": [1.0] * n_tile,
        "i_indg": i_indg,
        "j_indg": j_indg,
        "frac_cell": [1.0] * n_tile,
        "frac_pfaf": [1.0] * n_tile,
        "area": [1.0] * n_tile,
        "elev": [100.0] * n_tile,
    }
    _write_tilecoord_fields(path, fields)


def _write_tilecoord(path):
    fields = {
        "tile_id": [20, 10, 30, 5, 40],
        "typ": [1, 1, 1, 1, 1],
        "pfaf": [1, 1, 1, 1, 1],
        "com_lon": [0.0, 0.0, 0.0, 0.0, 0.0],
        "com_lat": [0.0, 0.0, 0.0, 0.0, 0.0],
        "min_lon": [0.0, 0.0, 0.0, 0.0, 0.0],
        "max_lon": [1.0, 1.0, 1.0, 1.0, 1.0],
        "min_lat": [0.0, 0.0, 0.0, 0.0, 0.0],
        "max_lat": [1.0, 1.0, 1.0, 1.0, 1.0],
        "i_indg": [2, 2, 3, 2, 4],
        "j_indg": [1, 1, 1, 1, 1],
        "frac_cell": [0.5, 0.7, 1.0, 0.7, 1.0],
        "frac_pfaf": [1.0, 1.0, 1.0, 1.0, 1.0],
        "area": [1.0, 1.0, 1.0, 1.0, 1.0],
        "elev": [100.0, 100.0, 100.0, 100.0, 100.0],
    }
    _write_tilecoord_fields(path, fields)


def _write_tilecoord_fields(path, fields):
    int_fields = {"tile_id", "typ", "pfaf", "i_indg", "j_indg"}
    n_tile = len(fields["tile_id"])

    with path.open("wb") as fp:
        _write_record(fp, np.asarray([n_tile], dtype="<i4"))
        for name, values in fields.items():
            dtype = "<i4" if name in int_fields else "<f4"
            _write_record(fp, np.asarray(values, dtype=dtype))


def _write_record(fp, values):
    data = values.tobytes()
    fp.write(struct.pack("<i", len(data)))
    fp.write(data)
    fp.write(struct.pack("<i", len(data)))
