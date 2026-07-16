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


def _write_model(path, values):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("tile", len(values))
        var = ds.createVariable("SFMC", "f4", ("time", "tile"))
        var.units = "m3 m-3"
        var[0, :] = np.asarray(values, dtype=np.float32)


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
    int_fields = {"tile_id", "typ", "pfaf", "i_indg", "j_indg"}

    with path.open("wb") as fp:
        _write_record(fp, np.asarray([5], dtype="<i4"))
        for name, values in fields.items():
            dtype = "<i4" if name in int_fields else "<f4"
            _write_record(fp, np.asarray(values, dtype=dtype))


def _write_record(fp, values):
    data = values.tobytes()
    fp.write(struct.pack("<i", len(data)))
    fp.write(data)
    fp.write(struct.pack("<i", len(data)))
