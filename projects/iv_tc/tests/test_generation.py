from datetime import date
from pathlib import Path
import struct
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

netCDF4 = pytest.importorskip("netCDF4")

from iv_tc.config import ProductRoots, RunConfig  # noqa: E402
from iv_tc.generation import (  # noqa: E402
    canonical_sensor,
    generate_daily_pairs,
    load_daily_pair_npz,
    pair_npz_is_valid,
    pair_output_path,
    save_daily_pair_npz,
)
from iv_tc.pairs import DailyPair  # noqa: E402


def test_save_load_and_validate_daily_pair_npz(tmp_path):
    pair = DailyPair(
        date=date(2020, 1, 2),
        sensor="SMOS-IC",
        run="OL",
        idx=np.array([12, 13], dtype=np.int64),
        obs=np.array([0.21, 0.32], dtype=np.float64),
        model=np.array([0.14, 0.33], dtype=np.float64),
        obs_units="m3 m-3",
        model_units="m3 m-3",
    )
    path = tmp_path / "pair.npz"

    save_daily_pair_npz(pair, path)

    assert pair_npz_is_valid(path)
    with np.load(path) as z:
        assert set(("idx0", "sm_obs", "sm_mod")).issubset(z.files)
        assert z["idx0"].dtype == np.int32
        assert z["sm_obs"].dtype == np.float32
        assert z["sm_mod"].dtype == np.float32

    loaded = load_daily_pair_npz(path)
    assert loaded.date == pair.date
    assert loaded.sensor == pair.sensor
    assert loaded.run == pair.run
    np.testing.assert_array_equal(loaded.idx, pair.idx)
    np.testing.assert_allclose(loaded.obs, pair.obs, rtol=1e-6)
    np.testing.assert_allclose(loaded.model, pair.model, rtol=1e-6)


def test_generate_daily_pairs_writes_smosic_and_skips_existing(tmp_path):
    day = date(2020, 1, 2)
    roots, run = _write_smosic_fixture_tree(tmp_path, day)
    output_root = tmp_path / "out"

    results = generate_daily_pairs(
        dates=[day],
        sensors=["smosic"],
        runs=[run],
        roots=roots,
        output_root=output_root,
        nx=10,
    )

    assert len(results) == 1
    result = results[0]
    assert result.status == "written"
    assert result.n_points == 2
    assert result.output_path == pair_output_path(output_root, "smosic", "OL", day)
    assert pair_npz_is_valid(result.output_path)

    loaded = load_daily_pair_npz(result.output_path)
    np.testing.assert_array_equal(loaded.idx, np.array([12, 13]))
    np.testing.assert_allclose(loaded.obs, np.array([0.21, 0.32]), rtol=1e-6)
    np.testing.assert_allclose(loaded.model, np.array([0.14, 0.33]), rtol=1e-6)

    again = generate_daily_pairs(
        dates=[day],
        sensors=["smosic"],
        runs=[run],
        roots=roots,
        output_root=output_root,
        nx=10,
    )
    assert again[0].status == "exists"


def test_generate_daily_pairs_accepts_label_run_root_with_single_fixture_file(tmp_path):
    day = date(2020, 1, 2)
    roots, run = _write_smosic_fixture_tree(tmp_path, day, run_dir_name="OL", exp_id="EXP_A")

    results = generate_daily_pairs(
        dates=[day],
        sensors=["smosic"],
        runs=[run],
        roots=roots,
        output_root=tmp_path / "out",
        nx=10,
    )

    assert results[0].status == "written"
    assert results[0].n_points == 2


def test_generate_daily_pairs_reports_missing_observation(tmp_path):
    day = date(2020, 1, 2)
    roots, run = _write_smosic_fixture_tree(tmp_path, day)
    roots = ProductRoots(
        ascat_h121_root=tmp_path / "ASCAT_SSM_CDR",
        ascat_h119_h120_root=tmp_path / "ASCAT_HSAF",
        smosic_root=roots.smosic_root,
        smap_l3_root=tmp_path / "SPL3SMP_v009",
    )

    results = generate_daily_pairs(
        dates=[day],
        sensors=["smap_l3"],
        runs=[run],
        roots=roots,
        output_root=tmp_path / "out",
        nx=10,
    )

    assert results[0].status == "missing"
    assert "SMAP_L3_SM_P_20200102" in results[0].message


def test_ascat_alias_requires_specific_product():
    with pytest.raises(ValueError, match="ambiguous"):
        canonical_sensor("ascat")


def _write_smosic_fixture_tree(tmp_path, day, run_dir_name="RUN_A", exp_id="RUN_A"):
    roots = ProductRoots(
        ascat_h121_root=tmp_path / "ASCAT_SSM_CDR",
        ascat_h119_h120_root=tmp_path / "ASCAT_HSAF",
        smosic_root=tmp_path / "SMOS_IC",
        smap_l3_root=tmp_path / "SPL3SMP_v009",
    )
    run_root = tmp_path / "runs" / run_dir_name
    run = RunConfig("OL", run_root)

    _write_smosic_sparse(
        roots.smosic_root / f"smos_ic_sm_m36_{day:%Y%m%d}.nc",
        idx=[12, 13, 99],
        values=[0.21, 0.32, 0.99],
    )
    _write_tilecoord(
        run_root
        / "output"
        / roots.domain
        / "rc_out"
        / f"{exp_id}.ldas_tilecoord.bin"
    )
    _write_model(
        run_root
        / "output"
        / roots.domain
        / "cat"
        / "ens_avg"
        / f"Y{day.year:04d}"
        / f"M{day.month:02d}"
        / f"{exp_id}{roots.collection}{day:%Y%m%d}_1200z.nc4",
        values=[0.01, 0.02, 0.33, 0.14, 0.55],
    )
    return roots, run


def _write_smosic_sparse(path, idx, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("n_points", len(idx))
        idx_var = ds.createVariable("idx_EASEv2_lonxlat", "i4", ("n_points",))
        sm_var = ds.createVariable("sm_obs", "f4", ("n_points",))
        idx_var[:] = np.asarray(idx, dtype=np.int32)
        sm_var[:] = np.asarray(values, dtype=np.float32)
        sm_var.units = "m3 m-3"
        ds.date = "20200102"


def _write_model(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("tile", len(values))
        var = ds.createVariable("SFMC", "f4", ("time", "tile"))
        var.units = "m3 m-3"
        var[0, :] = np.asarray(values, dtype=np.float32)


def _write_tilecoord(path):
    path.parent.mkdir(parents=True, exist_ok=True)
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
