from datetime import date
from pathlib import Path
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.climatology import (  # noqa: E402
    PENTAD_CENTER_DOYS,
    climatology_npz_is_valid,
    climatology_output_path,
    compute_pentad_climatology,
    compute_pentad_climatology_from_paths,
    daily_pair_paths,
    load_pentad_climatology_npz,
    matlab_default_min_count,
    normal_year_doy,
    pentad_index,
    save_pentad_climatology_npz,
)
from iv_tc.generation import save_daily_pair_npz  # noqa: E402
from iv_tc.pairs import DailyPair  # noqa: E402


def test_normal_year_calendar_and_pentads():
    assert normal_year_doy(date(2020, 2, 28)) == 59
    assert normal_year_doy(date(2020, 2, 29)) == 59
    assert normal_year_doy(date(2020, 3, 1)) == 60
    assert normal_year_doy(date(2020, 12, 31)) == 365
    assert pentad_index(date(2020, 1, 1)) == 0
    assert pentad_index(date(2020, 12, 31)) == 72
    np.testing.assert_array_equal(PENTAD_CENTER_DOYS[[0, -1]], np.array([3, 363]))


def test_matlab_default_min_count_uses_end_exclusive_year_span():
    assert matlab_default_min_count(date(2018, 8, 1), date(2024, 6, 30)) == 24
    assert matlab_default_min_count(date(2018, 1, 1), date(2018, 12, 31)) == 4
    assert matlab_default_min_count(date(2018, 10, 1), date(2018, 10, 31)) == 0


def test_climatology_wraps_31_day_window_across_year_boundary():
    pairs = [
        _pair(date(2019, 12, 31), idx=[7], obs=[10.0], model=[20.0]),
        _pair(date(2020, 1, 1), idx=[7], obs=[30.0], model=[40.0]),
    ]

    climatology = compute_pentad_climatology(pairs, min_count=1, grid_size=20)

    np.testing.assert_array_equal(climatology.idx, np.array([7]))
    assert climatology.obs.shape == (1, 73)
    assert climatology.obs[0, 0] == pytest.approx(20.0)
    assert climatology.model[0, 0] == pytest.approx(30.0)
    assert climatology.count[0, 0] == 2
    assert climatology.obs[0, -1] == pytest.approx(20.0)
    assert climatology.count[0, -1] == 2
    assert np.isnan(climatology.obs[0, 20])


def test_climatology_applies_count_threshold_per_sparse_cell():
    pairs = [
        _pair(date(2020, 1, 1), idx=[2, 9], obs=[1.0, 10.0], model=[2.0, 20.0]),
        _pair(date(2021, 1, 1), idx=[2], obs=[3.0], model=[4.0]),
    ]

    climatology = compute_pentad_climatology(pairs, min_count=2, grid_size=20)

    np.testing.assert_array_equal(climatology.idx, np.array([2, 9]))
    assert climatology.count[0, 0] == 2
    assert climatology.obs[0, 0] == pytest.approx(2.0)
    assert climatology.model[0, 0] == pytest.approx(3.0)
    assert climatology.count[1, 0] == 1
    assert np.isnan(climatology.obs[1, 0])
    assert np.isnan(climatology.model[1, 0])


def test_climatology_maps_leap_day_to_february_28():
    pairs = [
        _pair(date(2020, 2, 28), idx=[1], obs=[2.0], model=[3.0]),
        _pair(date(2020, 2, 29), idx=[1], obs=[4.0], model=[5.0]),
    ]

    climatology = compute_pentad_climatology(pairs, min_count=1, grid_size=10)
    column = int(np.where(PENTAD_CENTER_DOYS == 58)[0][0])

    assert climatology.count[0, column] == 2
    assert climatology.obs[0, column] == pytest.approx(3.0)


def test_climatology_save_load_round_trip(tmp_path):
    original = compute_pentad_climatology(
        [_pair(date(2020, 6, 15), idx=[2, 5], obs=[0.2, 0.5], model=[0.3, 0.6])],
        min_count=1,
        grid_size=10,
    )
    path = tmp_path / "climatology.npz"

    save_pentad_climatology_npz(original, path)
    loaded = load_pentad_climatology_npz(path)

    assert climatology_npz_is_valid(path)
    assert loaded.sensor == original.sensor
    assert loaded.run == original.run
    assert loaded.start_date == original.start_date
    assert loaded.end_date == original.end_date
    assert loaded.min_count == original.min_count
    np.testing.assert_array_equal(loaded.idx, original.idx)
    np.testing.assert_array_equal(loaded.count, original.count)
    np.testing.assert_allclose(loaded.obs, original.obs, equal_nan=True)
    np.testing.assert_allclose(loaded.model, original.model, equal_nan=True)


def test_climatology_streams_standard_daily_pair_paths(tmp_path):
    first_day = date(2020, 1, 1)
    second_day = date(2020, 1, 2)
    pair_root = tmp_path / "products"
    for pair in (
        _pair(first_day, idx=[2], obs=[1.0], model=[2.0]),
        _pair(second_day, idx=[2, 5], obs=[3.0, 4.0], model=[4.0, 5.0]),
    ):
        path = pair_root / "step2_pairs" / "smosic" / "OL" / f"{pair.date:%Y%m%d}.npz"
        save_daily_pair_npz(pair, path)

    paths, missing = daily_pair_paths(
        pair_root, "smosic", "OL", first_day, date(2020, 1, 3)
    )
    climatology = compute_pentad_climatology_from_paths(
        paths, min_count=1, grid_size=10
    )

    assert missing == [date(2020, 1, 3)]
    assert climatology.source_days == 2
    np.testing.assert_array_equal(climatology.idx, np.array([2, 5]))
    assert climatology.obs[0, 0] == pytest.approx(2.0)
    assert climatology.obs[1, 0] == pytest.approx(4.0)
    assert climatology_output_path(
        pair_root, "smosic", "OL", first_day, second_day
    ) == (
        pair_root
        / "step3_climatology"
        / "smosic"
        / "OL"
        / "20200101_20200102_w31.npz"
    )


def test_climatology_uses_only_jointly_finite_pair_values():
    pair = _pair(
        date(2020, 1, 1),
        idx=[1, 2, 3],
        obs=[1.0, np.nan, 3.0],
        model=[2.0, 4.0, np.nan],
    )

    climatology = compute_pentad_climatology([pair], min_count=1, grid_size=10)

    np.testing.assert_array_equal(climatology.idx, np.array([1, 2, 3]))
    assert climatology.count[0, 0] == 1
    assert climatology.count[1, 0] == 0
    assert climatology.count[2, 0] == 0
    assert np.isnan(climatology.obs[1, 0])
    assert np.isnan(climatology.model[2, 0])


def test_climatology_rejects_duplicate_dates_and_indices():
    pair = _pair(date(2020, 1, 1), idx=[1], obs=[1.0], model=[2.0])
    with pytest.raises(ValueError, match="Duplicate daily pair date"):
        compute_pentad_climatology([pair, pair], min_count=1, grid_size=10)

    duplicate_idx = _pair(
        date(2020, 1, 2), idx=[1, 1], obs=[1.0, 2.0], model=[3.0, 4.0]
    )
    with pytest.raises(ValueError, match="duplicate grid indices"):
        compute_pentad_climatology([duplicate_idx], min_count=1, grid_size=10)


def _pair(day, idx, obs, model):
    return DailyPair(
        date=day,
        sensor="test_sensor",
        run="OL",
        idx=np.asarray(idx, dtype=np.int64),
        obs=np.asarray(obs, dtype=np.float64),
        model=np.asarray(model, dtype=np.float64),
        obs_units="m3 m-3",
        model_units="m3 m-3",
    )
