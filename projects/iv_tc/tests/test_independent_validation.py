from datetime import date, timedelta
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.climatology import (  # noqa: E402
    PENTAD_CENTER_DOYS,
    PentadClimatology,
    climatology_output_path,
    pentad_index,
    save_pentad_climatology_npz,
)
from iv_tc.generation import save_daily_pair_npz  # noqa: E402
from iv_tc.independent_validation import (  # noqa: E402
    MATLAB_UNDEFINED_CORRELATION,
    compute_independent_validation,
    compute_independent_validation_from_paths,
    independent_validation_npz_is_valid,
    iv_output_path,
    load_independent_validation_npz,
    save_independent_validation_npz,
)
from iv_tc.pairs import DailyPair  # noqa: E402


def test_matlab_step4_formula_parity_for_synthetic_series():
    rng = np.random.default_rng(42)
    n_days = 240
    model_anom = np.zeros(n_days)
    obs_anom = np.zeros(n_days)
    for i in range(1, n_days):
        model_anom[i] = 0.92 * model_anom[i - 1] + rng.normal(scale=0.15)
        obs_anom[i] = (
            0.75 * model_anom[i]
            + 0.15 * obs_anom[i - 1]
            + rng.normal(scale=0.08)
        )
    climatology = _climatology([4])
    climatology.model[0] = np.linspace(0.1, 0.8, PENTAD_CENTER_DOYS.size)
    climatology.obs[0] = np.linspace(0.2, 1.2, PENTAD_CENTER_DOYS.size)
    days = [date(2019, 1, 1) + timedelta(days=i) for i in range(n_days)]
    pairs = [
        _pair(
            day,
            [4],
            [obs_anom[i] + climatology.obs[0, pentad_index(day)]],
            [model_anom[i] + climatology.model[0, pentad_index(day)]],
        )
        for i, day in enumerate(days)
    ]

    result = compute_independent_validation(
        pairs, climatology, lag_days=2, min_count=100
    )
    expected = _literal_matlab_step4(
        model_anom[2:], obs_anom[2:], model_anom[:-2], obs_anom[:-2]
    )

    assert result.count.tolist() == [n_days - 2]
    assert result.r_model_obs[0] == pytest.approx(expected[0])
    assert result.r2_ivd_model[0] == pytest.approx(expected[1])
    assert result.r2_ivd_obs[0] == pytest.approx(expected[2])
    assert result.r2_ivs_model[0] == pytest.approx(expected[3])
    assert result.r2_ivs_obs[0] == pytest.approx(expected[4])


def test_iv_uses_exact_calendar_lag_and_intersects_sparse_indices():
    start = date(2020, 1, 1)
    pairs = [
        _pair(start, [2, 5], [1.0, 2.0], [2.0, 4.0]),
        _pair(start + timedelta(days=2), [2, 7], [2.0, 8.0], [4.0, 9.0]),
        _pair(start + timedelta(days=4), [2, 5], [3.0, 5.0], [6.0, 8.0]),
    ]

    result = compute_independent_validation(
        pairs,
        _climatology([2, 5, 7]),
        lag_days=2,
        min_count=1,
        start_date=start,
        end_date=start + timedelta(days=4),
    )

    np.testing.assert_array_equal(result.idx, np.array([2]))
    np.testing.assert_array_equal(result.count, np.array([2]))
    assert result.candidate_days == 3
    assert result.paired_days == 2
    assert result.contributing_days == 2


def test_iv_matches_matlab_nan_propagation_from_lagged_climatology():
    start = date(2020, 1, 4)
    climatology = _climatology([3])
    climatology.model[0, 0] = np.nan
    pairs = [
        _pair(start, [3], [1.0], [2.0]),
        _pair(start + timedelta(days=2), [3], [2.0], [3.0]),
    ]

    result = compute_independent_validation(
        pairs, climatology, lag_days=2, min_count=1
    )

    assert result.count.tolist() == [1]
    assert result.r_model_obs.tolist() == [MATLAB_UNDEFINED_CORRELATION]
    assert np.isnan(result.r2_ivd_model[0])
    assert np.isnan(result.r2_ivd_obs[0])
    assert np.isnan(result.r2_ivs_model[0])
    assert np.isnan(result.r2_ivs_obs[0])


def test_iv_applies_minimum_count_and_matlab_undefined_correlation_fill():
    start = date(2020, 1, 1)
    pairs = [
        _pair(start + timedelta(days=i), [1], [float(i)], [float(2 * i)])
        for i in range(5)
    ]

    result = compute_independent_validation(
        pairs, _climatology([1]), lag_days=2, min_count=4
    )

    assert result.count.tolist() == [3]
    assert result.r_model_obs.tolist() == [MATLAB_UNDEFINED_CORRELATION]
    assert np.isnan(result.r2_ivd_model[0])


def test_iv_maps_infinite_correlation_from_zero_variance_model_to_matlab_undefined():
    # A model series that is bit-identical across all days has zero variance
    # (c_model_model == 0.0 exactly), but floating-point rounding in the
    # model*obs sum vs. mean*mean product can leave c_model_obs a tiny
    # nonzero value instead of exactly 0.0. That produces +/-inf, not NaN,
    # from c_model_obs / sqrt(c_model_model * c_obs_obs).
    start = date(2020, 1, 1)
    model_value = -2.302
    obs_values = [
        87.24998293,
        870.14484758,
        631.70710824,
        -994.52299966,
        714.80855318,
        -932.82884939,
        459.31089286,
    ]
    pairs = [
        _pair(start + timedelta(days=i), [1], [obs_values[i]], [model_value])
        for i in range(len(obs_values))
    ]

    result = compute_independent_validation(
        pairs, _climatology([1]), lag_days=2, min_count=3
    )

    assert result.count.tolist() == [5]
    assert np.isfinite(result.r_model_obs).all()
    assert result.r_model_obs.tolist() == [MATLAB_UNDEFINED_CORRELATION]


def test_iv_path_streaming_and_npz_round_trip(tmp_path):
    start = date(2020, 1, 1)
    pair_paths = []
    for i in range(6):
        pair = _pair(
            start + timedelta(days=i), [2], [float(i % 2)], [float(2 * (i % 2))]
        )
        path = tmp_path / "pairs" / f"{pair.date:%Y%m%d}.npz"
        save_daily_pair_npz(pair, path)
        pair_paths.append(path)

    result = compute_independent_validation_from_paths(
        pair_paths, _climatology([2]), lag_days=2, min_count=1
    )
    output = tmp_path / "iv.npz"
    save_independent_validation_npz(result, output)
    loaded = load_independent_validation_npz(output)

    assert independent_validation_npz_is_valid(output)
    assert loaded.sensor == result.sensor
    assert loaded.run == result.run
    assert loaded.start_date == result.start_date
    assert loaded.end_date == result.end_date
    assert loaded.lag_days == 2
    np.testing.assert_array_equal(loaded.idx, result.idx)
    np.testing.assert_array_equal(loaded.count, result.count)
    np.testing.assert_allclose(
        loaded.r2_ivd_model, result.r2_ivd_model, equal_nan=True
    )
    assert iv_output_path(
        tmp_path, "smosic", "OL", start, start + timedelta(days=5), 2, 1
    ) == (
        tmp_path / "step4_iv" / "smosic" / "OL" / "20200101_20200106_lag2_n1.npz"
    )


def test_iv_rejects_pair_metadata_that_does_not_match_climatology():
    pair = _pair(date(2020, 1, 1), [1], [1.0], [2.0], run="DA")
    with pytest.raises(ValueError, match="metadata does not match"):
        compute_independent_validation(
            [pair, _pair(date(2020, 1, 3), [1], [2.0], [3.0], run="DA")],
            _climatology([1]),
            min_count=1,
        )


def test_compute_iv_cli_end_to_end(tmp_path):
    start = date(2020, 1, 1)
    end = date(2020, 1, 6)
    for i in range(6):
        pair = _pair(
            start + timedelta(days=i), [2], [float(i % 2)], [float(2 * (i % 2))]
        )
        save_daily_pair_npz(
            pair,
            tmp_path
            / "step2_pairs"
            / "smosic"
            / "OL"
            / f"{pair.date:%Y%m%d}.npz",
        )
    save_pentad_climatology_npz(
        _climatology([2]),
        climatology_output_path(tmp_path, "smosic", "OL", start, end),
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "compute_iv.py"),
            "--start-date",
            str(start),
            "--end-date",
            str(end),
            "--sensor",
            "smosic",
            "--run",
            "OL",
            "--pair-root",
            str(tmp_path),
            "--min-count",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    output = iv_output_path(tmp_path, "smosic", "OL", start, end, 2, 1)
    assert independent_validation_npz_is_valid(output)
    assert "paired_days=4/4" in completed.stdout
    assert "Summary: written=1, exists=0, failed=0" in completed.stdout


def test_matlab_file_parity_utility(tmp_path):
    from scipy.io import savemat

    start = date(2020, 1, 1)
    pairs = [
        _pair(start + timedelta(days=i), [2], [float(i % 2)], [float(2 * (i % 2))])
        for i in range(8)
    ]
    result = compute_independent_validation(
        pairs, _climatology([2]), lag_days=2, min_count=1
    )
    python_path = tmp_path / "iv.npz"
    save_independent_validation_npz(result, python_path)

    matlab = {"N_sm": np.zeros((20, 1), dtype=np.int32)}
    matlab["N_sm"][result.idx, 0] = result.count
    result_fields = {
        "R2_ivd_mod": result.r2_ivd_model,
        "R2_ivd_obs": result.r2_ivd_obs,
        "R2_ivs_mod": result.r2_ivs_model,
        "R2_ivs_obs": result.r2_ivs_obs,
        "R_mod_obs": result.r_model_obs,
    }
    for field, values in result_fields.items():
        fill = MATLAB_UNDEFINED_CORRELATION if field == "R_mod_obs" else np.nan
        matlab[field] = np.full((20, 1), fill, dtype=np.float64)
        matlab[field][result.idx, 0] = values
    matlab_path = tmp_path / "iv.mat"
    savemat(matlab_path, matlab)

    completed = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "compare_iv_matlab.py"),
            "--python",
            str(python_path),
            "--matlab",
            str(matlab_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "N_sm:" in completed.stdout
    assert "Parity check passed." in completed.stdout


def _literal_matlab_step4(model, obs, model_lag, obs_lag):
    """Direct scalar translation of Compute_IVd_IVs.m lines 152-192."""

    n = float(len(model))
    model_mean = np.sum(model) / n
    obs_mean = np.sum(obs) / n
    c_model_model = np.sum(model**2) / n - model_mean**2
    c_obs_obs = np.sum(obs**2) / n - obs_mean**2
    c_model_obs = np.sum(model * obs) / n - model_mean * obs_mean
    r_model_obs = c_model_obs / np.sqrt(c_model_model * c_obs_obs)
    c_model_obs_positive = np.nan if c_model_obs < 0 else c_model_obs
    c_model_model_lag = np.sum(model * model_lag) / n - model_mean**2
    c_obs_obs_lag = np.sum(obs * obs_lag) / n - obs_mean**2
    c_model_model_lag = (
        np.nan if c_model_model_lag < 0 else c_model_model_lag
    )
    c_obs_obs_lag = np.nan if c_obs_obs_lag < 0 else c_obs_obs_lag
    c_model_obs_lag = np.sum(model * obs_lag) / n - model_mean * obs_mean
    scale_ivd = np.sqrt(c_model_model_lag / c_obs_obs_lag)
    scale_ivs_obs = c_model_obs_lag / c_obs_obs_lag
    values = [
        c_model_obs_positive * scale_ivd / c_model_model,
        c_model_obs_positive / c_obs_obs / scale_ivd,
        c_model_obs_positive * scale_ivs_obs / c_model_model,
        c_model_obs_positive / c_obs_obs / scale_ivs_obs,
    ]
    values = [np.nan if value < 0 else min(value, 1.0) for value in values]
    return (r_model_obs, *values)


def _climatology(idx):
    idx = np.asarray(idx, dtype=np.int64)
    shape = (idx.size, PENTAD_CENTER_DOYS.size)
    return PentadClimatology(
        idx=idx,
        obs=np.zeros(shape, dtype=np.float64),
        model=np.zeros(shape, dtype=np.float64),
        count=np.ones(shape, dtype=np.int64),
        center_doy=PENTAD_CENTER_DOYS.copy(),
        sensor="smosic",
        run="OL",
        obs_units="m3 m-3",
        model_units="m3 m-3",
        start_date=date(2019, 1, 1),
        end_date=date(2020, 12, 31),
        window_days=31,
        min_count=1,
        source_days=1,
        grid_size=20,
    )


def _pair(day, idx, obs, model, run="OL"):
    return DailyPair(
        date=day,
        sensor="smosic",
        run=run,
        idx=np.asarray(idx, dtype=np.int64),
        obs=np.asarray(obs, dtype=np.float64),
        model=np.asarray(model, dtype=np.float64),
        obs_units="m3 m-3",
        model_units="m3 m-3",
    )
