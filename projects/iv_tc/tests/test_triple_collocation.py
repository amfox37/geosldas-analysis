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
from iv_tc.pairs import DailyPair  # noqa: E402
from iv_tc.triple_collocation import (  # noqa: E402
    compute_triple_collocation,
    load_triple_collocation_npz,
    save_triple_collocation_npz,
    tc_output_path,
    triple_collocation_npz_is_valid,
)


def test_matlab_tc_formula_parity_for_synthetic_series():
    rng = np.random.default_rng(7)
    n_days = 120
    truth = rng.normal(size=n_days)
    primary_anom = 0.8 * truth + rng.normal(scale=0.35, size=n_days)
    model_anom = 1.1 * truth + rng.normal(scale=0.2, size=n_days)
    secondary_anom = 1.5 * truth + rng.normal(scale=0.5, size=n_days)
    days = [date(2019, 1, 1) + timedelta(days=i) for i in range(n_days)]
    primary_clim = _climatology("primary", [4])
    secondary_clim = _climatology(
        "secondary", [4], obs_units="percent saturation"
    )
    primary_clim.obs[0] = np.linspace(0.1, 0.7, 73)
    secondary_scaled_clim = np.linspace(0.1, 0.4, 73)
    secondary_clim.obs[0] = secondary_scaled_clim * 200.0
    secondary_clim.model[0] = np.linspace(0.2, 0.8, 73)
    primary_pairs = []
    secondary_pairs = []
    for i, day in enumerate(days):
        column = pentad_index(day)
        primary_pairs.append(
            _pair(
                day,
                "primary",
                [4],
                [primary_anom[i] + primary_clim.obs[0, column]],
                [model_anom[i] + secondary_clim.model[0, column]],
            )
        )
        secondary_pairs.append(
            _pair(
                day,
                "secondary",
                [4],
                [(secondary_anom[i] + secondary_scaled_clim[column]) * 200.0],
                [model_anom[i] + secondary_clim.model[0, column]],
                obs_units="percent saturation",
            )
        )

    result = compute_triple_collocation(
        primary_pairs,
        primary_clim,
        secondary_pairs,
        secondary_clim,
        primary_name="SMOS",
        secondary_name="ASCAT",
        secondary_scale=0.005,
        min_count=20,
    )
    expected = _literal_matlab_tc(primary_anom, model_anom, secondary_anom)

    assert result.source_names == ("SMOS", "model", "ASCAT")
    assert result.source_scales == (1.0, 1.0, 0.005)
    assert result.count.tolist() == [n_days]
    np.testing.assert_allclose(result.variance[0], expected[0])
    np.testing.assert_allclose(result.covariance[0], expected[1])
    np.testing.assert_allclose(result.correlation[0], expected[2])
    np.testing.assert_allclose(result.r2[0], expected[3])
    np.testing.assert_allclose(result.error_variance[0], expected[4])


def test_tc_uses_primary_model_values_and_secondary_model_climatology():
    start = date(2020, 1, 1)
    primary_clim = _climatology("primary", [2])
    secondary_clim = _climatology("secondary", [2])
    primary_clim.model[:] = 1000.0
    secondary_clim.model[:] = 10.0
    primary = []
    secondary = []
    model_anom = np.arange(6, dtype=float)
    for i in range(6):
        day = start + timedelta(days=i)
        primary.append(_pair(day, "primary", [2], [i + 1.0], [10.0 + model_anom[i]]))
        secondary.append(
            _pair(day, "secondary", [2], [2.0 * i + 1.0], [500.0 + i])
        )

    result = compute_triple_collocation(
        primary,
        primary_clim,
        secondary,
        secondary_clim,
        primary_name="one",
        secondary_name="three",
        min_count=1,
    )

    assert result.variance[0, 1] == pytest.approx(np.var(model_anom))
    assert result.model_max_abs_difference == pytest.approx(490.0)
    with pytest.raises(ValueError, match="exceeding tolerance"):
        compute_triple_collocation(
            primary,
            primary_clim,
            secondary,
            secondary_clim,
            primary_name="one",
            secondary_name="three",
            min_count=1,
            model_agreement_tolerance=1e-6,
        )


def test_tc_intersects_sparse_inputs_and_requires_all_three_finite():
    day = date(2020, 6, 1)
    primary = _pair(day, "primary", [1, 2, 3], [1.0, 2.0, np.nan], [2, 3, 4])
    secondary = _pair(day, "secondary", [2, 3, 4], [4.0, 5.0, 6.0], [3, 4, 5])

    result = compute_triple_collocation(
        [primary],
        _climatology("primary", [1, 2, 3]),
        [secondary],
        _climatology("secondary", [2, 3, 4]),
        primary_name="one",
        secondary_name="three",
        min_count=1,
    )

    np.testing.assert_array_equal(result.idx, np.array([2]))
    np.testing.assert_array_equal(result.count, np.array([1]))


def test_tc_seasons_missing_days_and_minimum_count():
    days = [date(2020, 1, 1), date(2020, 6, 1), date(2020, 6, 2)]
    primary = [_pair(day, "primary", [1], [i + 1], [i + 2]) for i, day in enumerate(days)]
    secondary = [
        _pair(day, "secondary", [1], [i + 3], [i + 2]) for i, day in enumerate(days)
    ]
    kwargs = dict(
        primary_name="one",
        secondary_name="three",
        start_date=days[0],
        end_date=days[-1],
        min_count=2,
        strict_missing=False,
    )

    summer = compute_triple_collocation(
        primary,
        _climatology("primary", [1]),
        secondary,
        _climatology("secondary", [1]),
        season="summer",
        **kwargs,
    )
    winter = compute_triple_collocation(
        primary,
        _climatology("primary", [1]),
        secondary,
        _climatology("secondary", [1]),
        season="winter",
        **kwargs,
    )

    assert summer.count.tolist() == [2]
    assert np.all(np.isfinite(summer.variance))
    assert winter.count.tolist() == [1]
    assert np.all(np.isnan(winter.variance))
    assert summer.missing_days == 151
    with pytest.raises(FileNotFoundError, match="Missing"):
        compute_triple_collocation(
            primary,
            _climatology("primary", [1]),
            secondary,
            _climatology("secondary", [1]),
            primary_name="one",
            secondary_name="three",
            start_date=days[0],
            end_date=days[-1],
            min_count=1,
        )


def test_tc_npz_round_trip_and_output_path(tmp_path):
    start = date(2020, 1, 1)
    days = [start + timedelta(days=i) for i in range(5)]
    result = compute_triple_collocation(
        [_pair(day, "primary", [2], [i], [2 * i]) for i, day in enumerate(days)],
        _climatology("primary", [2]),
        [_pair(day, "secondary", [2], [3 * i], [2 * i]) for i, day in enumerate(days)],
        _climatology("secondary", [2]),
        primary_name="SMOS IC",
        secondary_name="ASCAT/H119",
        min_count=1,
    )
    path = tmp_path / "tc.npz"
    save_triple_collocation_npz(result, path)
    loaded = load_triple_collocation_npz(path)

    assert triple_collocation_npz_is_valid(path)
    assert loaded.source_names == result.source_names
    np.testing.assert_allclose(loaded.r2, result.r2, equal_nan=True)
    assert tc_output_path(
        tmp_path,
        result.source_names,
        "OL",
        start,
        days[-1],
        min_count=1,
    ) == (
        tmp_path
        / "step4_tc"
        / "SMOS_IC__model__ASCAT_H119"
        / "OL"
        / "20200101_20200105_annual_n1.npz"
    )


def test_compute_tc_cli_end_to_end_with_ascat_default_scale(tmp_path):
    start = date(2020, 1, 1)
    end = date(2020, 1, 5)
    climatology_end = end + timedelta(days=1)
    primary_clim = _climatology("smosic", [2])
    secondary_clim = _climatology(
        "ascat_h119_h120", [2], obs_units="percent saturation"
    )
    for i in range(5):
        day = start + timedelta(days=i)
        for pair in (
            _pair(day, "smosic", [2], [i], [2 * i]),
            _pair(
                day,
                "ascat_h119_h120",
                [2],
                [600 * i],
                [2 * i],
                obs_units="percent saturation",
            ),
        ):
            save_daily_pair_npz(
                pair,
                tmp_path
                / "step2_pairs"
                / pair.sensor
                / "OL"
                / f"{day:%Y%m%d}.npz",
            )
    for climatology in (primary_clim, secondary_clim):
        save_pentad_climatology_npz(
            climatology,
            climatology_output_path(
                tmp_path, climatology.sensor, "OL", start, climatology_end
            ),
        )

    completed = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "compute_tc.py"),
            "--start-date",
            str(start),
            "--end-date",
            str(end),
            "--climatology-end-date",
            str(climatology_end),
            "--primary-sensor",
            "smosic",
            "--secondary-sensor",
            "ascat_h119_h120",
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

    output = tc_output_path(
        tmp_path,
        ("smosic", "model", "ascat_h119_h120"),
        "OL",
        start,
        end,
        min_count=1,
    )
    loaded = load_triple_collocation_npz(output)
    assert loaded.source_scales == (1.0, 1.0, 0.005)
    assert loaded.primary_climatology_end_date == primary_clim.end_date
    assert "paired_days=5/5" in completed.stdout
    assert "Summary: written=1, exists=0, failed=0" in completed.stdout


def test_matlab_tc_file_parity_utility(tmp_path):
    from scipy.io import savemat

    start = date(2020, 1, 1)
    days = [start + timedelta(days=i) for i in range(8)]
    result = compute_triple_collocation(
        [_pair(day, "primary", [2], [i], [2 * i]) for i, day in enumerate(days)],
        _climatology("primary", [2]),
        [_pair(day, "secondary", [2], [3 * i], [2 * i]) for i, day in enumerate(days)],
        _climatology("secondary", [2]),
        primary_name="SMOS",
        secondary_name="ASCAT",
        min_count=1,
    )
    python_path = tmp_path / "tc.npz"
    save_triple_collocation_npz(result, python_path)
    matlab = {"N_sm": _full(result, result.count, fill=0, dtype=np.int32)}
    mappings = {
        "R2_TC_SMOS": result.r2[:, 0],
        "R2_TC_mod": result.r2[:, 1],
        "R2_TC_ASC": result.r2[:, 2],
        "sigma2_SMOS": result.error_variance[:, 0],
        "sigma2_mod": result.error_variance[:, 1],
        "sigma2_ASC": result.error_variance[:, 2],
        "R_mod_SMOS": result.correlation[:, 0],
        "R_ASC_SMOS": result.correlation[:, 1],
        "R_mod_ASC": result.correlation[:, 2],
        "C_SMOS_mod": result.covariance[:, 0],
        "C_SMOS_ASC": result.covariance[:, 1],
        "C_mod_ASC": result.covariance[:, 2],
    }
    for field, values in mappings.items():
        matlab[field] = _full(result, values)
    matlab_path = tmp_path / "tc.mat"
    savemat(matlab_path, matlab)

    completed = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "compare_tc_matlab.py"),
            "--python",
            str(python_path),
            "--matlab",
            str(matlab_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Parity check passed." in completed.stdout


def _literal_matlab_tc(one, two, three):
    series = np.vstack((one, two, three))
    variance = np.mean(series**2, axis=1) - np.mean(series, axis=1) ** 2
    covariance = np.array(
        [
            np.mean(one * two) - np.mean(one) * np.mean(two),
            np.mean(one * three) - np.mean(one) * np.mean(three),
            np.mean(two * three) - np.mean(two) * np.mean(three),
        ]
    )
    correlation = np.array(
        [
            covariance[0] / np.sqrt(variance[0] * variance[1]),
            covariance[1] / np.sqrt(variance[0] * variance[2]),
            covariance[2] / np.sqrt(variance[1] * variance[2]),
        ]
    )
    r2 = np.array(
        [
            covariance[0] * covariance[1] / covariance[2] / variance[0],
            covariance[0] * covariance[2] / covariance[1] / variance[1],
            covariance[2] * covariance[1] / covariance[0] / variance[2],
        ]
    )
    r2[r2 < 0.001] = np.nan
    r2[r2 > 1] = 1
    error = np.array(
        [
            variance[0] - covariance[1] * covariance[0] / covariance[2],
            variance[1] - covariance[2] * covariance[0] / covariance[1],
            variance[2] - covariance[1] * covariance[2] / covariance[0],
        ]
    )
    return variance, covariance, correlation, r2, error


def _full(result, values, fill=np.nan, dtype=np.float64):
    output = np.full((result.grid_size, 1), fill, dtype=dtype)
    output[result.idx, 0] = values
    return output


def _pair(day, sensor, idx, obs, model, obs_units="m3 m-3"):
    return DailyPair(
        date=day,
        sensor=sensor,
        run="OL",
        idx=np.asarray(idx, dtype=np.int64),
        obs=np.asarray(obs, dtype=np.float64),
        model=np.asarray(model, dtype=np.float64),
        obs_units=obs_units,
        model_units="m3 m-3",
    )


def _climatology(sensor, idx, obs_units="m3 m-3"):
    idx = np.asarray(idx, dtype=np.int64)
    shape = (idx.size, PENTAD_CENTER_DOYS.size)
    return PentadClimatology(
        idx=idx,
        obs=np.zeros(shape, dtype=np.float64),
        model=np.zeros(shape, dtype=np.float64),
        count=np.ones(shape, dtype=np.int64),
        center_doy=PENTAD_CENTER_DOYS.copy(),
        sensor=sensor,
        run="OL",
        obs_units=obs_units,
        model_units="m3 m-3",
        start_date=date(2018, 1, 1),
        end_date=date(2024, 12, 31),
        window_days=31,
        min_count=1,
        source_days=1,
        grid_size=20,
    )
