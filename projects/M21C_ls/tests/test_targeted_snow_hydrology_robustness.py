from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import targeted_snow_hydrology_robustness as targeted  # noqa: E402


def test_oct_mar_dates_cross_year_and_do_not_overlap_responses():
    predictor = targeted.oct_mar_dates(2001)
    assert predictor.equals(pd.date_range("2000-10-01", "2001-03-01", freq="MS"))
    assert len(predictor.intersection(targeted.response_dates(2001, "AMJ"))) == 0
    assert len(predictor.intersection(targeted.response_dates(2001, "MJJ"))) == 0


def test_date_validation_reports_worked_example_and_six_seasons():
    config = {
        "modis_only_start": "2000-06-01",
        "microwave_soil_moisture_da_start": "2007-06-01",
    }
    message = targeted.validate_date_windows(list(range(2001, 2007)), config)
    assert "2000-10 through 2001-03" in message
    assert "overlap=0 months" in message


def test_monthly_increment_product_is_not_temporally_differenced():
    dataset = xr.Dataset(
        {"snow_net": (("time", "tile"), np.ones((288, 1)))},
        coords={"time": pd.date_range("2000-06-01", periods=288, freq="MS"), "tile": [0]},
        attrs={
            "note": (
                "Monthly cumulative sums from raw/submonthly catch_progn_incr files. "
                "Do not confuse with monthly files, which are time means."
            ),
            "source_root": "DA v3",
        },
    )
    dataset["snow_net"].attrs["units"] = "kg m-2"
    result = targeted.verify_monthly_increment_product(dataset)
    assert result["classification"] == "per-month increment sum"
    assert "do not difference" in result["action"]


def test_generic_models_recover_within_tile_slope():
    rows = []
    predictor_pattern = [0.0, 1.0, -1.0, 2.0, -2.0, 0.5]
    control_pattern = [1.0, -1.0, 0.5, 2.0, 0.0, -2.0]
    for tile in range(3):
        for year in range(2001, 2007):
            year_index = year - 2001
            predictor = float(year - 2000 + tile * predictor_pattern[year_index])
            control = float((year - 2000) ** 2 + tile * control_pattern[year_index])
            response = 2.0 * predictor + 0.5 * control + 10.0 * tile
            rows.append(
                {
                    "tile": tile,
                    "year": year,
                    "lat": 40.0 + tile,
                    "lon": -110.0 + tile,
                    "area": 1.0,
                    "predictor": predictor,
                    "control": control,
                    "response_name": response,
                }
            )
    sample = targeted.prepare_regression_sample(
        pd.DataFrame(rows), "response_name", "predictor", ["control"], 4
    ).rename(columns={"response_name": "response"})
    sample["response_anom"] = sample["response_name_anom"]
    assert np.isclose(targeted.fit_slope(sample, "predictor", "M3", "control"), 2.0)


def synthetic_baseline() -> xr.Dataset:
    years = [2001, 2002]
    tiles = [10, 20]
    data_vars = {}
    for name_index, name in enumerate(targeted.MONTHLY_BUDGET_NAMES):
        values = np.arange(2 * 12 * 2, dtype="float32").reshape(2, 12, 2) + name_index * 100
        data_vars[f"{name}_monthly"] = (("water_year", "water_month", "tile"), values)
    data_vars["dRunoff_total_monthly"] = (
        ("water_year", "water_month", "tile"),
        data_vars["dRunoff_surface_monthly"][1] + data_vars["dBaseflow_monthly"][1],
    )
    return xr.Dataset(
        data_vars,
        coords={
            "water_year": years,
            "water_month": np.arange(1, 13),
            "tile": tiles,
            "lat": ("tile", [45.0, 50.0]),
            "lon": ("tile", [-110.0, -100.0]),
            "area": ("tile", [1.0, 2.0]),
        },
    )


def test_september_august_reconstruction_uses_august_to_august(monkeypatch):
    baseline = synthetic_baseline()

    def fake_snapshot(_inputs, date, _tile_ids):
        month = pd.Timestamp(date).month
        result = {name: np.full(2, month + index * 100, dtype="float32") for index, name in enumerate(targeted.MONTHLY_BUDGET_NAMES)}
        result["dRunoff_total"] = result["dRunoff_surface"] + result["dBaseflow"]
        return result

    monkeypatch.setattr(targeted, "monthly_snapshot", fake_snapshot)
    result = targeted.build_september_august_dataset(baseline, {}, [2001, 2002])
    assert np.all(result["snow_net_monthly"].values[0, 0] == 9 + 2 * 100)
    assert np.all(
        result["snow_net_monthly"].values[1, 0]
        == baseline["snow_net_monthly"].values[0, 11]
    )
    expected_first_storage = (
        baseline["dTWLAND_monthly"].values[0, 10] - (8 + 10 * 100)
    )
    expected_second_storage = (
        baseline["dTWLAND_monthly"].values[1, 10]
        - baseline["dTWLAND_monthly"].values[0, 10]
    )
    assert np.allclose(result["dStorage"].values[0], expected_first_storage)
    assert np.allclose(result["dStorage"].values[1], expected_second_storage)
    closure = result["dRunoff_total"] + result["dET"] + result["dStorage"] + result["residual"]
    assert np.allclose(closure, result["I_snow"])


def test_september_diagnostic_uses_closing_september_and_august_change():
    baseline = synthetic_baseline()
    baseline["I_snow"] = (("water_year", "tile"), np.full((2, 2), 120.0))
    result = targeted.september_diagnostic(baseline).set_index("variable")
    september = baseline["snow_net_monthly"].values[:, 11].reshape(-1)
    weights = np.tile(baseline.area.values, 2)
    expected = np.average(september, weights=weights)
    assert np.isclose(result.loc["signed_snow_input", "all_tile_years_mean_kg_m2"], expected)
    expected_storage = (
        baseline["dTWLAND_monthly"].values[:, 11]
        - baseline["dTWLAND_monthly"].values[:, 10]
    ).reshape(-1)
    assert np.isclose(
        result.loc["august_to_september_dTWLAND", "all_tile_years_mean_kg_m2"],
        np.average(expected_storage, weights=weights),
    )
