from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from water_year_snow_da_budget import (  # noqa: E402
    BUDGET_TERMS,
    domain_budget_tables,
    monthly_flux_total,
    soil_metrics,
    water_year_dates,
    weighted_mean,
)


def test_water_year_dates_are_exactly_october_through_september():
    dates = water_year_dates(2004)
    assert len(dates) == 12
    assert dates[0] == pd.Timestamp("2003-10-01")
    assert dates[-1] == pd.Timestamp("2004-09-01")
    assert dates.is_unique


def test_monthly_flux_total_uses_actual_month_length():
    data = xr.DataArray(
        np.ones((2, 1)),
        dims=("time", "tile"),
        coords={"time": pd.to_datetime(["2004-02-01", "2004-03-01"]), "tile": [0]},
    )
    result = monthly_flux_total(data)
    assert result[:, 0].values.tolist() == [29 * 86400.0, 31 * 86400.0]


def test_weighted_mean_uses_requested_support_and_area():
    values = np.array([1.0, 3.0, 100.0])
    area = np.array([1.0, 3.0, 10.0])
    keep = np.array([True, True, False])
    assert weighted_mean(values, area, keep) == 2.5


def test_soil_metrics_use_water_year_months_and_measure_persistence():
    rzmc = np.array(
        [[
            [-0.010, -0.010],
            [0.000, -0.009],
            [0.002, -0.008],
            [0.005, -0.007],
            [0.004, -0.006],
            [0.0005, -0.005],
            [0.006, -0.004],
            [0.009, -0.003],
            [0.012, -0.002],
            [0.015, -0.001],
            [0.0004, -0.002],
            [0.0002, -0.003],
        ]],
        dtype="float64",
    )
    sfmc = 2.0 * rzmc
    result = soil_metrics(rzmc, sfmc, threshold=0.001)

    assert result["peak_dRZMC_water_month"][0, 0] == 10
    assert np.isnan(result["peak_dRZMC_water_month"][0, 1])
    assert result["peak_dSFMC_positive"][0, 1] == 0.0
    assert np.isclose(result["mjj_mean_dRZMC"][0, 0], np.mean([0.009, 0.012, 0.015]))
    assert np.isclose(result["amj_mean_dRZMC"][0, 0], np.mean([0.006, 0.009, 0.012]))
    assert result["rzmc_positive_month_count"][0, 0] == 10
    assert result["rzmc_months_from_peak_to_below_threshold"][0, 0] == 1
    assert not result["rzmc_persistence_censored"][0, 0]
    assert result["september_dRZMC"][0, 0] == 0.0002


def test_domain_budget_tables_preserve_closing_identity():
    years = [2001, 2002]
    tiles = [0, 1]
    area = np.array([1.0, 3.0])
    input_values = np.array([[10.0, 20.0], [20.0, 40.0]])
    runoff = 0.5 * input_values
    et = 0.25 * input_values
    storage = 0.10 * input_values
    residual = input_values - runoff - et - storage
    variables = {
        "I_snow": input_values,
        "dRunoff_total": runoff,
        "dET": et,
        "dStorage": storage,
        "residual": residual,
        "dRunoff_surface": 0.4 * input_values,
        "dBaseflow": 0.1 * input_values,
        "dSnowmelt": 0.9 * input_values,
        "dInfiltration": 0.7 * input_values,
        "dStorage_process_tendency": -0.75 * input_values,
        "process_balance_diagnostic": np.zeros_like(input_values),
        "dPrecipitation": np.zeros_like(input_values),
        "snow_abs_netpack": input_values,
    }
    dataset = xr.Dataset(
        {name: (("water_year", "tile"), value) for name, value in variables.items()},
        coords={
            "water_year": years,
            "tile": tiles,
            "area": ("tile", area),
            "lat": ("tile", [45.0, 50.0]),
            "lon": ("tile", [-110.0, -100.0]),
        },
    )
    config = {"minimum_abs_domain_input_for_fraction_kg_m2": 0.1}
    annual, partition, _ = domain_budget_tables(dataset, config)

    annual_only = annual.iloc[:2]
    closure = annual_only[["dRunoff_total", "dET", "dStorage", "residual"]].sum(axis=1)
    assert np.allclose(closure, annual_only["I_snow"])
    assert np.allclose(annual_only["fraction_dRunoff_total"], 0.5)
    assert np.isclose(partition.set_index("sample").loc["all", "fraction_residual"], 0.15)
    assert set(BUDGET_TERMS).issubset(annual.columns)
