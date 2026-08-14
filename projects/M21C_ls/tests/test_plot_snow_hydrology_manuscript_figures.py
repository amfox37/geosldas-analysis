from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from plot_snow_hydrology_manuscript_figures import (  # noqa: E402
    ANNUAL_TERMS,
    EXPECTED_ANNUAL,
    positive_partition_row,
    validate_annual_budget,
    validate_monthly_pathway,
    validate_modis_only_window,
    validate_octmar,
)


def test_annual_budget_validation_accepts_six_closing_years():
    exact = {
        2001: [32.567907, 19.605720, 12.601400, 2.821516, -2.460728],
        2002: [43.456407, 24.140219, 15.306204, 4.294713, -0.284729],
        2003: [64.504073, 40.982767, 24.152522, 3.117698, -3.748914],
        2004: [67.500405, 45.435482, 23.274197, 2.353629, -3.562903],
        2005: [72.968653, 51.142062, 24.527156, 0.233679, -2.934244],
        2006: [68.426690, 47.572181, 24.131001, -0.230220, -3.046272],
    }
    rows = [
        {"water_year": year, **dict(zip(ANNUAL_TERMS, values))}
        for year, values in exact.items()
    ]
    result = validate_annual_budget(pd.DataFrame(rows))
    assert result["water_year"].tolist() == list(range(2001, 2007))


def test_positive_partition_validation_requires_addition_sample_and_closure():
    values = {
        "I_snow": 71.320103,
        "dRunoff_total": 45.861031,
        "dRunoff_surface": 30.734008,
        "dBaseflow": 15.127023,
        "dET": 25.630015,
        "dStorage": 2.786458,
        "residual": -2.957401,
    }
    fractions = {f"fraction_{name}": value / values["I_snow"] for name, value in values.items() if name != "I_snow"}
    canonical = pd.DataFrame(
        [{"sample": "addition", "n_tile_years": 247545, **values, **fractions}]
    )
    boundary = pd.DataFrame(
        [
            {
                "boundary": "Oct-Sep",
                "scope": "positive_input_partition",
                "n_tile_years": 247545,
                **values,
                **fractions,
                "fraction_dRunoff_total_ci_low_5deg": 0.611077,
                "fraction_dRunoff_total_ci_high_5deg": 0.672333,
            }
        ]
    )
    row = positive_partition_row(canonical, boundary)
    assert np.isclose(row["fraction_dRunoff_total"], 0.643031, atol=1.0e-6)


def test_octmar_validation_retains_nonoverlap_and_infiltration_null():
    responses = ["snowmelt", "infiltration", "rzmc", "et", "total_runoff", "total_water"]
    table = pd.DataFrame(
        {
            "response": responses,
            "m3_march_snow_ci_low_5deg": [0.1, -0.01, 0.1, 0.1, 0.1, 0.1],
            "m3_march_snow_ci_high_5deg": [0.2, 0.02, 0.2, 0.2, 0.2, 0.2],
        }
    )
    validate_octmar(table)


def test_monthly_pathway_requires_finite_et_values():
    monthly = pd.DataFrame(
        {
            "water_month": range(1, 13),
            "month": ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"],
            "snow_net_monthly": np.arange(1.0, 13.0),
            "dSnowmelt_monthly": np.arange(1.0, 13.0),
            "dSFMC_monthly": np.full(12, 0.01),
            "dRZMC_monthly": np.full(12, 0.01),
            "dRunoff_total_monthly": np.arange(1.0, 13.0),
            "dET_monthly": np.arange(1.0, 13.0),
        }
    )
    partition = pd.DataFrame([{"sample": "addition", "n_tile_years": 247545}])
    soil_summary = pd.DataFrame(
        [
            {"metric": "peak_dRZMC", "n": 247545, "area_weighted_mean": 0.0189},
            {"metric": "mjj_mean_dRZMC", "n": 247545, "area_weighted_mean": 0.0108},
        ]
    )
    peak_timing = pd.DataFrame(
        [{"state": "RZMC", "month": "May", "area_weighted_fraction": 1.0}]
    )
    september = pd.DataFrame(
        [
            {
                "variable": "signed_snow_input",
                "snow_addition_tile_years_mean_kg_m2": 12.0,
            }
        ]
    )

    validate_monthly_pathway(monthly, partition, soil_summary, peak_timing, september)
    monthly.loc[0, "dET_monthly"] = np.nan
    with pytest.raises(AssertionError, match="non-finite"):
        validate_monthly_pathway(monthly, partition, soil_summary, peak_timing, september)


def test_modis_only_figure_windows_end_before_microwave_sm_da():
    validate_modis_only_window()
