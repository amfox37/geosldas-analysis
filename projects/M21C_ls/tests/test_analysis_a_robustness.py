from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analysis_a_robustness import (  # noqa: E402
    block_sufficient_statistics,
    restricted_sample,
    seasonal_aggregate,
    variation_diagnostics,
)


def test_seasonal_sum_preserves_original_partial_nan_rule():
    time = pd.date_range("2001-03-01", periods=3, freq="MS")
    data = xr.DataArray(
        [[1.0, np.nan], [2.0, np.nan], [np.nan, np.nan]],
        dims=("time", "tile"),
        coords={"time": time, "tile": [0, 1]},
    )
    result = seasonal_aggregate(data, [2001], "MAM", "sum")
    assert result.sel(year=2001, tile=0).item() == 3.0
    assert np.isnan(result.sel(year=2001, tile=1).item())


def test_restricted_sample_uses_one_four_year_complete_case_population():
    rows = []
    for tile, years in {1: range(2001, 2007), 2: range(2001, 2004)}.items():
        for year in years:
            rows.append(
                {
                    "tile": tile,
                    "year": year,
                    "lat": 45.0,
                    "lon": -100.0,
                    "snow_abs_netpack_mam": float(year - 2000),
                    "snow_net_mam": float(year - 2000),
                    "ol_snomasland_mam": 20.0 + year,
                    "et": 2.0 * (year - 2000),
                }
            )
    sample = restricted_sample(pd.DataFrame(rows), "et", minimum_years=4)
    assert sample["tile"].unique().tolist() == [1]
    assert len(sample) == 6
    assert np.isclose(sample.groupby("tile")["snow_net_mam_anom"].mean().item(), 0.0)
    assert np.isclose(sample.groupby("tile")["response_signed_anom"].mean().item(), 0.0)


def test_vectorized_block_sufficient_statistics_match_direct_products():
    x = np.array([[1.0, 2.0], [1.0, 3.0], [1.0, 4.0], [1.0, 8.0]])
    y = np.array([2.0, 4.0, 5.0, 11.0])
    codes = np.array([0, 0, 1, 1])
    xtx, xty = block_sufficient_statistics(x, y, codes)
    for block in (0, 1):
        keep = codes == block
        assert np.allclose(xtx[block], x[keep].T @ x[keep])
        assert np.allclose(xty[block], x[keep].T @ y[keep])


def test_variance_decomposition_recovers_total_and_flags_low_variation():
    table = pd.DataFrame(
        {
            "tile": np.repeat([0, 1], 4),
            "year": list(range(2001, 2005)) * 2,
            "snow_abs_netpack_mam": [1.0, 1.0, 1.0, 1.0, 4.0, 5.0, 6.0, 7.0],
        }
    )
    result = variation_diagnostics(
        table,
        "snow_abs_netpack_mam",
        minimum_years=4,
        minimum_sd=0.1,
    )
    assert result["variance_identity_error"] < 1.0e-12
    assert result["fraction_below_minimum_sd"] == 0.5
    assert np.isclose(
        result["within_variance_fraction"] + result["between_variance_fraction"],
        1.0,
    )
