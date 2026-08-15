from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from trend_breakpoint_series import MonthlySeriesLoader, load_variable_selection


TIME = pd.date_range("2020-06-01", "2020-08-01", freq="MS")
LAT = np.array([30.0, 40.0, 0.0])
LON = np.array([-110.0, -100.0, 10.0])
AREA = np.array([1.0, 3.0, 2.0])


def _land_dataset(experiment: str) -> xr.Dataset:
    is_da = experiment == "da"
    sfmc = np.array(
        [
            [2.0, 11.0, 7.0],
            [4.0, np.nan, 8.0],
            [6.0, 33.0, 9.0],
        ]
        if is_da
        else [
            [1.0, 10.0, 5.0],
            [2.0, 20.0, 6.0],
            [3.0, 30.0, np.nan],
        ],
        dtype="float32",
    )
    scf = np.array(
        [[0.10, 0.30, 0.0], [0.0, 0.30, 0.0], [0.0, 0.30, 0.0]],
        dtype="float32",
    )
    dataset = xr.Dataset(
        {
            "SFMC": (("time", "tile"), sfmc, {"units": "m3 m-3"}),
            "PRECTOTCORRLAND": (
                ("time", "tile"),
                np.full((3, 3), 2.0e-6 if is_da else 1.0e-6, dtype="float32"),
                {"units": "kg m-2 s-1"},
            ),
            "FRLANDSNO": (("time", "tile"), scf, {"units": "1"}),
            "SNOMASLAND": (
                ("time", "tile"),
                np.zeros((3, 3), dtype="float32"),
                {"units": "kg m-2"},
            ),
            "TSOIL1": (
                ("time", "tile"),
                np.tile(np.array([275.0, 270.0, 280.0], dtype="float32"), (3, 1)),
                {"units": "K"},
            ),
        },
        coords={"time": TIME, "lat": ("tile", LAT), "lon": ("tile", LON)},
        attrs={"source_root": f"/synthetic/{experiment.upper()}_SOURCE"},
    )
    return dataset


def _water_dataset(experiment: str) -> xr.Dataset:
    offset = 1.0e-6 if experiment == "da" else 0.0
    surface = np.full((3, 3), 2.0e-6 + offset, dtype="float32")
    baseflow = np.full((3, 3), 1.0e-6 + offset, dtype="float32")
    if experiment == "da":
        baseflow[1, 1] = np.nan
    return xr.Dataset(
        {
            "RUNSURFLAND": (
                ("time", "tile"), surface, {"units": "kg m-2 s-1"}
            ),
            "BASEFLOWLAND": (
                ("time", "tile"), baseflow, {"units": "kg m-2 s-1"}
            ),
        },
        coords={"time": TIME},
        attrs={"source_root": f"/synthetic/{experiment.upper()}_SOURCE"},
    )


@pytest.fixture()
def synthetic_inputs(tmp_path: Path) -> dict[str, Path]:
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    _land_dataset("ol").to_netcdf(data_dir / "ol_land.nc", engine="scipy")
    _land_dataset("da").to_netcdf(data_dir / "da_land.nc", engine="scipy")
    _water_dataset("ol").to_netcdf(data_dir / "ol_water.nc", engine="scipy")
    _water_dataset("da").to_netcdf(data_dir / "da_water.nc", engine="scipy")

    catch = xr.Dataset(
        {
            "snow_net": (
                ("time", "tile"),
                np.array([[-1.0, 2.0, 3.0], [4.0, -5.0, 6.0], [7.0, 8.0, -9.0]]),
                {"units": "kg m-2"},
            )
        },
        coords={"time": TIME, "lat": ("tile", LAT + 0.25)},
        attrs={"source_root": "/synthetic/DA_SOURCE"},
    )
    catch.to_netcdf(data_dir / "catch.nc", engine="scipy")

    inst3 = xr.Dataset(
        {
            "RZMC_INC_MEAN": (
                ("time", "tile"),
                np.array([[0.1, -0.2, 0.3], [0.4, -0.5, 0.6], [0.7, -0.8, 0.9]]),
                {"units": "m3 m-3"},
            )
        },
        coords={"time": TIME, "lat": ("tile", LAT - 0.25)},
        attrs={"source_root": "/synthetic/DA_SOURCE"},
    )
    inst3.to_netcdf(data_dir / "inst3.nc", engine="scipy")

    contract = {
        "schema_version": 1,
        "record_start": "2020-06-01",
        "record_end": "2020-08-01",
        "n_months": 3,
        "n_tiles": 3,
        "tilecoord_file": "unused.bin",
        "tile_area_units": "km2",
        "datasets": {
            "ol_land": {
                "filename": "ol_land.nc",
                "experiment": "ol",
                "family": "land",
                "source_token": "OL_SOURCE",
            },
            "da_land": {
                "filename": "da_land.nc",
                "experiment": "da",
                "family": "land",
                "source_token": "DA_SOURCE",
            },
            "catch_raw_cumulative": {
                "filename": "catch.nc",
                "experiment": "da",
                "family": "catch_raw_cumulative",
                "source_token": "DA_SOURCE",
            },
            "inst3_raw_diagnostic": {
                "filename": "inst3.nc",
                "experiment": "da",
                "family": "inst3_raw_diagnostic",
                "source_token": "DA_SOURCE",
            },
            "ol_water_budget": {
                "filename": "ol_water.nc",
                "experiment": "ol",
                "family": "water_budget",
                "source_token": "OL_SOURCE",
            },
            "da_water_budget": {
                "filename": "da_water.nc",
                "experiment": "da",
                "family": "water_budget",
                "source_token": "DA_SOURCE",
            },
        },
    }
    contract_path = tmp_path / "inputs.json"
    contract_path.write_text(json.dumps(contract))

    columns = [
        "analysis_variable",
        "ol_dataset",
        "da_dataset",
        "source_variable",
        "kind",
        "phase",
        "primary_estimand",
        "notes",
    ]
    selection = {
        "schema_version": 1,
        "derived_variables": {
            "TOTAL_RUNOFF": {
                "operation": "sum",
                "source_variables": ["RUNSURFLAND", "BASEFLOWLAND"],
                "long_name": "Total runoff",
            }
        },
        "columns": columns,
        "variables": [
            ["SFMC", "ol_land", "da_land", "SFMC", "paired_state", 1, "trend_of_DA_minus_OL", ""],
            ["PRECTOTCORRLAND", "ol_land", "da_land", "PRECTOTCORRLAND", "paired_flux", 1, "trend_of_DA_minus_OL", ""],
            ["snow_net", "", "catch_raw_cumulative", "snow_net", "cumulative_prognostic_increment", 1, "signed_trend", ""],
            ["RZMC_INC_MEAN", "", "inst3_raw_diagnostic", "RZMC_INC_MEAN", "diagnostic_state_correction", 1, "signed_trend", ""],
            ["TOTAL_RUNOFF", "ol_water_budget", "da_water_budget", "TOTAL_RUNOFF", "paired_flux", 2, "trend_of_DA_minus_OL", ""],
        ],
    }
    selection_path = tmp_path / "selection.json"
    selection_path.write_text(json.dumps(selection))

    return {
        "data_dir": data_dir,
        "contract": contract_path,
        "selection": selection_path,
    }


def _loader(paths: dict[str, Path], area: np.ndarray = AREA) -> MonthlySeriesLoader:
    return MonthlySeriesLoader(
        paths["data_dir"],
        input_contract=paths["contract"],
        variable_selection=paths["selection"],
        tile_area=area,
    )


def test_paired_delta_uses_identical_finite_sample(synthetic_inputs: dict[str, Path]) -> None:
    with _loader(synthetic_inputs) as loader:
        fields = loader.load_tile_series("SFMC").load()

    np.testing.assert_allclose(fields["delta"].isel(time=0), [1.0, 1.0, 2.0])
    assert np.isnan(fields["ol"].isel(time=1, tile=1))
    assert np.isnan(fields["da"].isel(time=1, tile=1))
    assert np.isnan(fields["ol"].isel(time=2, tile=2))
    assert np.isnan(fields["da"].isel(time=2, tile=2))
    assert fields["eligible"].all()
    assert fields.attrs["primary_estimand"] == "trend_of_DA_minus_OL"


def test_area_weighted_domain_mean_reports_support(synthetic_inputs: dict[str, Path]) -> None:
    with _loader(synthetic_inputs) as loader:
        series = loader.domain_mean("SFMC").load()

    np.testing.assert_allclose(series["delta"], [8.0 / 6.0, 2.0, 3.0])
    np.testing.assert_allclose(series["delta"], series["da"] - series["ol"])
    np.testing.assert_array_equal(series["n_tiles"], [3, 2, 2])
    np.testing.assert_allclose(series["area_sum"], [6.0, 3.0, 4.0])
    assert series["area_sum"].attrs["units"] == "km2"
    assert series.attrs["support_reference"] == "delta"


def test_water_rate_is_converted_to_monthly_total(synthetic_inputs: dict[str, Path]) -> None:
    with _loader(synthetic_inputs) as loader:
        fields = loader.load_tile_series("PRECTOTCORRLAND").load()

    june_seconds = 30 * 86400.0
    july_seconds = 31 * 86400.0
    assert fields["ol"].attrs["units"] == "kg m-2 month-1"
    assert fields["delta"].attrs["monthly_transform"].startswith("monthly mean rate")
    np.testing.assert_allclose(fields["delta"].isel(time=0), 1.0e-6 * june_seconds)
    np.testing.assert_allclose(fields["delta"].isel(time=1), 1.0e-6 * july_seconds)


def test_total_runoff_sums_components_before_strict_pairing(
    synthetic_inputs: dict[str, Path],
) -> None:
    with _loader(synthetic_inputs) as loader:
        fields = loader.load_tile_series("TOTAL_RUNOFF").load()

    june_seconds = 30 * 86400.0
    np.testing.assert_allclose(fields["ol"].isel(time=0), 3.0e-6 * june_seconds)
    np.testing.assert_allclose(fields["da"].isel(time=0), 5.0e-6 * june_seconds)
    np.testing.assert_allclose(fields["delta"].isel(time=0), 2.0e-6 * june_seconds)
    assert np.isnan(fields["ol"].isel(time=1, tile=1))
    assert np.isnan(fields["da"].isel(time=1, tile=1))
    assert fields["delta"].attrs["units"] == "kg m-2 month-1"
    assert fields["delta"].attrs["source_variable"] == "TOTAL_RUNOFF"


def test_existing_increment_and_diagnostic_monthly_values_are_not_reaggregated(
    synthetic_inputs: dict[str, Path],
) -> None:
    with _loader(synthetic_inputs) as loader:
        snow = loader.load_tile_series("snow_net").load()
        diagnostic = loader.load_tile_series("RZMC_INC_MEAN").load()

    np.testing.assert_allclose(snow["value"].isel(time=0), [-1.0, 2.0, 3.0])
    np.testing.assert_allclose(diagnostic["value"].isel(time=0), [0.1, -0.2, 0.3])
    assert snow["value"].attrs["monthly_transform"].startswith("none")
    assert diagnostic["value"].attrs["monthly_transform"].startswith("none")
    assert snow.attrs["kind"] == "cumulative_prognostic_increment"
    assert diagnostic.attrs["kind"] == "diagnostic_state_correction"
    np.testing.assert_allclose(snow["lat"], LAT)
    np.testing.assert_allclose(diagnostic["lat"], LAT)


def test_static_and_dynamic_masks_match_synthesis_definitions(
    synthetic_inputs: dict[str, Path],
) -> None:
    with _loader(synthetic_inputs) as loader:
        masks = loader.masks().load()

    np.testing.assert_array_equal(masks["snow_possible"], [True, True, False])
    np.testing.assert_array_equal(masks["seasonal_snow"], [True, False, False])
    np.testing.assert_array_equal(masks["warm_static"], [False, False, True])
    np.testing.assert_array_equal(
        masks["locally_snowy_monthly"],
        [[True, False, False], [False, False, False], [False, False, False]],
    )
    np.testing.assert_array_equal(
        masks["warm_snowfree_monthly"],
        [[False, False, True], [False, False, True], [False, False, True]],
    )

    with _loader(synthetic_inputs) as loader:
        fields = loader.load_tile_series("SFMC", mask="warm_snowfree_monthly").load()
    xr.testing.assert_equal(fields["eligible"], masks["warm_snowfree_monthly"])


def test_dynamic_mask_can_remove_all_support_for_a_month(
    synthetic_inputs: dict[str, Path],
) -> None:
    with _loader(synthetic_inputs) as loader:
        series = loader.domain_mean("SFMC", mask="warm_snowfree_monthly").load()

    np.testing.assert_allclose(series["delta"].isel(time=slice(0, 2)), [2.0, 2.0])
    assert np.isnan(series["delta"].isel(time=2))
    np.testing.assert_array_equal(series["n_tiles"], [1, 1, 0])


def test_invalid_area_weights_fail_loudly(synthetic_inputs: dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="no finite positive weights"):
        _loader(synthetic_inputs, area=np.array([np.nan, 0.0, -1.0]))


def test_contract_rejects_wrong_monthly_timeline(synthetic_inputs: dict[str, Path]) -> None:
    contract = json.loads(synthetic_inputs["contract"].read_text())
    contract["record_end"] = "2020-09-01"
    bad_contract = synthetic_inputs["contract"].with_name("bad_inputs.json")
    bad_contract.write_text(json.dumps(contract))

    with pytest.raises(ValueError, match="month count"):
        MonthlySeriesLoader(
            synthetic_inputs["data_dir"],
            input_contract=bad_contract,
            variable_selection=synthetic_inputs["selection"],
            tile_area=AREA,
        )


def test_contract_rejects_wrong_source_experiment(synthetic_inputs: dict[str, Path]) -> None:
    contract = json.loads(synthetic_inputs["contract"].read_text())
    contract["datasets"]["da_land"]["source_token"] = "WRONG_DA_VERSION"
    bad_contract = synthetic_inputs["contract"].with_name("bad_source_inputs.json")
    bad_contract.write_text(json.dumps(contract))

    with pytest.raises(ValueError, match="source_root does not contain"):
        MonthlySeriesLoader(
            synthetic_inputs["data_dir"],
            input_contract=bad_contract,
            variable_selection=synthetic_inputs["selection"],
            tile_area=AREA,
        )


def test_generic_legacy_increment_suffix_is_rejected(
    synthetic_inputs: dict[str, Path],
) -> None:
    selection = json.loads(synthetic_inputs["selection"].read_text())
    selection["variables"][0][3] = "CATDEF_INCR"
    bad_selection = synthetic_inputs["selection"].with_name("bad_selection.json")
    bad_selection.write_text(json.dumps(selection))

    with pytest.raises(ValueError, match="CATDEF_INCR"):
        load_variable_selection(bad_selection)


def test_contract_rejects_implausible_tile_area_sum(
    synthetic_inputs: dict[str, Path],
) -> None:
    contract = json.loads(synthetic_inputs["contract"].read_text())
    contract["tile_area_sum_min"] = 10.0
    contract["tile_area_sum_max"] = 20.0
    bad_contract = synthetic_inputs["contract"].with_name("bad_area_inputs.json")
    bad_contract.write_text(json.dumps(contract))

    with pytest.raises(ValueError, match="below contracted minimum"):
        MonthlySeriesLoader(
            synthetic_inputs["data_dir"],
            input_contract=bad_contract,
            variable_selection=synthetic_inputs["selection"],
            tile_area=AREA,
        )
