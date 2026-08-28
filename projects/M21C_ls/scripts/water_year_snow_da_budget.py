#!/usr/bin/env python3
"""Build the WY2001-WY2006 differential snow-DA water budget.

SFMC and RZMC are retained as diagnostic state responses. They are never
added to the closing equation because TWLAND already contains soil water.

The closing equation carries two terms that monthly TWLAND alone cannot
supply, both sourced from dedicated Discover deliveries:

- dPeatFreeStandingWater: PEATCLSM free-standing surface water. catch_calc_wtotl
  deliberately excludes this store from TWLAND, so on peat tiles
  (POROS >= 0.9) water entering or leaving it leaks out of the budget.
- dStorage: instantaneous 00Z Oct 1 TWLAND reconstructed from
  catch_internal_rst restarts, replacing the September monthly-mean
  endpoint proxy. The proxy is retained as dStorage_monthly_proxy.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from analysis_a_robustness import (
    block_codes,
    block_sufficient_statistics,
    confidence_interval,
)


MONTH_LABELS = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
BUDGET_TERMS = ["I_snow", "dRunoff_total", "dET", "dStorage", "dPeatFreeStandingWater", "residual"]
PARTITION_TERMS = [
    "dRunoff_surface",
    "dBaseflow",
    "dET",
    "dStorage",
    "dPeatFreeStandingWater",
    "residual",
]
# Terms whose partition fractions must sum to one against I_snow.
CLOSING_TERMS = ["dRunoff_total", "dET", "dStorage", "dPeatFreeStandingWater", "residual"]
SOIL_METRICS = [
    "peak_dRZMC_positive",
    "mjj_mean_dRZMC",
    "peak_dSFMC_positive",
    "rzmc_positive_month_count",
    "september_dRZMC",
]


def find_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    raise FileNotFoundError("Could not locate repository root")


def parse_args() -> argparse.Namespace:
    root = find_repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=root / "projects/M21C_ls/config/water_year_snow_da_budget.json",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=None,
        help="Testing override. Production uses the configured 1000 replicates.",
    )
    parser.add_argument(
        "--rebuild-tile-budget",
        action="store_true",
        help="Re-read the compressed source files instead of reusing the archived tile budget.",
    )
    return parser.parse_args()


def load_config(path: Path, bootstrap_override: int | None = None) -> dict:
    config = json.loads(path.read_text())
    if bootstrap_override is not None:
        config["bootstrap"]["replicates"] = int(bootstrap_override)
    elif config["bootstrap"]["replicates"] < 1000:
        raise ValueError("Production requires at least 1000 spatial bootstrap replicates")
    return config


def water_year_dates(water_year: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{water_year - 1}-10-01", f"{water_year}-09-01", freq="MS")


def weighted_mean(values: np.ndarray, weights: np.ndarray, keep: np.ndarray | None = None) -> float:
    values = np.asarray(values, dtype="float64")
    weights = np.asarray(weights, dtype="float64")
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if keep is not None:
        valid &= np.asarray(keep, dtype=bool)
    if not valid.any():
        return np.nan
    return float(np.sum(values[valid] * weights[valid]) / np.sum(weights[valid]))


def sum_monthly_flux(data: xr.DataArray) -> xr.DataArray:
    seconds = data.time.dt.days_in_month.astype("float64") * 86400.0
    return (data * seconds).sum("time", skipna=True)


def monthly_flux_total(data: xr.DataArray) -> xr.DataArray:
    seconds = data.time.dt.days_in_month.astype("float64") * 86400.0
    return data * seconds


def build_static_mask(
    ol_land: xr.Dataset,
    da_land: xr.Dataset,
    config: dict,
) -> xr.DataArray:
    scf = xr.concat([ol_land["FRLANDSNO"], da_land["FRLANDSNO"]], dim="experiment").max("experiment")
    swe = xr.concat([ol_land["SNOMASLAND"], da_land["SNOMASLAND"]], dim="experiment").max("experiment")
    mask_config = config["mask"]
    mask = (
        (ol_land["lat"] > mask_config["minimum_latitude_degrees_north"])
        & (
            (scf.max("time", skipna=True) > mask_config["snow_cover_possible_threshold"])
            | (swe.max("time", skipna=True) > mask_config["snow_mass_possible_threshold_kg_m2"])
        )
        & (
            scf.sel(time=scf.time.dt.month.isin([6, 7, 8])).mean("time", skipna=True)
            < mask_config["maximum_mean_jja_snow_cover_fraction"]
        )
        & np.isfinite(ol_land["lat"])
        & np.isfinite(ol_land["lon"])
    ).load()
    count = int(mask.sum().item())
    if count != mask_config["expected_tile_count"]:
        raise AssertionError(
            f"Seasonal-snow mask changed: expected {mask_config['expected_tile_count']}, found {count}"
        )
    return mask


def soil_metrics(rzmc: np.ndarray, sfmc: np.ndarray, threshold: float) -> dict[str, np.ndarray]:
    """Calculate peak, sustained, and persistence metrics on WY x month x tile arrays."""
    peak_rzmc = np.nanmax(rzmc, axis=1)
    peak_sfmc = np.nanmax(sfmc, axis=1)
    rzmc_argmax = np.nanargmax(np.where(np.isfinite(rzmc), rzmc, -np.inf), axis=1)
    sfmc_argmax = np.nanargmax(np.where(np.isfinite(sfmc), sfmc, -np.inf), axis=1)
    rzmc_month = np.where(peak_rzmc > 0, rzmc_argmax + 1, np.nan)
    sfmc_month = np.where(peak_sfmc > 0, sfmc_argmax + 1, np.nan)
    sfmc_at_rzmc_peak = np.take_along_axis(sfmc, rzmc_argmax[:, None, :], axis=1)[:, 0, :]
    months_to_below = np.full(peak_rzmc.shape, np.nan, dtype="float64")
    censored = np.zeros(peak_rzmc.shape, dtype=bool)
    for year_index in range(rzmc.shape[0]):
        for tile_index in range(rzmc.shape[2]):
            peak_index = rzmc_argmax[year_index, tile_index]
            peak_value = peak_rzmc[year_index, tile_index]
            if not np.isfinite(peak_value):
                continue
            if peak_value < threshold:
                months_to_below[year_index, tile_index] = 0.0
                continue
            later = rzmc[year_index, peak_index + 1 :, tile_index]
            hits = np.flatnonzero(np.isfinite(later) & (later < threshold))
            if hits.size:
                months_to_below[year_index, tile_index] = float(hits[0] + 1)
            else:
                censored[year_index, tile_index] = True
    return {
        "peak_dRZMC": peak_rzmc,
        "peak_dRZMC_positive": np.maximum(peak_rzmc, 0.0),
        "peak_dRZMC_water_month": rzmc_month,
        "peak_dSFMC": peak_sfmc,
        "peak_dSFMC_positive": np.maximum(peak_sfmc, 0.0),
        "peak_dSFMC_water_month": sfmc_month,
        "dSFMC_at_dRZMC_peak": sfmc_at_rzmc_peak,
        "mjj_mean_dRZMC": np.nanmean(rzmc[:, 7:10, :], axis=1),
        "amj_mean_dRZMC": np.nanmean(rzmc[:, 6:9, :], axis=1),
        "rzmc_positive_month_count": np.sum(rzmc > 0, axis=1).astype("float32"),
        "rzmc_months_from_peak_to_below_threshold": months_to_below,
        "rzmc_persistence_censored": censored,
        "september_dRZMC": rzmc[:, 11, :],
    }


def load_budget_dataset(config: dict) -> xr.Dataset:
    root = find_repo_root()
    data_dir = Path(config["data_directory"])
    paths = {
        "ol_land": data_dir / "OLv8_land_variables_2000_2024_compressed.nc",
        "da_land": data_dir / "DAv8_land_variables_2000_2024_compressed.nc",
        "ol_flux": data_dir / "OLv8_flux_core_2000_2024_compressed.nc",
        "da_flux": data_dir / "DAv8_flux_core_2000_2024_compressed.nc",
        "ol_water": data_dir / "OLv8_water_budget_2000_2024_compressed.nc",
        "da_water": data_dir / "DAv8_water_budget_2000_2024_compressed.nc",
        "catch": data_dir / "catch_progn_raw_monthly_cumulative_200006_202405.nc",
        "ol_fsw": data_dir / "OLv8_peat_fsw_2000_2024_compressed.nc",
        "da_fsw": data_dir / "DAv8_peat_fsw_2000_2024_compressed.nc",
        "ol_twland_rst": data_dir / "OLv8_twland_wy_endpoints_2000_2006.nc",
        "da_twland_rst": data_dir / "DAv8_twland_wy_endpoints_2000_2006.nc",
    }
    opened = {name: xr.open_dataset(path, decode_times=True) for name, path in paths.items()}
    try:
        mask = build_static_mask(opened["ol_land"], opened["da_land"], config)
        tile_indices = np.flatnonzero(mask.values)
        tile_ids = np.arange(mask.size, dtype="int32")[tile_indices]
        lat = opened["ol_land"]["lat"].values[tile_indices]
        lon = opened["ol_land"]["lon"].values[tile_indices]
        sys.path.insert(0, str(root / "common/python/io"))
        from read_GEOSldas import read_tilecoord

        tilecoord = read_tilecoord(str(data_dir / "LS_OLv8_M36.ldas_tilecoord.bin"))
        area = np.asarray(tilecoord["area"], dtype="float64")[tile_indices]
        if not np.isfinite(area).all() or np.any(area <= 0):
            raise AssertionError("Seasonal-snow tile areas must be finite and positive")

        n_all_tiles = int(mask.size)
        for name in ["ol_fsw", "da_fsw", "ol_twland_rst", "da_twland_rst"]:
            if opened[name].sizes["tile"] != n_all_tiles:
                raise AssertionError(
                    f"{name} carries {opened[name].sizes['tile']} tiles, expected {n_all_tiles}"
                )
            for axis in ["lat", "lon"]:
                if not np.allclose(
                    opened[name][axis].values[tile_indices],
                    opened["ol_land"][axis].values[tile_indices],
                    atol=1.0e-4,
                ):
                    raise AssertionError(f"{name} {axis} does not match the land-variable tile order")

        restart_time = pd.DatetimeIndex(opened["ol_twland_rst"]["time"].values)
        if not restart_time.equals(pd.DatetimeIndex(opened["da_twland_rst"]["time"].values)):
            raise AssertionError("OL and DA restart endpoint dates differ")
        if not ((restart_time.month == 10) & (restart_time.day == 1)).all():
            raise AssertionError("Restart TWLAND endpoints must all fall on October 1")
        restart_index = {int(stamp.year): position for position, stamp in enumerate(restart_time)}
        restart_delta = (
            opened["da_twland_rst"]["TWLAND"].values - opened["ol_twland_rst"]["TWLAND"].values
        )[:, tile_indices]

        water_years = config["water_years"]
        expected_dates = pd.DatetimeIndex(
            np.concatenate([water_year_dates(year).values for year in water_years])
        )
        if expected_dates.min() < pd.Timestamp(config["modis_only_start"]):
            raise AssertionError("Water-year sample starts before the MODIS-only interval")
        if expected_dates.max() >= pd.Timestamp(config["microwave_soil_moisture_da_start"]):
            raise AssertionError("Microwave soil-moisture DA entered the water-year sample")
        if expected_dates.duplicated().any() or len(expected_dates) != 72:
            raise AssertionError("Expected exactly 72 unique monthly samples")

        shape = (len(water_years), 12, len(tile_indices))
        monthly_names = [
            "ol_snow_mass",
            "dSnow_mass",
            "snow_net",
            "snow_abs_netpack",
            "dSnowmelt",
            "dInfiltration",
            "dET",
            "dRunoff_surface",
            "dBaseflow",
            "dRunoff_total",
            "dSFMC",
            "dRZMC",
            "dTWLAND",
            "dWCHANGELAND",
            "dPrecipitation",
            "dPeatFreeStandingWater",
        ]
        monthly = {name: np.full(shape, np.nan, dtype="float32") for name in monthly_names}
        annual_names = [
            "I_snow",
            "snow_abs_netpack",
            "dSnowmelt",
            "dInfiltration",
            "dET",
            "dRunoff_surface",
            "dBaseflow",
            "dRunoff_total",
            "dStorage",
            "dStorage_monthly_proxy",
            "dStorage_process_tendency",
            "dPeatFreeStandingWater",
            "residual",
            "process_balance_diagnostic",
            "dPrecipitation",
            "ol_mam_snow_mass",
        ]
        annual = {
            name: np.full((len(water_years), len(tile_indices)), np.nan, dtype="float32")
            for name in annual_names
        }

        for year_index, water_year in enumerate(water_years):
            dates = water_year_dates(water_year)
            date_slice = slice(dates[0], dates[-1])
            for name, dataset in opened.items():
                if name.endswith("_twland_rst"):
                    continue
                found = pd.DatetimeIndex(dataset.time.sel(time=date_slice).values)
                if not found.equals(dates):
                    raise AssertionError(f"{water_year}: monthly timestamps differ across inputs")

            ol_land = opened["ol_land"].sel(time=date_slice).isel(tile=tile_indices)
            da_land = opened["da_land"].sel(time=date_slice).isel(tile=tile_indices)
            ol_flux = opened["ol_flux"].sel(time=date_slice).isel(tile=tile_indices)
            da_flux = opened["da_flux"].sel(time=date_slice).isel(tile=tile_indices)
            ol_water = opened["ol_water"].sel(time=date_slice).isel(tile=tile_indices)
            da_water = opened["da_water"].sel(time=date_slice).isel(tile=tile_indices)
            catch = opened["catch"].sel(time=date_slice).isel(tile=tile_indices)
            ol_fsw = opened["ol_fsw"].sel(time=date_slice).isel(tile=tile_indices)
            da_fsw = opened["da_fsw"].sel(time=date_slice).isel(tile=tile_indices)

            monthly["ol_snow_mass"][year_index] = ol_land["SNOMASLAND"].values
            monthly["dSnow_mass"][year_index] = (da_land["SNOMASLAND"] - ol_land["SNOMASLAND"]).values
            monthly["snow_net"][year_index] = catch["snow_net"].values
            monthly["snow_abs_netpack"][year_index] = catch["snow_abs_netpack"].values
            monthly["dSnowmelt"][year_index] = monthly_flux_total(da_water["SMLAND"] - ol_water["SMLAND"]).values
            monthly["dInfiltration"][year_index] = monthly_flux_total(da_water["QINFILLAND"] - ol_water["QINFILLAND"]).values
            monthly["dET"][year_index] = monthly_flux_total(da_flux["EVLAND"] - ol_flux["EVLAND"]).values
            monthly["dRunoff_surface"][year_index] = monthly_flux_total(da_water["RUNSURFLAND"] - ol_water["RUNSURFLAND"]).values
            monthly["dBaseflow"][year_index] = monthly_flux_total(da_water["BASEFLOWLAND"] - ol_water["BASEFLOWLAND"]).values
            monthly["dRunoff_total"][year_index] = monthly["dRunoff_surface"][year_index] + monthly["dBaseflow"][year_index]
            monthly["dSFMC"][year_index] = (da_land["SFMC"] - ol_land["SFMC"]).values
            monthly["dRZMC"][year_index] = (da_land["RZMC"] - ol_land["RZMC"]).values
            monthly["dTWLAND"][year_index] = (da_water["TWLAND"] - ol_water["TWLAND"]).values
            monthly["dWCHANGELAND"][year_index] = monthly_flux_total(da_water["WCHANGELAND"] - ol_water["WCHANGELAND"]).values
            monthly["dPrecipitation"][year_index] = monthly_flux_total(da_land["PRECTOTCORRLAND"] - ol_land["PRECTOTCORRLAND"]).values
            monthly["dPeatFreeStandingWater"][year_index] = monthly_flux_total(
                da_fsw["PEATCLSM_FSWCHANGE"] - ol_fsw["PEATCLSM_FSWCHANGE"]
            ).values

            annual["I_snow"][year_index] = np.nansum(monthly["snow_net"][year_index], axis=0)
            annual["snow_abs_netpack"][year_index] = np.nansum(monthly["snow_abs_netpack"][year_index], axis=0)
            for name in ["dSnowmelt", "dInfiltration", "dET", "dRunoff_surface", "dBaseflow", "dRunoff_total", "dWCHANGELAND", "dPrecipitation", "dPeatFreeStandingWater"]:
                target = "dStorage_process_tendency" if name == "dWCHANGELAND" else name
                annual[target][year_index] = np.nansum(monthly[name][year_index], axis=0)

            tw_delta = opened["da_water"]["TWLAND"] - opened["ol_water"]["TWLAND"]
            prior_sep = tw_delta.sel(time=f"{water_year - 1}-09-01").values[tile_indices]
            current_sep = tw_delta.sel(time=f"{water_year}-09-01").values[tile_indices]
            annual["dStorage_monthly_proxy"][year_index] = current_sep - prior_sep
            annual["dStorage"][year_index] = (
                restart_delta[restart_index[water_year]]
                - restart_delta[restart_index[water_year - 1]]
            )
            annual["residual"][year_index] = (
                annual["I_snow"][year_index]
                - annual["dET"][year_index]
                - annual["dRunoff_total"][year_index]
                - annual["dStorage"][year_index]
                - annual["dPeatFreeStandingWater"][year_index]
            )
            annual["process_balance_diagnostic"][year_index] = (
                annual["dET"][year_index]
                + annual["dRunoff_total"][year_index]
                + annual["dStorage_process_tendency"][year_index]
            )
            mam = ol_land["SNOMASLAND"].isel(time=[5, 6, 7])
            mam_weights = mam.time.dt.days_in_month.astype("float64")
            annual["ol_mam_snow_mass"][year_index] = mam.weighted(mam_weights).mean("time").values

        finite_pair = np.isfinite(annual["I_snow"]) & np.isfinite(annual["snow_abs_netpack"])
        if np.any(annual["snow_abs_netpack"][finite_pair] + 1.0e-5 < np.abs(annual["I_snow"][finite_pair])):
            raise AssertionError("Annual snow_abs_netpack is smaller than abs(I_snow)")
        max_precip = float(np.nanmax(np.abs(annual["dPrecipitation"])))
        precip_limit = config["maximum_abs_annual_tile_precip_difference_kg_m2"]
        if max_precip > precip_limit:
            raise AssertionError(
                "Compressed DA and OL precipitation differ beyond the declared annual tile "
                f"limit: {max_precip} > {precip_limit} kg m-2"
            )

        soil = soil_metrics(
            monthly["dRZMC"],
            monthly["dSFMC"],
            config["rzmc_persistence_threshold_m3_m3"],
        )
        coords = {
            "water_year": np.asarray(water_years, dtype="int16"),
            "water_month": np.arange(1, 13, dtype="int8"),
            "tile": tile_ids,
            "lat": ("tile", lat),
            "lon": ("tile", lon),
            "area": ("tile", area),
        }
        data_vars = {}
        for name, values in monthly.items():
            data_vars[f"{name}_monthly"] = (("water_year", "water_month", "tile"), values)
        for name, values in {**annual, **soil}.items():
            data_vars[name] = (("water_year", "tile"), values)
        dataset = xr.Dataset(data_vars=data_vars, coords=coords)
        dataset.attrs.update(
            {
                "title": config["analysis_title"],
                "water_year_definition": "October through September",
                "domain": "Analysis A Northern Hemisphere seasonal-snow mask",
                "storage_closing_term": config["storage"]["primary_closing_term"],
                "storage_endpoint_source": opened["ol_twland_rst"].attrs.get("formula", ""),
                "peat_free_standing_water_role": config["peat_free_standing_water"]["role"],
                "peat_free_standing_water_source": opened["ol_fsw"].attrs.get("note", ""),
                "storage_limitation": config["storage"]["state_proxy_limitation"],
                "wchangeland_role": config["storage"]["wchangeland_role"],
                "soil_moisture_role": "SFMC/RZMC are diagnostic states and are excluded from the mass-closing equation",
                "maximum_abs_annual_tile_precip_difference_kg_m2": max_precip,
                "precipitation_check": config["precipitation_check"]["compressed_input_note"],
                "source_ol": opened["ol_land"].attrs.get("source_root", ""),
                "source_da": opened["da_land"].attrs.get("source_root", ""),
            }
        )
        units = {
            "I_snow": "kg m-2",
            "snow_abs_netpack": "kg m-2",
            "dSnowmelt": "kg m-2",
            "dInfiltration": "kg m-2",
            "dET": "kg m-2",
            "dRunoff_surface": "kg m-2",
            "dBaseflow": "kg m-2",
            "dRunoff_total": "kg m-2",
            "dStorage": "kg m-2",
            "dStorage_monthly_proxy": "kg m-2",
            "dStorage_process_tendency": "kg m-2",
            "dPeatFreeStandingWater": "kg m-2",
            "residual": "kg m-2",
            "process_balance_diagnostic": "kg m-2",
            "dPrecipitation": "kg m-2",
            "ol_mam_snow_mass": "kg m-2",
            "dSFMC_monthly": "m3 m-3",
            "dRZMC_monthly": "m3 m-3",
            "peak_dRZMC": "m3 m-3",
            "peak_dRZMC_positive": "m3 m-3",
            "peak_dSFMC": "m3 m-3",
            "peak_dSFMC_positive": "m3 m-3",
            "dSFMC_at_dRZMC_peak": "m3 m-3",
            "mjj_mean_dRZMC": "m3 m-3",
            "amj_mean_dRZMC": "m3 m-3",
            "september_dRZMC": "m3 m-3",
            "rzmc_positive_month_count": "months",
            "rzmc_months_from_peak_to_below_threshold": "months",
        }
        for name, unit in units.items():
            if name in dataset:
                dataset[name].attrs["units"] = unit
        return dataset
    finally:
        for dataset in opened.values():
            dataset.close()


def absolute_domain_budgets(config: dict) -> pd.DataFrame:
    """Area-weighted OL and DA annual water balances over the seasonal-snow mask.

    The differential budget cancels precipitation and every process common to
    both runs, which magnifies its residual: two per-run closure offsets of
    opposite sign add, and are then divided by the much smaller snow-DA input.
    These absolute balances supply the context that each run closes on its own.
    """
    data_dir = Path(config["data_directory"])
    root = find_repo_root()
    paths = {
        "ol_land": data_dir / "OLv8_land_variables_2000_2024_compressed.nc",
        "da_land": data_dir / "DAv8_land_variables_2000_2024_compressed.nc",
        "ol_flux": data_dir / "OLv8_flux_core_2000_2024_compressed.nc",
        "da_flux": data_dir / "DAv8_flux_core_2000_2024_compressed.nc",
        "ol_water": data_dir / "OLv8_water_budget_2000_2024_compressed.nc",
        "da_water": data_dir / "DAv8_water_budget_2000_2024_compressed.nc",
        "catch": data_dir / "catch_progn_raw_monthly_cumulative_200006_202405.nc",
        "ol_fsw": data_dir / "OLv8_peat_fsw_2000_2024_compressed.nc",
        "da_fsw": data_dir / "DAv8_peat_fsw_2000_2024_compressed.nc",
        "ol_twland_rst": data_dir / "OLv8_twland_wy_endpoints_2000_2006.nc",
        "da_twland_rst": data_dir / "DAv8_twland_wy_endpoints_2000_2006.nc",
    }
    opened = {name: xr.open_dataset(path, decode_times=True) for name, path in paths.items()}
    try:
        mask = build_static_mask(opened["ol_land"], opened["da_land"], config)
        tile_indices = np.flatnonzero(mask.values)
        sys.path.insert(0, str(root / "common/python/io"))
        from read_GEOSldas import read_tilecoord

        tilecoord = read_tilecoord(str(data_dir / "LS_OLv8_M36.ldas_tilecoord.bin"))
        area = np.asarray(tilecoord["area"], dtype="float64")[tile_indices]

        restart_time = pd.DatetimeIndex(opened["ol_twland_rst"]["time"].values)
        restart_index = {int(stamp.year): position for position, stamp in enumerate(restart_time)}

        rows = []
        for run in ["OL", "DA"]:
            prefix = run.lower()
            land = opened[f"{prefix}_land"]
            flux = opened[f"{prefix}_flux"]
            water = opened[f"{prefix}_water"]
            fsw = opened[f"{prefix}_fsw"]
            twland = opened[f"{prefix}_twland_rst"]["TWLAND"].values[:, tile_indices]
            for water_year in config["water_years"]:
                dates = water_year_dates(water_year)
                date_slice = slice(dates[0], dates[-1])
                select = dict(time=date_slice)
                row = {"run": run, "water_year": int(water_year)}
                row["precipitation"] = weighted_mean(
                    sum_monthly_flux(land["PRECTOTCORRLAND"].sel(**select)).values[tile_indices], area
                )
                # The snow analysis increment exists only in DA; OL is unassimilated.
                snow_net = np.nansum(
                    opened["catch"]["snow_net"].sel(**select).values[:, tile_indices], axis=0
                )
                row["I_snow"] = weighted_mean(snow_net, area) if run == "DA" else 0.0
                row["ET"] = weighted_mean(
                    sum_monthly_flux(flux["EVLAND"].sel(**select)).values[tile_indices], area
                )
                row["runoff_surface"] = weighted_mean(
                    sum_monthly_flux(water["RUNSURFLAND"].sel(**select)).values[tile_indices], area
                )
                row["baseflow"] = weighted_mean(
                    sum_monthly_flux(water["BASEFLOWLAND"].sel(**select)).values[tile_indices], area
                )
                row["storage"] = weighted_mean(
                    twland[restart_index[water_year]] - twland[restart_index[water_year - 1]], area
                )
                row["peat_free_standing_water"] = weighted_mean(
                    sum_monthly_flux(fsw["PEATCLSM_FSWCHANGE"].sel(**select)).values[tile_indices], area
                )
                row["runoff_total"] = row["runoff_surface"] + row["baseflow"]
                row["input_total"] = row["precipitation"] + row["I_snow"]
                row["residual"] = (
                    row["input_total"]
                    - row["ET"]
                    - row["runoff_total"]
                    - row["storage"]
                    - row["peat_free_standing_water"]
                )
                row["fraction_residual"] = row["residual"] / row["input_total"]
                rows.append(row)
        frame = pd.DataFrame(rows)
        means = []
        for run in ["OL", "DA"]:
            subset = frame[frame["run"] == run]
            mean_row = {"run": run, "water_year": "6-WY mean"}
            for column in subset.columns.drop(["run", "water_year"]):
                mean_row[column] = subset[column].mean()
            mean_row["fraction_residual"] = mean_row["residual"] / mean_row["input_total"]
            means.append(mean_row)
        return pd.concat([frame, pd.DataFrame(means)], ignore_index=True)
    finally:
        for dataset in opened.values():
            dataset.close()


def check_absolute_matches_differential(annual: pd.DataFrame, absolute: pd.DataFrame) -> None:
    """Assert the differential budget equals DA minus OL of the absolute budgets.

    The two are built by independent routes: the differential from per-tile DA-OL
    differences carried through the tile archive, the absolute from per-run domain
    means. Agreement is not automatic, and a mismatch would mean the mask, the tile
    ordering, or a term definition had diverged between them.
    """
    pairs = [
        ("dET", "ET"),
        ("dRunoff_total", "runoff_total"),
        ("dRunoff_surface", "runoff_surface"),
        ("dBaseflow", "baseflow"),
        ("dStorage", "storage"),
        ("dPeatFreeStandingWater", "peat_free_standing_water"),
        ("residual", "residual"),
        ("I_snow", "I_snow"),
    ]
    annual_only = annual[annual["water_year"].apply(lambda value: isinstance(value, (int, np.integer)))]
    ol = absolute[absolute["run"] == "OL"].set_index("water_year")
    da = absolute[absolute["run"] == "DA"].set_index("water_year")
    for water_year in annual_only["water_year"]:
        for differential_name, absolute_name in pairs:
            expected = da.loc[water_year, absolute_name] - ol.loc[water_year, absolute_name]
            found = float(annual_only.set_index("water_year").loc[water_year, differential_name])
            if not np.isclose(found, expected, atol=5.0e-3):
                raise AssertionError(
                    f"WY{water_year} {differential_name}: differential {found:.6f} does not match "
                    f"DA-OL absolute {expected:.6f}"
                )


def dataset_to_frame(dataset: xr.Dataset, variables: list[str]) -> pd.DataFrame:
    years = dataset.water_year.values
    tiles = dataset.tile.values
    frame = pd.DataFrame(
        {
            "water_year": np.repeat(years, len(tiles)),
            "tile": np.tile(tiles, len(years)),
            "lat": np.tile(dataset.lat.values, len(years)),
            "lon": np.tile(dataset.lon.values, len(years)),
            "area": np.tile(dataset.area.values, len(years)),
        }
    )
    for variable in variables:
        frame[variable] = dataset[variable].values.reshape(-1)
    return frame


def domain_budget_tables(dataset: xr.Dataset, config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    variables = BUDGET_TERMS + [
        "dRunoff_surface",
        "dBaseflow",
        "dSnowmelt",
        "dInfiltration",
        "dStorage_monthly_proxy",
        "dStorage_process_tendency",
        "process_balance_diagnostic",
        "dPrecipitation",
        "snow_abs_netpack",
    ]
    annual_rows = []
    weights = dataset.area.values
    for year_index, water_year in enumerate(dataset.water_year.values):
        row = {"water_year": int(water_year), "sample": "all"}
        for variable in variables:
            row[variable] = weighted_mean(dataset[variable].values[year_index], weights)
        denominator = row["I_snow"]
        valid_fraction = abs(denominator) >= config["minimum_abs_domain_input_for_fraction_kg_m2"]
        row["fraction_status"] = "reported" if valid_fraction else "net input too close to zero"
        for variable in CLOSING_TERMS:
            row[f"fraction_{variable}"] = row[variable] / denominator if valid_fraction else np.nan
        annual_rows.append(row)
    annual = pd.DataFrame(annual_rows)
    mean_row = {"water_year": "6-WY mean", "sample": "all"}
    for variable in variables:
        mean_row[variable] = annual[variable].mean()
        mean_row[f"sd_{variable}"] = annual[variable].std(ddof=1)
        mean_row[f"min_{variable}"] = annual[variable].min()
        mean_row[f"max_{variable}"] = annual[variable].max()
    mean_row["fraction_status"] = "reported"
    for variable in CLOSING_TERMS:
        mean_row[f"fraction_{variable}"] = mean_row[variable] / mean_row["I_snow"]
    annual_with_mean = pd.concat([annual, pd.DataFrame([mean_row])], ignore_index=True)

    frame = dataset_to_frame(dataset, variables)
    partition_rows = []
    sample_masks = {
        "all": np.ones(len(frame), dtype=bool),
        "addition": frame["I_snow"].to_numpy() > 0,
        "removal": frame["I_snow"].to_numpy() < 0,
    }
    for sample_name, keep in sample_masks.items():
        row = {
            "sample": sample_name,
            "n_tile_years": int(keep.sum()),
            "n_tiles": int(frame.loc[keep, "tile"].nunique()),
        }
        for variable in variables:
            row[variable] = weighted_mean(
                frame[variable].to_numpy(), frame["area"].to_numpy(), keep
            )
        for variable in ["dRunoff_total"] + PARTITION_TERMS:
            row[f"fraction_{variable}"] = row[variable] / row["I_snow"]
        partition_rows.append(row)
    partition = pd.DataFrame(partition_rows)

    stats_rows = []
    for variable in variables:
        stats_rows.append(
            {
                "variable": variable,
                "mean_across_water_years": annual[variable].mean(),
                "sd_across_water_years": annual[variable].std(ddof=1),
                "minimum_water_year": annual[variable].min(),
                "maximum_water_year": annual[variable].max(),
            }
        )
    return annual_with_mean, partition, pd.DataFrame(stats_rows)


def bootstrap_partition(
    frame: pd.DataFrame,
    config: dict,
    block_degrees: float,
    seed: int,
) -> pd.DataFrame:
    codes = block_codes(frame, block_degrees)
    n_blocks = int(codes.max()) + 1
    replicates = config["bootstrap"]["replicates"]
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=replicates)
    rows = []
    for sample_name, keep in {
        "all": np.ones(len(frame), dtype=bool),
        "addition": frame["I_snow"].to_numpy() > 0,
        "removal": frame["I_snow"].to_numpy() < 0,
    }.items():
        weighted_input = frame["I_snow"].to_numpy() * frame["area"].to_numpy() * keep
        input_by_block = np.bincount(codes, weights=weighted_input, minlength=n_blocks)
        denominator = draws @ input_by_block
        for variable in ["dRunoff_total"] + PARTITION_TERMS:
            weighted_response = frame[variable].to_numpy() * frame["area"].to_numpy() * keep
            response_by_block = np.bincount(codes, weights=weighted_response, minlength=n_blocks)
            ratio = np.divide(
                draws @ response_by_block,
                denominator,
                out=np.full(replicates, np.nan),
                where=np.abs(denominator) > 0,
            )
            low, high = confidence_interval(ratio)
            rows.append(
                {
                    "sample": sample_name,
                    "variable": variable,
                    "block_degrees": block_degrees,
                    "ci_low": low,
                    "ci_high": high,
                }
            )
    return pd.DataFrame(rows)


def m3_regressions(dataset: xr.Dataset, config: dict) -> pd.DataFrame:
    variables = ["I_snow", "ol_mam_snow_mass"] + [
        "dRunoff_total",
        "dRunoff_surface",
        "dBaseflow",
        "dET",
        "dStorage",
        "dPeatFreeStandingWater",
        "residual",
    ] + SOIL_METRICS
    frame = dataset_to_frame(dataset, variables).dropna().copy()
    demean = ["I_snow", "ol_mam_snow_mass"] + variables[2:]
    for variable in demean:
        frame[f"{variable}_anom"] = frame[variable] - frame.groupby("tile")[variable].transform("mean")
    year_dummies = pd.get_dummies(
        pd.Categorical(frame["water_year"], categories=sorted(frame["water_year"].unique())),
        drop_first=True,
        dtype=float,
    ).to_numpy()
    x = np.column_stack(
        [
            np.ones(len(frame)),
            frame["I_snow_anom"].to_numpy(),
            frame["ol_mam_snow_mass_anom"].to_numpy(),
            year_dummies,
        ]
    )
    rows = []
    for block_index, block_degrees in enumerate(
        [config["bootstrap"]["primary_block_degrees"], config["bootstrap"]["coarse_block_degrees"]]
    ):
        codes = block_codes(frame, block_degrees)
        n_blocks = int(codes.max()) + 1
        rng = np.random.default_rng(config["bootstrap"]["seed"] + 20 + block_index)
        draws = rng.multinomial(
            n_blocks,
            np.full(n_blocks, 1.0 / n_blocks),
            size=config["bootstrap"]["replicates"],
        )
        for variable in variables[2:]:
            y = frame[f"{variable}_anom"].to_numpy()
            beta, *_ = np.linalg.lstsq(x, y, rcond=None)
            xtx, xty = block_sufficient_statistics(x, y, codes)
            slopes = np.full(config["bootstrap"]["replicates"], np.nan)
            for draw_index, weights in enumerate(draws):
                lhs = np.tensordot(weights, xtx, axes=(0, 0))
                rhs = np.tensordot(weights, xty, axes=(0, 0))
                draw_beta, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
                slopes[draw_index] = draw_beta[1]
            low, high = confidence_interval(slopes)
            rows.append(
                {
                    "response": variable,
                    "block_degrees": block_degrees,
                    "m3_beta": float(beta[1]),
                    "ci_low": low,
                    "ci_high": high,
                    "n_tiles": int(frame["tile"].nunique()),
                    "n_tile_years": int(len(frame)),
                }
            )
    result = pd.DataFrame(rows)
    primary = result[result["block_degrees"] == config["bootstrap"]["primary_block_degrees"]].set_index("response")
    closure = primary.loc[CLOSING_TERMS, "m3_beta"].sum()
    if not np.isclose(closure, 1.0, atol=1.0e-8):
        raise AssertionError(f"M3 budget slopes do not close to one: {closure}")
    return result


def monthly_climatologies(dataset: xr.Dataset) -> tuple[pd.DataFrame, pd.DataFrame]:
    area = dataset.area.values
    all_rows = []
    positive_rows = []
    positive = dataset["I_snow"].values > 0
    variables = [
        "ol_snow_mass_monthly",
        "dSnow_mass_monthly",
        "snow_net_monthly",
        "snow_abs_netpack_monthly",
        "dSnowmelt_monthly",
        "dInfiltration_monthly",
        "dET_monthly",
        "dRunoff_total_monthly",
        "dPeatFreeStandingWater_monthly",
        "dSFMC_monthly",
        "dRZMC_monthly",
    ]
    repeated_weights = np.tile(area, len(dataset.water_year))
    for month_index, label in enumerate(MONTH_LABELS):
        all_row = {"water_month": month_index + 1, "month": label}
        positive_row = {"water_month": month_index + 1, "month": label}
        keep_positive = positive.reshape(-1)
        for variable in variables:
            values = dataset[variable].values[:, month_index, :].reshape(-1)
            all_row[variable] = weighted_mean(values, repeated_weights)
            positive_row[variable] = weighted_mean(values, repeated_weights, keep_positive)
        all_rows.append(all_row)
        positive_rows.append(positive_row)
    return pd.DataFrame(all_rows), pd.DataFrame(positive_rows)


def binned_soil_diagnostics(dataset: xr.Dataset, config: dict) -> pd.DataFrame:
    variables = ["I_snow"] + SOIL_METRICS
    frame = dataset_to_frame(dataset, variables).dropna().copy()
    frame["I_snow_anom"] = frame["I_snow"] - frame.groupby("tile")["I_snow"].transform("mean")
    for variable in SOIL_METRICS:
        frame[f"{variable}_anom"] = frame[variable] - frame.groupby("tile")[variable].transform("mean")
    frame["bin"] = pd.qcut(
        frame["I_snow_anom"],
        q=config["soil_moisture_binned_diagnostics"],
        labels=False,
        duplicates="drop",
    )
    codes = block_codes(frame, config["bootstrap"]["primary_block_degrees"])
    n_blocks = int(codes.max()) + 1
    n_bins = int(frame["bin"].max()) + 1
    combined = codes * n_bins + frame["bin"].to_numpy(dtype=int)
    length = n_blocks * n_bins
    counts = np.bincount(combined, minlength=length).reshape(n_blocks, n_bins)
    rng = np.random.default_rng(config["bootstrap"]["seed"] + 50)
    draws = rng.multinomial(
        n_blocks,
        np.full(n_blocks, 1.0 / n_blocks),
        size=config["bootstrap"]["replicates"],
    )
    draw_counts = draws @ counts
    rows = []
    for variable in SOIL_METRICS:
        y = frame[f"{variable}_anom"].to_numpy()
        ysum = np.bincount(combined, weights=y, minlength=length).reshape(n_blocks, n_bins)
        draw_mean = np.divide(
            draws @ ysum,
            draw_counts,
            out=np.full(draw_counts.shape, np.nan, dtype="float64"),
            where=draw_counts > 0,
        )
        for bin_index in range(n_bins):
            keep = frame["bin"].to_numpy() == bin_index
            low, high = confidence_interval(draw_mean[:, bin_index])
            rows.append(
                {
                    "response": variable,
                    "bin": bin_index + 1,
                    "n": int(keep.sum()),
                    "x_mean": float(frame.loc[keep, "I_snow_anom"].mean()),
                    "y_mean": float(frame.loc[keep, f"{variable}_anom"].mean()),
                    "ci_low": low,
                    "ci_high": high,
                }
            )
    return pd.DataFrame(rows)


def soil_summary(dataset: xr.Dataset, config: dict) -> pd.DataFrame:
    frame = dataset_to_frame(
        dataset,
        ["I_snow", "peak_dRZMC", "peak_dRZMC_water_month", "peak_dSFMC", "peak_dSFMC_water_month", "dSFMC_at_dRZMC_peak", "mjj_mean_dRZMC", "amj_mean_dRZMC", "rzmc_positive_month_count", "rzmc_months_from_peak_to_below_threshold", "rzmc_persistence_censored", "september_dRZMC"],
    )
    keep = frame["I_snow"] > 0
    positive = frame.loc[keep].copy()
    rows = []
    for variable in [
        "peak_dRZMC",
        "peak_dSFMC",
        "dSFMC_at_dRZMC_peak",
        "mjj_mean_dRZMC",
        "amj_mean_dRZMC",
        "rzmc_positive_month_count",
        "rzmc_months_from_peak_to_below_threshold",
        "september_dRZMC",
    ]:
        rows.append(
            {
                "sample": "I_snow > 0",
                "metric": variable,
                "n": int(positive[variable].notna().sum()),
                "area_weighted_mean": weighted_mean(
                    positive[variable].to_numpy(), positive["area"].to_numpy()
                ),
                "median": float(positive[variable].median()),
                "q25": float(positive[variable].quantile(0.25)),
                "q75": float(positive[variable].quantile(0.75)),
            }
        )
    censored = positive["rzmc_persistence_censored"].astype("float64")
    rows.append(
        {
            "sample": "I_snow > 0",
            "metric": "rzmc_persistence_censored_fraction",
            "n": int(len(positive)),
            "area_weighted_mean": weighted_mean(
                censored.to_numpy(), positive["area"].to_numpy()
            ),
            "median": float(censored.median()),
            "q25": float(censored.quantile(0.25)),
            "q75": float(censored.quantile(0.75)),
        }
    )
    peak_timing = []
    for state, column in [("RZMC", "peak_dRZMC_water_month"), ("SFMC", "peak_dSFMC_water_month")]:
        valid = positive[column].notna()
        for month in range(1, 13):
            month_keep = valid & positive[column].eq(month)
            peak_timing.append(
                {
                    "state": state,
                    "water_month": month,
                    "month": MONTH_LABELS[month - 1],
                    "n": int(month_keep.sum()),
                    "area_weighted_fraction": float(
                        positive.loc[month_keep, "area"].sum()
                        / positive.loc[valid, "area"].sum()
                    ),
                }
            )
    summary = pd.DataFrame(rows)
    summary.attrs["peak_timing"] = pd.DataFrame(peak_timing)
    summary.attrs["persistence_threshold"] = config["rzmc_persistence_threshold_m3_m3"]
    return summary


def plot_monthly_climatology(data: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.8), sharex=True)
    specs = [
        ("ol_snow_mass_monthly", "OL snow mass", "kg m-2", "#457b9d"),
        ("dSnow_mass_monthly", "DA - OL snow mass", "kg m-2", "#2a9d8f"),
        ("snow_net_monthly", "Signed snow-DA input", "kg m-2 month-1", "#d1495b"),
        ("snow_abs_netpack_monthly", "Absolute snow-DA activity", "kg m-2 month-1", "#6a4c93"),
    ]
    for ax, (variable, title, unit, color) in zip(axes.flat, specs):
        ax.plot(data["water_month"], data[variable], marker="o", color=color)
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(alpha=0.22)
        ax.set_xticks(range(1, 13), MONTH_LABELS, rotation=30)
    fig.suptitle("WY2001-WY2006 seasonal-snow monthly climatology")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_annual_budgets(annual: pd.DataFrame, path: Path) -> None:
    terms = ["I_snow"] + CLOSING_TERMS
    labels = ["Snow-DA input", "Runoff", "ET", "Storage change", "Peatland free-standing water", "Residual"]
    colors = ["#264653", "#2a9d8f", "#e9c46a", "#457b9d", "#8d6bb1", "#d1495b"]
    x = np.arange(len(annual))
    width = 0.14
    fig, ax = plt.subplots(figsize=(11, 5.2))
    for index, (term, label, color) in enumerate(zip(terms, labels, colors)):
        ax.bar(x + (index - 2.5) * width, annual[term], width, label=label, color=color)
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.set_xticks(x, [f"WY{value}" if isinstance(value, (int, np.integer)) else value for value in annual["water_year"]])
    ax.set_ylabel("Area-weighted mean (kg m-2 water year-1)")
    ax.set_title("Differential snow-DA water budget by water year")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncols=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_positive_partition(partition: pd.DataFrame, uncertainty: pd.DataFrame, path: Path) -> None:
    row = partition.set_index("sample").loc["addition"]
    uncertainty = uncertainty[uncertainty["sample"] == "addition"].set_index("variable")
    terms = PARTITION_TERMS
    labels = ["Surface runoff", "Baseflow", "ET", "Storage", "Peatland free-standing water", "Residual"]
    colors = ["#2a9d8f", "#57a773", "#e9c46a", "#457b9d", "#8d6bb1", "#d1495b"]
    values = np.array([row[f"fraction_{term}"] for term in terms])
    low = np.array([uncertainty.loc[term, "ci_low"] for term in terms])
    high = np.array([uncertainty.loc[term, "ci_high"] for term in terms])
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.bar(labels, 100 * values, color=colors, yerr=np.vstack([100 * (values - low), 100 * (high - values)]), capsize=3)
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.set_ylabel("Fraction of positive snow-DA input (%)")
    ax.set_title("Six-water-year partition for snow-addition tile-years")
    ax.grid(axis="y", alpha=0.22)
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_soil_pathway(data: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.8), sharex=True)
    specs = [
        (["snow_net_monthly"], ["Snow-DA input"], "kg m-2 month-1"),
        (["dSnowmelt_monthly"], ["DA - OL snowmelt"], "kg m-2 month-1"),
        (["dSFMC_monthly", "dRZMC_monthly"], ["SFMC", "RZMC"], "m3 m-3"),
        (["dRunoff_total_monthly"], ["DA - OL runoff"], "kg m-2 month-1"),
    ]
    titles = ["Snow increment", "Snowmelt response", "Soil-moisture response", "Runoff response"]
    colors = [["#d1495b"], ["#6a4c93"], ["#e76f51", "#2a9d8f"], ["#264653"]]
    for ax, variables, labels, unit, title, palette in zip(axes.flat, [s[0] for s in specs], [s[1] for s in specs], [s[2] for s in specs], titles, colors):
        for variable, label, color in zip(variables, labels, palette):
            ax.plot(data["water_month"], data[variable], marker="o", label=label, color=color)
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(alpha=0.22)
        ax.set_xticks(range(1, 13), MONTH_LABELS, rotation=30)
        if len(variables) > 1:
            ax.legend(fontsize=8)
    fig.suptitle("Snow-addition tile-years: snow to soil moisture to runoff")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_soil_bins(data: pd.DataFrame, path: Path) -> None:
    metrics = ["peak_dRZMC_positive", "mjj_mean_dRZMC", "rzmc_positive_month_count", "september_dRZMC"]
    titles = ["Peak positive RZMC", "MJJ mean RZMC", "Months with RZMC > 0", "September RZMC"]
    units = ["m3 m-3", "m3 m-3", "months", "m3 m-3"]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for ax, metric, title, unit in zip(axes.flat, metrics, titles, units):
        sub = data[data["response"] == metric].sort_values("x_mean")
        yerr = np.vstack([sub["y_mean"] - sub["ci_low"], sub["ci_high"] - sub["y_mean"]])
        ax.errorbar(sub["x_mean"], sub["y_mean"], yerr=yerr, marker="o", color="#15616d", capsize=2)
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.axvline(0, color="0.35", linewidth=0.8)
        ax.set_title(title)
        ax.set_xlabel("Within-tile I_snow anomaly (kg m-2)")
        ax.set_ylabel(f"Within-tile response anomaly ({unit})")
        ax.grid(alpha=0.22)
    fig.suptitle("Soil-moisture response versus signed water-year snow-DA input")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report_fragment(
    annual: pd.DataFrame,
    absolute: pd.DataFrame,
    partition: pd.DataFrame,
    uncertainty: pd.DataFrame,
    regressions: pd.DataFrame,
    soil: pd.DataFrame,
    config: dict,
    maximum_abs_annual_tile_precip_difference: float,
    path: Path,
) -> None:
    annual_only = annual[annual["water_year"].apply(lambda value: isinstance(value, (int, np.integer)))]
    addition = partition.set_index("sample").loc["addition"]
    ci = uncertainty[uncertainty["sample"] == "addition"].set_index("variable")
    primary_reg = regressions[regressions["block_degrees"] == config["bootstrap"]["primary_block_degrees"]].set_index("response")
    soil_index = soil.set_index("metric")
    rzmc_peak_timing = soil.attrs["peak_timing"].query("state == 'RZMC'")
    modal_rzmc_peak_month = rzmc_peak_timing.loc[
        rzmc_peak_timing["area_weighted_fraction"].idxmax(), "month"
    ]
    residual_fraction = annual_only["residual"] / annual_only["I_snow"]
    absolute_only = absolute[absolute["water_year"].apply(lambda value: isinstance(value, (int, np.integer)))]
    lines = [
        "## Water-year differential snow-DA budget and soil-moisture response",
        "",
        "This extension follows the native signed snow-water increments through six complete MODIS-only water years (October 2000 through September 2006). SFMC and RZMC diagnose the timing and persistence of the soil response but are excluded from the closing equation because total land-water storage already contains soil water.",
        "",
        "### Storage-definition audit",
        "",
        "Integrated `DA - OL WCHANGELAND` does not equal the change in the `DA - OL TWLAND` state anomaly. Instead, `dET + dRunoff + integrated dWCHANGELAND` is near zero, showing that WCHANGELAND is a model-process tendency that omits the discontinuous analysis mass injection. The closing storage term is therefore an endpoint difference in `DA - OL TWLAND`.",
        "",
        "That endpoint is now the instantaneous 00Z October 1 state reconstructed from `catch_internal_rst` restarts as the 24-member ensemble mean of `CDCR2/(1-WPWET) - CATDEF + RZEXC + SRFEXC + CAPAC + WESNN1-3`. It replaces the September monthly-mean proxy used previously, which is retained as `dStorage_monthly_proxy` for comparison. The proxy is close in six-year aggregate but differs by up to about 20% in individual water years, and gets the sign wrong where the true change is small.",
        "",
        "### Peat free-standing water",
        "",
        "`catch_calc_wtotl` builds TWLAND from soil, canopy, and snow stores only; PEATCLSM free-standing surface water is deliberately excluded. On peat tiles (`POROS >= 0.9`) water moving into or out of that store therefore leaves the TWLAND-based budget entirely and lands in the residual. `PEATCLSM_FSWCHANGE` closes that gap and enters the budget as `dPeatFreeStandingWater`. It is zero by construction on non-peat tiles, matching the model's own `FSW_CHANGE = 0.` initialization.",
        "",
        "### Annual domain budgets",
        "",
        "| Water year | Snow-DA input | Extra runoff | Extra ET | Storage change | Free-standing water | Residual | Residual / input |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in annual.iterrows():
        label = f"WY{row['water_year']}" if isinstance(row["water_year"], (int, np.integer)) else str(row["water_year"])
        fraction = row["residual"] / row["I_snow"]
        lines.append(
            f"| {label} | {row['I_snow']:.2f} | {row['dRunoff_total']:.2f} | {row['dET']:.2f} | {row['dStorage']:.2f} | {row['dPeatFreeStandingWater']:.2f} | {row['residual']:.2f} | {100 * fraction:.1f}% |"
        )
    lines.extend(
        [
            "",
            f"The annual budget residual ranges from {100 * residual_fraction.min():.1f}% to {100 * residual_fraction.max():.1f}% of domain snow input. Conceptually identical precipitation forcing differs slightly after independent float32 compression: the maximum absolute annual tile discrepancy is {maximum_abs_annual_tile_precip_difference:.3f} kg m-2 and the largest absolute annual area-weighted domain-mean discrepancy is {annual_only['dPrecipitation'].abs().max():.6f} kg m-2. Snowmelt and infiltration are retained as pathway diagnostics and are not added to the closing terms.",
            "",
            "![Water-year monthly climatology](monthly_synthesis_report_figures/water_year_budget_monthly_climatology.png)",
            "",
            "![Annual water-year budgets](monthly_synthesis_report_figures/water_year_budget_annual.png)",
            "",
            "### Absolute closure of each run",
            "",
            "The differential budget cancels precipitation and every process common to both runs, so its residual is magnified twice over: the two per-run closure offsets carry opposite signs and therefore add, and the sum is then divided by the much smaller snow-DA input rather than by total input. Each run closes far more tightly on its own.",
            "",
            "| Run | Total input | ET | Runoff | Storage | Free-standing water | Residual | Residual / input |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in absolute[absolute["water_year"] == "6-WY mean"].iterrows():
        lines.append(
            f"| {row['run']} | {row['input_total']:.2f} | {row['ET']:.2f} | {row['runoff_total']:.2f} | "
            f"{row['storage']:.2f} | {row['peat_free_standing_water']:.2f} | {row['residual']:+.2f} | "
            f"{100 * row['fraction_residual']:+.3f}% |"
        )
    lines.extend(
        [
            "",
            "All values are six-water-year means in kg m-2 water year-1 over the seasonal-snow mask.",
            "",
            "### Positive-input partition",
            "",
        ]
    )
    for variable, label in [
        ("dRunoff_total", "total runoff"),
        ("dRunoff_surface", "surface runoff"),
        ("dBaseflow", "baseflow"),
        ("dET", "ET"),
        ("dStorage", "storage"),
        ("dPeatFreeStandingWater", "peatland free-standing water"),
        ("residual", "residual"),
    ]:
        lines.append(
            f"- {label}: {100 * addition[f'fraction_{variable}']:.1f}% "
            f"[5-degree spatial-block 95% interval {100 * ci.loc[variable, 'ci_low']:.1f}%, {100 * ci.loc[variable, 'ci_high']:.1f}%]"
        )
    lines.extend(
        [
            "",
            "![Positive-input six-year partition](monthly_synthesis_report_figures/water_year_budget_positive_partition.png)",
            "",
            "| Sample | Native signed input | Runoff | ET | Storage | Free-standing water | Residual |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for sample, label in [("all", "All-sample net"), ("addition", "Snow addition"), ("removal", "Snow removal")]:
        row = partition.set_index("sample").loc[sample]
        lines.append(
            f"| {label} | {row['I_snow']:.2f} kg m-2 | {100 * row['fraction_dRunoff_total']:.1f}% | {100 * row['fraction_dET']:.1f}% | {100 * row['fraction_dStorage']:.1f}% | {100 * row['fraction_dPeatFreeStandingWater']:.1f}% | {100 * row['fraction_residual']:.1f}% |"
        )
    lines.extend(
        [
            "",
            "Fractions use native signed mass and are never based on absolute snow activity. The snow-removal row therefore reports signed response divided by signed input; its positive percentages describe same-direction water removal.",
            "",
            "### Controlled water-year response",
            "",
            "| Response | M3 beta | 5-degree block 95% CI |",
            "|---|---:|---:|",
        ]
    )
    for variable, label in [("dRunoff_total", "Runoff"), ("dET", "ET"), ("dStorage", "Storage"), ("dPeatFreeStandingWater", "Peatland free-standing water"), ("residual", "Residual")]:
        row = primary_reg.loc[variable]
        lines.append(f"| {label} | {row['m3_beta']:.3f} | [{row['ci_low']:.3f}, {row['ci_high']:.3f}] |")
    lines.extend(
        [
            "",
            "These dimensionless M3 slopes use within-tile signed snow input, year effects, and OL MAM snow amount. By construction, the runoff, ET, storage, peat free-water, and residual slopes sum to one; the direct domain accounting remains the primary budget result.",
            "",
            "### Soil-moisture consequence",
            "",
            f"For snow-addition tile-years, the area-weighted mean peak RZMC response is {soil_index.loc['peak_dRZMC', 'area_weighted_mean']:.4f} m3 m-3 and {modal_rzmc_peak_month} is the most common peak month, although peak timing is broad. The MJJ mean response is {soil_index.loc['mjj_mean_dRZMC', 'area_weighted_mean']:.4f} m3 m-3, RZMC is positive for {soil_index.loc['rzmc_positive_month_count', 'area_weighted_mean']:.1f} of 12 months on average, and the mean September response remains {soil_index.loc['september_dRZMC', 'area_weighted_mean']:.4f} m3 m-3.",
            "",
            f"Persistence is strongly right-censored: {100 * soil_index.loc['rzmc_persistence_censored_fraction', 'area_weighted_mean']:.1f}% of snow-addition tile-years never fall below {config['rzmc_persistence_threshold_m3_m3']:.3f} m3 m-3 after their within-year peak by September. Among the uncensored cases, the area-weighted mean time to the threshold is {soil_index.loc['rzmc_months_from_peak_to_below_threshold', 'area_weighted_mean']:.1f} months. Because the mean DA-minus-OL RZMC anomaly is already positive in October, these counts include inherited state differences from prior assimilation and should not be read as the residence time of only the current water-year increment. They are state diagnostics, not additional mass terms.",
            "",
            "![Snow-to-soil-moisture pathway](monthly_synthesis_report_figures/water_year_soil_moisture_pathway.png)",
            "",
            "![Soil-moisture binned diagnostics](monthly_synthesis_report_figures/water_year_soil_moisture_binned.png)",
            "",
            "The interpretation is deliberately model-internal: snow DA modifies snow water, the root zone becomes wetter during melt, and the model redistributes that perturbation through ET, runoff, and changing storage. RZMC and SFMC describe how strongly and for how long the soil state changes; they are not independent validation and do not enter the mass balance twice.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    config = load_config(args.config, args.bootstrap_replicates)
    root = find_repo_root()
    output_dir = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget"
    figure_dir = output_dir / "figures"
    report_figure_dir = root / "projects/M21C_ls/docs/monthly_synthesis_report_figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_figure_dir.mkdir(parents=True, exist_ok=True)

    archive = output_dir / "water_year_tile_budgets.nc"
    if archive.exists() and not args.rebuild_tile_budget:
        print(f"Reusing tile-level water-year archive: {archive}", flush=True)
        with xr.open_dataset(archive) as stored:
            dataset = stored.load()
    else:
        print("Building tile-level water-year budgets and monthly soil trajectories...", flush=True)
        dataset = load_budget_dataset(config)
        encoding = {
            name: {"zlib": True, "complevel": 3, "dtype": "float32"}
            for name, data in dataset.data_vars.items()
            if np.issubdtype(data.dtype, np.floating)
        }
        dataset.to_netcdf(archive, encoding=encoding)
        print(f"Archived tile-level water-year budgets: {archive}", flush=True)

    annual, partition, annual_stats = domain_budget_tables(dataset, config)
    annual.to_csv(output_dir / "annual_domain_budgets.csv", index=False)
    partition.to_csv(output_dir / "six_year_integrated_partitions.csv", index=False)
    annual_stats.to_csv(output_dir / "six_year_annual_statistics.csv", index=False)
    frame = dataset_to_frame(
        dataset,
        list(dict.fromkeys(BUDGET_TERMS + PARTITION_TERMS)),
    )
    uncertainty_primary = bootstrap_partition(
        frame,
        config,
        config["bootstrap"]["primary_block_degrees"],
        config["bootstrap"]["seed"],
    )
    uncertainty_coarse = bootstrap_partition(
        frame,
        config,
        config["bootstrap"]["coarse_block_degrees"],
        config["bootstrap"]["seed"] + 1,
    )
    uncertainty = pd.concat([uncertainty_primary, uncertainty_coarse], ignore_index=True)
    uncertainty.to_csv(output_dir / "partition_spatial_block_uncertainty.csv", index=False)

    absolute = absolute_domain_budgets(config)
    absolute.to_csv(output_dir / "annual_absolute_budgets.csv", index=False)
    check_absolute_matches_differential(annual, absolute)

    regressions = m3_regressions(dataset, config)
    regressions.to_csv(output_dir / "water_year_m3_regressions.csv", index=False)
    climatology, positive_climatology = monthly_climatologies(dataset)
    climatology.to_csv(output_dir / "monthly_climatology_all.csv", index=False)
    positive_climatology.to_csv(output_dir / "monthly_climatology_snow_addition.csv", index=False)
    bins = binned_soil_diagnostics(dataset, config)
    bins.to_csv(output_dir / "soil_moisture_binned_diagnostics.csv", index=False)
    soil = soil_summary(dataset, config)
    soil.to_csv(output_dir / "soil_moisture_summary_snow_addition.csv", index=False)
    soil.attrs["peak_timing"].to_csv(output_dir / "soil_moisture_peak_timing.csv", index=False)

    figures = {
        "water_year_budget_monthly_climatology.png": lambda path: plot_monthly_climatology(climatology, path),
        "water_year_budget_annual.png": lambda path: plot_annual_budgets(annual, path),
        "water_year_budget_positive_partition.png": lambda path: plot_positive_partition(
            partition,
            uncertainty_primary,
            path,
        ),
        "water_year_soil_moisture_pathway.png": lambda path: plot_soil_pathway(
            positive_climatology,
            path,
        ),
        "water_year_soil_moisture_binned.png": lambda path: plot_soil_bins(bins, path),
    }
    for name, create in figures.items():
        figure_path = figure_dir / name
        create(figure_path)
        shutil.copy2(figure_path, report_figure_dir / name)

    fragment = output_dir / "water_year_report_section.md"
    write_report_fragment(
        annual,
        absolute,
        partition,
        uncertainty_primary,
        regressions,
        soil,
        config,
        float(dataset.attrs["maximum_abs_annual_tile_precip_difference_kg_m2"]),
        fragment,
    )
    print("Water-year workflow complete")
    print(f"Outputs: {output_dir}")
    print(f"Report fragment: {fragment}")


if __name__ == "__main__":
    main()
