#!/usr/bin/env python3
"""Falsification analysis for the M21C monthly-synthesis Analysis A result.

The workflow intentionally rebuilds the original Analysis A sample from its
native monthly inputs before applying spatial and temporal controls. Derived
tables and figures are written below output/monthly_synthesis_diagnostics and
the compact report figures are copied to the tracked documentation directory.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


SEASONS = {
    "MAM": [3, 4, 5],
    "AMJ": [4, 5, 6],
    "MJJ": [5, 6, 7],
    "JJA": [6, 7, 8],
}
WATER_FLUX_VARS = {
    "PRECTOTCORRLAND",
    "EVLAND",
    "RUNSURFLAND",
    "BASEFLOWLAND",
    "QINFILLAND",
    "SMLAND",
    "WCHANGELAND",
}
RESPONSE_LABELS = {
    "snowmelt": "AMJ snowmelt",
    "infiltration": "AMJ infiltration",
    "rzmc": "MJJ RZMC",
    "et": "MJJ ET",
    "total_runoff": "MJJ total runoff",
    "total_water": "MJJ total water",
}
ORIGINAL_METRICS = {
    "snowmelt": "abs_dsmland_amj",
    "infiltration": "abs_dqinfilland_amj",
    "rzmc": "abs_drzmc_mjj",
    "et": "abs_devland_mjj",
    "total_runoff": "abs_dtotal_runoff_mjj",
    "total_water": "abs_dtwland_mjj",
}


@dataclass(frozen=True)
class ModelSpec:
    name: str
    x_columns: tuple[str, ...]
    coefficient_index: int
    intercept: bool


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / ".git").exists():
            return parent
    raise FileNotFoundError("Could not locate repository root")


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=root / "projects/M21C_ls/config/analysis_a_robustness.json",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=None,
        help="Testing override; production must use at least the configured 1000 replicates.",
    )
    return parser.parse_args()


def load_config(path: Path, override_replicates: int | None = None) -> dict:
    config = json.loads(path.read_text())
    if override_replicates is not None:
        config["bootstrap"]["replicates"] = int(override_replicates)
    if config["bootstrap"]["replicates"] < 1000 and override_replicates is None:
        raise ValueError("Production spatial bootstrap requires at least 1000 replicates")
    return config


def sum_preserve_nan(data: xr.DataArray, dim: str = "time") -> xr.DataArray:
    count = data.notnull().sum(dim)
    return data.fillna(0.0).sum(dim).where(count > 0)


def seasonal_aggregate(
    data: xr.DataArray,
    years: list[int],
    season: str,
    method: str,
) -> xr.DataArray:
    pieces = []
    for year in years:
        months = SEASONS[season]
        selected = data.sel(
            time=(data.time.dt.year == year) & data.time.dt.month.isin(months)
        )
        if selected.sizes.get("time", 0) != len(months):
            raise ValueError(
                f"{season} {year}: expected {len(months)} months, "
                f"found {selected.sizes.get('time', 0)}"
            )
        if method == "sum":
            value = sum_preserve_nan(selected)
        elif method == "mean":
            weights = selected.time.dt.days_in_month.astype("float64")
            value = selected.weighted(weights).mean("time", skipna=True)
        elif method == "flux_total":
            seconds = selected.time.dt.days_in_month.astype("float64") * 86400.0
            value = sum_preserve_nan(selected * seconds)
        else:
            raise ValueError(f"Unknown seasonal aggregation method: {method}")
        pieces.append(value)
    return xr.concat(pieces, dim=pd.Index(years, name="year"))


def model_season(
    datasets: dict[str, dict[str, xr.Dataset]],
    experiment: str,
    variable: str,
    years: list[int],
    season: str,
) -> xr.DataArray:
    if variable == "TOTAL_RUNOFF":
        return model_season(datasets, experiment, "RUNSURFLAND", years, season) + model_season(
            datasets, experiment, "BASEFLOWLAND", years, season
        )
    for dataset in datasets[experiment].values():
        if variable in dataset:
            method = "flux_total" if variable in WATER_FLUX_VARS else "mean"
            return seasonal_aggregate(dataset[variable], years, season, method)
    raise KeyError(f"Missing model variable {variable}")


def build_original_table(config: dict) -> tuple[pd.DataFrame, dict[str, int]]:
    data_dir = Path(config["data_directory"])
    years = config["years"]
    paths = {
        "land": ("OLv8_land_variables_2000_2024_compressed.nc", "DAv8_land_variables_2000_2024_compressed.nc"),
        "flux": ("OLv8_flux_core_2000_2024_compressed.nc", "DAv8_flux_core_2000_2024_compressed.nc"),
        "water": ("OLv8_water_budget_2000_2024_compressed.nc", "DAv8_water_budget_2000_2024_compressed.nc"),
    }
    datasets: dict[str, dict[str, xr.Dataset]] = {"ol": {}, "da": {}}
    for group, (ol_name, da_name) in paths.items():
        datasets["ol"][group] = xr.open_dataset(data_dir / ol_name, decode_times=True)
        datasets["da"][group] = xr.open_dataset(data_dir / da_name, decode_times=True)
    catch = xr.open_dataset(
        data_dir / "catch_progn_raw_monthly_cumulative_200006_202405.nc",
        decode_times=True,
    )
    try:
        land_ol = datasets["ol"]["land"]
        land_da = datasets["da"]["land"]
        all_times = pd.DatetimeIndex(land_ol.time.values)
        analysis_dates = all_times[
            all_times.year.isin(years)
            & all_times.month.isin(sorted({m for s in ("MAM", "AMJ", "MJJ") for m in SEASONS[s]}))
        ]
        start = pd.Timestamp(config["modis_only_start"])
        end = pd.Timestamp(config["modis_only_end"])
        if analysis_dates.empty or analysis_dates.min() < start or analysis_dates.max() > end:
            raise AssertionError("Analysis A dates escape the MODIS-only interval")
        if max(years) >= 2007:
            raise AssertionError("A microwave soil-moisture DA year entered Analysis A")

        scf_pair = xr.concat([land_ol["FRLANDSNO"], land_da["FRLANDSNO"]], dim="experiment")
        swe_pair = xr.concat([land_ol["SNOMASLAND"], land_da["SNOMASLAND"]], dim="experiment")
        scf_any = scf_pair.max("experiment").max("time", skipna=True).load()
        swe_any = swe_pair.max("experiment").max("time", skipna=True).load()
        jja_scf = (
            scf_pair.max("experiment")
            .sel(time=scf_pair.time.dt.month.isin(SEASONS["JJA"]))
            .mean("time", skipna=True)
            .load()
        )
        lat = land_ol["lat"].load()
        lon = land_ol["lon"].load()
        mask_cfg = config["mask"]
        mask = (
            (lat > mask_cfg["minimum_latitude_degrees_north"])
            & (
                (scf_any > mask_cfg["snow_cover_possible_threshold"])
                | (swe_any > mask_cfg["snow_mass_possible_threshold_kg_m2"])
            )
            & (jja_scf < mask_cfg["maximum_mean_jja_snow_cover_fraction"])
            & np.isfinite(lat)
            & np.isfinite(lon)
        ).load()

        snow_abs = seasonal_aggregate(catch["snow_abs_netpack"], years, "MAM", "sum").load()
        snow_net = seasonal_aggregate(catch["snow_net"], years, "MAM", "sum").load()
        ol_snow = model_season(datasets, "ol", "SNOMASLAND", years, "MAM").load()
        finite_pair = np.isfinite(snow_abs.values) & np.isfinite(snow_net.values)
        tolerance = 1.0e-5
        violations = finite_pair & (snow_abs.values + tolerance < np.abs(snow_net.values))
        if violations.any():
            year_idx, tile_idx = np.argwhere(violations)[0]
            raise AssertionError(
                "snow_abs_netpack < abs(snow_net) at "
                f"year={years[year_idx]}, tile={tile_idx}: "
                f"{snow_abs.values[year_idx, tile_idx]} < "
                f"{abs(snow_net.values[year_idx, tile_idx])}"
            )

        responses = {}
        for key, spec in config["responses"].items():
            da = model_season(datasets, "da", spec["variable"], years, spec["season"])
            ol = model_season(datasets, "ol", spec["variable"], years, spec["season"])
            responses[key] = (da - ol).load()

        valid_tiles = np.flatnonzero(mask.values)
        frames = []
        for year_index, year in enumerate(years):
            record = {
                "year": np.full(valid_tiles.size, year, dtype=np.int16),
                "tile": valid_tiles,
                "lat": lat.values[valid_tiles],
                "lon": lon.values[valid_tiles],
                "snow_abs_netpack_mam": snow_abs.values[year_index, valid_tiles],
                "snow_net_mam": snow_net.values[year_index, valid_tiles],
                "ol_snomasland_mam": ol_snow.values[year_index, valid_tiles],
            }
            for key, values in responses.items():
                record[key] = values.values[year_index, valid_tiles]
            frames.append(pd.DataFrame(record))
        table = pd.concat(frames, ignore_index=True)
        metadata = {
            "n_static_tiles": int(mask.sum().item()),
            "n_original_tile_years": len(table),
            "first_year": min(years),
            "last_year": max(years),
        }
        expected_rows = metadata["n_static_tiles"] * len(years)
        if len(table) != expected_rows:
            raise AssertionError(f"Expected {expected_rows} original rows, found {len(table)}")
        return table, metadata
    finally:
        catch.close()
        for experiment in datasets.values():
            for dataset in experiment.values():
                dataset.close()


def binned_summary(df: pd.DataFrame, x: str, y: str, bins: int) -> pd.DataFrame:
    sample = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    labels = pd.qcut(sample[x], q=bins, duplicates="drop")
    grouped = sample.assign(bin=labels).groupby("bin", observed=True)
    result = grouped.agg(
        n=(y, "size"),
        x_mean=(x, "mean"),
        x_median=(x, "median"),
        x_min=(x, "min"),
        x_max=(x, "max"),
        y_mean=(y, "mean"),
        y_median=(y, "median"),
        y_q25=(y, lambda values: np.nanpercentile(values, 25)),
        y_q75=(y, lambda values: np.nanpercentile(values, 75)),
        y_std=(y, "std"),
    ).reset_index()
    result["y_iqr"] = result["y_q75"] - result["y_q25"]
    result["y_se"] = result["y_std"] / np.sqrt(result["n"])
    result["bin_label"] = result["bin"].astype(str)
    return result.drop(columns="bin")


def reproduce_original_figure2(
    table: pd.DataFrame,
    config: dict,
    existing_path: Path,
) -> pd.DataFrame:
    frames = []
    for response in config["responses"]:
        work = table.assign(response_magnitude=np.abs(table[response]))
        summary = binned_summary(
            work,
            "snow_abs_netpack_mam",
            "response_magnitude",
            config["magnitude_bins"],
        )
        summary["response"] = response
        summary["y_metric"] = ORIGINAL_METRICS[response]
        frames.append(summary)
    reproduced = pd.concat(frames, ignore_index=True)
    if existing_path.exists():
        existing = pd.read_csv(existing_path)
        for response, metric in ORIGINAL_METRICS.items():
            old = existing.loc[existing["y_metric"] == metric].sort_values("x_mean")
            new = reproduced.loc[reproduced["response"] == response].sort_values("x_mean")
            if len(old) != len(new):
                raise AssertionError(f"Original Figure 2 bin count changed for {response}")
            for column in ("n", "x_mean", "y_mean", "y_std"):
                if not np.allclose(old[column], new[column], rtol=2.0e-6, atol=1.0e-8):
                    raise AssertionError(
                        f"Original Figure 2 reproduction failed for {response}/{column}"
                    )
    return reproduced


def restricted_sample(table: pd.DataFrame, response: str, minimum_years: int) -> pd.DataFrame:
    required = [
        "tile",
        "year",
        "lat",
        "lon",
        "snow_abs_netpack_mam",
        "snow_net_mam",
        "ol_snomasland_mam",
        response,
    ]
    sample = table[required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    counts = sample.groupby("tile")["year"].transform("size")
    sample = sample.loc[counts >= minimum_years].copy()
    if sample.empty:
        raise ValueError(f"No restricted sample remains for {response}")
    sample["response_signed"] = sample[response]
    sample["response_magnitude"] = sample[response].abs()
    for column in (
        "snow_abs_netpack_mam",
        "snow_net_mam",
        "ol_snomasland_mam",
        "response_signed",
        "response_magnitude",
    ):
        sample[f"{column}_anom"] = sample[column] - sample.groupby("tile")[column].transform("mean")
    return sample


def design_matrix(sample: pd.DataFrame, spec: ModelSpec) -> np.ndarray:
    columns = [sample[column].to_numpy(dtype="float64") for column in spec.x_columns]
    if any(column == "year_fe" for column in spec.x_columns):
        raise AssertionError("year_fe is expanded separately")
    matrix = np.column_stack(columns) if columns else np.empty((len(sample), 0))
    if spec.name in {"M2", "M3"}:
        dummies = pd.get_dummies(
            pd.Categorical(sample["year"], categories=sorted(sample["year"].unique())),
            drop_first=True,
            dtype=float,
        ).to_numpy()
        matrix = np.column_stack([matrix, dummies])
    if spec.intercept:
        matrix = np.column_stack([np.ones(len(sample)), matrix])
    return matrix


def model_specs(formulation: str) -> dict[str, ModelSpec]:
    if formulation == "magnitude":
        raw_x = "snow_abs_netpack_mam"
        anom_x = "snow_abs_netpack_mam_anom"
    elif formulation == "signed":
        raw_x = "snow_net_mam"
        anom_x = "snow_net_mam_anom"
    else:
        raise ValueError(formulation)
    return {
        "M0": ModelSpec("M0", (raw_x,), 1, True),
        "M1": ModelSpec("M1", (anom_x,), 0, False),
        "M2": ModelSpec("M2", (anom_x,), 1, True),
        "M3": ModelSpec("M3", (anom_x, "ol_snomasland_mam_anom"), 1, True),
    }


def response_column(formulation: str) -> str:
    return "response_magnitude_anom" if formulation == "magnitude" else "response_signed_anom"


def model_response_column(model: str, formulation: str) -> str:
    if model == "M0":
        return "response_magnitude" if formulation == "magnitude" else "response_signed"
    return response_column(formulation)


def fit_coefficient(sample: pd.DataFrame, spec: ModelSpec, formulation: str) -> float:
    x = design_matrix(sample, spec)
    y = sample[model_response_column(spec.name, formulation)].to_numpy(dtype="float64")
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return float(beta[spec.coefficient_index])


def between_coefficient(sample: pd.DataFrame, formulation: str) -> float:
    x_name = "snow_abs_netpack_mam" if formulation == "magnitude" else "snow_net_mam"
    y_name = "response_magnitude" if formulation == "magnitude" else "response_signed"
    means = sample.groupby("tile").agg(x=(x_name, "mean"), y=(y_name, "mean"), n=(y_name, "size"))
    x = np.column_stack([np.ones(len(means)), means["x"].to_numpy()])
    y = means["y"].to_numpy()
    root_weight = np.sqrt(means["n"].to_numpy(dtype="float64"))
    beta, *_ = np.linalg.lstsq(x * root_weight[:, None], y * root_weight, rcond=None)
    return float(beta[1])


def block_codes(sample: pd.DataFrame, block_degrees: float) -> np.ndarray:
    lon_bin = np.floor((sample["lon"].to_numpy() + 180.0) / block_degrees).astype(int)
    lat_bin = np.floor((sample["lat"].to_numpy() + 90.0) / block_degrees).astype(int)
    _, codes = np.unique(np.column_stack([lon_bin, lat_bin]), axis=0, return_inverse=True)
    return codes


def block_sufficient_statistics(
    x: np.ndarray,
    y: np.ndarray,
    codes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n_blocks = int(codes.max()) + 1
    n_columns = x.shape[1]
    xtx = np.empty((n_blocks, n_columns, n_columns), dtype="float64")
    xty = np.empty((n_blocks, n_columns), dtype="float64")
    for first in range(n_columns):
        xty[:, first] = np.bincount(
            codes, weights=x[:, first] * y, minlength=n_blocks
        )
        for second in range(first, n_columns):
            values = np.bincount(
                codes,
                weights=x[:, first] * x[:, second],
                minlength=n_blocks,
            )
            xtx[:, first, second] = values
            xtx[:, second, first] = values
    return xtx, xty


def bootstrap_models(
    sample: pd.DataFrame,
    formulation: str,
    block_degrees: float,
    replicates: int,
    seed: int,
) -> dict[str, np.ndarray]:
    specs = model_specs(formulation)
    codes = block_codes(sample, block_degrees)
    n_blocks = int(codes.max()) + 1
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=replicates)
    output = {}
    for name, spec in specs.items():
        x = design_matrix(sample, spec)
        y = sample[model_response_column(name, formulation)].to_numpy(dtype="float64")
        xtx, xty = block_sufficient_statistics(x, y, codes)
        coefficients = np.full(replicates, np.nan)
        for draw_index, weights in enumerate(draws):
            lhs = np.tensordot(weights, xtx, axes=(0, 0))
            rhs = np.tensordot(weights, xty, axes=(0, 0))
            beta, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
            coefficients[draw_index] = beta[spec.coefficient_index]
        output[name] = coefficients
    return output


def confidence_interval(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan, np.nan
    return tuple(np.percentile(finite, [2.5, 97.5]))


def tilewise_correlations(
    sample: pd.DataFrame,
    formulation: str,
    minimum_sd: float,
) -> np.ndarray:
    x = "snow_abs_netpack_mam" if formulation == "magnitude" else "snow_net_mam"
    y = "response_magnitude" if formulation == "magnitude" else "response_signed"
    work = sample[["tile", x, y]].copy()
    work["xy"] = work[x] * work[y]
    work["x2"] = work[x] ** 2
    work["y2"] = work[y] ** 2
    grouped = work.groupby("tile").agg(
        n=(x, "size"),
        sum_x=(x, "sum"),
        sum_y=(y, "sum"),
        sum_xy=("xy", "sum"),
        sum_x2=("x2", "sum"),
        sum_y2=("y2", "sum"),
        sd_x=(x, "std"),
    )
    numerator = grouped["n"] * grouped["sum_xy"] - grouped["sum_x"] * grouped["sum_y"]
    denominator = np.sqrt(
        (grouped["n"] * grouped["sum_x2"] - grouped["sum_x"] ** 2)
        * (grouped["n"] * grouped["sum_y2"] - grouped["sum_y"] ** 2)
    )
    corr = numerator / denominator.where(denominator > 0)
    corr = corr.where(grouped["sd_x"] >= minimum_sd).dropna()
    return corr.to_numpy(dtype="float64")


def variation_diagnostics(
    table: pd.DataFrame,
    predictor: str,
    minimum_years: int,
    minimum_sd: float,
) -> dict:
    sample = table[["tile", "year", predictor]].replace([np.inf, -np.inf], np.nan).dropna()
    counts = sample.groupby("tile")["year"].transform("size")
    sample = sample.loc[counts >= minimum_years].copy()
    grouped = sample.groupby("tile")[predictor]
    sd = grouped.std(ddof=1)
    means = grouped.mean()
    grand = sample[predictor].mean()
    residual = sample[predictor] - sample.groupby("tile")[predictor].transform("mean")
    total_ss = float(np.square(sample[predictor] - grand).sum())
    within_ss = float(np.square(residual).sum())
    tile_counts = grouped.size()
    between_ss = float((tile_counts * np.square(means - grand)).sum())
    return {
        "predictor": predictor,
        "n_tiles": int(sd.size),
        "n_tile_years": int(len(sample)),
        "within_sd_median": float(sd.median()),
        "within_sd_q25": float(sd.quantile(0.25)),
        "within_sd_q75": float(sd.quantile(0.75)),
        "fraction_below_minimum_sd": float((sd < minimum_sd).mean()),
        "fraction_at_or_above_minimum_sd": float((sd >= minimum_sd).mean()),
        "within_variance_fraction": within_ss / total_ss,
        "between_variance_fraction": between_ss / total_ss,
        "variance_identity_error": abs((within_ss + between_ss) / total_ss - 1.0),
    }


def bootstrap_binned_summary(
    sample: pd.DataFrame,
    formulation: str,
    bins: int,
    block_degrees: float,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    x = "snow_abs_netpack_mam_anom" if formulation == "magnitude" else "snow_net_mam_anom"
    y = response_column(formulation)
    work = sample[["lat", "lon", x, y]].dropna().copy()
    work["bin"] = pd.qcut(work[x], q=bins, labels=False, duplicates="drop")
    codes = block_codes(work, block_degrees)
    n_blocks = int(codes.max()) + 1
    n_bins_actual = int(work["bin"].max()) + 1
    combined = codes * n_bins_actual + work["bin"].to_numpy(dtype=int)
    length = n_blocks * n_bins_actual
    count = np.bincount(combined, minlength=length).reshape(n_blocks, n_bins_actual)
    xsum = np.bincount(
        combined, weights=work[x].to_numpy(), minlength=length
    ).reshape(n_blocks, n_bins_actual)
    ysum = np.bincount(
        combined, weights=work[y].to_numpy(), minlength=length
    ).reshape(n_blocks, n_bins_actual)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=replicates)
    draw_count = draws @ count
    draw_y = np.divide(
        draws @ ysum,
        draw_count,
        out=np.full(draw_count.shape, np.nan, dtype="float64"),
        where=draw_count > 0,
    )
    rows = []
    for bin_index in range(n_bins_actual):
        keep = work["bin"] == bin_index
        low, high = confidence_interval(draw_y[:, bin_index])
        rows.append(
            {
                "bin": bin_index + 1,
                "n": int(keep.sum()),
                "x_mean": float(work.loc[keep, x].mean()),
                "y_mean": float(work.loc[keep, y].mean()),
                "y_ci_low": low,
                "y_ci_high": high,
            }
        )
    return pd.DataFrame(rows)


def analyze_response(
    table: pd.DataFrame,
    response: str,
    config: dict,
    response_index: int,
) -> tuple[list[dict], dict[str, pd.DataFrame], dict[str, np.ndarray]]:
    sample = restricted_sample(table, response, config["minimum_valid_years"])
    boot = config["bootstrap"]
    min_sd = config["minimum_within_tile_sd_kg_m2"]
    rows = []
    bins = {}
    correlations = {}
    for form_index, formulation in enumerate(("magnitude", "signed")):
        specs = model_specs(formulation)
        points = {name: fit_coefficient(sample, spec, formulation) for name, spec in specs.items()}
        points["between"] = between_coefficient(sample, formulation)
        seed = boot["seed"] + 100 * response_index + 10 * form_index
        primary_draws = bootstrap_models(
            sample,
            formulation,
            boot["primary_block_degrees"],
            boot["replicates"],
            seed,
        )
        coarse_draws = bootstrap_models(
            sample,
            formulation,
            boot["coarse_block_degrees"],
            boot["replicates"],
            seed + 1,
        )
        x_anom = "snow_abs_netpack_mam_anom" if formulation == "magnitude" else "snow_net_mam_anom"
        lower = sample[x_anom].quantile(config["trim_lower_quantile"])
        upper = sample[x_anom].quantile(config["trim_upper_quantile"])
        trimmed = sample.loc[sample[x_anom].between(lower, upper)].copy()
        trimmed_points = {
            name: fit_coefficient(trimmed, spec, formulation) for name, spec in specs.items()
        }
        trimmed_draws = bootstrap_models(
            trimmed,
            formulation,
            boot["primary_block_degrees"],
            boot["replicates"],
            seed + 2,
        )
        pool_ci = confidence_interval(primary_draws["M0"])
        m3_ci = confidence_interval(primary_draws["M3"])
        coarse_ci = confidence_interval(coarse_draws["M3"])
        trim_ci = confidence_interval(trimmed_draws["M3"])
        ratio_draws = np.divide(
            primary_draws["M3"],
            primary_draws["M0"],
            out=np.full_like(primary_draws["M3"], np.nan),
            where=np.abs(primary_draws["M0"]) > 1.0e-12,
        )
        pooled_eligible = points["M0"] > 0 and pool_ci[0] > 0
        retained = points["M3"] / points["M0"] if pooled_eligible else np.nan
        retained_ci = confidence_interval(ratio_draws) if pooled_eligible else (np.nan, np.nan)
        corr = tilewise_correlations(sample, formulation, min_sd)
        correlations[formulation] = corr
        bins[formulation] = bootstrap_binned_summary(
            sample,
            formulation,
            config["magnitude_bins"] if formulation == "magnitude" else config["signed_bins"],
            boot["primary_block_degrees"],
            boot["replicates"],
            seed + 3,
        )
        rows.append(
            {
                "response": response,
                "formulation": formulation,
                "n_tiles": int(sample["tile"].nunique()),
                "n_tile_years": int(len(sample)),
                "pooled_beta_restricted": points["M0"],
                "pooled_ci_low": pool_ci[0],
                "pooled_ci_high": pool_ci[1],
                "between_beta_weighted_by_year_count": points["between"],
                "within_beta_m1": points["M1"],
                "within_year_fe_beta_m2": points["M2"],
                "within_year_fe_ol_snow_beta_m3": points["M3"],
                "m3_ci_low_5deg": m3_ci[0],
                "m3_ci_high_5deg": m3_ci[1],
                "m3_ci_low_10deg": coarse_ci[0],
                "m3_ci_high_10deg": coarse_ci[1],
                "retained_fraction": retained,
                "retained_fraction_ci_low": retained_ci[0],
                "retained_fraction_ci_high": retained_ci[1],
                "m3_beta_trim_1_99pct": trimmed_points["M3"],
                "m3_trim_ci_low": trim_ci[0],
                "m3_trim_ci_high": trim_ci[1],
                "tilewise_r_n": int(np.isfinite(corr).sum()),
                "tilewise_r_median": float(np.nanmedian(corr)),
                "tilewise_r_q25": float(np.nanpercentile(corr, 25)),
                "tilewise_r_q75": float(np.nanpercentile(corr, 75)),
                "tilewise_r_fraction_positive": float(np.mean(corr > 0)),
                "tilewise_r_fraction_negative": float(np.mean(corr < 0)),
                "outer_bin_center_x_range": float(bins[formulation]["x_mean"].iloc[-1] - bins[formulation]["x_mean"].iloc[0]),
                "ol_snow_change_m3_minus_m2": points["M3"] - points["M2"],
                "m4_status": "unavailable",
            }
        )
    return rows, bins, correlations


def original_sample_pooled(table: pd.DataFrame, response: str, formulation: str) -> float:
    x = "snow_abs_netpack_mam" if formulation == "magnitude" else "snow_net_mam"
    y_values = table[response].abs() if formulation == "magnitude" else table[response]
    sample = pd.DataFrame({"x": table[x], "y": y_values}).dropna()
    design = np.column_stack([np.ones(len(sample)), sample["x"]])
    beta, *_ = np.linalg.lstsq(design, sample["y"], rcond=None)
    return float(beta[1])


def classify(results: pd.DataFrame, variation: pd.DataFrame, config: dict) -> tuple[str, str]:
    mag_variation = variation.loc[
        variation["predictor"] == "snow_abs_netpack_mam"
    ].iloc[0]
    adequate = (
        mag_variation["fraction_at_or_above_minimum_sd"]
        >= config["adequate_variation_minimum_tile_fraction"]
        and mag_variation["within_variance_fraction"]
        >= config["adequate_variation_minimum_within_variance_fraction"]
    )
    if not adequate:
        return "D", "Within-tile snow-activity variation fails the predeclared adequacy thresholds."
    primary = config["classification"]["primary_responses"]
    mag = results[(results["formulation"] == "magnitude") & results["response"].isin(primary)]
    signed = results[(results["formulation"] == "signed") & results["response"].isin(primary)]
    needed = config["classification"]["minimum_responses"]
    positive = int((mag["within_year_fe_ol_snow_beta_m3"] > 0).sum())
    excludes_zero = int((mag["m3_ci_low_5deg"] > 0).sum())
    retains = int(
        (mag["retained_fraction"] >= config["classification"]["survives_retained_fraction"])
        .fillna(False)
        .sum()
    )
    signed_positive = int((signed["within_year_fe_ol_snow_beta_m3"] > 0).sum())
    if positive >= needed and excludes_zero >= needed and retains >= needed and signed_positive >= needed:
        return "A", "At least three primary responses satisfy every predeclared survival criterion."
    null_ci = mag["m3_ci_low_5deg"].le(0) & mag["m3_ci_high_5deg"].ge(0)
    weak_or_reverse = (
        mag["within_year_fe_ol_snow_beta_m3"].le(0)
        | mag["retained_fraction"].lt(config["classification"]["does_not_survive_retained_fraction"])
    )
    if int((null_ci & weak_or_reverse).sum()) >= needed:
        return "C", "At least three primary responses satisfy the predeclared non-survival rule."
    return "B", "The adequately identified result lies between the predeclared survival and non-survival rules."


def add_response_evidence(results: pd.DataFrame, config: dict) -> pd.DataFrame:
    output = results.copy()
    output["response_evidence"] = "not classified"
    primary = set(config["classification"]["primary_responses"])
    for response in primary:
        keep = (output["response"] == response) & (output["formulation"] == "magnitude")
        row = output.loc[keep].iloc[0]
        if (
            row["within_year_fe_ol_snow_beta_m3"] > 0
            and row["m3_ci_low_5deg"] > 0
            and row["retained_fraction"] >= config["classification"]["survives_retained_fraction"]
        ):
            label = "survives"
        elif (
            row["m3_ci_low_5deg"] <= 0 <= row["m3_ci_high_5deg"]
            and (
                row["within_year_fe_ol_snow_beta_m3"] <= 0
                or row["retained_fraction"] < config["classification"]["does_not_survive_retained_fraction"]
            )
        ):
            label = "does not survive"
        else:
            label = "mixed"
        output.loc[output["response"] == response, "response_evidence"] = label
    return output


def plot_binned(
    binned: dict[str, pd.DataFrame],
    formulation: str,
    config: dict,
    output_path: Path,
) -> None:
    responses = list(config["responses"])
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.2))
    for ax, response in zip(axes.flat, responses):
        data = binned[response]
        yerr = np.vstack([data["y_mean"] - data["y_ci_low"], data["y_ci_high"] - data["y_mean"]])
        ax.errorbar(data["x_mean"], data["y_mean"], yerr=yerr, marker="o", capsize=2, color="#15616d")
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.axvline(0, color="0.35", linewidth=0.8)
        ax.grid(alpha=0.22)
        ax.set_title(RESPONSE_LABELS[response])
        units = config["responses"][response]["units"]
        ylabel = f"Response anomaly\n[{units}]"
        if formulation == "magnitude":
            ylabel = f"Absolute-response anomaly\n[{units}]"
        ax.set_ylabel(ylabel, fontsize=9, labelpad=2)
        ax.set_xlabel(
            "MAM snow activity anomaly\n[kg m-2]"
            if formulation == "magnitude"
            else "MAM signed snow increment anomaly\n[kg m-2]"
        )
        x_range = data["x_mean"].iloc[-1] - data["x_mean"].iloc[0]
        ax.text(0.02, 0.98, f"outer-bin x range={x_range:.2f}\nN={data['n'].sum():,}", transform=ax.transAxes, va="top", fontsize=8)
    title = (
        "Analysis A robustness: within-tile snow-activity anomalies"
        if formulation == "magnitude"
        else "Analysis A robustness: signed snow-increment anomalies"
    )
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96), w_pad=2.4, h_pad=2.0)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_correlations(
    correlations: dict[str, dict[str, np.ndarray]],
    output_path: Path,
) -> None:
    responses = list(correlations)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, formulation in zip(axes, ("magnitude", "signed")):
        values = [correlations[response][formulation] for response in responses]
        boxes = ax.boxplot(values, tick_labels=[RESPONSE_LABELS[r] for r in responses], showfliers=False, patch_artist=True)
        for patch in boxes["boxes"]:
            patch.set_facecolor("#8ecae6" if formulation == "magnitude" else "#ffb703")
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.grid(axis="y", alpha=0.22)
        ax.tick_params(axis="x", rotation=35)
        ax.set_title(f"{formulation.capitalize()} formulation")
        ax.set_ylabel("Tile-wise Pearson r")
    fig.suptitle("Analysis A robustness: six-year tile-wise correlation distributions", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_control_sequence(results: pd.DataFrame, output_path: Path) -> None:
    primary = results[(results["formulation"] == "magnitude") & results["response"].isin(["et", "total_runoff", "total_water", "rzmc"])].copy()
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    models = [
        ("within_beta_m1", "M1: within"),
        ("within_year_fe_beta_m2", "M2: + year FE"),
        ("within_year_fe_ol_snow_beta_m3", "M3: + OL snow"),
    ]
    x = np.arange(len(primary))
    width = 0.22
    for index, (column, label) in enumerate(models):
        ratio = primary[column] / primary["pooled_beta_restricted"]
        ax.bar(x + (index - 1) * width, ratio, width=width, label=label)
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.axhline(0.4, color="0.35", linewidth=0.8, linestyle="--", label="40% survival rule")
    ax.axhline(0.15, color="0.5", linewidth=0.8, linestyle=":", label="15% non-survival rule")
    ax.set_xticks(x, [RESPONSE_LABELS[value] for value in primary["response"]])
    ax.set_ylabel("Coefficient / restricted pooled coefficient")
    ax.set_title("Analysis A robustness: effect retained through sequential controls")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncols=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def format_number(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "N/A"
    return f"{value:.{digits}g}"


def write_report_fragment(
    results: pd.DataFrame,
    variation: pd.DataFrame,
    classification: str,
    reason: str,
    config: dict,
    path: Path,
) -> None:
    mag = results[results["formulation"] == "magnitude"].set_index("response")
    signed = results[results["formulation"] == "signed"].set_index("response")
    lines = [
        "## Analysis A robustness: spatial and temporal snow-climatology controls",
        "",
        "This falsification analysis preserves the original 2001-2006 MODIS-only Analysis A and asks whether its pooled snow-activity relationship survives removal of persistent tile climatology, common year effects, and within-tile variation in open-loop MAM snow mass. It remains a model-internal `DA - OL` analysis, not causal proof or independent hydrologic validation.",
        "",
        "The original eight-bin Figure 2 sample and statistics were reproduced numerically before controls were applied. Every response then used one response-specific complete-case sample restricted to tiles with at least four of the six years. Between-tile fits are weighted by each tile's valid-year count. Confidence intervals use 1,000 spatial-block bootstrap replicates with approximately 5 by 5 degree blocks; the 10-degree results are retained as a sensitivity. The local inputs do not contain a clean pre-assimilation MODIS availability count, so M4 was not fitted.",
        "",
        "| Response | N tiles | N tile-years | Pooled beta | Between beta | Within M1 | + year FE M2 | + OL snow M3 [95% CI] | Retained | Tile-wise r median [IQR] | Signed M3 | 1-99% trim M3 | Evidence |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for response in config["responses"]:
        row = mag.loc[response]
        signed_row = signed.loc[response]
        retained = "N/A" if not np.isfinite(row["retained_fraction"]) else f"{100 * row['retained_fraction']:.0f}%"
        lines.append(
            "| " + " | ".join(
                [
                    RESPONSE_LABELS[response],
                    f"{int(row['n_tiles']):,}",
                    f"{int(row['n_tile_years']):,}",
                    format_number(row["pooled_beta_restricted"]),
                    format_number(row["between_beta_weighted_by_year_count"]),
                    format_number(row["within_beta_m1"]),
                    format_number(row["within_year_fe_beta_m2"]),
                    f"{format_number(row['within_year_fe_ol_snow_beta_m3'])} [{format_number(row['m3_ci_low_5deg'])}, {format_number(row['m3_ci_high_5deg'])}]",
                    retained,
                    f"{format_number(row['tilewise_r_median'])} [{format_number(row['tilewise_r_q25'])}, {format_number(row['tilewise_r_q75'])}]",
                    format_number(signed_row["within_year_fe_ol_snow_beta_m3"]),
                    format_number(row["m3_beta_trim_1_99pct"]),
                    str(row["response_evidence"]),
                ]
            ) + " |"
        )
    primary = config["classification"]["primary_responses"]
    primary_mag = mag.loc[primary]
    m2_to_m3_percent = 100.0 * (
        primary_mag["within_year_fe_ol_snow_beta_m3"]
        / primary_mag["within_year_fe_beta_m2"]
        - 1.0
    )
    positive_fraction_low = min(
        mag["tilewise_r_fraction_positive"].min(),
        signed["tilewise_r_fraction_positive"].min(),
    )
    positive_fraction_high = max(
        mag["tilewise_r_fraction_positive"].max(),
        signed["tilewise_r_fraction_positive"].max(),
    )
    lines.extend(
        [
            "",
            "Coefficient units are response units per `kg m-2` of the corresponding MAM snow-increment predictor. The signed M3 column uses signed `snow_net` and signed responses; the other coefficient columns use `snow_abs_netpack` and absolute response magnitude.",
            "",
            "All six selected fields were finite throughout the static mask, so the >=4-year restriction retained the original 48,067 tiles and 288,402 tile-years. Consequently, the restricted and original-sample pooled coefficients are identical here; both are retained separately in the output table so a future input change cannot silently conflate them.",
            "",
            "Adding OL MAM snow amount changes the primary M2 coefficients only modestly: "
            + ", ".join(
                f"{RESPONSE_LABELS[response]} {m2_to_m3_percent.loc[response]:+.1f}%"
                for response in primary
            )
            + ". Thus, year-to-year background snow amount does not explain away the within-tile relationship. All primary M3 intervals also remain above zero with 10-degree blocks and after trimming predictor anomalies below the 1st or above the 99th percentile.",
            "",
            f"Across the six pathways, {100 * positive_fraction_low:.1f}-{100 * positive_fraction_high:.1f}% of eligible tiles have positive tile-wise correlations, depending on response and formulation. The medians and IQRs in the table summarize their full distributions; no individual six-year tile correlation is assigned significance.",
            "",
            "### Identifying variation",
            "",
            "The practical near-zero threshold was fixed before response fitting at a within-tile predictor SD of `0.1 kg m-2`. Adequate identification required at least half the eligible tiles to exceed this threshold and within-tile variation to account for at least 5% of total predictor variance.",
            "",
            "| Predictor | Tiles | Tile-years | Within SD median [IQR] | Below threshold | Within / total variance | Between / total variance |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in variation.iterrows():
        label = "snow activity" if row["predictor"] == "snow_abs_netpack_mam" else "signed snow increment"
        lines.append(
            f"| {label} | {int(row['n_tiles']):,} | {int(row['n_tile_years']):,} | "
            f"{row['within_sd_median']:.2f} [{row['within_sd_q25']:.2f}, {row['within_sd_q75']:.2f}] | "
            f"{100 * row['fraction_below_minimum_sd']:.1f}% | "
            f"{100 * row['within_variance_fraction']:.1f}% | {100 * row['between_variance_fraction']:.1f}% |"
        )
    lines.extend(
        [
            "",
            "### Controlled diagnostics",
            "",
            "![Within-tile absolute-activity diagnostics](monthly_synthesis_report_figures/analysisA_robustness_binned_magnitude.png)",
            "",
            "The points are equal-count within-tile anomaly bins; bars are 95% spatial-block bootstrap intervals. Each panel reports the outer-bin center range so a visually flat relation can be interpreted against the available predictor contrast.",
            "",
            "![Within-tile signed diagnostics](monthly_synthesis_report_figures/analysisA_robustness_binned_signed.png)",
            "",
            "The signed version receives equal statistical treatment and greater interpretive weight where it conflicts with the absolute-magnitude result, because absolute responses can co-scale with background process variance.",
            "",
            "![Sequential controls](monthly_synthesis_report_figures/analysisA_robustness_control_sequence.png)",
            "",
            "This figure makes the M2-to-M3 change explicit: M2 removes persistent tile effects and common year effects; M3 additionally asks whether snow-DA activity predicts response after accounting for whether OL itself was snowier than usual in that tile-year.",
            "",
            "![Tile-wise correlation distributions](monthly_synthesis_report_figures/analysisA_robustness_tilewise_correlations.png)",
            "",
            "Individual six-year tile correlations are descriptive only; no tile-level significance is assigned. The 1st-99th percentile predictor-anomaly trim and the coarser 10-degree bootstrap are included in the machine-readable results table.",
            "",
            f"### Predeclared classification: {classification}",
            "",
            reason,
            "",
            "This classification is scoped to six MODIS-only seasons (2001-2006), the existing Northern Hemisphere seasonal-snow mask, and modeled `DA - OL` responses. Even a surviving relationship supports a physically coherent modeled propagation pathway under these controls; it does not establish a causal hydrologic effect through independent observations.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    config = load_config(args.config, args.bootstrap_replicates)
    root = repo_root()
    output_dir = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/analysisA_robustness"
    figure_dir = output_dir / "figures"
    report_figure_dir = root / "projects/M21C_ls/docs/monthly_synthesis_report_figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_figure_dir.mkdir(parents=True, exist_ok=True)

    print("Building the unchanged Analysis A tile-year sample...")
    table, metadata = build_original_table(config)
    original_bins_path = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/analysisA_binned_snow_activity_to_response_magnitude.csv"
    reproduced = reproduce_original_figure2(table, config, original_bins_path)
    reproduced.to_csv(output_dir / "original_figure2_reproduction.csv", index=False)
    (output_dir / "original_sample_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Original Figure 2 reproduced: {metadata['n_static_tiles']:,} tiles, {len(table):,} tile-years")

    variation = pd.DataFrame(
        [
            variation_diagnostics(
                table,
                predictor,
                config["minimum_valid_years"],
                config["minimum_within_tile_sd_kg_m2"],
            )
            for predictor in ("snow_abs_netpack_mam", "snow_net_mam")
        ]
    )
    identity_error = float(variation["variance_identity_error"].max())
    if identity_error > 1.0e-6:
        raise AssertionError(
            "Within/between variance decomposition does not recover total variance: "
            f"relative error={identity_error:.3e}"
        )
    variation.to_csv(output_dir / "identifying_variation.csv", index=False)

    rows = []
    all_bins: dict[str, dict[str, pd.DataFrame]] = {"magnitude": {}, "signed": {}}
    all_correlations: dict[str, dict[str, np.ndarray]] = {}
    for response_index, response in enumerate(config["responses"]):
        print(f"Analyzing {response} ({response_index + 1}/{len(config['responses'])})...")
        response_rows, response_bins, correlations = analyze_response(
            table, response, config, response_index
        )
        for row in response_rows:
            row["pooled_beta_original_sample"] = original_sample_pooled(
                table, response, row["formulation"]
            )
        rows.extend(response_rows)
        all_correlations[response] = correlations
        for formulation, summary in response_bins.items():
            all_bins[formulation][response] = summary
            out = summary.assign(response=response, formulation=formulation)
            out.to_csv(output_dir / f"binned_{formulation}_{response}.csv", index=False)

    results = pd.DataFrame(rows)
    classification, reason = classify(results, variation, config)
    results = add_response_evidence(results, config)
    results["overall_classification"] = classification
    results.to_csv(output_dir / "analysis_a_robustness_results.csv", index=False)
    pd.DataFrame(
        [{"classification": classification, "reason": reason}]
    ).to_csv(output_dir / "analysis_a_robustness_classification.csv", index=False)

    figures = {
        "analysisA_robustness_binned_magnitude.png": lambda path: plot_binned(
            all_bins["magnitude"], "magnitude", config, path
        ),
        "analysisA_robustness_binned_signed.png": lambda path: plot_binned(
            all_bins["signed"], "signed", config, path
        ),
        "analysisA_robustness_tilewise_correlations.png": lambda path: plot_correlations(
            all_correlations, path
        ),
        "analysisA_robustness_control_sequence.png": lambda path: plot_control_sequence(
            results, path
        ),
    }
    for name, create in figures.items():
        output_path = figure_dir / name
        create(output_path)
        shutil.copy2(output_path, report_figure_dir / name)

    fragment = output_dir / "analysis_a_robustness_report_section.md"
    write_report_fragment(results, variation, classification, reason, config, fragment)
    print(f"Classification: {classification} - {reason}")
    print(f"Results: {output_dir}")
    print(f"Report fragment: {fragment}")


if __name__ == "__main__":
    main()
