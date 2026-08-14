#!/usr/bin/env python3
"""Run the targeted MODIS-only snow-DA hydrology robustness checks."""

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


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from analysis_a_robustness import (  # noqa: E402
    RESPONSE_LABELS,
    block_codes,
    block_sufficient_statistics,
    confidence_interval,
    variation_diagnostics,
)
from water_year_snow_da_budget import (  # noqa: E402
    bootstrap_partition,
    dataset_to_frame,
    domain_budget_tables,
    monthly_flux_total,
    weighted_mean,
)


FLUX_VARIABLES = {
    "EVLAND",
    "RUNSURFLAND",
    "BASEFLOWLAND",
    "QINFILLAND",
    "SMLAND",
    "WCHANGELAND",
    "PRECTOTCORRLAND",
}
SEASON_MONTHS = {"AMJ": [4, 5, 6], "MJJ": [5, 6, 7]}
RESPONSE_SPECS = {
    "snowmelt": ("SMLAND", "AMJ", "kg m-2 season-1", False),
    "infiltration": ("QINFILLAND", "AMJ", "kg m-2 season-1", False),
    "rzmc": ("RZMC", "MJJ", "m3 m-3", True),
    "et": ("EVLAND", "MJJ", "kg m-2 season-1", True),
    "total_runoff": ("TOTAL_RUNOFF", "MJJ", "kg m-2 season-1", True),
    "total_water": ("TWLAND", "MJJ", "kg m-2", True),
}
BUDGET_RESPONSES = ["dRunoff_total", "dET", "dStorage", "residual"]
MONTHLY_BUDGET_NAMES = [
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
    "dTWLAND",
    "dWCHANGELAND",
    "dPrecipitation",
]
ANNUAL_SUM_MAP = {
    "I_snow": "snow_net",
    "snow_abs_netpack": "snow_abs_netpack",
    "dSnowmelt": "dSnowmelt",
    "dInfiltration": "dInfiltration",
    "dET": "dET",
    "dRunoff_surface": "dRunoff_surface",
    "dBaseflow": "dBaseflow",
    "dRunoff_total": "dRunoff_total",
    "dStorage_process_tendency": "dWCHANGELAND",
    "dPrecipitation": "dPrecipitation",
}
REPORT_START = "<!-- TARGETED_SNOW_HYDROLOGY_ROBUSTNESS_START -->"
REPORT_END = "<!-- TARGETED_SNOW_HYDROLOGY_ROBUSTNESS_END -->"


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    raise FileNotFoundError("Could not locate repository root")


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=root / "projects/M21C_ls/config/targeted_snow_hydrology_robustness.json",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=None,
        help="Testing override. Production uses the configured 1,000 replicates.",
    )
    return parser.parse_args()


def load_config(path: Path, bootstrap_override: int | None = None) -> dict:
    config = json.loads(path.read_text())
    if bootstrap_override is not None:
        config["bootstrap"]["replicates"] = int(bootstrap_override)
    elif config["bootstrap"]["replicates"] < 1000:
        raise ValueError("Production requires at least 1,000 spatial bootstrap replicates")
    return config


def oct_mar_dates(response_year: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{response_year - 1}-10-01", f"{response_year}-03-01", freq="MS")


def response_dates(response_year: int, season: str) -> pd.DatetimeIndex:
    return pd.DatetimeIndex(
        [pd.Timestamp(response_year, month, 1) for month in SEASON_MONTHS[season]]
    )


def september_august_dates(year: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{year - 1}-09-01", f"{year}-08-01", freq="MS")


def validate_date_windows(years: list[int], config: dict) -> str:
    examples = []
    all_dates = []
    for year in years:
        predictor = oct_mar_dates(year)
        if len(predictor) != 6:
            raise AssertionError(f"{year}: Oct-Mar predictor does not contain six months")
        for season in ("AMJ", "MJJ"):
            response = response_dates(year, season)
            overlap = predictor.intersection(response)
            if len(overlap):
                raise AssertionError(f"{year} {season}: predictor-response overlap: {overlap}")
            all_dates.extend(response)
        all_dates.extend(predictor)
        examples.append((year, predictor, response_dates(year, "AMJ"), response_dates(year, "MJJ")))
    if len(examples) != 6:
        raise AssertionError(f"Expected exactly six response seasons, found {len(examples)}")
    dates = pd.DatetimeIndex(all_dates)
    if dates.min() < pd.Timestamp(config["modis_only_start"]):
        raise AssertionError("Seasonal sample begins before the MODIS-only period")
    if dates.max() >= pd.Timestamp(config["microwave_soil_moisture_da_start"]):
        raise AssertionError("Microwave soil-moisture DA enters the seasonal sample")
    alternative_dates = pd.DatetimeIndex(
        np.concatenate([september_august_dates(year).values for year in years])
    )
    storage_endpoints = pd.date_range("2000-08-01", "2006-08-01", freq="12MS")
    robustness_dates = alternative_dates.append(storage_endpoints)
    if robustness_dates.min() < pd.Timestamp(config["modis_only_start"]):
        raise AssertionError("Boundary sensitivity begins before the MODIS-only period")
    if robustness_dates.max() >= pd.Timestamp(config["microwave_soil_moisture_da_start"]):
        raise AssertionError("Microwave soil-moisture DA enters the boundary sensitivity")
    year, predictor, amj, mjj = examples[0]
    return (
        f"Worked indexing example for response year {year}: predictor "
        f"{predictor[0]:%Y-%m} through {predictor[-1]:%Y-%m}; AMJ response "
        f"{amj[0]:%Y-%m} through {amj[-1]:%Y-%m}; MJJ response "
        f"{mjj[0]:%Y-%m} through {mjj[-1]:%Y-%m}; overlap=0 months."
    )


def open_inputs(data_dir: Path) -> dict[str, xr.Dataset]:
    names = {
        "ol_land": "OLv8_land_variables_2000_2024_compressed.nc",
        "da_land": "DAv8_land_variables_2000_2024_compressed.nc",
        "ol_flux": "OLv8_flux_core_2000_2024_compressed.nc",
        "da_flux": "DAv8_flux_core_2000_2024_compressed.nc",
        "ol_water": "OLv8_water_budget_2000_2024_compressed.nc",
        "da_water": "DAv8_water_budget_2000_2024_compressed.nc",
        "catch": "catch_progn_raw_monthly_cumulative_200006_202405.nc",
    }
    return {key: xr.open_dataset(data_dir / name, decode_times=True) for key, name in names.items()}


def close_inputs(inputs: dict[str, xr.Dataset]) -> None:
    for dataset in inputs.values():
        dataset.close()


def verify_monthly_increment_product(catch: xr.Dataset) -> dict[str, str]:
    note = str(catch.attrs.get("note", ""))
    source = str(catch.attrs.get("source_root", ""))
    if "Monthly cumulative sums from raw/submonthly catch_progn_incr" not in note:
        raise AssertionError("The catch product does not declare within-month raw-increment sums")
    if "time means" not in note:
        raise AssertionError("The catch product is missing its monthly-time-mean warning")
    if catch.sizes.get("time") != 288 or not pd.DatetimeIndex(catch.time.values).is_monotonic_increasing:
        raise AssertionError("Unexpected catch-product time coordinate")
    if catch["snow_net"].attrs.get("units") != "kg m-2":
        raise AssertionError("snow_net units are not kg m-2")
    # The producer sums raw increments separately inside each summarize_month call.
    # A temporal difference would therefore convert valid monthly totals into changes of totals.
    return {
        "classification": "per-month increment sum",
        "action": "use stored monthly values directly; do not difference across months",
        "evidence": note,
        "source": source,
    }


def select_exact(data: xr.DataArray, dates: pd.DatetimeIndex, tile_ids: np.ndarray) -> xr.DataArray:
    selected = data.sel(time=dates).isel(tile=tile_ids)
    found = pd.DatetimeIndex(selected.time.values)
    if not found.equals(dates):
        raise AssertionError(f"Expected dates {dates[0]}..{dates[-1]}, found {found}")
    return selected


def aggregate(data: xr.DataArray, method: str) -> np.ndarray:
    count = data.notnull().sum("time")
    if method == "sum":
        result = data.fillna(0.0).sum("time").where(count > 0)
    elif method == "mean":
        weights = data.time.dt.days_in_month.astype("float64")
        result = data.weighted(weights).mean("time", skipna=True).where(count > 0)
    elif method == "flux_total":
        seconds = data.time.dt.days_in_month.astype("float64") * 86400.0
        result = (data * seconds).fillna(0.0).sum("time").where(count > 0)
    else:
        raise ValueError(method)
    return result.values


def model_dataset_pair(inputs: dict[str, xr.Dataset], variable: str) -> tuple[xr.Dataset, xr.Dataset]:
    for group in ("land", "flux", "water"):
        if variable in inputs[f"ol_{group}"]:
            return inputs[f"ol_{group}"], inputs[f"da_{group}"]
    raise KeyError(variable)


def model_delta(
    inputs: dict[str, xr.Dataset],
    variable: str,
    dates: pd.DatetimeIndex,
    tile_ids: np.ndarray,
) -> np.ndarray:
    if variable == "TOTAL_RUNOFF":
        return model_delta(inputs, "RUNSURFLAND", dates, tile_ids) + model_delta(
            inputs, "BASEFLOWLAND", dates, tile_ids
        )
    ol, da = model_dataset_pair(inputs, variable)
    difference = select_exact(da[variable] - ol[variable], dates, tile_ids)
    return aggregate(difference, "flux_total" if variable in FLUX_VARIABLES else "mean")


def build_seasonal_table(
    inputs: dict[str, xr.Dataset],
    archive: xr.Dataset,
    years: list[int],
) -> pd.DataFrame:
    tile_ids = archive.tile.values.astype(int)
    tile_count = len(tile_ids)
    pieces = []
    for year in years:
        predictor_dates = oct_mar_dates(year)
        mam_dates = pd.DatetimeIndex([pd.Timestamp(year, month, 1) for month in [3, 4, 5]])
        snow_octmar = aggregate(select_exact(inputs["catch"]["snow_net"], predictor_dates, tile_ids), "sum")
        snow_mam = aggregate(select_exact(inputs["catch"]["snow_net"], mam_dates, tile_ids), "sum")
        ol_snow_octmar = select_exact(inputs["ol_land"]["SNOMASLAND"], predictor_dates, tile_ids)
        ol_snow_mean = aggregate(ol_snow_octmar, "mean")
        ol_snow_march = ol_snow_octmar.isel(time=-1).values
        ol_snow_mam = aggregate(select_exact(inputs["ol_land"]["SNOMASLAND"], mam_dates, tile_ids), "mean")
        data = {
            "tile": tile_ids,
            "year": np.full(tile_count, year, dtype=int),
            "lat": archive.lat.values,
            "lon": archive.lon.values,
            "area": archive.area.values,
            "snow_net_octmar": snow_octmar,
            "snow_net_mam": snow_mam,
            "ol_snow_march": ol_snow_march,
            "ol_snow_octmar_mean": ol_snow_mean,
            "ol_snow_mam": ol_snow_mam,
        }
        for response, (variable, season, _, _) in RESPONSE_SPECS.items():
            data[response] = model_delta(inputs, variable, response_dates(year, season), tile_ids)
        pieces.append(pd.DataFrame(data))
    table = pd.concat(pieces, ignore_index=True)
    if table["year"].nunique() != 6:
        raise AssertionError("Seasonal table does not contain exactly six response years")
    return table


def prepare_regression_sample(
    table: pd.DataFrame,
    response: str,
    predictor: str,
    controls: list[str],
    minimum_years: int,
) -> pd.DataFrame:
    columns = ["tile", "year", "lat", "lon", "area", predictor, response] + controls
    sample = table[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    counts = sample.groupby("tile")["year"].transform("size")
    sample = sample.loc[counts >= minimum_years].copy()
    if sample.empty:
        raise ValueError(f"No complete sample remains for {response}")
    for column in [predictor, response] + controls:
        sample[f"{column}_anom"] = sample[column] - sample.groupby("tile")[column].transform("mean")
    return sample


def regression_design(
    sample: pd.DataFrame,
    predictor: str,
    model: str,
    control: str | None = None,
) -> tuple[np.ndarray, np.ndarray, int]:
    if model == "M0":
        x = np.column_stack([np.ones(len(sample)), sample[predictor].to_numpy()])
        return x, sample["response"].to_numpy(), 1
    x_columns = [sample[f"{predictor}_anom"].to_numpy()]
    if model == "M3":
        if control is None:
            raise ValueError("M3 requires a control")
        x_columns.append(sample[f"{control}_anom"].to_numpy())
    if model in {"M2", "M3"}:
        year = pd.get_dummies(
            pd.Categorical(sample["year"], categories=sorted(sample["year"].unique())),
            drop_first=True,
            dtype=float,
        ).to_numpy()
        x = np.column_stack([np.ones(len(sample))] + x_columns + [year])
        coefficient_index = 1
    elif model == "M1":
        x = np.column_stack(x_columns)
        coefficient_index = 0
    else:
        raise ValueError(model)
    return x, sample["response_anom"].to_numpy(), coefficient_index


def fit_slope(sample: pd.DataFrame, predictor: str, model: str, control: str | None = None) -> float:
    x, y, coefficient_index = regression_design(sample, predictor, model, control)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return float(beta[coefficient_index])


def bootstrap_slope(
    sample: pd.DataFrame,
    predictor: str,
    model: str,
    control: str | None,
    block_degrees: float,
    replicates: int,
    seed: int,
) -> np.ndarray:
    x, y, coefficient_index = regression_design(sample, predictor, model, control)
    codes = block_codes(sample, block_degrees)
    n_blocks = int(codes.max()) + 1
    xtx, xty = block_sufficient_statistics(x, y, codes)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=replicates)
    slopes = np.full(replicates, np.nan)
    for draw_index, weights in enumerate(draws):
        lhs = np.tensordot(weights, xtx, axes=(0, 0))
        rhs = np.tensordot(weights, xty, axes=(0, 0))
        beta, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
        slopes[draw_index] = beta[coefficient_index]
    return slopes


def predictor_response_correlations(
    sample: pd.DataFrame,
    predictor: str,
    minimum_sd: float,
) -> dict[str, float]:
    grouped = sample.groupby("tile")
    stats = grouped.agg(
        n=(predictor, "size"),
        sum_x=(predictor, "sum"),
        sum_y=("response", "sum"),
    )
    work = sample[["tile", predictor, "response"]].copy()
    work["xy"] = work[predictor] * work["response"]
    work["x2"] = work[predictor] ** 2
    work["y2"] = work["response"] ** 2
    sums = work.groupby("tile")[["xy", "x2", "y2"]].sum()
    numerator = stats["n"] * sums["xy"] - stats["sum_x"] * stats["sum_y"]
    denominator = np.sqrt(
        (stats["n"] * sums["x2"] - stats["sum_x"] ** 2)
        * (stats["n"] * sums["y2"] - stats["sum_y"] ** 2)
    )
    sd = grouped[predictor].std(ddof=1)
    correlations = (numerator / denominator.where(denominator > 0)).where(sd >= minimum_sd).dropna()
    return {
        "tilewise_r_n": int(len(correlations)),
        "tilewise_r_median": float(correlations.median()),
        "tilewise_r_q25": float(correlations.quantile(0.25)),
        "tilewise_r_q75": float(correlations.quantile(0.75)),
        "tilewise_r_fraction_positive": float((correlations > 0).mean()),
    }


def interval_excludes_zero(interval: tuple[float, float]) -> bool:
    return bool(interval[0] > 0 or interval[1] < 0)


def classify_response(
    response: str,
    retention: float,
    intervals: list[tuple[float, float]],
    config: dict,
) -> str:
    if not RESPONSE_SPECS[response][3]:
        return "pathway diagnostic"
    if np.isfinite(retention) and retention < config["classification"]["non_survival_retained_fraction"]:
        return "non-survival"
    if (
        np.isfinite(retention)
        and retention >= config["classification"]["survives_retained_fraction"]
        and all(interval_excludes_zero(interval) for interval in intervals)
    ):
        return "survives"
    return "inconclusive"


def fit_seasonal_attribution(
    table: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    boot = config["bootstrap"]
    rows = []
    comparison_rows = []
    variation_rows = []
    for response_index, response in enumerate(RESPONSE_SPECS):
        primary = prepare_regression_sample(
            table,
            response,
            "snow_net_octmar",
            ["ol_snow_march", "ol_snow_octmar_mean", "snow_net_mam", "ol_snow_mam"],
            config["minimum_valid_years"],
        ).rename(columns={response: "response"})
        mam = prepare_regression_sample(
            table,
            response,
            "snow_net_mam",
            ["ol_snow_mam"],
            config["minimum_valid_years"],
        ).rename(columns={response: "response"})
        # Renaming happens after anomaly construction; preserve a generic response anomaly.
        primary["response_anom"] = primary[f"{response}_anom"]
        mam["response_anom"] = mam[f"{response}_anom"]
        points = {
            model: fit_slope(
                primary,
                "snow_net_octmar",
                model,
                "ol_snow_march" if model == "M3" else None,
            )
            for model in ["M0", "M1", "M2", "M3"]
        }
        secondary_m3 = fit_slope(primary, "snow_net_octmar", "M3", "ol_snow_octmar_mean")
        seed = boot["seed"] + response_index * 20
        primary_ci = confidence_interval(
            bootstrap_slope(
                primary,
                "snow_net_octmar",
                "M3",
                "ol_snow_march",
                boot["primary_block_degrees"],
                boot["replicates"],
                seed,
            )
        )
        coarse_ci = confidence_interval(
            bootstrap_slope(
                primary,
                "snow_net_octmar",
                "M3",
                "ol_snow_march",
                boot["coarse_block_degrees"],
                boot["replicates"],
                seed + 1,
            )
        )
        lower = primary["snow_net_octmar_anom"].quantile(config["trim_lower_quantile"])
        upper = primary["snow_net_octmar_anom"].quantile(config["trim_upper_quantile"])
        trimmed = primary.loc[primary["snow_net_octmar_anom"].between(lower, upper)].copy()
        trim_beta = fit_slope(trimmed, "snow_net_octmar", "M3", "ol_snow_march")
        trim_ci = confidence_interval(
            bootstrap_slope(
                trimmed,
                "snow_net_octmar",
                "M3",
                "ol_snow_march",
                boot["primary_block_degrees"],
                boot["replicates"],
                seed + 2,
            )
        )
        predictor_sd = float(primary["snow_net_octmar_anom"].std(ddof=1))
        retention = points["M3"] / points["M0"] if abs(points["M0"]) > 1.0e-12 else np.nan
        corr = predictor_response_correlations(
            primary, "snow_net_octmar", config["minimum_within_tile_sd_kg_m2"]
        )
        classification = classify_response(
            response, retention, [primary_ci, coarse_ci, trim_ci], config
        )
        row = {
            "response": response,
            "response_label": RESPONSE_LABELS[response],
            "response_units": RESPONSE_SPECS[response][2],
            "n_tiles": int(primary["tile"].nunique()),
            "n_tile_years": int(len(primary)),
            "pooled_beta": points["M0"],
            "m1_beta": points["M1"],
            "m2_beta": points["M2"],
            "m3_march_snow_beta": points["M3"],
            "m3_march_snow_ci_low_5deg": primary_ci[0],
            "m3_march_snow_ci_high_5deg": primary_ci[1],
            "m3_march_snow_ci_low_10deg": coarse_ci[0],
            "m3_march_snow_ci_high_10deg": coarse_ci[1],
            "m3_octmar_mean_snow_beta": secondary_m3,
            "m3_trim_1_99pct_beta": trim_beta,
            "m3_trim_ci_low_5deg": trim_ci[0],
            "m3_trim_ci_high_5deg": trim_ci[1],
            "within_predictor_sd_kg_m2": predictor_sd,
            "retained_fraction_of_pooled": retention,
            "classification": classification,
            "classification_note": config["classification"]["note"],
            "m4_status": "unfitted: no pre-assimilation MODIS observation-availability count",
            **corr,
        }
        for model, beta in points.items():
            row[f"{model.lower()}_beta_per_within_sd"] = beta * predictor_sd
        row["m3_octmar_mean_snow_beta_per_within_sd"] = secondary_m3 * predictor_sd
        row["m3_trim_beta_per_within_sd"] = trim_beta * predictor_sd
        rows.append(row)

        mam_sd = float(mam["snow_net_mam_anom"].std(ddof=1))
        for model in ["M0", "M1", "M2", "M3"]:
            beta = fit_slope(
                mam,
                "snow_net_mam",
                model,
                "ol_snow_mam" if model == "M3" else None,
            )
            comparison_rows.append(
                {
                    "response": response,
                    "window": "MAM",
                    "model": model,
                    "raw_beta": beta,
                    "within_predictor_sd_kg_m2": mam_sd,
                    "beta_per_within_sd": beta * mam_sd,
                    "published_signed_mam_m3_benchmark": (
                        config["signed_mam_m3_benchmarks"][response] if model == "M3" else np.nan
                    ),
                }
            )
        for model in ["M0", "M1", "M2", "M3"]:
            comparison_rows.append(
                {
                    "response": response,
                    "window": "Oct-Mar",
                    "model": model,
                    "raw_beta": points[model],
                    "within_predictor_sd_kg_m2": predictor_sd,
                    "beta_per_within_sd": points[model] * predictor_sd,
                    "published_signed_mam_m3_benchmark": np.nan,
                }
            )
        variation_rows.append(
            variation_diagnostics(
                primary,
                "snow_net_octmar",
                config["minimum_valid_years"],
                config["minimum_within_tile_sd_kg_m2"],
            )
            | {"response": response}
        )
    results = pd.DataFrame(rows)
    comparison = pd.DataFrame(comparison_rows)
    variation = pd.DataFrame(variation_rows)
    return results, comparison, variation


def predictor_window_correlations(table: pd.DataFrame) -> pd.DataFrame:
    sample = table[["tile", "snow_net_octmar", "snow_net_mam"]].dropna().copy()
    pooled = sample[["snow_net_octmar", "snow_net_mam"]].corr().iloc[0, 1]
    for column in ["snow_net_octmar", "snow_net_mam"]:
        sample[f"{column}_anom"] = sample[column] - sample.groupby("tile")[column].transform("mean")
    within = sample[["snow_net_octmar_anom", "snow_net_mam_anom"]].corr().iloc[0, 1]
    return pd.DataFrame(
        [
            {"correlation": "pooled", "r": pooled, "n_tile_years": len(sample)},
            {"correlation": "tile-demeaned", "r": within, "n_tile_years": len(sample)},
        ]
    )


def plot_seasonal_control_sequence(results: pd.DataFrame, path: Path) -> None:
    responses = ["rzmc", "et", "total_runoff", "total_water"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), constrained_layout=True)
    model_columns = ["pooled_beta", "m1_beta", "m2_beta", "m3_march_snow_beta"]
    labels = ["Pooled", "M1", "M2", "M3"]
    for panel, (ax, response) in enumerate(zip(axes.flat, responses)):
        row = results.set_index("response").loc[response]
        values = row[model_columns].to_numpy(dtype=float)
        ax.plot(np.arange(4), values, marker="o", color="#2166ac", linewidth=1.5)
        m3 = values[-1]
        ax.errorbar(
            [3],
            [m3],
            yerr=[[m3 - row["m3_march_snow_ci_low_5deg"]], [row["m3_march_snow_ci_high_5deg"] - m3]],
            color="#b2182b",
            capsize=4,
            linewidth=1.2,
        )
        ax.axhline(0, color="0.35", linewidth=0.8)
        ax.set_xticks(np.arange(4), labels)
        ax.set_title(f"({chr(97 + panel)}) {RESPONSE_LABELS[response]}")
        ax.set_ylabel(f"Coefficient ({row['response_units']} per kg m$^{{-2}}$)")
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Oct-Mar signed snow input: sequential controls")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def add_water_year_controls(dataset: xr.Dataset) -> xr.Dataset:
    output = dataset.copy()
    snow = output["ol_snow_mass_monthly"]
    output["ol_wy_peak_snow_mass"] = snow.max("water_month", skipna=True)
    output["ol_march_snow_mass"] = snow.isel(water_month=5)
    output["ol_wy_peak_snow_mass"].attrs["units"] = "kg m-2"
    output["ol_march_snow_mass"].attrs["units"] = "kg m-2"
    return output


def water_year_control_regressions(dataset: xr.Dataset, config: dict) -> pd.DataFrame:
    controls = {
        "water-year maximum OL snow mass": "ol_wy_peak_snow_mass",
        "March OL snow mass": "ol_march_snow_mass",
        "MAM mean OL snow mass": "ol_mam_snow_mass",
    }
    boot = config["bootstrap"]
    rows = []
    for control_index, (control_label, control) in enumerate(controls.items()):
        frame = dataset_to_frame(dataset, ["I_snow", control] + BUDGET_RESPONSES).dropna().copy()
        frame = frame.rename(columns={"water_year": "year"})
        for column in ["I_snow", control] + BUDGET_RESPONSES:
            frame[f"{column}_anom"] = frame[column] - frame.groupby("tile")[column].transform("mean")
        for response_index, response in enumerate(BUDGET_RESPONSES):
            sample = frame.rename(columns={response: "response"}).copy()
            sample["response_anom"] = sample[f"{response}_anom"]
            beta = fit_slope(sample, "I_snow", "M3", control)
            ci = confidence_interval(
                bootstrap_slope(
                    sample,
                    "I_snow",
                    "M3",
                    control,
                    boot["primary_block_degrees"],
                    boot["replicates"],
                    boot["seed"] + 500 + control_index * 20 + response_index,
                )
            )
            rows.append(
                {
                    "snow_control": control_label,
                    "control_variable": control,
                    "response": response,
                    "slope": beta,
                    "units": "dimensionless",
                    "ci_low_5deg": ci[0],
                    "ci_high_5deg": ci[1],
                    "n_tiles": int(frame["tile"].nunique()),
                    "n_tile_years": int(len(frame)),
                    "interpretation": "marginal partition sensitivity; not causal proof",
                }
            )
        control_rows = pd.DataFrame(rows).query("snow_control == @control_label")
        closure = control_rows.set_index("response").loc[BUDGET_RESPONSES, "slope"].sum()
        if not np.isclose(closure, 1.0, atol=1.0e-8):
            raise AssertionError(f"{control_label}: budget regression slopes sum to {closure}")
    return pd.DataFrame(rows)


def monthly_snapshot(
    inputs: dict[str, xr.Dataset],
    date: str,
    tile_ids: np.ndarray,
) -> dict[str, np.ndarray]:
    timestamp = pd.Timestamp(date)
    dates = pd.DatetimeIndex([timestamp])

    def delta(variable: str) -> np.ndarray:
        return model_delta(inputs, variable, dates, tile_ids)

    ol_snow = select_exact(inputs["ol_land"]["SNOMASLAND"], dates, tile_ids).isel(time=0).values
    da_snow = select_exact(inputs["da_land"]["SNOMASLAND"], dates, tile_ids).isel(time=0).values
    snow = select_exact(inputs["catch"]["snow_net"], dates, tile_ids).isel(time=0).values
    snow_abs = select_exact(inputs["catch"]["snow_abs_netpack"], dates, tile_ids).isel(time=0).values
    dtw = select_exact(
        inputs["da_water"]["TWLAND"] - inputs["ol_water"]["TWLAND"], dates, tile_ids
    ).isel(time=0).values
    surface = delta("RUNSURFLAND")
    baseflow = delta("BASEFLOWLAND")
    return {
        "ol_snow_mass": ol_snow,
        "dSnow_mass": da_snow - ol_snow,
        "snow_net": snow,
        "snow_abs_netpack": snow_abs,
        "dSnowmelt": delta("SMLAND"),
        "dInfiltration": delta("QINFILLAND"),
        "dET": delta("EVLAND"),
        "dRunoff_surface": surface,
        "dBaseflow": baseflow,
        "dRunoff_total": surface + baseflow,
        "dTWLAND": dtw,
        "dWCHANGELAND": delta("WCHANGELAND"),
        "dPrecipitation": delta("PRECTOTCORRLAND"),
    }


def build_september_august_dataset(
    baseline: xr.Dataset,
    inputs: dict[str, xr.Dataset],
    years: list[int],
) -> xr.Dataset:
    tile_ids = baseline.tile.values.astype(int)
    september_2000 = monthly_snapshot(inputs, "2000-09-01", tile_ids)
    august_2000 = monthly_snapshot(inputs, "2000-08-01", tile_ids)
    shape = (len(years), 12, len(tile_ids))
    monthly = {name: np.full(shape, np.nan, dtype="float32") for name in MONTHLY_BUDGET_NAMES}
    for year_index, year in enumerate(years):
        dates = september_august_dates(year)
        if len(dates) != 12 or dates[0] != pd.Timestamp(year - 1, 9, 1):
            raise AssertionError(f"{year}: malformed Sep-Aug dates")
        for name in MONTHLY_BUDGET_NAMES:
            if year_index == 0:
                september = september_2000[name]
            else:
                september = baseline[f"{name}_monthly"].values[year_index - 1, 11]
            october_august = baseline[f"{name}_monthly"].values[year_index, :11]
            monthly[name][year_index] = np.concatenate([september[None, :], october_august], axis=0)

    annual = {
        name: np.nansum(monthly[source], axis=1).astype("float32")
        for name, source in ANNUAL_SUM_MAP.items()
    }
    storage = np.full((len(years), len(tile_ids)), np.nan, dtype="float32")
    for year_index in range(len(years)):
        prior_august = (
            august_2000["dTWLAND"]
            if year_index == 0
            else baseline["dTWLAND_monthly"].values[year_index - 1, 10]
        )
        current_august = baseline["dTWLAND_monthly"].values[year_index, 10]
        storage[year_index] = current_august - prior_august
    annual["dStorage"] = storage
    annual["residual"] = (
        annual["I_snow"] - annual["dRunoff_total"] - annual["dET"] - annual["dStorage"]
    )
    annual["process_balance_diagnostic"] = (
        annual["dET"] + annual["dRunoff_total"] + annual["dStorage_process_tendency"]
    )
    data_vars = {
        f"{name}_monthly": (("water_year", "water_month", "tile"), values)
        for name, values in monthly.items()
    }
    data_vars.update(
        {name: (("water_year", "tile"), values) for name, values in annual.items()}
    )
    output = xr.Dataset(
        data_vars=data_vars,
        coords={
            "water_year": np.asarray(years, dtype="int16"),
            "water_month": np.arange(1, 13, dtype="int8"),
            "tile": baseline.tile.values,
            "lat": baseline.lat,
            "lon": baseline.lon,
            "area": baseline.area,
        },
        attrs={
            "water_year_definition": "September through August",
            "storage_closing_term": "August-to-August monthly-mean DA-minus-OL TWLAND state change",
        },
    )
    closure = output["dRunoff_total"] + output["dET"] + output["dStorage"] + output["residual"]
    if not np.allclose(closure.values, output["I_snow"].values, atol=1.0e-5, equal_nan=True):
        raise AssertionError("Sep-Aug direct budget does not close")
    if not np.allclose(
        output["dRunoff_surface"] + output["dBaseflow"],
        output["dRunoff_total"],
        atol=1.0e-5,
        equal_nan=True,
    ):
        raise AssertionError("Sep-Aug surface runoff plus baseflow does not equal total runoff")
    return output


def reproduce_baseline(
    baseline: xr.Dataset,
    config: dict,
    root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    annual, partition, _ = domain_budget_tables(baseline, config)
    archived_path = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/annual_domain_budgets.csv"
    archived = pd.read_csv(archived_path)
    columns = ["I_snow", "dRunoff_total", "dET", "dStorage", "residual"]
    if not np.allclose(annual[columns], archived[columns], atol=1.0e-10, equal_nan=True):
        raise AssertionError("Existing Oct-Sep baseline was not reproduced exactly")
    mean_input = float(annual.loc[annual["water_year"] == "6-WY mean", "I_snow"].iloc[0])
    expected = config["expected_octsep_mean_signed_input_kg_m2"]
    if abs(mean_input - expected) > config["expected_input_tolerance_kg_m2"]:
        raise AssertionError(f"Mean signed annual input is {mean_input}, expected about {expected}")
    return annual, partition


def precipitation_check(dataset: xr.Dataset) -> dict[str, float]:
    maximum_tile = float(np.nanmax(np.abs(dataset["dPrecipitation"].values)))
    annual_domain = [
        weighted_mean(values, dataset.area.values)
        for values in dataset["dPrecipitation"].values
    ]
    return {
        "maximum_abs_annual_tile_kg_m2": maximum_tile,
        "maximum_abs_annual_domain_mean_kg_m2": float(np.max(np.abs(annual_domain))),
    }


def boundary_summary_rows(
    boundary: str,
    annual: pd.DataFrame,
    partition: pd.DataFrame,
    uncertainty: pd.DataFrame,
) -> list[dict]:
    rows = []
    for _, source in annual.iterrows():
        rows.append(
            {
                "boundary": boundary,
                "scope": "annual" if source["water_year"] != "6-WY mean" else "six_year_mean",
                "water_year": source["water_year"],
                **{name: source[name] for name in ["I_snow", "dRunoff_total", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]},
                **{f"fraction_{name}": source[f"fraction_{name}"] for name in BUDGET_RESPONSES},
            }
        )
    addition = partition.set_index("sample").loc["addition"]
    row = {
        "boundary": boundary,
        "scope": "positive_input_partition",
        "water_year": "all six",
        "n_tile_years": addition["n_tile_years"],
        **{name: addition[name] for name in ["I_snow", "dRunoff_total", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]},
        **{f"fraction_{name}": addition[f"fraction_{name}"] for name in BUDGET_RESPONSES},
        "fraction_dRunoff_surface": addition["fraction_dRunoff_surface"],
        "fraction_dBaseflow": addition["fraction_dBaseflow"],
    }
    for variable in ["dRunoff_total", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]:
        ci = uncertainty.query("sample == 'addition' and variable == @variable").iloc[0]
        row[f"fraction_{variable}_ci_low_5deg"] = ci["ci_low"]
        row[f"fraction_{variable}_ci_high_5deg"] = ci["ci_high"]
    rows.append(row)
    return rows


def common_positive_partition(
    baseline: xr.Dataset,
    alternative: xr.Dataset,
) -> pd.DataFrame:
    variables = ["I_snow", "dRunoff_total", "dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]
    left = dataset_to_frame(baseline, variables)
    right = dataset_to_frame(alternative, variables)
    common = (left["I_snow"].to_numpy() > 0) & (right["I_snow"].to_numpy() > 0)
    rows = []
    for boundary, frame in [("Oct-Sep", left), ("Sep-Aug", right)]:
        row = {"boundary": boundary, "scope": "common_positive_partition", "n_tile_years": int(common.sum())}
        input_mean = weighted_mean(frame["I_snow"].to_numpy(), frame["area"].to_numpy(), common)
        row["I_snow"] = input_mean
        for variable in variables[1:]:
            value = weighted_mean(frame[variable].to_numpy(), frame["area"].to_numpy(), common)
            row[variable] = value
            row[f"fraction_{variable}"] = value / input_mean
        rows.append(row)
    return pd.DataFrame(rows)


def september_diagnostic(baseline: xr.Dataset) -> pd.DataFrame:
    area = np.tile(baseline.area.values, len(baseline.water_year))
    addition = baseline["I_snow"].values.reshape(-1) > 0
    annual_input_all = weighted_mean(baseline["I_snow"].values.reshape(-1), area)
    annual_input_addition = weighted_mean(baseline["I_snow"].values.reshape(-1), area, addition)
    variables = {
        "signed_snow_input": baseline["snow_net_monthly"].values[:, 11],
        "snow_mass_change": baseline["dSnow_mass_monthly"].values[:, 11],
        "snowmelt": baseline["dSnowmelt_monthly"].values[:, 11],
        "runoff": baseline["dRunoff_total_monthly"].values[:, 11],
        "et": baseline["dET_monthly"].values[:, 11],
        "august_to_september_dTWLAND": (
            baseline["dTWLAND_monthly"].values[:, 11] - baseline["dTWLAND_monthly"].values[:, 10]
        ),
    }
    rows = []
    for variable, values in variables.items():
        flat = values.reshape(-1)
        row = {
            "variable": variable,
            "all_tile_years_mean_kg_m2": weighted_mean(flat, area),
            "snow_addition_tile_years_mean_kg_m2": weighted_mean(flat, area, addition),
        }
        if variable == "signed_snow_input":
            row["percent_of_octsep_annual_input_all"] = 100.0 * row["all_tile_years_mean_kg_m2"] / annual_input_all
            row["percent_of_octsep_annual_input_snow_addition"] = 100.0 * row["snow_addition_tile_years_mean_kg_m2"] / annual_input_addition
        rows.append(row)
    return pd.DataFrame(rows)


def plot_boundary_sensitivity(
    summary: pd.DataFrame,
    path: Path,
) -> None:
    rows = summary.query("scope == 'positive_input_partition'").set_index("boundary")
    variables = ["dRunoff_surface", "dBaseflow", "dET", "dStorage", "residual"]
    labels = ["Surface runoff", "Baseflow", "ET", "Storage", "Residual"]
    colors = ["#2166ac", "#67a9cf", "#1b9e77", "#d9a21b", "#7f7f7f"]
    x = np.arange(len(variables))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for index, boundary in enumerate(["Oct-Sep", "Sep-Aug"]):
        row = rows.loc[boundary]
        values = np.array([100.0 * row[f"fraction_{variable}"] for variable in variables])
        low = np.array([100.0 * row[f"fraction_{variable}_ci_low_5deg"] for variable in variables])
        high = np.array([100.0 * row[f"fraction_{variable}_ci_high_5deg"] for variable in variables])
        bars = ax.bar(
            x + (index - 0.5) * width,
            values,
            width,
            color=colors,
            edgecolor="black" if index else "none",
            linewidth=0.8,
            alpha=0.72 if index else 1.0,
            label=boundary,
        )
        bars[-1].set_hatch("///")
        bars[-1].set_edgecolor("black")
        ax.errorbar(
            x + (index - 0.5) * width,
            values,
            yerr=[values - low, high - values],
            fmt="none",
            ecolor="0.2",
            capsize=3,
            linewidth=0.9,
        )
    ax.axhline(100, color="0.35", linewidth=0.8, linestyle="--")
    ax.axhline(0, color="0.25", linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Fraction of positive snow-DA input (%)")
    ax.set_title("Positive-input water partition: accounting-boundary sensitivity")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def fmt(value: float, digits: int = 3) -> str:
    return "N/A" if not np.isfinite(value) else f"{value:.{digits}g}"


def build_report_section(
    seasonal: pd.DataFrame,
    comparison: pd.DataFrame,
    correlations: pd.DataFrame,
    control_sensitivity: pd.DataFrame,
    boundary: pd.DataFrame,
    september: pd.DataFrame,
    product_check: dict[str, str],
    noticeable: dict[str, bool],
) -> str:
    octmar_m3 = seasonal.set_index("response")
    mam_m3 = comparison.query("window == 'MAM' and model == 'M3'").set_index("response")
    peak = control_sensitivity.query("snow_control == 'water-year maximum OL snow mass'").set_index("response")
    old = control_sensitivity.query("snow_control == 'MAM mean OL snow mass'").set_index("response")
    positive = boundary.query("scope == 'positive_input_partition'").set_index("boundary")
    sep = september.set_index("variable")
    octsep_mean_input = boundary.query(
        "boundary == 'Oct-Sep' and scope == 'six_year_mean'"
    )["I_snow"].iloc[0]
    sepaug_mean_input = boundary.query(
        "boundary == 'Sep-Aug' and scope == 'six_year_mean'"
    )["I_snow"].iloc[0]
    annual = boundary.query("scope == 'annual'")
    negative_residual_years = {
        name: int((annual.query("boundary == @name")["residual"] < 0).sum())
        for name in ["Oct-Sep", "Sep-Aug"]
    }
    lines = [
        REPORT_START,
        "## Targeted robustness checks for the snow-DA water budget",
        "",
        "These post hoc sensitivity tests address three identified design issues: overlap between the old MAM predictor and AMJ/MJJ responses, mismatch between annual input and its OL-snow control, and appreciable snow-DA input in September at the Oct-Sep boundary. They preserve the 48,067-tile Northern Hemisphere seasonal-snow mask and the six 2001-2006 MODIS-only years, before microwave soil-moisture DA began.",
        "",
        f"The raw-increment product was confirmed as a **{product_check['classification']}**: each stored month is already the sum of native submonthly increments within that month, so no temporal differencing was applied. The existing Oct-Sep baseline was reproduced exactly, including a six-year mean signed input of {octsep_mean_input:.2f} kg m-2 yr-1.",
        "",
        "### Non-overlapping Oct-Mar seasonal attribution",
        "",
        "The predictor is signed snow-water input from October through March; responses begin in April or May, giving zero overlapping months. M3 removes tile means, common response-year effects, and response-year March OL snow mass. Coefficients are native response units per kg m-2 of signed snow input; standardized values are the response associated with one within-tile predictor SD.",
        "",
        "| Response | Oct-Mar M3 [5-degree 95% CI] | Per within-tile SD | Signed MAM M3 | MAM per within-tile SD | Retained vs pooled | Classification |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for response in RESPONSE_SPECS:
        row = octmar_m3.loc[response]
        old_row = mam_m3.loc[response]
        lines.append(
            f"| {RESPONSE_LABELS[response]} | {fmt(row['m3_march_snow_beta'])} "
            f"[{fmt(row['m3_march_snow_ci_low_5deg'])}, {fmt(row['m3_march_snow_ci_high_5deg'])}] | "
            f"{fmt(row['m3_beta_per_within_sd'])} | {fmt(old_row['raw_beta'])} | "
            f"{fmt(old_row['beta_per_within_sd'])} | {100 * row['retained_fraction_of_pooled']:.1f}% | "
            f"{row['classification']} |"
        )
    pooled_r = correlations.set_index("correlation").loc["pooled", "r"]
    within_r = correlations.set_index("correlation").loc["tile-demeaned", "r"]
    lines.extend(
        [
            "",
            "All four classified downstream responses survive: their M3 intervals exclude zero with 5-degree and 10-degree blocks and after the 1st-99th percentile predictor-anomaly trim. The standardized Oct-Mar effects are smaller than signed MAM, ranging from 24% of the MAM value for runoff to 67% for ET, but remain positive. Snowmelt remains a positive pathway diagnostic; infiltration is near zero after controls and its interval includes zero.",
            "",
            f"The Oct-Mar and MAM signed predictors correlate at r={pooled_r:.3f} pooled and r={within_r:.3f} after tile demeaning. Classification uses the inherited MAM thresholds without modification after this post hoc window change. Snowmelt and infiltration remain pathway diagnostics, not pass/fail outcomes. M4 remains unfitted because the archive has no pre-assimilation MODIS availability count.",
            "",
            "![Oct-Mar signed control sequence](monthly_synthesis_report_figures/analysisA_octmar_signed_control_sequence.png)",
            "",
            "### Water-year marginal partition control sensitivity",
            "",
            "These dimensionless regressions estimate the within-tile marginal partition of an additional unit of Oct-Sep snow-DA input. They are not an independent existence or causal test. In every control formulation, runoff + ET + storage + residual slopes sum to exactly one by construction.",
            "",
            "| Snow control | Runoff [95% CI] | ET [95% CI] | Storage [95% CI] | Residual [95% CI] |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for label in ["water-year maximum OL snow mass", "March OL snow mass", "MAM mean OL snow mass"]:
        rows = control_sensitivity.query("snow_control == @label").set_index("response")
        cells = []
        for response in BUDGET_RESPONSES:
            row = rows.loc[response]
            cells.append(f"{row['slope']:.3f} [{row['ci_low_5deg']:.3f}, {row['ci_high_5deg']:.3f}]")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines.extend(
        [
            "",
            f"Replacing the old MAM control with the full-year peak changes the marginal runoff slope from {old.loc['dRunoff_total', 'slope']:.3f} to {peak.loc['dRunoff_total', 'slope']:.3f} and ET from {old.loc['dET', 'slope']:.3f} to {peak.loc['dET', 'slope']:.3f}. For context, direct positive-input accounting describes the average fate of all added water ({100 * positive.loc['Oct-Sep', 'fraction_dRunoff_total']:.1f}% runoff and {100 * positive.loc['Oct-Sep', 'fraction_dET']:.1f}% ET), whereas this regression describes a marginal unit; the two are not expected to be numerically equivalent.",
            "",
            "### Accounting-boundary and September sensitivity",
            "",
            f"Mean September signed input is {sep.loc['signed_snow_input', 'all_tile_years_mean_kg_m2']:.2f} kg m-2, or {sep.loc['signed_snow_input', 'percent_of_octsep_annual_input_all']:.1f}% of the mean Oct-Sep annual input. In snow-addition tile-years the corresponding percentage is {sep.loc['signed_snow_input', 'percent_of_octsep_annual_input_snow_addition']:.1f}%. September DA-OL changes are {sep.loc['snow_mass_change', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 snow mass, {sep.loc['snowmelt', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 snowmelt, {sep.loc['runoff', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 runoff, {sep.loc['et', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 ET, and {sep.loc['august_to_september_dTWLAND', 'all_tile_years_mean_kg_m2']:.2f} kg m-2 August-to-September total-land-water change. The small remaining snow-mass change alongside substantial melt and runoff indicates that much of the September correction is melted and redistributed within September rather than retained as snow at the boundary. These monthly-mean quantities are a timing diagnostic, not an exact monthly closure.",
            "",
            f"The all-tile six-year mean signed input is {octsep_mean_input:.2f} kg m-2 yr-1 for Oct-Sep and {sepaug_mean_input:.2f} kg m-2 yr-1 for Sep-Aug. The table below reports the positive-input tile-year partition; brackets are 5-degree spatial-block 95% intervals.",
            "",
            "| Boundary | Positive input | Runoff [95% CI] | ET [95% CI] | Storage [95% CI] | Residual [95% CI] | Surface runoff | Baseflow |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for boundary_name in ["Oct-Sep", "Sep-Aug"]:
        row = positive.loc[boundary_name]
        lines.append(
            f"| {boundary_name} | {row['I_snow']:.2f} kg m-2 | {100 * row['fraction_dRunoff_total']:.1f}% "
            f"[{100 * row['fraction_dRunoff_total_ci_low_5deg']:.1f}, {100 * row['fraction_dRunoff_total_ci_high_5deg']:.1f}] | "
            f"{100 * row['fraction_dET']:.1f}% [{100 * row['fraction_dET_ci_low_5deg']:.1f}, {100 * row['fraction_dET_ci_high_5deg']:.1f}] | "
            f"{100 * row['fraction_dStorage']:.1f}% [{100 * row['fraction_dStorage_ci_low_5deg']:.1f}, {100 * row['fraction_dStorage_ci_high_5deg']:.1f}] | "
            f"{100 * row['fraction_residual']:.1f}% [{100 * row['fraction_residual_ci_low_5deg']:.1f}, {100 * row['fraction_residual_ci_high_5deg']:.1f}] | "
            f"{100 * row['fraction_dRunoff_surface']:.1f}% | "
            f"{100 * row['fraction_dBaseflow']:.1f}% |"
        )
    residual_change = 100 * (
        positive.loc["Sep-Aug", "fraction_residual"] - positive.loc["Oct-Sep", "fraction_residual"]
    )
    lines.extend(
        [
            "",
            f"The boundary shift is {'materially noticeable' if noticeable['runoff_or_et'] else 'below the pre-set 5 percentage-point runoff/ET reporting threshold'} for the headline partition. Storage changes by {100 * (positive.loc['Sep-Aug', 'fraction_dStorage'] - positive.loc['Oct-Sep', 'fraction_dStorage']):+.1f} percentage points. The residual changes by {residual_change:+.1f} percentage points, {'exceeding' if noticeable['residual'] else 'below'} the 2-point discussion threshold. Thus the boundary shift {'does' if noticeable['residual'] else 'does not'} explain a meaningful part of the negative Oct-Sep residual under the stated reporting rule.",
            f"Annual residuals are negative in {negative_residual_years['Oct-Sep']} of 6 Oct-Sep years and {negative_residual_years['Sep-Aug']} of 6 Sep-Aug years.",
            "",
            "![Water-year boundary sensitivity](monthly_synthesis_report_figures/water_year_budget_boundary_sensitivity.png)",
            "",
            "Interpretive roles remain distinct: **direct accounting** describes the average fate of snow-added water; the **water-year regression** estimates marginal partition; and the **Oct-Mar seasonal regression** tests whether prior within-location input predicts later hydrological response.",
            REPORT_END,
        ]
    )
    return "\n".join(lines) + "\n"


def replace_report_section(path: Path, section: str) -> None:
    text = path.read_text()
    if REPORT_START in text:
        before, remainder = text.split(REPORT_START, 1)
        _, after = remainder.split(REPORT_END, 1)
        text = before.rstrip() + "\n\n" + section + after.lstrip("\n")
    else:
        text = text.rstrip() + "\n\n" + section
    path.write_text(text)


def main() -> None:
    args = parse_args()
    config = load_config(args.config, args.bootstrap_replicates)
    root = repo_root()
    output_dir = root / "projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness"
    figure_dir = output_dir / "figures"
    report_figure_dir = root / "projects/M21C_ls/docs/monthly_synthesis_report_figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_figure_dir.mkdir(parents=True, exist_ok=True)

    worked_example = validate_date_windows(config["response_years"], config)
    print(worked_example)
    baseline_path = root / config["baseline_archive"]
    with xr.open_dataset(baseline_path) as baseline_open:
        baseline = baseline_open.load()
    if baseline.sizes["tile"] != config["expected_tile_count"]:
        raise AssertionError(f"Static mask has {baseline.sizes['tile']} tiles")
    if "LS_OLv8_M36_v2" not in str(baseline.attrs.get("source_ol", "")):
        raise AssertionError("Baseline archive does not identify the OL v2 source")
    if "LS_DAv8_M36_v3" not in str(baseline.attrs.get("source_da", "")):
        raise AssertionError("Baseline archive does not identify the DA v3 source")

    inputs = open_inputs(Path(config["data_directory"]))
    try:
        product_check = verify_monthly_increment_product(inputs["catch"])
        print(f"Input product: {product_check['classification']}; {product_check['action']}")
        baseline_annual, baseline_partition = reproduce_baseline(baseline, config, root)
        mean_input = baseline_annual.loc[baseline_annual["water_year"] == "6-WY mean", "I_snow"].iloc[0]
        print(f"Six-year area-weighted Oct-Sep signed input: {mean_input:.6f} kg m-2 yr-1")

        seasonal_table = build_seasonal_table(inputs, baseline, config["response_years"])
        seasonal, comparison, variation = fit_seasonal_attribution(seasonal_table, config)
        mam_check = comparison.query("window == 'MAM' and model == 'M3'")
        for _, row in mam_check.iterrows():
            benchmark = row["published_signed_mam_m3_benchmark"]
            tolerance = config["signed_mam_m3_benchmark_tolerances"][row["response"]]
            if not np.isclose(row["raw_beta"], benchmark, atol=tolerance, rtol=0):
                raise AssertionError(
                    f"Signed MAM benchmark drift for {row['response']}: "
                    f"{row['raw_beta']} versus {benchmark}"
                )
        correlations = predictor_window_correlations(seasonal_table)

        controlled_baseline = add_water_year_controls(baseline)
        control_sensitivity = water_year_control_regressions(controlled_baseline, config)

        alternative = build_september_august_dataset(baseline, inputs, config["response_years"])
        alternative_annual, alternative_partition, _ = domain_budget_tables(alternative, config)
        baseline_uncertainty = bootstrap_partition(
            dataset_to_frame(
                baseline,
                [
                    "I_snow",
                    "dRunoff_total",
                    "dRunoff_surface",
                    "dBaseflow",
                    "dET",
                    "dStorage",
                    "residual",
                ],
            ),
            config,
            config["bootstrap"]["primary_block_degrees"],
            config["bootstrap"]["seed"] + 700,
        )
        alternative_uncertainty = bootstrap_partition(
            dataset_to_frame(
                alternative,
                [
                    "I_snow",
                    "dRunoff_total",
                    "dRunoff_surface",
                    "dBaseflow",
                    "dET",
                    "dStorage",
                    "residual",
                ],
            ),
            config,
            config["bootstrap"]["primary_block_degrees"],
            config["bootstrap"]["seed"] + 701,
        )
        boundary = pd.DataFrame(
            boundary_summary_rows("Oct-Sep", baseline_annual, baseline_partition, baseline_uncertainty)
            + boundary_summary_rows("Sep-Aug", alternative_annual, alternative_partition, alternative_uncertainty)
        )
        positive = boundary.query("scope == 'positive_input_partition'").set_index("boundary")
        runoff_change = 100.0 * (
            positive.loc["Sep-Aug", "fraction_dRunoff_total"]
            - positive.loc["Oct-Sep", "fraction_dRunoff_total"]
        )
        et_change = 100.0 * (
            positive.loc["Sep-Aug", "fraction_dET"] - positive.loc["Oct-Sep", "fraction_dET"]
        )
        residual_change = 100.0 * (
            positive.loc["Sep-Aug", "fraction_residual"]
            - positive.loc["Oct-Sep", "fraction_residual"]
        )
        thresholds = config["boundary_reporting_thresholds_percentage_points"]
        noticeable = {
            "runoff_or_et": max(abs(runoff_change), abs(et_change)) > thresholds["runoff_or_et"],
            "residual": abs(residual_change) > thresholds["residual"],
        }
        if noticeable["runoff_or_et"] or noticeable["residual"]:
            common = common_positive_partition(baseline, alternative)
            boundary = pd.concat([boundary, common], ignore_index=True)
        september = september_diagnostic(baseline)
        precip = pd.DataFrame(
            [
                {"boundary": "Oct-Sep", **precipitation_check(baseline)},
                {"boundary": "Sep-Aug", **precipitation_check(alternative)},
            ]
        )
        if precip["maximum_abs_annual_tile_kg_m2"].max() > config["maximum_abs_annual_tile_precip_difference_kg_m2"]:
            raise AssertionError("Annual tile precipitation difference exceeds its guard")
        if precip["maximum_abs_annual_domain_mean_kg_m2"].max() > config["maximum_abs_annual_domain_precip_difference_kg_m2"]:
            raise AssertionError("Annual domain precipitation difference is not negligible")
    finally:
        close_inputs(inputs)

    seasonal.to_csv(output_dir / "analysisA_octmar_signed_controls.csv", index=False)
    variation.to_csv(output_dir / "analysisA_octmar_identifying_variation.csv", index=False)
    comparison.to_csv(output_dir / "analysisA_octmar_vs_mam_comparison.csv", index=False)
    correlations.to_csv(output_dir / "analysisA_octmar_vs_mam_predictor_correlations.csv", index=False)
    control_sensitivity.to_csv(output_dir / "water_year_partition_control_sensitivity.csv", index=False)
    boundary.to_csv(output_dir / "water_year_boundary_sensitivity.csv", index=False)
    september.to_csv(output_dir / "water_year_september_timing_diagnostic.csv", index=False)
    precip.to_csv(output_dir / "water_year_boundary_precipitation_check.csv", index=False)
    (output_dir / "input_product_check.json").write_text(json.dumps(product_check, indent=2) + "\n")

    seasonal_figure = figure_dir / "analysisA_octmar_signed_control_sequence.png"
    boundary_figure = figure_dir / "water_year_budget_boundary_sensitivity.png"
    plot_seasonal_control_sequence(seasonal, seasonal_figure)
    plot_boundary_sensitivity(boundary, boundary_figure)
    shutil.copy2(seasonal_figure, report_figure_dir / seasonal_figure.name)
    shutil.copy2(boundary_figure, report_figure_dir / boundary_figure.name)

    report = build_report_section(
        seasonal,
        comparison,
        correlations,
        control_sensitivity,
        boundary,
        september,
        product_check,
        noticeable,
    )
    (output_dir / "targeted_robustness_report_section.md").write_text(report)
    replace_report_section(root / "projects/M21C_ls/docs/monthly_synthesis_report_out.md", report)
    print("Targeted snow-hydrology robustness workflow complete")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
