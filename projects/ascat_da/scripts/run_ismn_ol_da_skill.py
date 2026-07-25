#!/usr/bin/env python
"""ISMN in-situ soil-moisture skill for OL vs DA GEOSldas experiments.

Single analysis window (the full experiment span), every ISMN network in the
local archive -- no predetermined network list. Surface skill is computed
wherever a station has a sensor at or above `--surface-max-depth-m`; root-zone
skill is computed only where a station's profile supports the generic
0-1 m weighted composite, so root-zone site counts are reported separately.

Outputs (under --output-dir):
    cache_obs_daily.nc            daily ISMN surface/root-zone per station
    cache_model_daily_<run>.nc    daily model surface/root-zone per station
    ismn_station_inventory.csv    per-station metadata + obs-day counts
    ismn_skill_stations.csv       per station/domain/run skill statistics
    ismn_skill_network_summary.csv  per network/domain deltas vs the reference run
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
sys.path.insert(0, str(REPO_ROOT / "projects/matlab2python/scripts"))
sys.path.insert(0, str(REPO_ROOT / "common/python/io"))

from ismn_io import (  # noqa: E402
    build_station_series,
    filter_inventory_by_window,
    load_sm_sensor_inventory,
)
from read_GEOSldas import read_tilecoord  # noqa: E402
from sm_skill_vs_insitu import compute_anom, get_validation_stats  # noqa: E402


DEFAULT_ISMN_ROOT = Path("/discover/nobackup/projects/land_da/ISMN_data")
DEFAULT_OUTPUT_DIR = Path(
    "/discover/nobackup/projects/land_da/Evaluation/ISMN/ascat_da_ol_da_skill"
)
DEFAULT_DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
DEFAULT_COLLECTION = ".SMAP_L4_SM_gph."
DEFAULT_RUNS = (
    "OL=/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121",
    "DA_H121=/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_H121",
    "DA_legacy=/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_legacy",
    "DA_SMAP_comb_fp_scaled=/discover/nobackup/amfox/Experiments/"
    "DAv7_M36_SMAP_type_13_comb_fp_scaled/DAv7_M36_SMAP_type_13_comb_fp_scaled",
)

DOMAINS = ("surface", "rz")
METRICS = ("R", "anomR", "bias", "RMSE", "ubRMSE")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    start = parse_date(args.start_date)
    end = parse_date(args.end_date)
    runs = [parse_run_spec(spec) for spec in args.run]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Window: {start:%Y-%m-%d} .. {end:%Y-%m-%d}")
    print(f"Runs: {', '.join(name for name, _ in runs)}")
    print(f"Output: {output_dir}")

    obs = load_or_build_observations(args, start, end, output_dir)
    station_df = obs["stations"]
    print(f"Stations with usable ISMN observations: {len(station_df)}")

    model = {}
    for name, root in runs:
        model[name] = load_or_build_model(
            args, name, root, station_df, start, end, output_dir
        )

    stats_df = compute_skill(args, obs, model, station_df)
    stats_path = output_dir / "ismn_skill_stations.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"Wrote per-station skill: {stats_path} ({len(stats_df)} rows)")

    summary_df = summarize_networks(stats_df, args.reference_run, args.min_sites_per_network)
    summary_path = output_dir / "ismn_skill_network_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Wrote network summary: {summary_path} ({len(summary_df)} rows)")

    report(stats_df, summary_df, args.reference_run)


# ----------------------
# Observations
# ----------------------
def load_or_build_observations(args, start, end, output_dir: Path) -> dict:
    cache_path = output_dir / "cache_obs_daily.nc"
    inventory_path = output_dir / "ismn_station_inventory.csv"

    if cache_path.exists() and not args.overwrite_obs:
        print(f"Loading ISMN observation cache: {cache_path}")
        with xr.open_dataset(cache_path) as ds:
            ds.load()
        station_df = pd.read_csv(inventory_path)
        return {"dataset": ds, "stations": station_df}

    print(f"Building ISMN observations from archive: {args.archive_root}")
    inventory = load_sm_sensor_inventory(args.metadata_csv, args.archive_root)
    print(f"Soil-moisture sensors in archive: {len(inventory)}")

    inventory = filter_inventory_by_window(inventory, start, end)
    print(f"Sensors overlapping the window: {len(inventory)}")

    if args.network:
        wanted = {n.strip().upper() for n in args.network}
        inventory = inventory[inventory["network"].str.upper().isin(wanted)].reset_index(drop=True)
        print(f"Sensors after --network filter: {len(inventory)}")

    days = pd.date_range(start, end, freq="D") + pd.Timedelta(hours=args.daily_shift_hours)
    groups = list(inventory.groupby("station_key", sort=True))
    if args.max_stations:
        groups = groups[: args.max_stations]
    print(f"Candidate stations: {len(groups)}")

    surface_cols: dict[str, np.ndarray] = {}
    rz_cols: dict[str, np.ndarray] = {}
    records = []

    for i, (key, rows) in enumerate(groups, start=1):
        pack = build_station_series(
            rows,
            required_flag=args.required_flag,
            target_surface_depth_m=args.target_surface_depth_m,
            surface_max_depth_m=args.surface_max_depth_m,
            daily_shift_hours=args.daily_shift_hours,
            rz_top_m=args.rz_top_m,
            rz_bottom_m=args.rz_bottom_m,
            rz_min_sensors=args.rz_min_sensors,
            rz_shallow_max_m=args.rz_shallow_max_m,
            rz_deep_min_m=args.rz_deep_min_m,
            start=start,
            end=end,
        )

        if i % 250 == 0 or i == len(groups):
            print(f"  processed {i}/{len(groups)} stations")

        if not pack.get("ok"):
            continue

        surface = pack["surface_daily"].reindex(days)
        rootzone = pack["rz_daily"].reindex(days)

        n_surface = int(np.isfinite(surface.to_numpy(dtype=float)).sum())
        n_rz = int(np.isfinite(rootzone.to_numpy(dtype=float)).sum())
        if n_surface == 0 and n_rz == 0:
            continue

        surface_cols[key] = surface.to_numpy(dtype=np.float32)
        rz_cols[key] = rootzone.to_numpy(dtype=np.float32)

        first = rows.iloc[0]
        records.append(
            {
                "station_key": key,
                "network": first["network"],
                "station": first["station"],
                "lat": float(first["lat"]),
                "lon": float(first["lon"]),
                "surface_depth_m": pack["surface_depth_m"],
                "n_depths": pack["n_depths"],
                "depths_m": ",".join(f"{d:.3f}" for d in pack["depths_m"]),
                "surface_obs_days": n_surface,
                "rz_obs_days": n_rz,
            }
        )

    if not records:
        raise RuntimeError("No stations produced usable ISMN observation series")

    station_df = pd.DataFrame(records).sort_values("station_key").reset_index(drop=True)
    stations = station_df["station_key"].tolist()

    ds = xr.Dataset(
        data_vars={
            "obs_surface": (
                ("time", "station"),
                np.column_stack([surface_cols[s] for s in stations]),
            ),
            "obs_rz": (
                ("time", "station"),
                np.column_stack([rz_cols[s] for s in stations]),
            ),
            "station_lat": (("station",), station_df["lat"].to_numpy(dtype=np.float32)),
            "station_lon": (("station",), station_df["lon"].to_numpy(dtype=np.float32)),
        },
        coords={"time": days, "station": np.array(stations, dtype=object)},
        attrs={
            "required_flag": args.required_flag,
            "target_surface_depth_m": args.target_surface_depth_m,
            "surface_max_depth_m": args.surface_max_depth_m,
            "daily_shift_hours": args.daily_shift_hours,
            "rz_min_sensors": args.rz_min_sensors,
            "rz_shallow_max_m": args.rz_shallow_max_m,
            "rz_deep_min_m": args.rz_deep_min_m,
            "created_utc": pd.Timestamp.utcnow().isoformat(),
        },
    )

    ds.to_netcdf(cache_path)
    station_df.to_csv(inventory_path, index=False)
    print(f"Wrote ISMN observation cache: {cache_path}")
    print(f"Wrote station inventory: {inventory_path}")
    return {"dataset": ds, "stations": station_df}


# ----------------------
# Model
# ----------------------
def load_or_build_model(args, name, root: Path, station_df, start, end, output_dir: Path):
    cache_path = output_dir / f"cache_model_daily_{name}.nc"
    if cache_path.exists() and not args.overwrite_model:
        print(f"Loading model cache for {name}: {cache_path}")
        with xr.open_dataset(cache_path) as ds:
            ds.load()
        return ds

    exp_id = root.name
    tilecoord_path = root / "output" / args.domain / "rc_out" / f"{exp_id}.ldas_tilecoord.bin"
    if not tilecoord_path.exists():
        matches = sorted((root / "output" / args.domain / "rc_out").glob("*.ldas_tilecoord.bin"))
        if not matches:
            raise FileNotFoundError(f"No tilecoord found for {name} under {root}")
        tilecoord_path = matches[0]

    tc = read_tilecoord(str(tilecoord_path))
    tile_lat = np.asarray(tc["com_lat"], dtype=float)
    tile_lon = np.asarray(tc["com_lon"], dtype=float)

    stations = station_df["station_key"].tolist()
    lat = station_df["lat"].to_numpy(dtype=float)
    lon = station_df["lon"].to_numpy(dtype=float)

    tile_index = np.full(len(stations), -1, dtype=np.int64)
    tile_distance = np.full(len(stations), np.nan, dtype=float)
    for i in range(len(stations)):
        if not (np.isfinite(lat[i]) and np.isfinite(lon[i])):
            continue
        metric = (tile_lat - lat[i]) ** 2 + (tile_lon - lon[i]) ** 2
        j = int(np.argmin(metric))
        if metric[j] > float(args.max_distance_deg2):
            continue
        tile_index[i] = j
        tile_distance[i] = float(metric[j])

    mapped = tile_index >= 0
    print(f"{name}: stations mapped to tiles: {int(mapped.sum())}/{len(stations)}")
    if not mapped.any():
        raise RuntimeError(f"No stations mapped to tiles for run {name}")

    read_tiles = tile_index[mapped]
    days = pd.date_range(start, end, freq="D")
    surface = np.full((len(days), len(stations)), np.nan, dtype=np.float32)
    rootzone = np.full((len(days), len(stations)), np.nan, dtype=np.float32)

    n_found = 0
    for i, day in enumerate(days):
        path = resolve_model_path(root, exp_id, args.domain, args.collection, day)
        if path is None:
            continue
        sf, rz = read_daily_tiles(path, args.model_surface_variable, args.model_rz_variable, read_tiles)
        surface[i, mapped] = sf
        rootzone[i, mapped] = rz
        n_found += 1
        if (i + 1) % 365 == 0 or (i + 1) == len(days):
            print(f"  {name}: {i + 1}/{len(days)} days, files found={n_found}")

    print(f"{name}: model files found {n_found}/{len(days)}")

    ds = xr.Dataset(
        data_vars={
            "model_surface": (("time", "station"), surface),
            "model_rz": (("time", "station"), rootzone),
            "tile_index": (("station",), tile_index),
            "tile_distance_deg2": (("station",), tile_distance.astype(np.float32)),
        },
        coords={
            "time": days + pd.Timedelta(hours=args.daily_shift_hours),
            "station": np.array(stations, dtype=object),
        },
        attrs={
            "run": name,
            "run_root": str(root),
            "collection": args.collection,
            "surface_variable": args.model_surface_variable,
            "rz_variable": args.model_rz_variable,
            "days_found": n_found,
            "created_utc": pd.Timestamp.utcnow().isoformat(),
        },
    )
    ds.to_netcdf(cache_path)
    print(f"Wrote model cache: {cache_path}")
    return ds


def resolve_model_path(root: Path, exp_id: str, domain: str, collection: str, day) -> Path | None:
    directory = root / "output" / domain / "cat" / "ens_avg" / f"Y{day.year:04d}" / f"M{day.month:02d}"
    stamp = f"{day:%Y%m%d}"
    for name in (f"{exp_id}{collection}{stamp}.nc4", f"{exp_id}{collection}{stamp}_1200z.nc4"):
        candidate = directory / name
        if candidate.exists():
            return candidate
    return None


def read_daily_tiles(path: Path, surface_variable: str, rz_variable: str, tiles: np.ndarray):
    """Daily-mean model surface/root-zone soil moisture for selected tiles.

    SMAP_L4_SM_gph bundles a day's eight 3-hourly instants into one file, so
    the time axis is averaged; single-record collections pass straight through.
    """

    from netCDF4 import Dataset

    with Dataset(path, "r") as ds:
        ds.set_auto_mask(False)
        surface = np.asarray(ds[surface_variable][:], dtype=np.float64)
        rootzone = np.asarray(ds[rz_variable][:], dtype=np.float64)

    def reduce(values: np.ndarray) -> np.ndarray:
        values = np.where(values > 1e14, np.nan, values)
        if values.ndim == 2:
            values = np.nanmean(values, axis=0)
        subset = values[tiles]
        # Model zeros are treated as missing, matching the legacy skill scripts.
        return np.where(subset == 0.0, np.nan, subset)

    return reduce(surface), reduce(rootzone)


# ----------------------
# Statistics
# ----------------------
def compute_skill(args, obs, model, station_df) -> pd.DataFrame:
    obs_ds = obs["dataset"]
    obs_time = pd.DatetimeIndex(obs_ds["time"].values)
    stations = [str(s) for s in obs_ds["station"].values]
    station_meta = station_df.set_index("station_key")

    obs_arrays = {
        "surface": np.asarray(obs_ds["obs_surface"].values, dtype=float),
        "rz": np.asarray(obs_ds["obs_rz"].values, dtype=float),
    }

    records = []
    for run_name, run_ds in model.items():
        run_time = pd.DatetimeIndex(run_ds["time"].values)
        if not run_time.equals(obs_time):
            raise RuntimeError(f"Time axis mismatch between obs cache and run {run_name}")
        # A stale model cache built against a different station set would
        # silently pair the wrong station with the wrong tile.
        if [str(s) for s in run_ds["station"].values] != stations:
            raise RuntimeError(
                f"Station axis mismatch between obs cache and run {run_name}; "
                "rerun with --overwrite-model"
            )

        model_arrays = {
            "surface": np.asarray(run_ds["model_surface"].values, dtype=float),
            "rz": np.asarray(run_ds["model_rz"].values, dtype=float),
        }
        tile_index = np.asarray(run_ds["tile_index"].values, dtype=np.int64)
        tile_distance = np.asarray(run_ds["tile_distance_deg2"].values, dtype=float)

        for j, key in enumerate(stations):
            if tile_index[j] < 0:
                continue
            meta = station_meta.loc[key]

            for domain in DOMAINS:
                obs_series = obs_arrays[domain][:, j]
                mod_series = model_arrays[domain][:, j]
                paired = np.isfinite(obs_series) & np.isfinite(mod_series)
                n_pairs = int(paired.sum())
                if n_pairs <= args.nmin:
                    continue

                stats = skill_stats(
                    obs_series[paired],
                    mod_series[paired],
                    obs_time[paired],
                    nmin=args.nmin,
                    nmin_day=args.nmin_day,
                )
                records.append(
                    {
                        "station_key": key,
                        "network": meta["network"],
                        "station": meta["station"],
                        "lat": meta["lat"],
                        "lon": meta["lon"],
                        "domain": domain,
                        "run": run_name,
                        "surface_depth_m": meta["surface_depth_m"],
                        "tile_index": int(tile_index[j]),
                        "tile_distance_deg2": tile_distance[j],
                        **stats,
                    }
                )

    if not records:
        raise RuntimeError("No station/run combinations met the pairing threshold")

    stats_df = pd.DataFrame(records)

    # A station is kept only if every run cleared the threshold for that domain,
    # so run-to-run differences never reflect a changing site population.
    counts = stats_df.groupby(["station_key", "domain"])["run"].nunique()
    complete = counts[counts == len(model)].index
    if complete.empty:
        raise RuntimeError(
            "No station/domain cleared the pairing threshold in every run; "
            "lower --nmin or check model file coverage"
        )
    stats_df = stats_df.set_index(["station_key", "domain"]).loc[complete].reset_index()
    return stats_df.sort_values(["network", "station", "domain", "run"]).reset_index(drop=True)


def skill_stats(obs_values, mod_values, times, nmin: int, nmin_day: int) -> dict:
    data = np.column_stack([obs_values, mod_values])
    stats = get_validation_stats(data, AC=True, complete=True, ref_col=1, select_cols=[1, 2], Nmin=nmin)

    out = {
        "N_pairs": int(stats["N_pairs"]),
        "R": stats["R"],
        "bias": stats["bias"],
        "MSE": stats["MSE"],
        "ubMSE": stats["ubMSE"],
        "RMSE": float(np.sqrt(stats["MSE"])) if np.isfinite(stats["MSE"]) else np.nan,
        "ubRMSE": float(np.sqrt(max(stats["ubMSE"], 0.0))) if np.isfinite(stats["ubMSE"]) else np.nan,
        "anomN_pairs": 0,
        "anomR": np.nan,
    }

    doy = times.dayofyear.to_numpy(dtype=int)
    obs_anom = compute_anom(np.asarray(obs_values, dtype=float), doy, Nmin_day=nmin_day)
    mod_anom = compute_anom(np.asarray(mod_values, dtype=float), doy, Nmin_day=nmin_day)
    anom_ok = np.isfinite(obs_anom) & np.isfinite(mod_anom)
    out["anomN_pairs"] = int(anom_ok.sum())

    if anom_ok.sum() > nmin:
        anom_stats = get_validation_stats(
            np.column_stack([obs_anom[anom_ok], mod_anom[anom_ok]]),
            AC=True,
            complete=True,
            ref_col=1,
            select_cols=[1, 2],
            Nmin=nmin,
        )
        out["anomR"] = anom_stats["R"]

    return out


def summarize_networks(stats_df: pd.DataFrame, reference_run: str, min_sites: int) -> pd.DataFrame:
    if reference_run not in set(stats_df["run"]):
        raise ValueError(f"Reference run {reference_run!r} not present in results")

    rows = []
    for (network, domain), group in stats_df.groupby(["network", "domain"]):
        n_sites = group["station_key"].nunique()
        if n_sites < min_sites:
            continue

        for run_name in sorted(set(group["run"])):
            record = {
                "network": network,
                "domain": domain,
                "run": run_name,
                "n_sites": n_sites,
                "reference_run": reference_run,
            }
            for metric in METRICS:
                pivot = group.pivot_table(index="station_key", columns="run", values=metric)
                if run_name not in pivot.columns:
                    continue
                values = pivot[run_name].to_numpy(dtype=float)
                record[f"{metric}_mean"] = float(np.nanmean(values)) if np.isfinite(values).any() else np.nan

                if run_name == reference_run or reference_run not in pivot.columns:
                    continue
                delta = paired_delta(metric, pivot[reference_run].to_numpy(dtype=float), values)
                finite = np.isfinite(delta)
                record[f"{metric}_delta_mean"] = float(np.mean(delta[finite])) if finite.any() else np.nan
                record[f"{metric}_delta_sem"] = (
                    float(np.std(delta[finite], ddof=1) / np.sqrt(finite.sum())) if finite.sum() > 1 else np.nan
                )
            rows.append(record)

    return pd.DataFrame(rows).sort_values(["domain", "network", "run"]).reset_index(drop=True)


def paired_delta(metric: str, reference: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Per-station improvement over the reference run.

    Signed so that positive always means the run beats the reference: skill
    metrics gain, error metrics shrink, and bias is compared in magnitude
    because either sign of bias is a miss.
    """

    if metric == "bias":
        return np.abs(reference) - np.abs(values)
    if metric in ("RMSE", "ubRMSE", "MSE", "ubMSE"):
        return reference - values
    return values - reference


def report(stats_df: pd.DataFrame, summary_df: pd.DataFrame, reference_run: str) -> None:
    print("\n=== Site counts by domain ===")
    for domain in DOMAINS:
        subset = stats_df[stats_df["domain"] == domain]
        n_sites = subset["station_key"].nunique()
        n_networks = subset["network"].nunique()
        print(f"  {domain:8s}: {n_sites} sites across {n_networks} networks")

    print(f"\n=== Global mean skill by run (reference: {reference_run}) ===")
    for domain in DOMAINS:
        subset = stats_df[stats_df["domain"] == domain]
        if subset.empty:
            continue
        print(f"  -- {domain} --")
        for run_name in sorted(set(subset["run"])):
            run_rows = subset[subset["run"] == run_name]
            print(
                f"     {run_name:26s} R={run_rows['R'].mean():.4f} "
                f"anomR={run_rows['anomR'].mean():.4f} "
                f"ubRMSE={run_rows['ubRMSE'].mean():.4f} "
                f"n={len(run_rows)}"
            )


# ----------------------
# CLI
# ----------------------
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start-date", default="2015-04-01")
    parser.add_argument("--end-date", default="2021-03-31")
    parser.add_argument(
        "--run",
        action="append",
        default=None,
        help="NAME=/path/to/run; repeatable. Defaults to OL + the three DA runs.",
    )
    parser.add_argument("--reference-run", default="OL", help="Run that deltas are measured against.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)

    parser.add_argument("--archive-root", type=Path, default=DEFAULT_ISMN_ROOT)
    parser.add_argument("--metadata-csv", type=Path, default=None)
    parser.add_argument("--required-flag", default="G", help="ISMN quality flag to retain.")
    parser.add_argument(
        "--network",
        action="append",
        default=None,
        help="Restrict to these ISMN networks; repeatable. Default is every network in the archive.",
    )
    parser.add_argument("--max-stations", type=int, default=0, help="Debug cap on stations processed.")

    parser.add_argument("--target-surface-depth-m", type=float, default=0.05)
    parser.add_argument(
        "--surface-max-depth-m",
        type=float,
        default=0.10,
        help="Deepest sensor still eligible to represent surface soil moisture.",
    )
    parser.add_argument("--rz-top-m", type=float, default=0.0)
    parser.add_argument("--rz-bottom-m", type=float, default=1.0)
    parser.add_argument("--rz-min-sensors", type=int, default=3)
    parser.add_argument("--rz-shallow-max-m", type=float, default=0.20)
    parser.add_argument("--rz-deep-min-m", type=float, default=0.50)
    parser.add_argument("--daily-shift-hours", type=int, default=12)

    parser.add_argument("--domain", default=DEFAULT_DOMAIN)
    parser.add_argument("--collection", default=DEFAULT_COLLECTION)
    parser.add_argument("--model-surface-variable", default="sm_surface")
    parser.add_argument("--model-rz-variable", default="sm_rootzone")
    parser.add_argument("--max-distance-deg2", type=float, default=0.1)

    parser.add_argument(
        "--nmin",
        type=int,
        default=1000,
        help="Minimum paired obs/model days for a station/domain to be scored.",
    )
    parser.add_argument("--nmin-day", type=int, default=30, help="Minimum samples per day-of-year climatology.")
    parser.add_argument("--min-sites-per-network", type=int, default=1)

    parser.add_argument("--overwrite-obs", action="store_true")
    parser.add_argument("--overwrite-model", action="store_true")

    args = parser.parse_args(argv)
    if args.run is None:
        args.run = list(DEFAULT_RUNS)
    if args.metadata_csv is None:
        args.metadata_csv = args.archive_root / "python_metadata" / "ISMN_data.csv"
    return args


def parse_run_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"Run must be NAME=PATH, got {spec!r}")
    name, root = spec.split("=", 1)
    name = name.strip()
    root = root.strip()
    if not name or not root:
        raise ValueError(f"Run must be NAME=PATH, got {spec!r}")
    return name, Path(root)


def parse_date(value: str):
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            continue
    raise ValueError(f"Unsupported date {value!r}")


if __name__ == "__main__":
    main()
