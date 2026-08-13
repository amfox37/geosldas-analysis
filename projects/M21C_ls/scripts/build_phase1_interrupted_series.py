#!/usr/bin/env python3
"""Build area-weighted P1-P9 interrupted-series results for Phase 1."""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

from interrupted_time_series import DEFAULT_CONFIG, fit_interrupted_series, load_interrupted_config
from m21c_periods import DEFAULT_REGISTRY
from phase1_interrupted_workflow import (
    DEFAULT_AUDIT,
    DEFAULT_COEFFICIENTS,
    DEFAULT_MODELS,
    DEFAULT_MONTHLY,
    assign_transition_fdr,
    audit_interrupted_outputs,
    phase1_series_catalogue,
)
from phase1_trend_workflow import DEFAULT_RUN_MATRIX
from trend_breakpoint_series import (
    DEFAULT_DATA_DIR,
    DEFAULT_INPUT_CONTRACT,
    DEFAULT_VARIABLE_SELECTION,
    MonthlySeriesLoader,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument("--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION)
    parser.add_argument("--run-matrix", type=Path, default=DEFAULT_RUN_MATRIX)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--coefficients", type=Path, default=DEFAULT_COEFFICIENTS)
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--monthly", type=Path, default=DEFAULT_MONTHLY)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--parallel-prefer", choices=["threads", "processes"], default="processes")
    return parser.parse_args()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.incomplete")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _atomic_netcdf(dataset: xr.Dataset, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.incomplete{path.suffix}")
    temporary.unlink(missing_ok=True)
    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in dataset.data_vars
        if dataset[name].dtype.kind in "fiu"
    }
    try:
        dataset.to_netcdf(temporary, encoding=encoding)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _fit_one(
    row: dict[str, Any],
    values: np.ndarray,
    time: np.ndarray,
    settings: dict[str, Any],
    registry_path: Path,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any], np.ndarray, np.ndarray]:
    result = fit_interrupted_series(
        values,
        time,
        config=settings,
        registry_path=registry_path,
        random_seed=seed,
    )
    coefficient = pd.concat(
        [result.coefficients, result.period_slopes], ignore_index=True, sort=False
    )
    for name, value in row.items():
        coefficient[name] = value
    return coefficient, {**row, **result.summary}, result.fitted, result.residual


def main() -> int:
    args = parse_args()
    if args.n_jobs == 0:
        raise ValueError("--n-jobs cannot be zero")
    settings = load_interrupted_config(args.config)
    catalogue = phase1_series_catalogue(
        args.run_matrix, variable_selection=args.variable_selection
    )
    series_rows: list[dict[str, Any]] = []
    observed: list[np.ndarray] = []
    n_tiles: list[np.ndarray] = []
    area_sum: list[np.ndarray] = []
    time_values: np.ndarray | None = None

    print(f"Loading {len(catalogue)} area-weighted Phase 1 domain series", flush=True)
    loaded_count = 0
    with MonthlySeriesLoader(
        args.data_dir,
        input_contract=args.input_contract,
        variable_selection=args.variable_selection,
    ) as loader:
        for _, group in catalogue.groupby("run_id", sort=False):
            first = group.iloc[0]
            domain = loader.domain_mean(first["variable"], mask=first["mask"]).load()
            for _, row in group.iterrows():
                source = row["source_series"]
                values = domain[source]
                if time_values is None:
                    time_values = np.asarray(values.time.values)
                elif not np.array_equal(time_values, values.time.values):
                    raise ValueError("Domain series have inconsistent time coordinates")
                current = row.to_dict()
                current["units"] = str(values.attrs.get("units", ""))
                series_rows.append(current)
                observed.append(np.asarray(values.values, dtype="float64"))
                n_tiles.append(np.asarray(domain["n_tiles"].values, dtype="int32"))
                area_sum.append(np.asarray(domain["area_sum"].values, dtype="float64"))
                loaded_count += 1
                print(
                    f"  [{loaded_count:02d}/{len(catalogue)}] {row['series_id']}",
                    flush=True,
                )
    if time_values is None:
        raise RuntimeError("No interrupted series were loaded")

    base_seed = int(settings["bootstrap"]["random_seed"])
    print(
        f"Fitting {len(series_rows)} models with "
        f"{settings['bootstrap']['replicates']} bootstrap replicates each",
        flush=True,
    )
    fitted_results = Parallel(n_jobs=args.n_jobs, prefer=args.parallel_prefer)(
        delayed(_fit_one)(
            row,
            values,
            time_values,
            settings,
            args.registry,
            base_seed + index,
        )
        for index, (row, values) in enumerate(zip(series_rows, observed))
    )
    coefficient_frames, model_rows, fitted, residual = zip(*fitted_results)
    coefficients = pd.concat(coefficient_frames, ignore_index=True, sort=False)
    coefficients = assign_transition_fdr(
        coefficients, alpha=float(settings["fdr"]["alpha"])
    )
    models = pd.DataFrame(model_rows)
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    commit = git_commit()
    models["generated_utc"] = generated
    models["source_git_commit"] = commit
    coefficients["generated_utc"] = generated
    coefficients["source_git_commit"] = commit

    metadata = pd.DataFrame(series_rows)
    series_ids = metadata["series_id"].to_numpy(dtype=str)
    monthly = xr.Dataset(
        {
            "observed": (("series", "time"), np.stack(observed)),
            "fitted": (("series", "time"), np.stack(fitted)),
            "residual": (("series", "time"), np.stack(residual)),
            "n_tiles": (("series", "time"), np.stack(n_tiles)),
            "area_sum": (("series", "time"), np.stack(area_sum)),
        },
        coords={
            "series": series_ids,
            "time": time_values,
            **{
                column: ("series", metadata[column].fillna("").astype(str).to_numpy())
                for column in metadata.columns
                if column != "series_id"
            },
            "series_id": ("series", series_ids),
        },
        attrs={
            "title": "M21C Phase 1 area-weighted P1-P9 interrupted time series",
            "configuration": json.dumps(settings, sort_keys=True),
            "config_path": str(args.config),
            "observing_system_registry": str(args.registry),
            "run_matrix": str(args.run_matrix),
            "input_contract": str(args.input_contract),
            "variable_selection": str(args.variable_selection),
            "data_directory": str(args.data_dir.resolve()),
            "generated_utc": generated,
            "source_git_commit": commit,
            "fdr_scope": "separate BH family for each transition coefficient across all domain series",
            "p7_constraint": "level change estimated; slope constrained to P6 until P8",
        },
    )
    monthly["n_tiles"].attrs = {"long_name": "finite contributing land tiles", "units": "1"}
    monthly["area_sum"].attrs = {
        "long_name": "area represented by finite contributing tiles",
        "units": "km2",
    }

    _atomic_csv(coefficients, args.coefficients)
    _atomic_csv(models, args.models)
    _atomic_netcdf(monthly, args.monthly)
    audit = audit_interrupted_outputs(
        coefficients_path=args.coefficients,
        models_path=args.models,
        monthly_path=args.monthly,
        run_matrix=args.run_matrix,
        variable_selection=args.variable_selection,
        config_path=args.config,
    )
    _atomic_csv(audit, args.audit)
    print(audit.to_string(index=False), flush=True)
    if audit.iloc[0]["status"] != "PASS":
        return 1
    print(f"Wrote {args.coefficients}")
    print(f"Wrote {args.models}")
    print(f"Wrote {args.monthly}")
    print(f"Wrote {args.audit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
