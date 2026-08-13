#!/usr/bin/env python3
"""Detect and compare independent changepoints in all 43 Phase 1 domain series."""

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

from changepoint_detection import (
    DEFAULT_CONFIG,
    compare_known_boundaries,
    detect_changepoints,
    known_boundary_table,
    load_changepoint_config,
)
from m21c_periods import DEFAULT_REGISTRY
from phase1_changepoint_workflow import (
    DEFAULT_AUDIT,
    DEFAULT_COMPARISON,
    DEFAULT_DETECTIONS,
    DEFAULT_MODELS,
    DEFAULT_MONTHLY_OUTPUT,
    DEFAULT_PENALTY_GRID,
    audit_changepoint_outputs,
)
from phase1_interrupted_workflow import DEFAULT_MONTHLY
from phase1_trend_workflow import DEFAULT_RUN_MATRIX
from trend_breakpoint_series import DEFAULT_VARIABLE_SELECTION


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-monthly", type=Path, default=DEFAULT_MONTHLY)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--run-matrix", type=Path, default=DEFAULT_RUN_MATRIX)
    parser.add_argument("--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION)
    parser.add_argument("--detections", type=Path, default=DEFAULT_DETECTIONS)
    parser.add_argument("--penalty-grid", type=Path, default=DEFAULT_PENALTY_GRID)
    parser.add_argument("--comparison", type=Path, default=DEFAULT_COMPARISON)
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--monthly", type=Path, default=DEFAULT_MONTHLY_OUTPUT)
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
    try:
        dataset.to_netcdf(
            temporary,
            encoding={name: {"zlib": True, "complevel": 4} for name in dataset.data_vars},
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _concat_frames(frames: tuple[pd.DataFrame, ...]) -> pd.DataFrame:
    """Concatenate nonempty frames while preserving an all-empty schema."""

    nonempty = [frame for frame in frames if not frame.empty]
    return pd.concat(nonempty, ignore_index=True, sort=False) if nonempty else frames[0].copy()


def _fit_one(
    series_id: str,
    values: np.ndarray,
    time: np.ndarray,
    metadata: dict[str, Any],
    settings: dict[str, Any],
    registry: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], np.ndarray, np.ndarray]:
    result = detect_changepoints(values, time, config=settings)
    detections = result.detections.copy()
    penalty_grid = result.penalty_detections.copy()
    comparison = compare_known_boundaries(
        detections,
        time,
        config=settings,
        registry_path=registry,
    )
    boundaries = known_boundary_table(time, registry_path=registry)
    if not detections.empty:
        nearest_rows = []
        for detection in detections.itertuples():
            distance = np.abs(boundaries["boundary_index"] - int(detection.break_index))
            nearest = boundaries.loc[distance.idxmin()]
            nearest_rows.append(
                {
                    "nearest_period_id": nearest["period_id"],
                    "nearest_boundary_date": nearest["boundary_date"],
                    "nearest_boundary_distance_months": int(distance.min()),
                    "nearest_boundary_detection_exempt": bool(nearest["detection_exempt"]),
                }
            )
        detections = pd.concat([detections.reset_index(drop=True), pd.DataFrame(nearest_rows)], axis=1)
    for frame in (detections, penalty_grid, comparison):
        for key, value in metadata.items():
            frame[key] = value
    return (
        detections,
        penalty_grid,
        comparison,
        {**metadata, **result.summary},
        result.seasonal_adjusted,
        result.global_residual,
    )


def main() -> int:
    args = parse_args()
    if args.n_jobs == 0:
        raise ValueError("--n-jobs cannot be zero")
    settings = load_changepoint_config(args.config)
    with xr.open_dataset(args.input_monthly) as source:
        if source.sizes.get("series") != 43 or source.sizes.get("time") != 288:
            raise ValueError("Input interrupted-series monthly file is not the 43 x 288 product")
        source = source.load()
    series_ids = source["series_id"].values.astype(str)
    time = source.time.values
    observed = source["observed"].values.astype("float64")
    metadata_rows = [
        {
            "series_id": series_id,
            **{
                name: str(source[name].sel(series=series_id).item())
                for name in (
                    "run_id",
                    "variable",
                    "source_series",
                    "mask",
                    "role",
                    "series_role",
                    "units",
                    "rationale",
                )
            },
        }
        for series_id in series_ids
    ]
    print(f"Detecting changepoints in {len(series_ids)} domain series", flush=True)
    results = Parallel(n_jobs=args.n_jobs, prefer=args.parallel_prefer)(
        delayed(_fit_one)(series_id, values, time, metadata, settings, args.registry)
        for series_id, values, metadata in zip(series_ids, observed, metadata_rows)
    )
    detections, penalty_grid, comparison, models, adjusted, residual = zip(*results)
    detections_frame = _concat_frames(detections)
    penalty_frame = _concat_frames(penalty_grid)
    comparison_frame = pd.concat(comparison, ignore_index=True, sort=False)
    models_frame = pd.DataFrame(models)
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    commit = git_commit()
    for frame in (detections_frame, penalty_frame, comparison_frame, models_frame):
        frame["generated_utc"] = generated
        frame["source_git_commit"] = commit
        frame["input_source_git_commit"] = str(source.attrs.get("source_git_commit", ""))

    monthly = xr.Dataset(
        {
            "observed": (("series", "time"), observed),
            "seasonal_adjusted": (("series", "time"), np.stack(adjusted)),
            "global_linear_residual": (("series", "time"), np.stack(residual)),
        },
        coords={
            "series": series_ids,
            "time": time,
            **{
                name: ("series", source[name].values.astype(str))
                for name in (
                    "series_id",
                    "run_id",
                    "variable",
                    "source_series",
                    "mask",
                    "role",
                    "series_role",
                    "units",
                    "rationale",
                )
            },
        },
        attrs={
            "title": "M21C Phase 1 independent changepoint preprocessing",
            "configuration": json.dumps(settings, sort_keys=True),
            "config_path": str(args.config),
            "observing_system_registry": str(args.registry),
            "input_monthly": str(args.input_monthly),
            "input_source_git_commit": str(source.attrs.get("source_git_commit", "")),
            "generated_utc": generated,
            "source_git_commit": commit,
        },
    )
    _atomic_csv(detections_frame, args.detections)
    _atomic_csv(penalty_frame, args.penalty_grid)
    _atomic_csv(comparison_frame, args.comparison)
    _atomic_csv(models_frame, args.models)
    _atomic_netcdf(monthly, args.monthly)
    audit = audit_changepoint_outputs(
        input_monthly=args.input_monthly,
        detections_path=args.detections,
        penalty_grid_path=args.penalty_grid,
        comparison_path=args.comparison,
        models_path=args.models,
        monthly_path=args.monthly,
        run_matrix=args.run_matrix,
        variable_selection=args.variable_selection,
        config_path=args.config,
    )
    _atomic_csv(audit, args.audit)
    print(audit.to_string(index=False), flush=True)
    return 0 if audit.iloc[0]["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
