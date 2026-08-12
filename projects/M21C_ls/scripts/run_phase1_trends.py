#!/usr/bin/env python3
"""Run the restartable M21C Phase 1 production trend matrix."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from phase1_trend_workflow import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_RUN_MATRIX,
    audit_phase1_outputs,
    audit_trend_output,
    load_phase1_runs,
    output_filename,
)
from trend_breakpoint_series import (
    DEFAULT_DATA_DIR,
    DEFAULT_INPUT_CONTRACT,
    DEFAULT_VARIABLE_SELECTION,
)
from trend_statistics import DEFAULT_CONFIG


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUILDER = PROJECT_ROOT / "scripts" / "build_trend_statistics.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-matrix", type=Path, default=DEFAULT_RUN_MATRIX)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument(
        "--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION
    )
    parser.add_argument("--trend-config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument(
        "--parallel-prefer", choices=["threads", "processes"], default="processes"
    )
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        metavar="RUN_ID",
        help="Run only a named matrix row; may be repeated",
    )
    parser.add_argument("--role", choices=["primary", "sensitivity"])
    parser.add_argument("--force", action="store_true", help="Rebuild outputs that pass audit")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--diagnostic-tile-stop",
        type=int,
        help="Run tiles [0, N) with explicitly non-global FDR outputs",
    )
    return parser.parse_args()


def _write_manifest(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.incomplete")
    pd.DataFrame(rows).to_csv(temporary, index=False)
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    _, runs = load_phase1_runs(
        args.run_matrix, variable_selection=args.variable_selection
    )
    if args.only:
        unknown = set(args.only).difference(run["run_id"] for run in runs)
        if unknown:
            raise ValueError(f"Unknown --only run IDs: {sorted(unknown)}")
        runs = [run for run in runs if run["run_id"] in set(args.only)]
    if args.role:
        runs = [run for run in runs if run["role"] == args.role]
    if not runs:
        raise ValueError("No Phase 1 runs selected")
    if args.n_jobs == 0:
        raise ValueError("--n-jobs cannot be zero")

    contract = json.loads(args.input_contract.read_text())
    total_tiles = int(contract["n_tiles"])
    diagnostic = args.diagnostic_tile_stop is not None
    expected_tiles = int(args.diagnostic_tile_stop or total_tiles)
    if expected_tiles <= 0 or expected_tiles > total_tiles:
        raise ValueError(
            f"Invalid --diagnostic-tile-stop {expected_tiles} for {total_tiles} tiles"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / (
        "phase1_trend_diagnostic_manifest.csv"
        if diagnostic
        else "phase1_trend_batch_manifest.csv"
    )
    rows: list[dict[str, Any]] = []
    failures = 0

    print(
        f"Phase 1 trend batch: {len(runs)} runs, {expected_tiles} tiles each, "
        f"n_jobs={args.n_jobs}, parallel_prefer={args.parallel_prefer}, "
        f"diagnostic={diagnostic}",
        flush=True,
    )
    for index, run in enumerate(runs, start=1):
        filename = output_filename(run)
        if diagnostic:
            filename = filename.replace(
                "_trend_statistics.nc",
                f"_diagnostic_tiles_0_{expected_tiles}_trend_statistics.nc",
            )
        output_path = args.output_dir / filename
        prior = audit_trend_output(
            output_path,
            run,
            expected_tiles=expected_tiles,
            expected_total_input_tiles=total_tiles,
            require_global_fdr=not diagnostic,
        )
        if prior["status"] == "PASS" and not args.force:
            print(f"[{index}/{len(runs)}] SKIP {run['run_id']} (existing output passes)", flush=True)
            rows.append(
                {
                    **prior,
                    "action": "skipped_existing",
                    "started_utc": "",
                    "duration_seconds": 0.0,
                    "returncode": 0,
                }
            )
            _write_manifest(rows, manifest_path)
            continue

        command = [
            sys.executable,
            str(BUILDER),
            "--variable",
            run["variable"],
            "--series",
            run["series"],
            "--mask",
            run["mask"],
            "--data-dir",
            str(args.data_dir),
            "--input-contract",
            str(args.input_contract),
            "--variable-selection",
            str(args.variable_selection),
            "--config",
            str(args.trend_config),
            "--output",
            str(output_path),
            "--n-jobs",
            str(args.n_jobs),
            "--parallel-prefer",
            args.parallel_prefer,
            "--run-id",
            run["run_id"],
            "--run-role",
            run["role"],
            "--run-matrix",
            str(args.run_matrix),
        ]
        if diagnostic:
            command.extend(["--tile-stop", str(expected_tiles)])
        print(f"[{index}/{len(runs)}] RUN  {run['run_id']}", flush=True)
        if args.dry_run:
            print("  " + " ".join(command), flush=True)
            rows.append(
                {
                    **prior,
                    "action": "dry_run",
                    "started_utc": "",
                    "duration_seconds": 0.0,
                    "returncode": 0,
                }
            )
            continue

        started_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
        started = time.monotonic()
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        duration = time.monotonic() - started
        log_path = output_path.with_suffix(".log")
        log_path.write_text(
            "$ " + " ".join(command) + "\n\nSTDOUT\n" + completed.stdout
            + "\nSTDERR\n" + completed.stderr
        )
        current = audit_trend_output(
            output_path,
            run,
            expected_tiles=expected_tiles,
            expected_total_input_tiles=total_tiles,
            require_global_fdr=not diagnostic,
        )
        action = "completed" if completed.returncode == 0 and current["status"] == "PASS" else "failed"
        rows.append(
            {
                **current,
                "action": action,
                "started_utc": started_utc,
                "duration_seconds": round(duration, 3),
                "returncode": completed.returncode,
                "log_path": str(log_path),
            }
        )
        _write_manifest(rows, manifest_path)
        if action == "failed":
            failures += 1
            print(
                f"  FAIL after {duration / 60.0:.1f} min: {current['detail']} "
                f"(see {log_path})",
                flush=True,
            )
        else:
            print(
                f"  PASS in {duration / 60.0:.1f} min; "
                f"success={current['n_success']:,}, "
                f"FDR significant={current['n_significant_fdr']:,}",
                flush=True,
            )

    if args.dry_run:
        print(f"Dry run complete; {len(runs)} commands checked")
        return 0
    if failures:
        print(f"Batch finished with {failures} failed runs; rerun to resume")
        return 1
    if not diagnostic and not args.only and args.role is None:
        final_report = audit_phase1_outputs(
            run_matrix=args.run_matrix,
            variable_selection=args.variable_selection,
            input_contract=args.input_contract,
            output_dir=args.output_dir,
        )
        final_report_path = args.output_dir / "phase1_trend_output_audit.csv"
        final_report.to_csv(final_report_path, index=False)
        n_pass = int((final_report["status"] == "PASS").sum())
        print(f"Final production audit: {n_pass}/{len(final_report)} PASS")
        if n_pass != len(final_report):
            return 1
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
