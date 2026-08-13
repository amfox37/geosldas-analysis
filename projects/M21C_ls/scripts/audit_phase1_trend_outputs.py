#!/usr/bin/env python3
"""Audit all production NetCDFs declared by the M21C Phase 1 trend matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

from phase1_trend_workflow import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_RUN_MATRIX,
    audit_phase1_outputs,
)
from trend_breakpoint_series import DEFAULT_INPUT_CONTRACT, DEFAULT_VARIABLE_SELECTION
from trend_statistics import DEFAULT_CONFIG


DEFAULT_REPORT = DEFAULT_OUTPUT_DIR / "phase1_trend_output_audit.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-matrix", type=Path, default=DEFAULT_RUN_MATRIX)
    parser.add_argument(
        "--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION
    )
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument("--trend-config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = audit_phase1_outputs(
        run_matrix=args.run_matrix,
        variable_selection=args.variable_selection,
        input_contract=args.input_contract,
        output_dir=args.output_dir,
        trend_config=args.trend_config,
    )
    if not args.no_write:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(args.report, index=False)
        print(f"Wrote {args.report}")
    passed = int((report["status"] == "PASS").sum())
    print(f"Phase 1 output audit: {passed}/{len(report)} PASS")
    failed = report.loc[report["status"] != "PASS", ["run_id", "detail"]]
    if not failed.empty:
        print(failed.to_string(index=False))
        return 1
    print(
        report[
            [
                "run_id",
                "n_success",
                "n_significant_fdr",
                "fraction_significant_fdr",
                "fraction_area_significant_fdr",
                "n_significant_positive_slope",
                "n_significant_negative_slope",
                "n_significant_zero_slope",
                "median_slope_success",
                "n_fdr_ci_disagreement",
                "fraction_significant_with_ci_disagreement",
            ]
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
