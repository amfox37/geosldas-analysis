#!/usr/bin/env python3
"""Audit M21C Phase 1 area-weighted interrupted-series outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from interrupted_time_series import DEFAULT_CONFIG
from phase1_interrupted_workflow import (
    DEFAULT_AUDIT,
    DEFAULT_COEFFICIENTS,
    DEFAULT_MODELS,
    DEFAULT_MONTHLY,
    audit_interrupted_outputs,
)
from phase1_trend_workflow import DEFAULT_RUN_MATRIX
from trend_breakpoint_series import DEFAULT_VARIABLE_SELECTION


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", type=Path, default=DEFAULT_COEFFICIENTS)
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--monthly", type=Path, default=DEFAULT_MONTHLY)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--run-matrix", type=Path, default=DEFAULT_RUN_MATRIX)
    parser.add_argument("--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = audit_interrupted_outputs(
        coefficients_path=args.coefficients,
        models_path=args.models,
        monthly_path=args.monthly,
        run_matrix=args.run_matrix,
        variable_selection=args.variable_selection,
        config_path=args.config,
    )
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.audit, index=False)
    print(report.to_string(index=False))
    return 0 if report.iloc[0]["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
