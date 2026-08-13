#!/usr/bin/env python3
"""Summarize matched OL, DA, and DA-minus-OL state trends."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from phase1_trend_workflow import (
    DEFAULT_OUTPUT_DIR,
    audit_phase1_outputs,
    load_phase1_runs,
    output_filename,
)
from trend_breakpoint_series import (
    DEFAULT_DATA_DIR,
    DEFAULT_INPUT_CONTRACT,
    DEFAULT_VARIABLE_SELECTION,
    MonthlySeriesLoader,
)
from trend_statistics import (
    DEFAULT_CONFIG,
    analyze_monthly_series,
    benjamini_hochberg,
    load_trend_config,
    monthly_time_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTEXT_MATRIX = PROJECT_ROOT / "config" / "phase1_state_trend_runs.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-matrix", type=Path, default=DEFAULT_CONTEXT_MATRIX)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument(
        "--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION
    )
    parser.add_argument("--trend-config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def _load_result(path: Path) -> xr.Dataset:
    if not path.exists():
        raise FileNotFoundError(f"Missing production trend field: {path}")
    with xr.open_dataset(path) as dataset:
        return dataset.load()


def main() -> int:
    args = parse_args()
    _, context_runs = load_phase1_runs(
        args.context_matrix, variable_selection=args.variable_selection
    )
    audit = audit_phase1_outputs(
        run_matrix=args.context_matrix,
        variable_selection=args.variable_selection,
        input_contract=args.input_contract,
        output_dir=args.output_dir,
        trend_config=args.trend_config,
    )
    if not (audit["status"] == "PASS").all():
        failed = audit.loc[audit["status"] != "PASS", ["run_id", "detail"]]
        raise RuntimeError(f"Context trend audit failed:\n{failed.to_string(index=False)}")

    settings = load_trend_config(args.trend_config)
    run_lookup = {
        (run["variable"], run["series"]): run for run in context_runs
    }
    variable_masks = {run["variable"]: run["mask"] for run in context_runs}
    variables = list(dict.fromkeys(run["variable"] for run in context_runs))
    tile_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    domain_rows: list[dict[str, object]] = []

    with MonthlySeriesLoader(
        args.data_dir,
        input_contract=args.input_contract,
        variable_selection=args.variable_selection,
    ) as loader:
        for variable in variables:
            mask = variable_masks[variable]
            results: dict[str, xr.Dataset] = {}
            for series in ("ol", "da"):
                run = run_lookup[(variable, series)]
                results[series] = _load_result(args.output_dir / output_filename(run))
            results["delta"] = _load_result(
                args.output_dir / f"{variable}_delta_{mask}_trend_statistics.nc"
            )

            for series, result in results.items():
                success = result["status"].values == 0
                significant = result["significant_fdr"].values.astype(bool)
                slopes = result["slope"].values
                area = result["tile_area"].values
                tile_rows.append(
                    {
                        "variable": variable,
                        "series": series,
                        "mask": mask,
                        "n_supported": int(success.sum()),
                        "n_significant_fdr": int(significant.sum()),
                        "fraction_supported_significant": float(
                            significant.sum() / success.sum()
                        ),
                        "fraction_supported_area_significant": float(
                            area[significant].sum() / area[success].sum()
                        ),
                        "n_significant_positive": int(
                            (significant & (slopes > 0)).sum()
                        ),
                        "n_significant_negative": int(
                            (significant & (slopes < 0)).sum()
                        ),
                    }
                )

            ol_sig = results["ol"]["significant_fdr"].values.astype(bool)
            da_sig = results["da"]["significant_fdr"].values.astype(bool)
            ol_slope = results["ol"]["slope"].values
            da_slope = results["da"]["slope"].values
            common = ol_sig & da_sig
            finite = np.isfinite(ol_slope) & np.isfinite(da_slope)
            comparison_rows.append(
                {
                    "variable": variable,
                    "mask": mask,
                    "n_ol_da_significant_overlap": int(common.sum()),
                    "n_ol_da_significant_union": int((ol_sig | da_sig).sum()),
                    "n_overlap_same_direction": int(
                        (common & (np.sign(ol_slope) == np.sign(da_slope))).sum()
                    ),
                    "n_overlap_opposite_direction": int(
                        (common & (np.sign(ol_slope) != np.sign(da_slope))).sum()
                    ),
                    "ol_da_slope_correlation": float(
                        np.corrcoef(ol_slope[finite], da_slope[finite])[0, 1]
                    ),
                }
            )

            monthly = loader.domain_mean(variable, mask=mask).load()
            time_years, months = monthly_time_axis(monthly.time.values)
            for series in ("ol", "da", "delta"):
                values = monthly[series].values
                estimate = analyze_monthly_series(
                    values,
                    np.isfinite(values),
                    time_years,
                    months,
                    settings,
                )
                domain_rows.append(
                    {
                        "variable": variable,
                        "series": series,
                        "mask": mask,
                        "units": monthly[series].attrs.get("units", ""),
                        "slope_per_year": estimate["slope"],
                        "slope_ci_low": estimate["slope_ci_low"],
                        "slope_ci_high": estimate["slope_ci_high"],
                        "p_value": estimate["p_value"],
                        "mk_variance_factor": estimate["mk_variance_factor"],
                    }
                )

    domain = pd.DataFrame(domain_rows)
    q_value, significant = benjamini_hochberg(
        domain["p_value"].to_numpy(), alpha=float(settings["fdr"]["alpha"])
    )
    domain["q_value"] = q_value
    domain["significant_fdr"] = significant

    outputs = {
        "phase1_state_trend_tile_summary.csv": pd.DataFrame(tile_rows),
        "phase1_state_trend_pair_comparison.csv": pd.DataFrame(comparison_rows),
        "phase1_state_trend_domain_summary.csv": domain,
    }
    for filename, frame in outputs.items():
        path = args.output_dir / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.incomplete")
        frame.to_csv(temporary, index=False)
        temporary.replace(path)
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
