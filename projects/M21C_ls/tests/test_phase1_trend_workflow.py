from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from phase1_trend_workflow import audit_trend_output, load_phase1_runs, load_trend_runs
from trend_statistics import compute_trend_statistics


STATE_CONTEXT_MATRIX = PROJECT_ROOT / "config" / "phase1_state_trend_runs.json"
PHASE2_PRIMARY_MATRIX = (
    PROJECT_ROOT / "config" / "phase2_flux_storage_trend_runs.json"
)
PHASE2_CONTEXT_MATRIX = (
    PROJECT_ROOT / "config" / "phase2_flux_storage_state_trend_runs.json"
)


def test_phase1_matrix_has_one_primary_run_per_selected_variable() -> None:
    payload, runs = load_phase1_runs()

    primary = [run for run in runs if run["role"] == "primary"]
    sensitivity = [run for run in runs if run["role"] == "sensitivity"]
    assert payload["expected_run_count"] == 21
    assert len(primary) == 17
    assert len(sensitivity) == 4
    assert {run["variable"] for run in sensitivity} == {
        "SFMC",
        "RZMC",
        "SNOMASLAND",
        "SNODPLAND",
    }


def test_state_context_matrix_has_matched_ol_and_da_runs() -> None:
    payload, runs = load_phase1_runs(STATE_CONTEXT_MATRIX)

    assert payload["matrix_kind"] == "state_context"
    assert len(runs) == 12
    assert {run["role"] for run in runs} == {"context"}
    signatures = {(run["variable"], run["series"], run["mask"]) for run in runs}
    for variable, mask in {
        "PRECTOTCORRLAND": "valid_land",
        "SFMC": "valid_land",
        "RZMC": "valid_land",
        "SNOMASLAND": "seasonal_snow",
        "SNODPLAND": "seasonal_snow",
        "FRLANDSNO": "seasonal_snow",
    }.items():
        assert (variable, "ol", mask) in signatures
        assert (variable, "da", mask) in signatures


def test_focused_phase2_matrices_cover_delta_ol_and_da() -> None:
    primary_payload, primary = load_trend_runs(PHASE2_PRIMARY_MATRIX)
    context_payload, context = load_trend_runs(PHASE2_CONTEXT_MATRIX)

    variables = {"EVLAND", "TOTAL_RUNOFF", "TWLAND"}
    assert primary_payload["phase"] == context_payload["phase"] == 2
    assert len(primary) == 3
    assert len(context) == 6
    assert {(run["variable"], run["series"]) for run in primary} == {
        (variable, "delta") for variable in variables
    }
    assert {(run["variable"], run["series"]) for run in context} == {
        (variable, series) for variable in variables for series in ("ol", "da")
    }


def _trend_output(run: dict[str, str]) -> xr.Dataset:
    time = pd.date_range("2000-06-01", periods=288, freq="MS")
    years = np.arange(time.size) / 12.0
    values = xr.DataArray(
        np.column_stack([0.01 * years, -0.01 * years]),
        dims=("time", "tile"),
        coords={
            "time": time,
            "tile": [0, 1],
            "lat": ("tile", [40.0, 50.0]),
            "lon": ("tile", [-110.0, -100.0]),
            "tile_area": ("tile", [1.0, 2.0]),
        },
        name=run["series"],
        attrs={
            "analysis_variable": run["variable"],
            "source_variable": run["variable"],
            "units": "test-unit",
        },
    )
    output = compute_trend_statistics(values)
    output.attrs.update(
        mask=run["mask"],
        selected_series=run["series"],
        tile_start=0,
        tile_stop=2,
        total_input_tiles=2,
        fdr_scope="all finite tiles eligible under the selected mask",
        run_id=run["run_id"],
        run_role=run["role"],
        source_git_commit="synthetic-test-commit",
    )
    return output


def test_output_audit_accepts_complete_global_field(tmp_path: Path) -> None:
    _, runs = load_phase1_runs()
    run = runs[0]
    path = tmp_path / "trend.nc"
    _trend_output(run).to_netcdf(path)

    summary = audit_trend_output(
        path,
        run,
        expected_tiles=2,
        expected_total_input_tiles=2,
        require_global_fdr=True,
    )

    assert summary["status"] == "PASS"
    assert summary["n_success"] == 2
    assert summary["n_significant_fdr"] == 2
    assert summary["fraction_area_significant_fdr"] == 1.0
    assert summary["n_significant_positive_slope"] == 1
    assert summary["n_significant_negative_slope"] == 1
    assert summary["source_git_commit"] == "synthetic-test-commit"
    assert summary["fraction_significant_with_ci_disagreement"] == 0.0


def test_output_audit_rejects_wrong_embedded_configuration(tmp_path: Path) -> None:
    _, runs = load_phase1_runs()
    run = runs[0]
    path = tmp_path / "trend.nc"
    _trend_output(run).to_netcdf(path)

    summary = audit_trend_output(
        path,
        run,
        expected_tiles=2,
        expected_total_input_tiles=2,
        require_global_fdr=True,
        expected_configuration="not-the-embedded-configuration",
    )

    assert summary["status"] == "FAIL"
    assert "embedded trend configuration differs" in summary["detail"]


def test_output_audit_rejects_subset_as_production(tmp_path: Path) -> None:
    _, runs = load_phase1_runs()
    run = runs[0]
    output = _trend_output(run)
    output.attrs["total_input_tiles"] = 112573
    output.attrs["fdr_scope"] = "diagnostic selected tile subset; not global FDR"
    path = tmp_path / "subset.nc"
    output.to_netcdf(path)

    summary = audit_trend_output(
        path,
        run,
        expected_tiles=2,
        expected_total_input_tiles=112573,
        require_global_fdr=True,
    )

    assert summary["status"] == "FAIL"
    assert "non-global FDR scope" in summary["detail"]
