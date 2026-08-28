from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from phase1_interrupted_workflow import assign_transition_fdr, phase1_series_catalogue


PHASE2_MATRIX = PROJECT_ROOT / "config" / "phase2_flux_storage_trend_runs.json"


def test_phase1_catalogue_expands_to_43_domain_series() -> None:
    catalogue = phase1_series_catalogue()

    assert len(catalogue) == 43
    assert catalogue["series_id"].is_unique
    assert (catalogue["series_role"] == "primary_estimand").sum() == 21
    assert (catalogue["series_role"] == "paired_control").sum() == 22
    assert set(catalogue.groupby("run_id").size()) == {1, 3}


def test_phase2_catalogue_expands_three_paired_variables_to_nine_series() -> None:
    catalogue = phase1_series_catalogue(PHASE2_MATRIX)

    assert len(catalogue) == 9
    assert catalogue["series_id"].is_unique
    assert set(catalogue["variable"]) == {"EVLAND", "TOTAL_RUNOFF", "TWLAND"}
    assert set(catalogue.groupby("run_id").size()) == {3}


def test_fdr_is_separate_by_boundary_and_excludes_p7_slope() -> None:
    catalogue = phase1_series_catalogue()
    rows: list[dict[str, object]] = []
    for coefficient, coefficient_type, independent in (
        ("level_change_P7", "level_change", True),
        ("slope_change_P8", "slope_change", True),
        ("period_slope_P7", "period_slope", False),
        ("intercept", "nuisance", True),
    ):
        for index, series_id in enumerate(catalogue["series_id"]):
            rows.append(
                {
                    "series_id": series_id,
                    "coefficient": coefficient,
                    "coefficient_type": coefficient_type,
                    "independently_estimated": independent,
                    "p_value": 0.001 if index == 0 else 0.5,
                }
            )
    output = assign_transition_fdr(pd.DataFrame(rows), alpha=0.05)

    families = output.loc[output["fdr_family"] != ""].groupby("fdr_family").size()
    assert families.to_dict() == {"level_change_P7": 43, "slope_change_P8": 43}
    assert not output.loc[
        output["coefficient"] == "period_slope_P7", "q_value"
    ].notna().any()
    assert not output.loc[
        output["coefficient"] == "intercept", "q_value"
    ].notna().any()
    assert np.isclose(
        output.loc[
            (output["coefficient"] == "level_change_P7")
            & (output["p_value"] == 0.001),
            "q_value",
        ].iloc[0],
        0.043,
    )
