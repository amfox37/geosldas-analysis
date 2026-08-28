from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import plot_trends_breakpoints_manuscript_figures as figures  # noqa: E402


def test_periods_come_from_registry_and_preserve_p6_p9_dates():
    periods = figures.load_periods()
    assert periods["period_id"].tolist() == [f"P{i}" for i in range(1, 10)]
    assert periods.set_index("period_id").loc["P6", "start"] == np.datetime64("2015-04-01")
    assert periods.set_index("period_id").loc["P9", "start"] == np.datetime64("2021-12-01")


def test_segmented_diverging_scale_is_symmetric_with_white_zero_bin():
    cmap, norm, bounds = figures.segmented_diverging_scale(2.0)
    assert np.allclose(bounds, -bounds[::-1])
    assert np.allclose(cmap.colors[cmap.N // 2], (1, 1, 1, 1))
    assert norm(-2.0) == 0
    assert norm(2.0) == cmap.N


def test_main_figure_uses_variable_specific_delta_scales():
    assert figures.MAIN_TREND_DELTA_SCALES == {"RZMC": "separate", "FRLANDSNO": "state"}


def test_trend_rows_rejects_unknown_delta_scale():
    with pytest.raises(ValueError, match="delta_scales"):
        figures.plot_trend_rows(["RZMC"], "unused", delta_scales={"RZMC": "magnify"})
