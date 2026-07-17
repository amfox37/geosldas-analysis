"""Regression tests for the M36 grid-origin constants in superob.py.

The M36 grid is a fixed (964, 406) cell-center lattice symmetric about its
own center; it does not extend all the way to the geographic pole. An
earlier version derived X_ORIGIN/Y_ORIGIN by reprojecting the pole itself
(EASE2(-180, 90)), which put Y_ORIGIN 0.768 cells short of the grid's true
top-row edge -- close enough to be "correct by accident" for cell-center
queries, but wrong for arbitrary scattered points (see
project_ascat_da_superob_latlon_bug in the auto-memory for the real-tile and
real-granule evidence that motivated this fix).
"""

import numpy as np
import pytest

pytest.importorskip("pyproj")

from superob import EASE2, M36_CELL, X_ORIGIN, Y_ORIGIN, latlon_to_ij  # noqa: E402


def test_grid_origin_is_exact_cell_edge_not_reprojected_pole():
    # The true left/top edges are exactly -482 and 203 whole cells from the
    # grid center -- integers, unlike the reprojected pole's fractional
    # 203.768...-cell position.
    assert X_ORIGIN == pytest.approx(-482 * M36_CELL, abs=1e-6)
    assert Y_ORIGIN == pytest.approx(203 * M36_CELL, abs=1e-6)

    pole_x, pole_y = EASE2(-180, 90)
    assert pole_x == pytest.approx(X_ORIGIN, abs=1e-6)  # x is longitude-only, unaffected
    assert (pole_y - Y_ORIGIN) / M36_CELL == pytest.approx(0.7684595957222012, rel=1e-9)


def test_latlon_to_ij_round_trips_arbitrary_cell_centers():
    # Cell centers scattered across the full grid (not just near the
    # equator, where a small origin error is least likely to flip a row).
    nx, ny = 964, 406
    r0, s0 = (nx - 1) / 2.0, (ny - 1) / 2.0
    i_true = np.array([0, 1, 480, 481, 500, 700, 962, 963])
    j_true = np.array([0, 1, 200, 202, 205, 300, 404, 405])

    x = (i_true - r0) * M36_CELL
    y = (s0 - j_true) * M36_CELL
    lon, lat = EASE2(x, y, inverse=True)

    i_back, j_back = latlon_to_ij(lat, lon)

    np.testing.assert_array_equal(i_back, i_true)
    np.testing.assert_array_equal(j_back, j_true)


def test_latlon_to_ij_round_trips_off_center_points_within_a_cell():
    # Points away from cell centers (e.g. near an edge) are the case the
    # reprojected-pole origin got wrong for ~77% of a real H121 granule --
    # a fractional-cell origin error doesn't only show up at centers.
    nx, ny = 964, 406
    r0, s0 = (nx - 1) / 2.0, (ny - 1) / 2.0
    i_true = np.full(9, 500, dtype=np.int64)
    j_true = np.full(9, 200, dtype=np.int64)
    frac = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.95])

    x = (i_true - r0 + (frac - 0.5)) * M36_CELL
    y = (s0 - j_true - (frac - 0.5)) * M36_CELL
    lon, lat = EASE2(x, y, inverse=True)

    i_back, j_back = latlon_to_ij(lat, lon)

    np.testing.assert_array_equal(i_back, i_true)
    np.testing.assert_array_equal(j_back, j_true)
