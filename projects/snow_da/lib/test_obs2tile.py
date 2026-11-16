"""
Tests for lib/obs2tile.py (cube-sphere cell lookup of GEOSldas_GridComp PR #252).

run:  /Users/amfox/mamba/envs/regrid/bin/python -m pytest projects/snow_da/lib/test_obs2tile.py
The full check against the PR's Fortran is c3_offline_harness/offline_map.x (see the report).
"""

import numpy as np

from .obs2tile import cs_ij


def test_cs_ij_range_and_faces():
    rng = np.random.default_rng(0)
    lon, lat = rng.uniform(-180, 180, 100000), rng.uniform(-90, 90, 100000)
    i, j = cs_ij(lon, lat, 360)
    assert i.min() >= 1 and i.max() <= 360 and j.min() >= 1 and j.max() <= 2160
    assert set(np.unique((j - 1) // 360)) == set(range(6))     # all six faces used


def test_cs_ij_poles_on_faces_3_and_6():
    i, j = cs_ij(np.array([0., 0.]), np.array([89.99, -89.99]), 360)
    assert (j[0] - 1) // 360 == 2 and (j[1] - 1) // 360 == 5
