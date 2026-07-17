"""Regression tests for the vectorized H121 time-axis handling in readers.py.

The original implementation converted every raw observation's CF time value
via cftime.num2date(...) followed by a per-observation Python loop pulling
.hour/.minute off each resulting object. At H121's ~1M raw obs/day, that was
the dominant cost of reading a day of files (~34s of ~65s in a profiled
sample day, almost all of it numpy.ma.core.__getitem__ overhead from
num2date's element-wise handling of a masked input array) -- far more than
the actual QC or tile-binning work downstream. _hour_of_day_from_time_var
replaces it with plain vectorized arithmetic on the numeric time axis.
"""

from pathlib import Path
import sys

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

ASCAT_ROOT = Path(__file__).resolve().parents[1]
if str(ASCAT_ROOT) not in sys.path:
    sys.path.insert(0, str(ASCAT_ROOT))

from lib.readers import _hour_of_day_from_time_var  # noqa: E402


def test_hour_of_day_matches_num2date_to_the_minute():
    # Real H121 units: "days since 1970-01-01 00:00:00" (a midnight reference).
    units = "days since 1970-01-01 00:00:00"
    # A handful of fractional-day offsets, including one just past midnight
    # and one just before it, to exercise the modulo wraparound.
    t_vals = np.array([17809.0, 17809.25, 17809.5, 17809.99998, 17809.00000669])

    import netCDF4 as nc

    times = nc.num2date(t_vals, units)
    expected = np.array([t.hour + t.minute / 60.0 for t in times])

    actual = _hour_of_day_from_time_var(t_vals, units)

    # The original code only ever used whole .hour/.minute (no .second), so
    # the vectorized version -- which keeps full sub-minute precision -- can
    # differ by up to a minute from that truncated reference.
    assert np.max(np.abs(actual - expected)) < 1.0 / 60.0 + 1e-9


def test_hour_of_day_handles_non_midnight_reference():
    # A reference time-of-day that isn't 00:00:00 must be added in, not
    # dropped -- otherwise every hour-of-day would be shifted by a constant.
    units = "seconds since 2000-01-01 12:00:00"
    t_vals = np.array([0.0, 3600.0, 12 * 3600.0])

    hour_frac = _hour_of_day_from_time_var(t_vals, units)

    np.testing.assert_allclose(hour_frac, [12.0, 13.0, 0.0], atol=1e-9)


def test_hour_of_day_wraps_multi_day_offsets_to_time_of_day_only():
    units = "hours since 2015-06-01 00:00:00"
    t_vals = np.array([0.0, 25.0, 49.5, -1.0])

    hour_frac = _hour_of_day_from_time_var(t_vals, units)

    np.testing.assert_allclose(hour_frac, [0.0, 1.0, 1.5, 23.0], atol=1e-9)


def test_hour_of_day_rejects_unrecognized_units():
    with pytest.raises(ValueError):
        _hour_of_day_from_time_var(np.array([0.0]), "furlongs since yesterday")
