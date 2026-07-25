"""Tests for the direct ISMN .stm archive readers in ismn_io.py.

The `ismn` package is unavailable on Discover, so these readers parse the
archive format themselves; these tests pin the QC, depth-selection, and
root-zone rules that the skill workflow depends on.
"""

import numpy as np
import pandas as pd
import pytest

from ismn_io import (
    build_station_series,
    depth_layer_weights,
    read_stm_series,
    station_key,
    weighted_rootzone,
)


HEADER = "SCAN       SCAN       Testville       34.78534 -86.55906                 241.0 {df} {dt} Probe\n"


def write_stm(path, depth_from, depth_to, rows):
    lines = [HEADER.format(df=depth_from, dt=depth_to)]
    lines.extend(f"{stamp} {value} {flag} V\n" for stamp, value, flag in rows)
    path.write_text("".join(lines))
    return path


def hourly_rows(start, hours, value, flag="G"):
    stamps = pd.date_range(start, periods=hours, freq="h")
    return [(s.strftime("%Y/%m/%d %H:%M"), value, flag) for s in stamps]


def test_read_stm_series_keeps_only_good_flags(tmp_path):
    path = write_stm(
        tmp_path / "a.stm",
        "0.0500",
        "0.0500",
        [
            ("2018/01/01 00:00", 0.30, "G"),
            ("2018/01/01 01:00", 0.99, "C03"),
            ("2018/01/01 02:00", 0.32, "G"),
            ("2018/01/01 03:00", -0.70, "C01,D01"),
        ],
    )

    series = read_stm_series(path)

    assert list(series.to_numpy()) == [0.30, 0.32]
    assert series.index[0] == pd.Timestamp("2018-01-01 00:00")


def test_read_stm_series_averages_duplicate_timestamps(tmp_path):
    path = write_stm(
        tmp_path / "a.stm",
        "0.0500",
        "0.0500",
        [("2018/01/01 00:00", 0.20, "G"), ("2018/01/01 00:00", 0.40, "G")],
    )

    series = read_stm_series(path)

    assert len(series) == 1
    assert series.iloc[0] == pytest.approx(0.30)


def test_depth_layer_weights_span_the_full_profile():
    weights = depth_layer_weights([0.05, 0.20, 0.50], z_top=0.0, z_bottom=1.0)

    # Layer edges fall midway between neighbours; the ends clamp to the profile.
    assert weights[0.05] == pytest.approx(0.125)
    assert weights[0.20] == pytest.approx(0.225)
    assert weights[0.50] == pytest.approx(0.65)
    assert sum(weights.values()) == pytest.approx(1.0)


def test_weighted_rootzone_requires_shallow_and_deep_coverage():
    index = pd.date_range("2018-01-01", periods=3, freq="D")
    frame = pd.DataFrame(
        {
            0.05: [0.30, 0.30, np.nan],
            0.20: [0.25, 0.25, 0.25],
            0.60: [0.20, np.nan, 0.20],
        },
        index=index,
    )

    series = weighted_rootzone(frame, min_sensors=3, shallow_max=0.20, deep_min=0.50)

    # Row 0 has all three layers; rows 1 and 2 each drop below three finite
    # layers, so no partial profile is passed off as root zone.
    assert np.isfinite(series.iloc[0])
    assert np.isnan(series.iloc[1])
    assert np.isnan(series.iloc[2])


def test_build_station_series_picks_nearest_shallow_sensor(tmp_path):
    rows = pd.DataFrame(
        [
            {"depth_center": 0.05, "path": write_stm(tmp_path / "s05.stm", "0.05", "0.05", hourly_rows("2018-01-01", 48, 0.30))},
            {"depth_center": 0.10, "path": write_stm(tmp_path / "s10.stm", "0.10", "0.10", hourly_rows("2018-01-01", 48, 0.40))},
            {"depth_center": 0.60, "path": write_stm(tmp_path / "s60.stm", "0.60", "0.60", hourly_rows("2018-01-01", 48, 0.50))},
        ]
    )

    pack = build_station_series(rows, start="2018-01-01", end="2018-01-02")

    assert pack["ok"]
    assert pack["surface_depth_m"] == pytest.approx(0.05)
    assert pack["surface_daily"].iloc[0] == pytest.approx(0.30)
    # The 12-hour shift aligns the in-situ day with the model day.
    assert pack["surface_daily"].index[0] == pd.Timestamp("2018-01-01 12:00")


def test_build_station_series_rejects_deep_only_station(tmp_path):
    rows = pd.DataFrame(
        [
            {"depth_center": 0.50, "path": write_stm(tmp_path / "d50.stm", "0.50", "0.50", hourly_rows("2018-01-01", 48, 0.28))},
            {"depth_center": 1.00, "path": write_stm(tmp_path / "d100.stm", "1.00", "1.00", hourly_rows("2018-01-01", 48, 0.31))},
        ]
    )

    pack = build_station_series(rows, start="2018-01-01", end="2018-01-02")

    # A 0.5 m sensor must never stand in for surface soil moisture.
    assert pack["ok"]
    assert np.isnan(pack["surface_depth_m"])
    assert pack["surface_daily"].empty


def test_build_station_series_trims_to_window(tmp_path):
    rows = pd.DataFrame(
        [
            {
                "depth_center": 0.05,
                "path": write_stm(
                    tmp_path / "s.stm",
                    "0.05",
                    "0.05",
                    hourly_rows("2017-12-30", 24, 0.10) + hourly_rows("2018-01-01", 24, 0.30),
                ),
            }
        ]
    )

    pack = build_station_series(rows, start="2018-01-01", end="2018-01-01")

    assert len(pack["surface_daily"]) == 1
    assert pack["surface_daily"].iloc[0] == pytest.approx(0.30)


def test_station_key_disambiguates_networks():
    assert station_key("SCAN", "Alpha") != station_key("SNOTEL", "Alpha")
