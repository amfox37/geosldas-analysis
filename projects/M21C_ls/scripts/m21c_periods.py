#!/usr/bin/env python3
"""Load and validate the shared M21C observing-system period registry."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


DEFAULT_REGISTRY = (
    Path(__file__).resolve().parents[1] / "config" / "observing_system_registry.json"
)


def _month_count(start: pd.Timestamp, end: pd.Timestamp) -> int:
    return len(pd.period_range(start=start, end=end, freq="M"))


def _validate_contiguous(df: pd.DataFrame, name: str) -> None:
    for previous, current in zip(df.iloc[:-1].itertuples(), df.iloc[1:].itertuples()):
        expected = previous.end + pd.Timedelta(days=1)
        if current.start != expected:
            raise ValueError(
                f"{name} is not contiguous between {previous.period_id} and "
                f"{current.period_id}: expected {expected.date()}, got {current.start.date()}"
            )


def load_period_frames(
    registry_path: str | Path = DEFAULT_REGISTRY,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return validated registry metadata and fine, validation, and sensor frames."""

    path = Path(registry_path)
    registry = json.loads(path.read_text())
    fine = pd.DataFrame(registry["fine_periods"])
    validation = pd.DataFrame(registry["validation_periods"])
    sensors = pd.DataFrame(registry["sensors"])

    for frame in (fine, validation, sensors):
        frame["start"] = pd.to_datetime(frame["start"])
        frame["end"] = pd.to_datetime(frame["end"])

    for frame, name in ((fine, "fine periods"), (validation, "validation periods")):
        if frame["period_id"].duplicated().any():
            raise ValueError(f"Duplicate IDs in {name}")
        computed_months = frame.apply(lambda row: _month_count(row["start"], row["end"]), axis=1)
        if not computed_months.equals(frame["n_months"]):
            bad = frame.loc[computed_months != frame["n_months"], ["period_id", "n_months"]]
            raise ValueError(f"Stored month counts do not match dates in {name}:\n{bad}")
        frame["n_days_inclusive"] = (frame["end"] - frame["start"]).dt.days + 1
        _validate_contiguous(frame, name)

    paper_start = pd.Timestamp(registry["paper_start"])
    paper_end = pd.Timestamp(registry["paper_end"])
    for frame, name in ((fine, "fine periods"), (validation, "validation periods")):
        if frame["start"].iloc[0] != paper_start or frame["end"].iloc[-1] != paper_end:
            raise ValueError(f"{name} does not cover the complete paper period")

    slope_floor = int(registry["slope_min_segment_months"])
    expected_reliability = fine["n_months"] >= slope_floor
    if not expected_reliability.equals(fine["reliable_for_slope"]):
        raise ValueError("reliable_for_slope flags do not match the configured slope floor")

    changepoint_floor = int(registry["changepoint_min_segment_months"])
    too_short = fine["n_months"] < changepoint_floor
    if not fine.loc[too_short, "changepoint_detection_exempt"].all():
        raise ValueError("Every sub-floor period must be exempt from changepoint scoring")

    fine_lookup = fine.set_index("period_id")
    listed_members: list[str] = []
    for row in validation.itertuples():
        members = list(row.fine_periods)
        listed_members.extend(members)
        member_rows = fine_lookup.loc[members]
        if member_rows["start"].iloc[0] != row.start or member_rows["end"].iloc[-1] != row.end:
            raise ValueError(f"{row.period_id} boundaries do not match its fine periods")
        if not (member_rows["validation_id"] == row.period_id).all():
            raise ValueError(f"{row.period_id} membership conflicts with fine-period validation_id")
    if listed_members != fine["period_id"].tolist():
        raise ValueError("Validation-period members do not partition P1-P9 in order")

    validation["fine_periods"] = validation["fine_periods"].map(lambda values: " + ".join(values))
    if (sensors["start"] < paper_start).any() or (sensors["end"] > paper_end).any():
        raise ValueError("Sensor availability falls outside the paper period")

    return registry, fine, validation, sensors


if __name__ == "__main__":
    metadata, fine_periods, validation_periods, sensor_rows = load_period_frames()
    print(f"Validated {len(fine_periods)} fine periods and {len(validation_periods)} validation periods")
    print(fine_periods[["period_id", "n_months", "reliable_for_slope", "changepoint_detection_exempt"]].to_string(index=False))
