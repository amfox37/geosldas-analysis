#!/usr/bin/env python3
"""Build raw and tile-seasonally-adjusted regional DA-OL monthly series.

The paired DA-OL series at every tile receives the same trend-preserving
calendar-month adjustment used by the production tile-trend workflow. The raw
and adjusted tile series are then area-aggregated over static regional support.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from trend_breakpoint_series import MonthlySeriesLoader  # noqa: E402
from trend_statistics import (  # noqa: E402
    DEFAULT_CONFIG,
    calendar_month_anomalies,
    load_trend_config,
    monthly_time_axis,
)

ROOT = HERE.parent
REGIONS = ROOT / "config" / "regional_rzmc_regions.json"
OUT = ROOT / "output" / "regional_rzmc_transitions"


def region_tile_mask(lat, lon, bounds):
    if bounds is None:
        return np.ones(lat.size, dtype=bool)
    return (
        (lat >= bounds["lat_min"]) & (lat <= bounds["lat_max"])
        & (lon >= bounds["lon_min"]) & (lon <= bounds["lon_max"])
    )


def _adjust_tile_chunk(
    start: int,
    values: np.ndarray,
    time_years: np.ndarray,
    months: np.ndarray,
    minimum_samples: int,
) -> tuple[int, np.ndarray]:
    """Apply the production tile-level seasonal adjustment to one tile chunk."""

    adjusted = np.empty_like(values, dtype="float64")
    valid = np.ones(values.shape[1], dtype=bool)
    for index, tile_values in enumerate(values):
        adjusted[index], month_mask = calendar_month_anomalies(
            tile_values,
            time_years,
            months,
            valid,
            minimum_samples,
        )
        if month_mask != (1 << 12) - 1:
            raise RuntimeError(
                f"Tile {start + index} did not retain all calendar months: {month_mask}"
            )
    return start, adjusted


def aggregate_adjusted_tile_series(
    values: np.ndarray,
    region_weights: np.ndarray,
    time: pd.DatetimeIndex,
    *,
    minimum_samples: int,
    n_jobs: int = 1,
    chunk_size: int = 256,
    progress_every: int = 25,
) -> tuple[np.ndarray, np.ndarray]:
    """Return regional raw and tile-adjusted area-weighted monthly series."""

    values = np.asarray(values, dtype="float64")
    region_weights = np.asarray(region_weights, dtype="float64")
    if values.ndim != 2 or region_weights.ndim != 2:
        raise ValueError("values and region_weights must both be two-dimensional")
    if region_weights.shape[1] != values.shape[0]:
        raise ValueError("region_weights must have one column per input tile")
    if not np.all(np.isfinite(values)):
        raise ValueError("tile adjustment requires complete static-support series")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")

    time_years, months = monthly_time_axis(time)
    raw = region_weights @ values
    adjusted_regional = np.zeros_like(raw)
    tasks = (
        delayed(_adjust_tile_chunk)(
            start,
            values[start : start + chunk_size],
            time_years,
            months,
            minimum_samples,
        )
        for start in range(0, values.shape[0], chunk_size)
    )
    results = Parallel(
        n_jobs=n_jobs,
        prefer="threads" if n_jobs != 1 else None,
        return_as="generator",
        batch_size=1,
    )(tasks)
    n_chunks = (values.shape[0] + chunk_size - 1) // chunk_size
    for chunk_number, (start, adjusted_chunk) in enumerate(results, start=1):
        stop = start + adjusted_chunk.shape[0]
        adjusted_regional += region_weights[:, start:stop] @ adjusted_chunk
        if progress_every > 0 and (
            chunk_number % progress_every == 0 or chunk_number == n_chunks
        ):
            print(
                f"  seasonally adjusted {stop:,}/{values.shape[0]:,} tiles",
                flush=True,
            )
    return raw, adjusted_regional


def main(
    var: str = "RZMC",
    *,
    n_jobs: int = 1,
    chunk_size: int = 256,
    trend_config: Path = DEFAULT_CONFIG,
) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    registry = json.loads(REGIONS.read_text())
    settings = load_trend_config(trend_config)

    with MonthlySeriesLoader() as loader:
        lat = np.asarray(loader.lat.values, dtype="float64")
        lon = np.asarray(loader.lon.values, dtype="float64")
        fields = loader.load_tile_series(var, mask="valid_land")
        time = pd.DatetimeIndex(fields.time.values)
        area = np.asarray(fields["tile_area"].fillna(0.0).values, dtype="float64")
        ol = fields["ol"].transpose("tile", "time").values.astype("float64")
        delta = fields["delta"].transpose("tile", "time").values.astype("float64")

    n_tile, n_month = delta.shape
    print(f"loaded {var}: {n_tile} tiles x {n_month} months", flush=True)

    # tiles finite in BOTH runs in EVERY month -> support cannot change in time
    always_finite = np.all(np.isfinite(delta), axis=1)
    print(f"tiles finite in all months: {int(always_finite.sum())}", flush=True)

    rows = []
    region_ids = []
    region_weights = []
    for reg in registry["regions"]:
        rid, label = reg["region_id"], reg["label"]
        sel = region_tile_mask(lat, lon, reg["bounds"]) & always_finite & (area > 0)
        n = int(sel.sum())
        if n == 0:
            print(f"  !! {rid}: no tiles"); continue
        w = area[sel]
        wsum = w.sum()
        normalized_weights = np.zeros(n_tile, dtype="float64")
        normalized_weights[sel] = w / wsum
        region_ids.append(rid)
        region_weights.append(normalized_weights)
        # static-support audit: recompute per-month contributing count
        per_month = np.isfinite(delta[sel, :]).sum(axis=0)
        constant = bool(per_month.min() == per_month.max() == n)
        rows.append({
            "region_id": rid, "label": label, "n_tiles": n,
            "area_km2": float(wsum),
            "support_constant": constant,
            "n_tiles_min": int(per_month.min()), "n_tiles_max": int(per_month.max()),
        })
        print(f"  {rid:9s} n={n:6d}  support_constant={constant}", flush=True)

    weights = np.asarray(region_weights)
    raw_delta, adjusted_delta = aggregate_adjusted_tile_series(
        delta,
        weights,
        time,
        minimum_samples=int(settings["minimum_climatology_samples_per_month"]),
        n_jobs=n_jobs,
        chunk_size=chunk_size,
    )
    raw_ol = weights @ ol
    for index, row in enumerate(rows):
        row.update(
            delta_raw_mean=float(raw_delta[index].mean()),
            delta_raw_sd=float(raw_delta[index].std(ddof=1)),
            delta_adjusted_mean=float(adjusted_delta[index].mean()),
            delta_adjusted_sd=float(adjusted_delta[index].std(ddof=1)),
        )
        print(
            f"  {row['region_id']:9s} mean raw d{var}={raw_delta[index].mean():+.5f}  "
            f"mean adjusted={adjusted_delta[index].mean():+.5f}",
            flush=True,
        )

    support = pd.DataFrame(rows)
    support.to_csv(OUT / f"regional_{var.lower()}_support_audit.csv", index=False)
    assert support["support_constant"].all(), "non-static support detected"

    ds = xr.Dataset(
        {
            "delta_raw": (("region", "time"), raw_delta),
            "delta_adjusted": (("region", "time"), adjusted_delta),
            # Backward-compatible raw alias for superseded exploratory scripts.
            "delta": (("region", "time"), raw_delta),
            "ol": (("region", "time"), raw_ol),
        },
        coords={"region": region_ids, "time": time},
    )
    ds["delta_raw"].attrs = {
        "units": "m3 m-3",
        "long_name": f"raw area-weighted {var} DA-OL",
    }
    ds["delta_adjusted"].attrs = {
        "units": "m3 m-3",
        "long_name": f"area-weighted tile-seasonally-adjusted {var} DA-OL",
        "seasonal_adjustment": settings["seasonal_adjustment"],
        "adjustment_order": "paired tile DA-OL adjustment, then regional area weighting",
    }
    ds["delta"].attrs = {
        "units": "m3 m-3",
        "long_name": f"deprecated alias of raw area-weighted {var} DA-OL",
    }
    ds["ol"].attrs = {"units": "m3 m-3", "long_name": f"area-weighted OL {var}"}
    ds.attrs.update(
        seasonal_adjustment=settings["seasonal_adjustment"],
        adjustment_source="trend_statistics.calendar_month_anomalies",
        adjustment_order="paired tile DA-OL adjustment, then regional area weighting",
        trend_config=str(trend_config),
    )
    ds.to_netcdf(OUT / f"regional_{var.lower()}_monthly.nc")
    print(f"\nwrote {OUT}/regional_{var.lower()}_monthly.nc")
    return 0


if __name__ == "__main__":
    _ap = argparse.ArgumentParser()
    _ap.add_argument("--variable", default="RZMC")
    _ap.add_argument("--n-jobs", type=int, default=1)
    _ap.add_argument("--chunk-size", type=int, default=256)
    _ap.add_argument("--trend-config", type=Path, default=DEFAULT_CONFIG)
    _args = _ap.parse_args()
    raise SystemExit(
        main(
            _args.variable,
            n_jobs=_args.n_jobs,
            chunk_size=_args.chunk_size,
            trend_config=_args.trend_config,
        )
    )
