#!/usr/bin/env python3
"""Build area-weighted regional RZMC DA-OL monthly series.

Reuses the audited MonthlySeriesLoader without modifying it. Regional masks are static: a lat/lon box intersected with
valid_land and with the set of tiles finite in OL and DA in every month, so the
contributing sample cannot change through the record.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from trend_breakpoint_series import MonthlySeriesLoader  # noqa: E402
from m21c_periods import load_period_frames  # noqa: E402

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


def main(var: str = "RZMC") -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    registry = json.loads(REGIONS.read_text())
    _, fine, _, _ = load_period_frames()

    with MonthlySeriesLoader() as loader:
        lat = np.asarray(loader.lat.values, dtype="float64")
        lon = np.asarray(loader.lon.values, dtype="float64")
        fields = loader.load_tile_series(var, mask="valid_land")
        time = pd.DatetimeIndex(fields.time.values)
        area = np.asarray(fields["tile_area"].fillna(0.0).values, dtype="float64")
        ol = fields["ol"].transpose("tile", "time").values.astype("float64")
        da = fields["da"].transpose("tile", "time").values.astype("float64")

    n_tile, n_month = ol.shape
    print(f"loaded {var}: {n_tile} tiles x {n_month} months", flush=True)

    # tiles finite in BOTH runs in EVERY month -> support cannot change in time
    always_finite = np.all(np.isfinite(ol) & np.isfinite(da), axis=1)
    print(f"tiles finite in all months: {int(always_finite.sum())}", flush=True)

    delta = da - ol
    rows, series = [], {}
    for reg in registry["regions"]:
        rid, label = reg["region_id"], reg["label"]
        sel = region_tile_mask(lat, lon, reg["bounds"]) & always_finite & (area > 0)
        n = int(sel.sum())
        if n == 0:
            print(f"  !! {rid}: no tiles"); continue
        w = area[sel]
        wsum = w.sum()
        d_series = (delta[sel, :] * w[:, None]).sum(axis=0) / wsum
        o_series = (ol[sel, :] * w[:, None]).sum(axis=0) / wsum
        # static-support audit: recompute per-month contributing count
        per_month = np.isfinite(delta[sel, :]).sum(axis=0)
        constant = bool(per_month.min() == per_month.max() == n)
        rows.append({
            "region_id": rid, "label": label, "n_tiles": n,
            "area_km2": float(wsum),
            "support_constant": constant,
            "n_tiles_min": int(per_month.min()), "n_tiles_max": int(per_month.max()),
            "delta_mean": float(d_series.mean()), "delta_sd": float(d_series.std(ddof=1)),
        })
        series[rid] = {"delta": d_series, "ol": o_series}
        print(f"  {rid:9s} n={n:6d}  support_constant={constant}  "
              f"mean d{var}={d_series.mean():+.5f}", flush=True)

    support = pd.DataFrame(rows)
    support.to_csv(OUT / f"regional_{var.lower()}_support_audit.csv", index=False)
    assert support["support_constant"].all(), "non-static support detected"

    ds = xr.Dataset(
        {
            "delta": (("region", "time"), np.array([series[r]["delta"] for r in series])),
            "ol":    (("region", "time"), np.array([series[r]["ol"] for r in series])),
        },
        coords={"region": list(series), "time": time},
    )
    ds["delta"].attrs = {"units": "m3 m-3", "long_name": f"area-weighted {var} DA-OL"}
    ds["ol"].attrs = {"units": "m3 m-3", "long_name": f"area-weighted OL {var}"}
    ds.to_netcdf(OUT / f"regional_{var.lower()}_monthly.nc")
    print(f"\nwrote {OUT}/regional_{var.lower()}_monthly.nc")
    return 0


if __name__ == "__main__":
    _ap = argparse.ArgumentParser(); _ap.add_argument("--variable", default="RZMC")
    raise SystemExit(main(_ap.parse_args().variable))
