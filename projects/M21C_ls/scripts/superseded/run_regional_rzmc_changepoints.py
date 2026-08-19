#!/usr/bin/env python3
"""Blind changepoint detection on regional RZMC DA-OL series, plus OL controls.

Includes the methodological sanity check that PELT on the standardized
seasonally adjusted series and on the same series in native units accepts
identical break dates.
"""
from __future__ import annotations

import json, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import changepoint_detection as cpd  # noqa: E402
from m21c_periods import load_period_frames  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "output" / "regional_rzmc_transitions"


def accepted_dates(result) -> list[pd.Timestamp]:
    det = result.detections
    if det.empty or "accepted_detection" not in det.columns:
        return []
    acc = det.loc[det["accepted_detection"].astype(bool), "break_date"]
    return sorted(pd.Timestamp(v) for v in acc.tolist())


def main() -> int:
    cfg = cpd.load_changepoint_config()
    ds = xr.open_dataset(OUT / "regional_rzmc_monthly.nc")
    time = pd.DatetimeIndex(ds.time.values)
    regions = [str(r) for r in ds.region.values]
    labels = json.loads((ROOT/"config"/"regional_rzmc_regions.json").read_text())["regions"]
    label_of = {r["region_id"]: r["label"] for r in labels}

    # ---- sanity check: standardized vs native ----
    print("=== sanity check: standardized vs native-unit detection ===")
    original = cpd._standardize
    def identity(values):
        scale = float(np.std(values, ddof=1))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("no variation")
        return values.copy(), 0.0, 1.0
    check = []
    for rid in regions:
        v = ds["delta"].sel(region=rid).values.astype("float64")
        std_dates = accepted_dates(cpd.detect_changepoints(v, time, config=cfg))
        cpd._standardize = identity
        try:
            nat_dates = accepted_dates(cpd.detect_changepoints(v, time, config=cfg))
        finally:
            cpd._standardize = original
        same = std_dates == nat_dates
        check.append({"region_id": rid, "standardized": [str(d.date()) for d in std_dates],
                      "native": [str(d.date()) for d in nat_dates], "identical": same})
        print(f"  {rid:9s} identical={same}  std={[str(d.date()) for d in std_dates]}")
    pd.DataFrame(check).to_csv(OUT/"standardization_sanity_check.csv", index=False)
    all_same = all(c["identical"] for c in check)
    print(f"  -> all regions identical: {all_same}\n")

    # ---- production detection on delta and OL control ----
    _, fine, _, _ = load_period_frames()
    bounds = [(row.period_id, pd.Timestamp(row.start)) for row in fine.itertuples()][1:]  # P2..P9
    rows = []
    for rid in regions:
        for kind in ("delta", "ol"):
            v = ds[kind].sel(region=rid).values.astype("float64")
            dates = accepted_dates(cpd.detect_changepoints(v, time, config=cfg))
            if not dates:
                rows.append({"region_id": rid, "label": label_of[rid], "series": kind,
                             "detected": None, "nearest_boundary": None, "offset_months": None,
                             "match_3mo": None, "match_6mo": None})
                continue
            for d in dates:
                offs = [(pid, (d.year-b.year)*12 + (d.month-b.month)) for pid, b in bounds]
                pid, off = min(offs, key=lambda t: abs(t[1]))
                rows.append({"region_id": rid, "label": label_of[rid], "series": kind,
                             "detected": str(d.date()), "nearest_boundary": pid,
                             "offset_months": off, "match_3mo": abs(off) <= 3,
                             "match_6mo": abs(off) <= 6})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT/"regional_breakpoint_table.csv", index=False)
    print("=== regional breakpoint table ===")
    print(tab.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
