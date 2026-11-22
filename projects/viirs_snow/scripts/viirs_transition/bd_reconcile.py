#!/usr/bin/env python3
"""Reconcile analysis B (incr_summary.py) and D (wb_analyze.py) snow-increment totals.

Both sum WESNN1-3_INCR (kg m-2) over the 8 daily steps of cat/ens0000/*catch_progn_incr*. Differences are only in
(1) tile set, (2) weighting, (3) which days, (4) the 0.01 kg m-2 threshold, (5) the 365/n_days scaling.
This script recomputes the Terra+Aqua / Terra-only DA-arm totals step by step from one definition to the other and
writes docs/viirs_transition_report_figures/BD_reconcile_table.csv (+ prints a markdown table).
Run with GEOSpyD python3 -I (no SLURM needed, ~1 min).
"""
import csv
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

PROJ = Path(__file__).resolve().parents[2]
REPO = PROJ.parents[1]
FIG = PROJ / "docs" / "viirs_transition_report_figures"
CATCHMENT_DEF = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/land/EASEv2_M36/clsm/catchment.def"
ARMS = ["T_to_J2", "TA_to_J1J2", "T_control", "TA_control"]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import incr_summary as B  # noqa: E402  (excluded days, NH groups; main() is not run on import)


def rd(path, name):
    with Dataset(path) as d:
        d.set_auto_mask(False)
        return np.array(d[name][:], float)


def main():
    cfg = json.load(open(PROJ / "config" / "arms.json"))
    sy_ids = [s["id"] for s in cfg["snow_years"]]
    area = np.loadtxt(CATCHMENT_DEF, usecols=(6,), skiprows=1)
    ol = Dataset(PROJ / "output/water_balance/wb_OL.nc")
    lat = np.array(ol["lat"][:], float)
    swe_max = np.nanmax(np.ma.filled(ol["SWE_max"][:], np.nan), 0)
    scf = np.ma.filled(ol["m_SCF"][:], np.nan)
    scf_max = np.nanmax(scf.reshape(-1, lat.size), 0)
    jja = np.nanmean(scf[:, [10, 11, 0], :].reshape(-1, lat.size), 0)
    mask = (lat > 20) & ((scf_max > 0.05) | (swe_max > 5)) & (jja < 0.2) & np.isfinite(lat)
    n30 = lat >= 30
    print("tiles lat>=30:", int(n30.sum()), " mask:", int(mask.sum()), " mask & lat>=30:", int((mask & n30).sum()))
    rows = []
    for arm in ARMS:
        st = Dataset(PROJ / f"output/incr_stats/incr_stats_{arm}.nc"); st.set_auto_mask(False)
        wb = Dataset(PROJ / f"output/water_balance/wb_{arm}.nc"); wb.set_auto_mask(False)
        wsy = [str(x) for x in wb["sy_id"][:]]
        dd = B.dates(st["date"][:])
        msk = B.sy_masks(dd)
        step_net = np.array(st["step_net"][:], float)
        tile_net = np.array(st["tile_net"][:], float)
        slat = np.array(st["lat"][:], float)
        assert np.allclose(slat, lat, atol=1e-4)
        for sid in wsy:
            mk = msk[sid]
            if not mk.any():
                continue
            i = sy_ids.index(sid)
            j = wsy.index(sid)
            Id = np.nan_to_num(np.array(wb["I_snow"][j], float))
            nd = mk.sum()
            ndays_full = int(wb["n_days"][j])
            q = {}
            q["1_B_headline (sum step_net, N30, excl days, x365/nd)"] = B.grp_nh(step_net[mk]).sum() / n30.sum() * 365.0 / nd
            q["2_B no 365/nd scaling"] = B.grp_nh(step_net[mk]).sum() / n30.sum()
            tn = tile_net[i]
            q["3_B tile_net (all days incl. gap days), N30 mean"] = tn[n30].mean()
            q["4_D I_snow (no 0.01 threshold), N30 mean"] = Id[n30].mean()
            q["5_D I_snow, N30 area-weighted"] = np.sum(Id[n30] * area[n30]) / area[n30].sum()
            q["6_D I_snow, mask tiles, unweighted"] = Id[mask].mean()
            q["7_D I_snow, mask tiles, area-weighted (= D table, C0 vs OL)"] = np.sum(Id[mask] * area[mask]) / area[mask].sum()
            # extras
            q["x_frac N30 tiles with any |incr|>=0.01"] = float((tn[n30] != 0).mean())
            q["x_frac mask tiles within N30"] = float(n30[mask].mean())
            q["x_mean over N30 tiles with increments only"] = tn[n30][tn[n30] != 0].mean()
            q["x_n_days B / D"] = f"{nd}/{ndays_full}"
            for k, v in q.items():
                rows.append((arm, sid, k, v))
    out = FIG / "BD_reconcile_table.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["arm", "sy", "quantity", "value"]); w.writerows(rows)
    # print wide
    keys = sorted({r[2] for r in rows})
    print("arm sy | " + " | ".join(k[:2] for k in keys))
    for arm in ARMS:
        for sid in sy_ids:
            vals = {r[2]: r[3] for r in rows if r[0] == arm and r[1] == sid}
            if vals:
                print(arm, sid, " | ".join(f"{vals[k]:.3f}" if not isinstance(vals[k], str) else vals[k] for k in keys))
    print("wrote", out)


if __name__ == "__main__":
    main()
