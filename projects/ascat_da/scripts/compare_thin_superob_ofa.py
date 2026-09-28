#!/usr/bin/env python
"""Compare GEOSldas ObsFcstAna diagnostics across an OL and several DA runs.

Built for the CF0360 H SAF ASCAT thinning vs super-obbing tests (M21C_testing),
but works for any set of runs that write ens_avg ldas_ObsFcstAna files.

Two kinds of diagnostics:

1. Monitor-obs skill (independent obs). For every obs type, keep only the obs
   present in *all* runs (matched on cycle, species, tile, lat, lon), restricted
   to cycles that all runs have reached, so partially finished runs can be
   compared. Reports O-F mean/std per run and day, the change in O-F std vs the
   reference (OL) run, and a latitude-band breakdown.

2. Filter consistency for the species each run assimilates (each run's own obs):
   - normalized innovations (O-F)/sqrt(obsvar+fcstvar): std should be ~1
   - Desroziers ratios: mean((O-A)(O-F))/mean(obsvar) should be ~1 if R is right,
     mean((A-F)(O-F))/mean(fcstvar) should be ~1 if HPH^T is right
   - O-A std vs O-F std (analysis should move toward the obs)

Outputs CSV tables in --out-dir and prints them.
"""

from __future__ import annotations

import argparse
import glob
import os
import time
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


M21C = Path("/gpfsm/dnb06/projects/p284/M21C_testing")
DOMAIN = "CF0360x6C_GLOBAL"

DEFAULT_RUNS = {
    "OL": "M21C_test_CF0360_OL_monitor_HSAF",
    "baseline": "M21C_test_CF0360_baseline_HSAF_xc03125",
    "superob025": "M21C_test_CF0360_superob025_HSAF_xc03125",
    "thin025": "M21C_test_CF0360_thin025_HSAF_xc03125",
}

# obs_param descr prefix -> reporting group
GROUPS = [
    ("SMOS_", "SMOS Tb"),
    ("SMAP_L1C_Tb", "SMAP Tb"),
    ("CYGNSS_SM", "CYGNSS SM"),
    ("MYD10C1", "MODIS SCF"),
    ("MOD10C1", "MODIS SCF"),
    ("ASCAT_HSAF", "HSAF SM"),
]

LAT_BANDS = [-90, -30, 0, 30, 45, 55, 90]


def species_groups(run_dir: Path) -> dict[int, str]:
    """Map ObsFcstAna species index (1-based order of species with innov/assim) to group.

    The ObsFcstAna 'species' index counts the obs_param entries that were read,
    in obs_param order; take them from the run's special nml (descr with
    assim or getinnov true).
    """
    nml = run_dir / "run" / "LDASsa_SPECIAL_inputs_ensupd.nml"
    vals: dict[int, dict[str, str]] = {}
    for line in nml.read_text().splitlines():
        line = line.split("!")[0].strip()
        if not line.startswith("obs_param_nml(") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        idx = int(key[key.index("(") + 1:key.index(")")])
        field = key.split("%")[1].strip()
        vals.setdefault(idx, {})[field] = val.strip().strip("'").strip()
    used = [i for i in sorted(vals)
            if vals[i].get("assim") == ".true." or vals[i].get("getinnov") == ".true."]
    out = {}
    for k, i in enumerate(used, 1):
        descr = vals[i].get("descr", "")
        out[k] = next((g for p, g in GROUPS if descr.startswith(p)), descr)
    return out


def list_ofa_files(run_dir: Path) -> list[str]:
    """ObsFcstAna files moved to output/, plus those of a running segment in scratch/.

    Scratch files modified in the last 2 minutes may still be open and are skipped.
    """
    done = glob.glob(str(run_dir / "output" / DOMAIN / "ana" / "ens_avg" / "Y*" / "M*"
                         / "*ldas_ObsFcstAna.*.nc4"))
    names = {os.path.basename(f) for f in done}
    now = time.time()
    running = [f for f in glob.glob(str(run_dir / "scratch" / "*ldas_ObsFcstAna.*.nc4"))
               if os.path.basename(f) not in names and now - os.path.getmtime(f) > 120]
    return sorted(done + running)


def read_run(name: str, run_dir: Path) -> pd.DataFrame:
    files = list_ofa_files(run_dir)
    groups = species_groups(run_dir)
    frames = []
    for f in files:
        cycle = os.path.basename(f).split(".")[-2]
        with nc.Dataset(f) as d:
            if d.dimensions["n_obs"].size == 0:
                continue
            g = lambda v: np.ma.filled(d[v][:].astype(np.float64), np.nan)
            frames.append(pd.DataFrame({
                "cycle": cycle,
                "species": d["species"][:].astype(np.int16),
                "tile": d["tilenum"][:].astype(np.int32),
                "lat": np.round(g("lat"), 4),
                "lon": np.round(g("lon"), 4),
                "obs": g("obs"), "obsvar": g("obsvar"),
                "fcst": g("fcst"), "fcstvar": g("fcstvar"),
                "ana": g("ana"), "assim": d["assim_flag"][:].astype(np.int8),
            }))
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["group"] = df["species"].map(groups)
    df["run"] = name
    df["day"] = df["cycle"].str[:8]
    return df


def monitor_skill(data: dict[str, pd.DataFrame], ref: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs = list(data)
    cycles = set.intersection(*(set(df["cycle"]) for df in data.values()))
    keys = ["cycle", "species", "tile", "lat", "lon"]
    common = None
    for r, df in data.items():
        # monitor obs only: H SAF is assimilated in the DA runs, and its obs sets
        # differ between them (raw, super-obbed, thinned)
        x = df[df["cycle"].isin(cycles) & (df["group"] != "HSAF SM")]
        x = x[keys + ["group", "day", "obs", "fcst"]].copy()
        x[f"omf_{r}"] = x["obs"] - x["fcst"]
        x = x.drop(columns=["obs", "fcst"]).drop_duplicates(keys)
        common = x if common is None else common.merge(x.drop(columns=["group", "day"]), on=keys)
    common = common.dropna()

    rows = []
    for (grp, day), x in common.groupby(["group", "day"]):
        row = {"group": grp, "day": day, "N": len(x)}
        for r in runs:
            row[f"{r}_mean"] = x[f"omf_{r}"].mean()
            row[f"{r}_std"] = x[f"omf_{r}"].std()
        for r in runs:
            if r != ref:
                row[f"dstd_{r}_vs_{ref}"] = row[f"{r}_std"] - row[f"{ref}_std"]
        rows.append(row)
    by_day = pd.DataFrame(rows)

    common["band"] = pd.cut(common["lat"], LAT_BANDS)
    rows = []
    for (grp, band), x in common.groupby(["group", "band"], observed=True):
        row = {"group": grp, "lat_band": str(band), "N": len(x)}
        for r in runs:
            row[f"{r}_std"] = x[f"omf_{r}"].std()
        rows.append(row)
    by_band = pd.DataFrame(rows)
    return by_day, by_band


def filter_consistency(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for r, df in data.items():
        for (grp, assim), x in df.groupby(["group", "assim"]):
            if grp != "HSAF SM":
                continue
            omf = x["obs"] - x["fcst"]
            oma = x["obs"] - x["ana"]
            amf = x["ana"] - x["fcst"]
            tot = np.sqrt(x["obsvar"] + x["fcstvar"])
            rows.append({
                "run": r, "group": grp, "assimilated": bool(assim), "N": len(x),
                "omf_std": omf.std(), "oma_std": oma.std(),
                "norm_innov_std": (omf / tot).std(),
                # Desroziers ratios need an analysis; for monitor-only obs A == F
                "desroziers_R_ratio": (oma * omf).mean() / x["obsvar"].mean() if assim else np.nan,
                "desroziers_HBH_ratio": (amf * omf).mean() / x["fcstvar"].mean() if assim else np.nan,
                "sqrt_mean_obsvar": np.sqrt(x["obsvar"].mean()),
                "sqrt_mean_fcstvar": np.sqrt(x["fcstvar"].mean()),
            })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=M21C)
    parser.add_argument("--runs", nargs="+", default=[f"{k}={v}" for k, v in DEFAULT_RUNS.items()],
                        help="name=experiment_dir pairs (the first is the reference)")
    parser.add_argument("--out-dir", type=Path, default=M21C / "analysis_thin_superob")
    args = parser.parse_args()

    runs = dict(s.split("=", 1) for s in args.runs)
    data = {}
    for name, exp in runs.items():
        df = read_run(name, args.root / exp)
        print(f"{name:11s} {exp}: {df['cycle'].nunique() if len(df) else 0} cycles, "
              f"{len(df):,} obs", flush=True)
        if len(df):
            data[name] = df
    ref = next(iter(data))
    print(f"Reference run for O-F std differences: {ref}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.float_format", lambda v: f"{v:.4g}")

    by_day, by_band = monitor_skill(data, ref)
    cons = filter_consistency(data)
    for label, tab in [("monitor_omf_by_day", by_day), ("monitor_omf_by_latband", by_band),
                       ("hsaf_filter_consistency", cons)]:
        tab.to_csv(args.out_dir / f"{label}.csv", index=False)
        print(f"\n=== {label} ===")
        print(tab.to_string(index=False))
    print(f"\nWrote CSVs to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
