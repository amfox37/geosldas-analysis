#!/usr/bin/env python
"""Validate cached global ASCAT superobs against GEOSldas ObsFcstAna.

This is the global companion to check_raw_superobs_vs_ofa.py. It does not
reread raw BUFR/H121 files. Instead it uses the derived global super-ob cache
created by build_global_superob_cache.py and compares those tile/cycle means
with ASCAT observations stored in GEOSldas ObsFcstAna files.

Raw-only cache rows are expected because OFA has additional model-based QC.
The primary validation metric is therefore agreement for OFA tile/cycles that
also exist in the cached raw-derived superobs.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


PRODUCTS = {
    "legacy": {"Metop-A": 9, "Metop-B": 10, "Metop-C": 11},
    "h121": {"Metop-A": 14, "Metop-B": 15, "Metop-C": 16},
}


def parse_date(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%d")


def iter_dates(start: datetime, end: datetime):
    date = start
    while date <= end:
        yield date
        date += timedelta(days=1)


def filename_window(path: Path) -> int:
    """Return GEOSldas 3-hour cycle label from an OFA filename timestamp."""
    token = [p for p in path.name.replace(".", "_").split("_") if p.endswith("z")][0]
    hour = int(token[:2])
    return hour // 3


def cache_file(cache_dir: Path, date: datetime, version: str) -> Path:
    return cache_dir / f"ascat_global_superobs_{date:%Y%m%d}_{version}.pkl"


def read_cache(cache_dir: Path, date: datetime, version: str) -> pd.DataFrame:
    path = cache_file(cache_dir, date, version)
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_pickle(path)


def read_ofa_day(ofa_dir: Path, date: datetime) -> pd.DataFrame:
    rows = []
    species_keep = {sp for mapping in PRODUCTS.values() for sp in mapping.values()}
    for fpath in sorted(ofa_dir.glob(f"*{date:%Y%m%d}*.nc4")):
        cycle = filename_window(fpath)
        with nc.Dataset(fpath) as ds:
            species = ds.variables["species"][:]
            mask = np.isin(species, list(species_keep))
            if not np.any(mask):
                continue
            rows.append(
                pd.DataFrame(
                    {
                        "date": date.strftime("%Y-%m-%d"),
                        "species": species[mask].astype("int16"),
                        "tilenum": ds.variables["tilenum"][:][mask].astype("int64"),
                        "cycle": np.int16(cycle),
                        "ofa_obs_pct": (ds.variables["obs"][:][mask] * 100.0).astype("float32"),
                        "ofa_assim_flag": ds.variables["assim_flag"][:][mask].astype("int16"),
                    }
                )
            )
    if not rows:
        return pd.DataFrame(columns=["date", "species", "tilenum", "cycle", "ofa_obs_pct", "ofa_assim_flag"])

    # Duplicates should be rare, but collapse defensively to the same key as
    # the cached raw-derived superobs.
    ofa = pd.concat(rows, ignore_index=True)
    return (
        ofa.groupby(["date", "species", "tilenum", "cycle"], as_index=False)
        .agg(ofa_obs_pct=("ofa_obs_pct", "mean"), ofa_count=("ofa_obs_pct", "size"))
    )


def annotate_species(cache: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for product, platforms in PRODUCTS.items():
        for platform, species in platforms.items():
            sel = cache[(cache["product"] == product) & (cache["platform"] == platform)].copy()
            sel["species"] = np.int16(species)
            frames.append(sel)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def normalize_cache_to_analysis_cycles(cache: pd.DataFrame) -> pd.DataFrame:
    """Collapse source-day superobs onto GEOS analysis-date/cycle keys."""
    cache = annotate_species(cache)
    source_date = pd.to_datetime(cache["date"])
    cycle = cache["cycle"].astype("int16")
    cache["date"] = (source_date + pd.to_timedelta((cycle // 8).astype(int), unit="D")).dt.strftime("%Y-%m-%d")
    cache["cycle"] = (cycle % 8).astype("int16")
    cache["window"] = cache["cycle"].astype("int8")

    weighted = cache.copy()
    n_obs = weighted["n_obs"].astype(float)
    weighted["_ssm_sum"] = weighted["ssm_pct"] * n_obs
    weighted["_lat_sum"] = weighted["lat"] * n_obs
    weighted["_lon_sum"] = weighted["lon"] * n_obs
    weighted["_ssm2_sum"] = (weighted["ssm_std_pct"] ** 2 + weighted["ssm_pct"] ** 2) * n_obs

    grouped = (
        weighted.groupby(["date", "product", "platform", "species", "tilenum", "cycle"], as_index=False)
        .agg(
            window=("window", "first"),
            i_indg=("i_indg", "first"),
            j_indg=("j_indg", "first"),
            n_obs=("n_obs", "sum"),
            _ssm_sum=("_ssm_sum", "sum"),
            _lat_sum=("_lat_sum", "sum"),
            _lon_sum=("_lon_sum", "sum"),
            _ssm2_sum=("_ssm2_sum", "sum"),
            ssm_min_pct=("ssm_min_pct", "min"),
            ssm_max_pct=("ssm_max_pct", "max"),
        )
    )
    n = grouped["n_obs"].astype(float)
    grouped["ssm_pct"] = grouped["_ssm_sum"] / n
    grouped["lat"] = grouped["_lat_sum"] / n
    grouped["lon"] = grouped["_lon_sum"] / n
    var = grouped["_ssm2_sum"] / n - grouped["ssm_pct"] ** 2
    grouped["ssm_std_pct"] = np.sqrt(np.maximum(var, 0.0))
    return grouped.drop(columns=["_ssm_sum", "_lat_sum", "_lon_sum", "_ssm2_sum"])


def compare_tables(cache: pd.DataFrame, ofa: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = cache.merge(
        ofa,
        on=["date", "species", "tilenum", "cycle"],
        how="outer",
        indicator=True,
    )
    matched = merged[merged["_merge"] == "both"].copy()
    matched["diff_pct"] = matched["ofa_obs_pct"] - matched["ssm_pct"]

    rows = []
    for product, platforms in PRODUCTS.items():
        for platform, species in platforms.items():
            sel = merged[(merged["species"] == species)]
            pair = matched[(matched["species"] == species)]
            raw_total = int(((sel["_merge"] == "both") | (sel["_merge"] == "left_only")).sum())
            ofa_total = int(((sel["_merge"] == "both") | (sel["_merge"] == "right_only")).sum())
            diff = pair["diff_pct"].to_numpy(float)
            rows.append(
                {
                    "date": "all",
                    "product": product,
                    "platform": platform,
                    "species": species,
                    "raw_total": raw_total,
                    "ofa_total": ofa_total,
                    "matched": int(len(pair)),
                    "matched_frac_ofa": len(pair) / ofa_total if ofa_total else np.nan,
                    "raw_only": int((sel["_merge"] == "left_only").sum()),
                    "ofa_only": int((sel["_merge"] == "right_only").sum()),
                    "bias_pct": float(np.mean(diff)) if len(diff) else np.nan,
                    "rmse_pct": float(np.sqrt(np.mean(diff**2))) if len(diff) else np.nan,
                    "mae_pct": float(np.mean(np.abs(diff))) if len(diff) else np.nan,
                    "median_abs_pct": float(np.median(np.abs(diff))) if len(diff) else np.nan,
                    "max_abs_pct": float(np.max(np.abs(diff))) if len(diff) else np.nan,
                }
            )
    return matched, pd.DataFrame(rows)


def aggregate_summary(daily_summary: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for product, platforms in PRODUCTS.items():
        for platform, species in platforms.items():
            s = daily_summary[daily_summary["species"] == species]
            p = pairs[pairs["species"] == species]
            diff = p["diff_pct"].to_numpy(float)
            rows.append(
                {
                    "product": product,
                    "platform": platform,
                    "species": species,
                    "raw_total": int(s["raw_total"].sum()),
                    "ofa_total": int(s["ofa_total"].sum()),
                    "matched": int(len(p)),
                    "matched_frac_ofa": len(p) / s["ofa_total"].sum() if s["ofa_total"].sum() else np.nan,
                    "raw_only": int(s["raw_only"].sum()),
                    "ofa_only": int(s["ofa_only"].sum()),
                    "bias_pct": float(np.mean(diff)) if len(diff) else np.nan,
                    "rmse_pct": float(np.sqrt(np.mean(diff**2))) if len(diff) else np.nan,
                    "mae_pct": float(np.mean(np.abs(diff))) if len(diff) else np.nan,
                    "median_abs_pct": float(np.median(np.abs(diff))) if len(diff) else np.nan,
                    "max_abs_pct": float(np.max(np.abs(diff))) if len(diff) else np.nan,
                    "raw_range_pct": f"{p['ssm_pct'].min():.1f}-{p['ssm_pct'].max():.1f}" if len(p) else "",
                    "ofa_range_pct": f"{p['ofa_obs_pct'].min():.1f}-{p['ofa_obs_pct'].max():.1f}" if len(p) else "",
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=parse_date, default=parse_date("2020-06-01"))
    parser.add_argument("--end-date", type=parse_date, default=parse_date("2020-06-10"))
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("projects/ascat_da/.cache/global_superobs"),
        help="Directory containing daily global super-ob cache pickle files.",
    )
    parser.add_argument("--version", default="geos_cycle_global_v1")
    parser.add_argument(
        "--ofa-dir",
        type=Path,
        default=Path(
            "/Users/amfox/Desktop/ASCAT_SSM_CDR/"
            "hsaf_cdr_test_DAv8_M36_202006_innov/output/SMAP_EASEv2_M36_GLOBAL/"
            "ana/ens_avg/Y2020/M06"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("projects/ascat_da/.cache/global_superobs"),
        help="Directory for validation CSV outputs.",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cache_frames = []
    ofa_frames = []

    for date in iter_dates(args.start_date, args.end_date):
        print(f"{date:%Y-%m-%d}: reading cache and OFA", flush=True)
        cache_frames.append(read_cache(args.cache_dir, date, args.version))
        ofa_frames.append(read_ofa_day(args.ofa_dir, date))
        print(
            f"  cache rows={len(cache_frames[-1]):,} OFA rows={len(ofa_frames[-1]):,}",
            flush=True,
        )

    print("Normalizing source dates to GEOS analysis cycles and comparing", flush=True)
    cache = pd.concat(cache_frames, ignore_index=True)
    cache = normalize_cache_to_analysis_cycles(cache)
    date_min = args.start_date.strftime("%Y-%m-%d")
    date_max = args.end_date.strftime("%Y-%m-%d")
    cache = cache[(cache["date"] >= date_min) & (cache["date"] <= date_max)].copy()
    ofa = pd.concat(ofa_frames, ignore_index=True)
    pairs, summary = compare_tables(cache, ofa)
    aggregate = aggregate_summary(summary, pairs)

    tag = f"{args.start_date:%Y%m%d}_{args.end_date:%Y%m%d}_{args.version}"
    pair_path = args.out_dir / f"global_superobs_vs_ofa_pairs_{tag}.csv"
    daily_path = args.out_dir / f"global_superobs_vs_ofa_daily_{tag}.csv"
    agg_path = args.out_dir / f"global_superobs_vs_ofa_summary_{tag}.csv"

    keep_pairs = [
        "date",
        "product",
        "platform",
        "species",
        "tilenum",
        "cycle",
        "window",
        "ssm_pct",
        "ofa_obs_pct",
        "diff_pct",
        "n_obs",
        "ofa_count",
        "ssm_std_pct",
        "ssm_min_pct",
        "ssm_max_pct",
    ]
    pairs[keep_pairs].to_csv(pair_path, index=False)
    summary.to_csv(daily_path, index=False)
    aggregate.to_csv(agg_path, index=False)

    print("\nAggregate global validation:")
    print(aggregate.round(3).to_string(index=False))
    print(f"\nWrote matched pairs: {pair_path}")
    print(f"Wrote daily summary: {daily_path}")
    print(f"Wrote aggregate summary: {agg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
