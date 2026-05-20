#!/usr/bin/env python3
"""
Thoroughly compare every matched .bin / .nc4 ObsFcstAna pair in a directory.

Default target: /Users/amfox/Desktop/ens_avg/Y2020/M05
Override with a positional argument or --dir.

Outputs:
  - Per-pair console summary (PASS / issues list)
  - CSV report written next to this script (or --out)
  - Final totals

Usage:
  python diff_obsfcstana_dir.py
  python diff_obsfcstana_dir.py /some/other/dir --out /tmp/report.csv
  python diff_obsfcstana_dir.py --atol 1e-5
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------------------
# Repo path wiring
# ---------------------------------------------------------------------------
THIS_DIR = Path(__file__).resolve().parent
for _candidate in [
    THIS_DIR / "../../../common/python/io",
    THIS_DIR / "../../shared/python",
]:
    _p = _candidate.resolve()
    if _p.exists():
        sys.path.insert(0, str(_p))
        break

from read_GEOSldas import read_ObsFcstAna  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_DIR = Path("/Users/amfox/Desktop/ens_avg/Y2020/M05")
MAPL_UNDEF  = 9.999999870e+14   # MAPL fill value used by nc4 writer
BIN_MISSING = -9999.0           # legacy sentinel used by read_ObsFcstAna

# Fields present in both bin and nc4
INT_FIELDS = [
    ("obs_assim",   "assim_flag"),
    ("obs_species", "species"),
    ("obs_tilenum", "tilenum"),
]
FLOAT_FIELDS = [
    ("obs_lon",     "lon"),
    ("obs_lat",     "lat"),
    ("obs_obs",     "obs"),
    ("obs_obsvar",  "obsvar"),
    ("obs_fcst",    "fcst"),
    ("obs_fcstvar", "fcstvar"),
    ("obs_ana",     "ana"),
    ("obs_anavar",  "anavar"),
]
FLOAT_BIN_KEYS = {bk for bk, _ in FLOAT_FIELDS}

# Time key mapping  bin_key -> nc4_attr
TIME_MAP = {
    "year": "year", "month": "month", "day": "day",
    "hour": "hour", "min": "minute", "sec": "second",
    "dofyr": "dofyr", "pentad": "pentad",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _to_f64(arr):
    """Return plain float64 ndarray (handles masked arrays)."""
    if np.ma.isMaskedArray(arr):
        return arr.filled(np.nan).astype(np.float64)
    return np.asarray(arr, dtype=np.float64)


def _normalise(arr, sentinel):
    """Replace sentinel fill value with NaN for uniform comparison."""
    out = _to_f64(arr).copy()
    out[np.isclose(out, sentinel, rtol=0, atol=abs(sentinel) * 1e-6 + 1.0)] = np.nan
    return out


def _missing_summary(arr_norm):
    """Count NaNs in a normalised array."""
    return int(np.sum(np.isnan(arr_norm)))


def _diff_stats(a, b):
    """
    Given two NaN-normalised float64 arrays of the same shape, return a dict:
      n_both_nan   : positions where both are NaN (expected, not a problem)
      n_only_a_nan : positions where only a is NaN
      n_only_b_nan : positions where only b is NaN
      n_value_diff : positions where both are finite but differ > atol
      max_abs_diff : max |a-b| over value-diff positions
      max_abs_diff_all : max |a-b| ignoring NaNs entirely
    """
    nan_a = np.isnan(a)
    nan_b = np.isnan(b)
    both_nan   = nan_a & nan_b
    only_a_nan = nan_a & ~nan_b
    only_b_nan = ~nan_a & nan_b
    both_finite = ~nan_a & ~nan_b

    abs_diff = np.where(both_finite, np.abs(a - b), np.nan)
    return {
        "n_both_nan":        int(np.sum(both_nan)),
        "n_only_a_nan":      int(np.sum(only_a_nan)),
        "n_only_b_nan":      int(np.sum(only_b_nan)),
        "n_value_diff":      0,          # filled below after atol applied
        "max_abs_diff":      0.0,
        "max_abs_diff_all":  float(np.nanmax(abs_diff)) if np.any(both_finite) else np.nan,
        "_abs_diff":         abs_diff,   # kept for atol check
    }


# ---------------------------------------------------------------------------
# NC4 reader
# ---------------------------------------------------------------------------

def read_nc4(path):
    """Return (time_attrs, data_dict, nc4_fill_values).

    data_dict keys use the *bin* naming convention for easy pairing.
    nc4_fill_values maps bin_key -> fill_value used in nc4.
    """
    data   = {}
    fills  = {}
    tattrs = {}

    with Dataset(path, mode="r") as nc:
        # Global time attributes
        for bk, nk in TIME_MAP.items():
            tattrs[bk] = int(getattr(nc, nk, -9999))

        # Integer fields
        for bk, nk in INT_FIELDS:
            data[bk] = np.asarray(nc.variables[nk][:], dtype=np.int32)

        # Float fields
        for bk, nk in FLOAT_FIELDS:
            var  = nc.variables[nk]
            fv   = float(getattr(var, "_FillValue",
                         getattr(var, "missing_value", BIN_MISSING)))
            raw  = np.asarray(var[:], dtype=np.float64)
            fills[bk] = fv
            data[bk]  = raw   # keep raw; normalise during comparison

    return tattrs, data, fills


# ---------------------------------------------------------------------------
# Pair comparison
# ---------------------------------------------------------------------------

def compare_pair(bin_path, nc4_path, atol=1e-6):
    """Compare one pair.  Returns a list of issue strings (empty = PASS)."""

    issues = []

    # --- read ---
    b = read_ObsFcstAna(str(bin_path))
    t_nc4, n, nc4_fills = read_nc4(str(nc4_path))

    # ---- N_obs ----
    n_obs_bin = int(np.asarray(b["obs_species"]).size)
    n_obs_nc4 = int(np.asarray(n["obs_species"]).size)
    if n_obs_bin != n_obs_nc4:
        issues.append(f"N_obs MISMATCH  bin={n_obs_bin}  nc4={n_obs_nc4}")
        # Can't do element-wise below if sizes differ
        return issues, n_obs_bin, {}

    N = n_obs_bin

    # ---- time / header ----
    t_bin = b["date_time"]
    for bk, nk in TIME_MAP.items():
        vb = int(t_bin.get(bk, -9999))
        vn = int(t_nc4.get(bk, -9999))
        if vb != vn and not (vb == -9999 or vn == -9999):
            issues.append(f"time[{bk}]  bin={vb}  nc4={vn}")

    # ---- integer fields ----
    field_details = {}
    for bk, _ in INT_FIELDS:
        ab = np.asarray(b[bk], dtype=np.int32)
        an = np.asarray(n[bk], dtype=np.int32)
        if ab.shape != an.shape:
            issues.append(f"{bk}  shape mismatch  bin={ab.shape}  nc4={an.shape}")
            continue
        bad = ab != an
        n_bad = int(np.sum(bad))
        max_diff = int(np.max(np.abs(ab.astype(np.int64) - an.astype(np.int64)))) if n_bad else 0
        field_details[bk] = {"type": "int", "n_bad": n_bad, "max_diff": max_diff}
        if n_bad:
            issues.append(f"{bk}  {n_bad}/{N} values differ  max_abs={max_diff}")

    # ---- float fields ----
    for bk, _ in FLOAT_FIELDS:
        ab_raw  = _to_f64(b[bk])
        an_raw  = _to_f64(n[bk])
        nc4_fv  = nc4_fills.get(bk, MAPL_UNDEF)

        # Normalise both sides to NaN-based representation
        ab_norm = _normalise(ab_raw, BIN_MISSING)
        an_norm = _normalise(an_raw, nc4_fv)

        n_miss_bin = _missing_summary(ab_norm)
        n_miss_nc4 = _missing_summary(an_norm)

        stats = _diff_stats(ab_norm, an_norm)

        # Apply atol to the finite-vs-finite diffs
        value_diff_mask = stats["_abs_diff"] > atol
        n_val_diff  = int(np.sum(value_diff_mask))
        max_val_diff = float(np.nanmax(stats["_abs_diff"][value_diff_mask])) if n_val_diff else 0.0

        stats["n_value_diff"] = n_val_diff
        stats["max_abs_diff"] = max_val_diff
        del stats["_abs_diff"]

        field_details[bk] = {
            "type":           "float",
            "N":              N,
            "n_miss_bin":     n_miss_bin,
            "n_miss_nc4":     n_miss_nc4,
            "nc4_fill_value": nc4_fv,
            **stats,
        }

        # Report any problems
        if n_miss_bin != n_miss_nc4:
            issues.append(
                f"{bk}  missing count mismatch  "
                f"bin(NaN/-9999)={n_miss_bin}  nc4(NaN/MAPL_UNDEF)={n_miss_nc4}"
            )
        if stats["n_only_a_nan"]:
            idx = np.where(np.isnan(ab_norm) & ~np.isnan(an_norm))[0]
            lines = [f"{bk}  {stats['n_only_a_nan']} positions NaN in bin only:"]
            for pos in idx:
                lines.append(
                    f"    idx={pos:6d}  bin=NaN  nc4={an_norm[pos]:.8g}  "
                    f"(nc4_raw={an_raw[pos]:.8g}  bin_raw={ab_raw[pos]:.8g})"
                )
            issues.append("\n".join(lines))
        if stats["n_only_b_nan"]:
            idx = np.where(~np.isnan(ab_norm) & np.isnan(an_norm))[0]
            lines = [f"{bk}  {stats['n_only_b_nan']} positions NaN in nc4 only:"]
            for pos in idx:
                lines.append(
                    f"    idx={pos:6d}  bin={ab_norm[pos]:.8g}  nc4=NaN  "
                    f"(bin_raw={ab_raw[pos]:.8g}  nc4_raw={an_raw[pos]:.8g})"
                )
            issues.append("\n".join(lines))
        if n_val_diff:
            abs_diff_vals = np.where(value_diff_mask, np.abs(ab_norm - an_norm), 0.0)
            worst_idx = np.argsort(abs_diff_vals)[-min(n_val_diff, 10):][::-1]
            lines = [f"{bk}  {n_val_diff}/{N} values differ (atol={atol:.2g})  max_abs={max_val_diff:.6g}:"]
            for pos in worst_idx:
                lines.append(
                    f"    idx={pos:6d}  bin={ab_norm[pos]:.8g}  nc4={an_norm[pos]:.8g}  "
                    f"abs_diff={abs(ab_norm[pos]-an_norm[pos]):.4g}"
                )
            issues.append("\n".join(lines))

    return issues, N, field_details


# ---------------------------------------------------------------------------
# Directory scan
# ---------------------------------------------------------------------------

def find_pairs(directory):
    d = Path(directory)
    bin_files = sorted(d.glob("*.bin"))
    nc4_set   = {f.stem: f for f in d.glob("*.nc4")}

    pairs         = []
    missing_nc4   = []
    missing_bin   = []

    for bf in bin_files:
        stem = bf.stem  # e.g. LS_DA...20200501_0300z
        if stem in nc4_set:
            pairs.append((bf, nc4_set[stem]))
        else:
            missing_nc4.append(bf)

    paired_stems = {bf.stem for bf, _ in pairs}
    for stem, nf in nc4_set.items():
        if stem not in paired_stems:
            missing_bin.append(nf)

    return pairs, missing_nc4, missing_bin


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compare ObsFcstAna .bin vs .nc4 pairs in a flat directory.")
    parser.add_argument("dir", nargs="?", default=str(DEFAULT_DIR),
                        help=f"Directory to scan (default: {DEFAULT_DIR})")
    parser.add_argument("--atol", type=float, default=1e-6,
                        help="Absolute tolerance for float value comparison (default 1e-6)")
    parser.add_argument("--out", default=None,
                        help="Path for CSV report (default: <dir>/bin_nc4_diff_report.csv)")
    args = parser.parse_args()

    target_dir = Path(args.dir).expanduser().resolve()
    if not target_dir.is_dir():
        print(f"ERROR: not a directory: {target_dir}", file=sys.stderr)
        return 1

    out_csv = Path(args.out) if args.out else target_dir / "bin_nc4_diff_report.csv"

    pairs, missing_nc4, missing_bin = find_pairs(target_dir)

    print(f"\n{'='*70}")
    print(f"Directory : {target_dir}")
    print(f"Pairs     : {len(pairs)}")
    print(f"No nc4    : {len(missing_nc4)}")
    print(f"No bin    : {len(missing_bin)}")
    print(f"atol      : {args.atol:.2g}")
    print(f"{'='*70}\n")

    if missing_nc4:
        print("Files with no matching .nc4:")
        for p in missing_nc4:
            print(f"  {p.name}")
    if missing_bin:
        print("Files with no matching .bin:")
        for p in missing_bin:
            print(f"  {p.name}")

    # CSV accumulator
    csv_rows = []

    n_pass = n_fail = 0
    all_field_diffs = {}   # field -> total n_bad across all files

    for i, (bf, nf) in enumerate(pairs, 1):
        label = bf.stem
        print(f"[{i:02d}/{len(pairs)}]  {label}")

        issues, N, details = compare_pair(bf, nf, atol=args.atol)

        # Accumulate per-field totals
        for fld, d in details.items():
            if fld not in all_field_diffs:
                all_field_diffs[fld] = {"n_bad": 0, "n_miss_mismatch": 0, "files": []}
            if d.get("type") == "int":
                if d["n_bad"]:
                    all_field_diffs[fld]["n_bad"] += d["n_bad"]
                    all_field_diffs[fld]["files"].append(label)
            else:
                bad = d.get("n_value_diff", 0) + d.get("n_only_a_nan", 0) + d.get("n_only_b_nan", 0)
                miss_mm = abs(d.get("n_miss_bin", 0) - d.get("n_miss_nc4", 0))
                all_field_diffs[fld]["n_bad"] += bad
                all_field_diffs[fld]["n_miss_mismatch"] += miss_mm
                if bad or miss_mm:
                    all_field_diffs[fld]["files"].append(label)

        # CSV row per file per field
        for fld, d in details.items():
            row = {"file": label, "field": fld, "N_obs": N}
            if d.get("type") == "int":
                row.update({"n_bad": d["n_bad"], "max_diff": d.get("max_diff", ""),
                             "n_miss_bin": "", "n_miss_nc4": "",
                             "n_only_bin_nan": "", "n_only_nc4_nan": "",
                             "n_value_diff": "", "max_abs_diff": ""})
            else:
                row.update({
                    "n_bad":          d.get("n_value_diff", 0),
                    "max_diff":       d.get("max_abs_diff", ""),
                    "n_miss_bin":     d.get("n_miss_bin", ""),
                    "n_miss_nc4":     d.get("n_miss_nc4", ""),
                    "n_only_bin_nan": d.get("n_only_a_nan", ""),
                    "n_only_nc4_nan": d.get("n_only_b_nan", ""),
                    "n_value_diff":   d.get("n_value_diff", ""),
                    "max_abs_diff":   d.get("max_abs_diff", ""),
                })
            csv_rows.append(row)

        if issues:
            n_fail += 1
            for msg in issues:
                print(f"        !! {msg}")
        else:
            n_pass += 1
            print(f"         → PASS  (N_obs={N})")

    # ---- summary table ----
    print(f"\n{'='*70}")
    print(f"SUMMARY   pass={n_pass}  fail={n_fail}  total={len(pairs)}")
    print(f"{'='*70}")

    any_field_issue = False
    for fld, d in all_field_diffs.items():
        total_bad = d["n_bad"] + d["n_miss_mismatch"]
        if total_bad:
            any_field_issue = True
            print(f"  {fld:20s}  total_bad_elements={d['n_bad']:6d}  "
                  f"missing_count_mismatches={d['n_miss_mismatch']:4d}  "
                  f"affected_files={len(d['files'])}")
    if not any_field_issue:
        print("  All fields match across all pairs.")

    # ---- write CSV ----
    fieldnames = ["file", "field", "N_obs", "n_bad", "max_diff",
                  "n_miss_bin", "n_miss_nc4", "n_only_bin_nan", "n_only_nc4_nan",
                  "n_value_diff", "max_abs_diff"]
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    print(f"\nCSV report : {out_csv}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
