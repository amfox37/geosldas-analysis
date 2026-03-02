#!/usr/bin/env python3
"""
Build IMS vs GEOSldas OL/DA verification outputs with map-ready cell metrics.

What this script does:
1) Reads regridded IMS-on-M36 daily category files by year.
2) Reads daily GEOSldas OL and DA cat files.
3) Uses a fair paired mask (same valid collocated cells for OL and DA each day).
4) Accumulates per-cell contingency counts (A/B/C/D/N) for map-making.
5) Writes map-ready NetCDF + scope metadata CSV.
6) Writes the same table products used in the notebook workflow:
   - daily extraction table
   - paired daily OL-vs-DA table
   - bootstrap comparison table with CI

Speed-oriented design choices:
- Uses netCDF4 readers for low overhead.
- Works in representative-cell 1D vectors (no per-day full-grid rebuild).
- Updates per-cell counts directly for requested scopes.

Example run (Discover, your requested config):
  python run_ims_ol_da_cell_metrics.py \
    --domain SMAP_EASEv2_M36_GLOBAL \
    --year-start 2000 \
    --year-end 2024 \
    --ims-regrid-dir /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output \
    --ims-regrid-template 'ims_snowcover_24km_{year}_on_m36_nearest.nc4' \
    --ims-var ims_category \
    --ol-exp-name LS_OLv8_M36 \
    --ol-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36 \
    --da-exp-name LS_DAv8_M36 \
    --da-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36 \
    --output-dir /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from netCDF4 import Dataset


SEASON_ORDER = ["DJF", "MAM", "JJA", "SON"]
SEASON_TO_CODE = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}
CODE_TO_SEASON = {v: k for k, v in SEASON_TO_CODE.items()}


DEFAULT_EXPERIMENTS = {
    "OL": {
        "exp_name": "LS_OLv8_M36",
        "run_root": Path("/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36"),
    },
    "DA": {
        "exp_name": "LS_DAv8_M36",
        "run_root": Path("/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36"),
    },
}

DEFAULT_DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
DEFAULT_YEAR_START = 2000
DEFAULT_YEAR_END = 2024

DEFAULT_IMS_REGRID_DIR = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output")
DEFAULT_IMS_REGRID_TEMPLATE = "ims_snowcover_24km_{year}_on_m36_nearest.nc4"
DEFAULT_IMS_VAR = "ims_category"

DEFAULT_MODEL_SCF_VAR_CANDIDATES = ("FRLANDSNO", "FRSNO", "SNCOVFR", "SNOWCOVERFR", "SCF")
DEFAULT_SCF_THRESHOLD = 0.5

DEFAULT_AUTO_INFER_IMS_CODES = True
DEFAULT_IMS_SNOW_CODES = {4}
DEFAULT_IMS_NO_SNOW_CODES = {0, 1, 2, 3}
DEFAULT_IMS_FILL_VALUES = {-32768}

DEFAULT_N_BOOTSTRAP = 2000
DEFAULT_BOOTSTRAP_SEED = 42
DEFAULT_CI_LOW = 2.5
DEFAULT_CI_HIGH = 97.5
DEFAULT_MIN_IMS_SNOW_DAYS = 10

METRIC_ORDER = [
    "accuracy",
    "hit_rate",
    "miss_rate",
    "false_alarm_ratio",
    "correct_rejection_rate",
]

COUNT_ORDER = ["A", "B", "C", "D", "N"]


def parse_int_set(text: str) -> set[int]:
    """Parse comma-separated ints into a Python set."""
    out: set[int] = set()
    for tok in str(text).split(","):
        s = tok.strip()
        if not s:
            continue
        out.add(int(s))
    return out


def parse_str_tuple(text: str) -> tuple[str, ...]:
    """Parse comma-separated strings into a tuple, preserving order."""
    vals = [x.strip() for x in str(text).split(",") if x.strip()]
    if not vals:
        raise ValueError("At least one model variable candidate is required")
    return tuple(vals)


def season_name(ts: pd.Timestamp) -> str:
    """Return climatological season name for a timestamp."""
    m = int(ts.month)
    if m in (12, 1, 2):
        return "DJF"
    if m in (3, 4, 5):
        return "MAM"
    if m in (6, 7, 8):
        return "JJA"
    return "SON"


def find_repo_root(start: Path) -> Path:
    """Find repo root by locating common/python/io/read_GEOSldas.py upward."""
    p = start.resolve()
    for c in [p] + list(p.parents):
        if (c / "common/python/io/read_GEOSldas.py").exists():
            return c
    raise FileNotFoundError("Could not locate repo root containing common/python/io/read_GEOSldas.py")


def locate_tilecoord_file(run_root: Path, exp_name: str, domain: str, explicit: Path | None = None) -> Path:
    """Locate a tilecoord file using common GEOSldas directory conventions."""
    if explicit is not None:
        if explicit.exists():
            return explicit
        raise FileNotFoundError(f"--tilecoord does not exist: {explicit}")

    candidates = [
        run_root / "output" / domain / "rc_out" / f"{exp_name}.ldas_tilecoord.bin",
        run_root / exp_name / "output" / domain / "rc_out" / f"{exp_name}.ldas_tilecoord.bin",
        run_root / f"{exp_name}.ldas_tilecoord.bin",
    ]
    for p in candidates:
        if p.exists():
            return p

    raise FileNotFoundError(
        "Could not find tilecoord file. Checked: " + ", ".join(str(x) for x in candidates)
    )


def locate_daily_cat_file(run_root: Path, exp_name: str, domain: str, day: pd.Timestamp) -> Path | None:
    """Locate GEOSldas daily cat file for one date."""
    y = f"Y{day.year:04d}"
    m = f"M{day.month:02d}"
    stamp = day.strftime("%Y%m%d")
    fname = f"{exp_name}.tavg24_1d_lnd_Nt.{stamp}_1200z.nc4"
    candidates = [
        run_root / "output" / domain / "cat" / "ens_avg" / y / m / fname,
        run_root / exp_name / "output" / domain / "cat" / "ens_avg" / y / m / fname,
        run_root / "cat" / "ens_avg" / y / m / fname,
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def choose_representative_tile_per_cell(tc: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Deterministically pick one tile per EASE cell.

    Rule:
    - Primary: highest frac_cell in the cell.
    - Tie-break: smallest tile_id.
    """
    i_indg = np.asarray(tc["i_indg"], dtype=np.int64)
    j_indg = np.asarray(tc["j_indg"], dtype=np.int64)
    tile_id = np.asarray(tc["tile_id"], dtype=np.int64)
    frac_cell = np.asarray(tc["frac_cell"], dtype=np.float64)

    frac_sort = np.where(np.isfinite(frac_cell), frac_cell, -np.inf)

    nx = int(i_indg.max()) + 1
    ny = int(j_indg.max()) + 1
    cell_code = j_indg * nx + i_indg

    order = np.lexsort((tile_id, -frac_sort, cell_code))
    code_sorted = cell_code[order]

    first = np.empty(code_sorted.size, dtype=bool)
    first[0] = True
    first[1:] = code_sorted[1:] != code_sorted[:-1]
    rep_idx = order[first]

    return {
        "rep_idx": rep_idx.astype(np.int64),
        "rep_i": i_indg[rep_idx].astype(np.int32),
        "rep_j": j_indg[rep_idx].astype(np.int32),
        "rep_tile_id": tile_id[rep_idx].astype(np.int32),
        "rep_frac_cell": frac_cell[rep_idx].astype(np.float32),
        "rep_lat": np.asarray(tc["com_lat"][rep_idx], dtype=np.float64),
        "rep_lon": np.asarray(tc["com_lon"][rep_idx], dtype=np.float64),
        "rep_elev_m": np.asarray(tc["elev"][rep_idx], dtype=np.float64),
        "nx": np.int32(nx),
        "ny": np.int32(ny),
    }


def find_model_var_in_dataset(ds: Dataset, candidates: tuple[str, ...]) -> str:
    """Find the first model SCF variable present in the file."""
    for name in candidates:
        if name in ds.variables:
            return name
    raise KeyError(f"None of model SCF variable candidates found: {candidates}")


def _read_model_var_first_time(var) -> np.ndarray:
    """
    Read model variable as a 1D tile vector.

    Expected common layouts:
    - (tile,)
    - (time, tile), use time index 0
    """
    if getattr(var, "ndim", None) == 1:
        return np.asarray(var[:], dtype=np.float32).reshape(-1)

    if getattr(var, "ndim", None) == 2:
        # Daily files usually have one time record.
        return np.asarray(var[0, :], dtype=np.float32).reshape(-1)

    raise ValueError(
        f"Unsupported model variable dimensions: ndim={getattr(var, 'ndim', None)}, dims={getattr(var, 'dimensions', None)}"
    )


def read_model_scf_rep_values(
    nc_path: Path,
    rep_idx: np.ndarray,
    candidates: tuple[str, ...],
    forced_var_name: str | None = None,
) -> tuple[np.ndarray, str]:
    """
    Read one model file and return representative-cell values.

    Returns
    -------
    values : np.ndarray shape (n_rep,)
    used_var : str
    """
    with Dataset(nc_path, "r") as ds:
        used_var = forced_var_name if forced_var_name is not None else find_model_var_in_dataset(ds, candidates)
        if used_var not in ds.variables:
            raise KeyError(f"Model SCF variable '{used_var}' not found in {nc_path}")

        vals = _read_model_var_first_time(ds.variables[used_var])

    if vals.size <= int(rep_idx.max()):
        raise ValueError(
            f"Model tile vector length {vals.size} smaller than representative index max {int(rep_idx.max())}"
        )

    sel = np.asarray(vals[rep_idx], dtype=np.float32)
    # GEOS fill and physically invalid values.
    sel[(sel > 1e14) | (sel < 0.0)] = np.nan
    return sel, used_var


def _decode_time_like(values: np.ndarray, units: str) -> pd.DatetimeIndex | None:
    """Decode common numeric time encodings if possible."""
    if np.issubdtype(values.dtype, np.datetime64):
        return pd.DatetimeIndex(pd.to_datetime(values))

    if "since" in str(units):
        try:
            base_txt = str(units).split("since", 1)[1].strip().split()[0]
            base = pd.Timestamp(base_txt)
            return pd.DatetimeIndex(base + pd.to_timedelta(values.astype(float), unit="D"))
        except Exception:
            return None

    return None


def get_ims_time_dim_and_dates(ds: Dataset, ims_var: str, year: int) -> tuple[str, pd.DatetimeIndex]:
    """Resolve IMS time dimension and per-slice dates robustly across variants."""
    v = ds.variables[ims_var]
    if getattr(v, "ndim", None) != 3:
        raise ValueError(f"IMS var {ims_var} must be 3D, found ndim={getattr(v, 'ndim', None)}")

    time_dim = str(v.dimensions[0])
    n_time = int(v.shape[0])

    if "doy" in ds.variables:
        doy_var = ds.variables["doy"]
        if getattr(doy_var, "ndim", None) == 1 and len(doy_var) == n_time:
            doy = np.asarray(doy_var[:], dtype=float)
            base = pd.Timestamp(f"{year}-01-01")
            dates = pd.DatetimeIndex(base + pd.to_timedelta(doy - 1.0, unit="D"))
            return time_dim, dates

    for name in (time_dim, "time", "day_of_year"):
        if name not in ds.variables:
            continue
        tv = ds.variables[name]
        if getattr(tv, "ndim", None) != 1 or len(tv) != n_time:
            continue
        vals = np.asarray(tv[:])
        units = getattr(tv, "units", "")
        decoded = _decode_time_like(vals, str(units))
        if decoded is not None:
            return time_dim, decoded

    base = pd.Timestamp(f"{year}-01-01")
    dates = pd.DatetimeIndex(base + pd.to_timedelta(np.arange(n_time), unit="D"))
    return time_dim, dates


def infer_ims_code_sets(unique_codes: set[int], fill_values: set[int]) -> tuple[set[int], set[int]]:
    """Infer snow/no-snow code sets from observed IMS categories."""
    code_set = {int(c) for c in unique_codes if int(c) not in fill_values}
    if not code_set:
        raise ValueError("No usable IMS category codes available for inference")

    if code_set.issubset({0, 1}):
        return {1}, {0}
    if 4 in code_set:
        return {4}, {c for c in code_set if c != 4}
    if 1 in code_set:
        return {1}, {c for c in code_set if c != 1}

    raise ValueError(
        "Could not infer IMS snow/no-snow code sets from observed codes: "
        f"{sorted(code_set)}. Set --no-auto-infer-ims-codes and provide explicit sets."
    )


def ims_category_to_binary_scf(
    cat_grid: np.ndarray,
    snow_codes: set[int],
    no_snow_codes: set[int],
    fill_values: set[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert IMS category grid into binary snow cover fraction (1/0/NaN).

    Returns
    -------
    obs_scf : float32 array with {1.0, 0.0, NaN}
    valid_mask : bool array where category is known snow/no-snow and not fill
    unknown_mask : bool array where category is finite/not-fill but not in known sets
    arr_i : int32 rounded category codes
    """
    arr = np.asarray(cat_grid, dtype=np.float32)
    finite = np.isfinite(arr)

    arr_i = np.full(arr.shape, -99999, dtype=np.int32)
    arr_i[finite] = np.rint(arr[finite]).astype(np.int32)

    valid_not_fill = finite.copy()
    for fv in fill_values:
        valid_not_fill &= arr_i != int(fv)

    known_codes = set(snow_codes) | set(no_snow_codes)
    known_mask = np.isin(arr_i, np.array(sorted(known_codes), dtype=np.int32))

    valid_mask = valid_not_fill & known_mask
    unknown_mask = valid_not_fill & (~known_mask)

    obs_scf = np.full(arr.shape, np.nan, dtype=np.float32)
    snow_mask = valid_mask & np.isin(arr_i, np.array(sorted(snow_codes), dtype=np.int32))
    nosnow_mask = valid_mask & np.isin(arr_i, np.array(sorted(no_snow_codes), dtype=np.int32))
    obs_scf[snow_mask] = np.float32(1.0)
    obs_scf[nosnow_mask] = np.float32(0.0)

    return obs_scf, valid_mask, unknown_mask, arr_i


def contingency_masks_from_vectors(
    obs_bin_rep: np.ndarray,
    mod_scf_rep: np.ndarray,
    valid_mask: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build A/B/C/D boolean masks and return also final valid mask."""
    valid = np.asarray(valid_mask, dtype=bool)
    valid &= np.isfinite(obs_bin_rep)
    valid &= np.isfinite(mod_scf_rep)

    snow_obs = obs_bin_rep >= np.float32(0.5)
    snow_mod = mod_scf_rep >= np.float32(threshold)

    A = valid & snow_mod & snow_obs
    B = valid & snow_mod & (~snow_obs)
    C = valid & (~snow_mod) & snow_obs
    D = valid & (~snow_mod) & (~snow_obs)
    return A, B, C, D, valid


def counts_from_masks(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray) -> dict[str, int]:
    """Convert contingency masks to integer counts."""
    a = int(np.sum(A))
    b = int(np.sum(B))
    c = int(np.sum(C))
    d = int(np.sum(D))
    return {"A": a, "B": b, "C": c, "D": d, "N": int(a + b + c + d)}


def snow_scores_from_counts(A: int, B: int, C: int, D: int) -> dict[str, float]:
    """Compute contingency-table verification metrics with NaN-safe zero-denominator handling."""
    N = A + B + C + D
    accuracy = (A + D) / N if N else np.nan

    den_snow_obs = A + C
    hit_rate = A / den_snow_obs if den_snow_obs else np.nan
    miss_rate = C / den_snow_obs if den_snow_obs else np.nan

    den_snow_mod = A + B
    false_alarm_ratio = B / den_snow_mod if den_snow_mod else np.nan

    den_nosnow_obs = B + D
    correct_rejection_rate = D / den_nosnow_obs if den_nosnow_obs else np.nan

    return {
        "accuracy": float(accuracy) if np.isfinite(accuracy) else np.nan,
        "hit_rate": float(hit_rate) if np.isfinite(hit_rate) else np.nan,
        "miss_rate": float(miss_rate) if np.isfinite(miss_rate) else np.nan,
        "false_alarm_ratio": float(false_alarm_ratio) if np.isfinite(false_alarm_ratio) else np.nan,
        "correct_rejection_rate": float(correct_rejection_rate) if np.isfinite(correct_rejection_rate) else np.nan,
    }


def _pct_finite(x, p) -> float:
    a = np.asarray(x, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    return float(np.percentile(a, p))


def bootstrap_compare_table_from_pair_days(
    days_df: pd.DataFrame,
    exp_ol: str,
    exp_da: str,
    n_boot: int,
    seed: int,
    ci_low: float,
    ci_high: float,
) -> pd.DataFrame:
    """
    Paired day-block bootstrap for OL-vs-DA contingency metrics.

    Each replicate resamples day rows with replacement and keeps OL/DA paired
    counts from the same selected day.
    """
    req_cols = [
        f"A_{exp_ol}",
        f"B_{exp_ol}",
        f"C_{exp_ol}",
        f"D_{exp_ol}",
        f"A_{exp_da}",
        f"B_{exp_da}",
        f"C_{exp_da}",
        f"D_{exp_da}",
    ]
    missing = [c for c in req_cols if c not in days_df.columns]
    if missing:
        raise KeyError(f"Missing required paired columns: {missing}")

    n_days = len(days_df)
    if n_days == 0:
        return pd.DataFrame()

    ol_counts = days_df[[f"A_{exp_ol}", f"B_{exp_ol}", f"C_{exp_ol}", f"D_{exp_ol}"]].to_numpy(dtype=np.int64)
    da_counts = days_df[[f"A_{exp_da}", f"B_{exp_da}", f"C_{exp_da}", f"D_{exp_da}"]].to_numpy(dtype=np.int64)

    A_ol, B_ol, C_ol, D_ol = np.sum(ol_counts, axis=0)
    A_da, B_da, C_da, D_da = np.sum(da_counts, axis=0)

    point_ol = snow_scores_from_counts(int(A_ol), int(B_ol), int(C_ol), int(D_ol))
    point_da = snow_scores_from_counts(int(A_da), int(B_da), int(C_da), int(D_da))

    rng = np.random.default_rng(int(seed))

    boot_ol = {m: [] for m in METRIC_ORDER}
    boot_da = {m: [] for m in METRIC_ORDER}
    boot_delta = {m: [] for m in METRIC_ORDER}

    for _ in range(int(n_boot)):
        idx = rng.integers(0, n_days, size=n_days)

        b_ol = np.sum(ol_counts[idx, :], axis=0)
        b_da = np.sum(da_counts[idx, :], axis=0)

        m_ol = snow_scores_from_counts(int(b_ol[0]), int(b_ol[1]), int(b_ol[2]), int(b_ol[3]))
        m_da = snow_scores_from_counts(int(b_da[0]), int(b_da[1]), int(b_da[2]), int(b_da[3]))

        for m in METRIC_ORDER:
            v_ol = m_ol[m]
            v_da = m_da[m]
            boot_ol[m].append(v_ol)
            boot_da[m].append(v_da)
            if np.isfinite(v_ol) and np.isfinite(v_da):
                boot_delta[m].append(v_da - v_ol)
            else:
                boot_delta[m].append(np.nan)

    rows = []
    for m in METRIC_ORDER:
        ol = point_ol[m]
        da = point_da[m]
        delta = (da - ol) if (np.isfinite(ol) and np.isfinite(da)) else np.nan
        rows.append(
            {
                "metric": m,
                "ol": float(ol) if np.isfinite(ol) else np.nan,
                "ol_ci_lo": _pct_finite(boot_ol[m], ci_low),
                "ol_ci_hi": _pct_finite(boot_ol[m], ci_high),
                "da": float(da) if np.isfinite(da) else np.nan,
                "da_ci_lo": _pct_finite(boot_da[m], ci_low),
                "da_ci_hi": _pct_finite(boot_da[m], ci_high),
                "delta_da_minus_ol": float(delta) if np.isfinite(delta) else np.nan,
                "delta_ci_lo": _pct_finite(boot_delta[m], ci_low),
                "delta_ci_hi": _pct_finite(boot_delta[m], ci_high),
                "n_days": int(n_days),
                "A_ol": int(A_ol),
                "B_ol": int(B_ol),
                "C_ol": int(C_ol),
                "D_ol": int(D_ol),
                "N_ol": int(A_ol + B_ol + C_ol + D_ol),
                "A_da": int(A_da),
                "B_da": int(B_da),
                "C_da": int(C_da),
                "D_da": int(D_da),
                "N_da": int(A_da + B_da + C_da + D_da),
            }
        )

    return pd.DataFrame(rows)


def build_scope_metadata(year_start: int, year_end: int) -> tuple[pd.DataFrame, dict[str, dict]]:
    """Build scope table and index mappings used for per-cell accumulators."""
    rows = []
    sid = 0

    idx_all = sid
    rows.append({"scope_id": sid, "scope": "ALL_PERIOD", "year": -1, "season": "ALL", "scope_type_code": 0, "season_code": -1})
    sid += 1

    season_map: dict[str, int] = {}
    for s in SEASON_ORDER:
        season_map[s] = sid
        rows.append({"scope_id": sid, "scope": "SEASON_ALL_YEARS", "year": -1, "season": s, "scope_type_code": 1, "season_code": SEASON_TO_CODE[s]})
        sid += 1

    year_map: dict[int, int] = {}
    for y in range(year_start, year_end + 1):
        year_map[y] = sid
        rows.append({"scope_id": sid, "scope": "YEAR", "year": int(y), "season": "ALL", "scope_type_code": 2, "season_code": -1})
        sid += 1

    year_season_map: dict[tuple[int, str], int] = {}
    for y in range(year_start, year_end + 1):
        for s in SEASON_ORDER:
            year_season_map[(int(y), s)] = sid
            rows.append(
                {
                    "scope_id": sid,
                    "scope": "YEAR_SEASON",
                    "year": int(y),
                    "season": s,
                    "scope_type_code": 3,
                    "season_code": SEASON_TO_CODE[s],
                }
            )
            sid += 1

    meta = pd.DataFrame(rows)
    maps = {
        "all": {"ALL_PERIOD": idx_all},
        "season": season_map,
        "year": year_map,
        "year_season": year_season_map,
    }
    return meta, maps


def init_cell_count_arrays(n_exp: int, n_scope: int, n_cell: int) -> dict[str, np.ndarray]:
    """Initialize per-cell count accumulators."""
    shape = (n_exp, n_scope, n_cell)
    return {k: np.zeros(shape, dtype=np.int32) for k in COUNT_ORDER}


def update_cell_counts_for_day(
    cell_counts: dict[str, np.ndarray],
    exp_idx: int,
    scope_ids: list[int],
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    N_valid: np.ndarray,
) -> None:
    """
    Increment per-cell accumulators for one experiment/day for all requested scopes.

    `A/B/C/D/N_valid` are boolean vectors over representative cells.
    """
    idx_a = np.flatnonzero(A)
    idx_b = np.flatnonzero(B)
    idx_c = np.flatnonzero(C)
    idx_d = np.flatnonzero(D)
    idx_n = np.flatnonzero(N_valid)

    for sid in scope_ids:
        if idx_a.size:
            arr = cell_counts["A"][exp_idx, sid, :]
            arr[idx_a] += 1
        if idx_b.size:
            arr = cell_counts["B"][exp_idx, sid, :]
            arr[idx_b] += 1
        if idx_c.size:
            arr = cell_counts["C"][exp_idx, sid, :]
            arr[idx_c] += 1
        if idx_d.size:
            arr = cell_counts["D"][exp_idx, sid, :]
            arr[idx_d] += 1
        if idx_n.size:
            arr = cell_counts["N"][exp_idx, sid, :]
            arr[idx_n] += 1


def counts_to_metrics_arrays(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    N: np.ndarray,
) -> dict[str, np.ndarray]:
    """Convert per-cell contingency counts to metric arrays."""
    A_f = A.astype(np.float64)
    B_f = B.astype(np.float64)
    C_f = C.astype(np.float64)
    D_f = D.astype(np.float64)
    N_f = N.astype(np.float64)

    out = {m: np.full(A.shape, np.nan, dtype=np.float32) for m in METRIC_ORDER}

    # Accuracy
    m = N_f > 0.0
    out["accuracy"][m] = ((A_f + D_f) / N_f)[m].astype(np.float32)

    # Hit / miss rates
    den = A_f + C_f
    m = den > 0.0
    out["hit_rate"][m] = (A_f / den)[m].astype(np.float32)
    out["miss_rate"][m] = (C_f / den)[m].astype(np.float32)

    # False alarm ratio
    den = A_f + B_f
    m = den > 0.0
    out["false_alarm_ratio"][m] = (B_f / den)[m].astype(np.float32)

    # Correct rejection rate
    den = B_f + D_f
    m = den > 0.0
    out["correct_rejection_rate"][m] = (D_f / den)[m].astype(np.float32)

    return out


def write_cell_counts_netcdf(
    out_nc: Path,
    rep: dict[str, np.ndarray],
    exp_keys: list[str],
    scope_meta: pd.DataFrame,
    cell_counts: dict[str, np.ndarray],
    threshold: float,
    ims_var: str,
    snow_codes: set[int],
    no_snow_codes: set[int],
    fill_values: set[int],
    eligible_mask: np.ndarray | None = None,
    extra_note: str = "",
) -> None:
    """Write map-ready per-cell counts + metrics to NetCDF."""
    ny = int(rep["ny"])
    nx = int(rep["nx"])
    n_cell = int(rep["rep_idx"].size)
    n_scope = int(scope_meta.shape[0])
    n_exp = len(exp_keys)

    metrics = counts_to_metrics_arrays(
        A=cell_counts["A"],
        B=cell_counts["B"],
        C=cell_counts["C"],
        D=cell_counts["D"],
        N=cell_counts["N"],
    )

    ds = Dataset(out_nc, "w", format="NETCDF4")
    try:
        ds.createDimension("experiment", n_exp)
        ds.createDimension("scope", n_scope)
        ds.createDimension("cell", n_cell)

        exp_idx = ds.createVariable("experiment_index", "i4", ("experiment",))
        exp_idx[:] = np.arange(n_exp, dtype=np.int32)
        exp_idx.long_name = "Experiment index"
        exp_idx.codes = ", ".join([f"{i}:{k}" for i, k in enumerate(exp_keys)])

        scope_id = ds.createVariable("scope_id", "i4", ("scope",))
        scope_id[:] = scope_meta["scope_id"].to_numpy(dtype=np.int32)

        scope_type_code = ds.createVariable("scope_type_code", "i4", ("scope",))
        scope_type_code[:] = scope_meta["scope_type_code"].to_numpy(dtype=np.int32)
        scope_type_code.codes = "0:ALL_PERIOD,1:SEASON_ALL_YEARS,2:YEAR,3:YEAR_SEASON"

        scope_year = ds.createVariable("scope_year", "i4", ("scope",), fill_value=-1)
        scope_year[:] = scope_meta["year"].to_numpy(dtype=np.int32)

        scope_season_code = ds.createVariable("scope_season_code", "i4", ("scope",), fill_value=-1)
        scope_season_code[:] = scope_meta["season_code"].to_numpy(dtype=np.int32)
        scope_season_code.codes = "-1:ALL,0:DJF,1:MAM,2:JJA,3:SON"

        rep_i = ds.createVariable("cell_i", "i4", ("cell",))
        rep_j = ds.createVariable("cell_j", "i4", ("cell",))
        rep_lat = ds.createVariable("cell_lat", "f4", ("cell",))
        rep_lon = ds.createVariable("cell_lon", "f4", ("cell",))
        rep_elev = ds.createVariable("cell_elev_m", "f4", ("cell",))
        rep_tile = ds.createVariable("cell_rep_tile_id", "i4", ("cell",))
        rep_frac = ds.createVariable("cell_rep_frac_cell", "f4", ("cell",))
        if eligible_mask is not None:
            rep_eligible = ds.createVariable("cell_eligible", "i1", ("cell",))
        else:
            rep_eligible = None

        rep_i[:] = rep["rep_i"].astype(np.int32)
        rep_j[:] = rep["rep_j"].astype(np.int32)
        rep_lat[:] = rep["rep_lat"].astype(np.float32)
        rep_lon[:] = rep["rep_lon"].astype(np.float32)
        rep_elev[:] = rep["rep_elev_m"].astype(np.float32)
        rep_tile[:] = rep["rep_tile_id"].astype(np.int32)
        rep_frac[:] = rep["rep_frac_cell"].astype(np.float32)
        if rep_eligible is not None:
            rep_eligible[:] = np.asarray(eligible_mask, dtype=np.int8)
            rep_eligible.long_name = "Eligibility flag after IMS snow-day filter (1=kept,0=excluded)"

        rep_i.long_name = "M36 x index of representative land cell"
        rep_j.long_name = "M36 y index of representative land cell"
        rep_lat.units = "degrees_north"
        rep_lon.units = "degrees_east"
        rep_elev.units = "m"

        # Counts
        for name in COUNT_ORDER:
            v = ds.createVariable(name, "i4", ("experiment", "scope", "cell"), zlib=True, complevel=4)
            v[:] = cell_counts[name].astype(np.int32)
            v.long_name = f"Per-cell contingency count {name}"

        # Metrics derived from counts.
        for name in METRIC_ORDER:
            v = ds.createVariable(name, "f4", ("experiment", "scope", "cell"), zlib=True, complevel=4)
            v[:] = metrics[name].astype(np.float32)
            v.long_name = f"Per-cell metric {name}"

        ds.title = "IMS vs GEOSldas OL/DA per-cell contingency counts and metrics"
        ds.history = f"{datetime.utcnow().isoformat()}Z created by run_ims_ol_da_cell_metrics.py"
        ds.grid_nx = int(nx)
        ds.grid_ny = int(ny)
        ds.model_threshold_scf = float(threshold)
        ds.ims_var = str(ims_var)
        ds.ims_snow_codes = ",".join(str(x) for x in sorted(snow_codes))
        ds.ims_no_snow_codes = ",".join(str(x) for x in sorted(no_snow_codes))
        ds.ims_fill_values = ",".join(str(x) for x in sorted(fill_values))
        ds.note = "Counts were accumulated using fair daily common masks where both OL and DA were finite."
        if str(extra_note).strip():
            ds.note = f"{ds.note} {str(extra_note).strip()}"

    finally:
        ds.close()


def maybe_write_parquet(df: pd.DataFrame, path: Path, label: str) -> None:
    """Write parquet if possible; otherwise keep going with CSV outputs."""
    try:
        df.to_parquet(path, index=False)
        print(f"Wrote {label} parquet: {path}")
    except Exception as exc:
        warnings.warn(f"Could not write parquet for {label} ({path}): {exc}")


def read_scope_meta_from_cell_counts_nc(nc_path: Path) -> pd.DataFrame:
    """Read scope metadata table from an existing cell-counts NetCDF file."""
    with Dataset(nc_path, "r") as ds:
        scope_id = np.asarray(ds.variables["scope_id"][:], dtype=np.int32)
        scope_type_code = np.asarray(ds.variables["scope_type_code"][:], dtype=np.int32)
        scope_year = np.asarray(ds.variables["scope_year"][:], dtype=np.int32)
        scope_season_code = np.asarray(ds.variables["scope_season_code"][:], dtype=np.int32)

    rows = []
    for sid, stc, yr, ssc in zip(scope_id, scope_type_code, scope_year, scope_season_code):
        if int(stc) == 0:
            scope = "ALL_PERIOD"
        elif int(stc) == 1:
            scope = "SEASON_ALL_YEARS"
        elif int(stc) == 2:
            scope = "YEAR"
        elif int(stc) == 3:
            scope = "YEAR_SEASON"
        else:
            scope = "UNKNOWN"

        season = "ALL" if int(ssc) < 0 else CODE_TO_SEASON.get(int(ssc), "ALL")
        rows.append(
            {
                "scope_id": int(sid),
                "scope": str(scope),
                "year": int(yr),
                "season": str(season),
                "scope_type_code": int(stc),
                "season_code": int(ssc),
            }
        )
    return pd.DataFrame(rows)


def read_counts_and_rep_from_cell_counts_nc(nc_path: Path) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, str]]:
    """Load per-cell count arrays and representative-cell metadata from an existing NetCDF."""
    with Dataset(nc_path, "r") as ds:
        cell_counts = {}
        for name in COUNT_ORDER:
            if name not in ds.variables:
                raise KeyError(f"{nc_path} missing required variable '{name}'")
            cell_counts[name] = np.asarray(ds.variables[name][:], dtype=np.int32)

        n_cell = int(ds.dimensions["cell"].size)
        rep = {
            "rep_idx": np.arange(n_cell, dtype=np.int64),
            "rep_i": np.asarray(ds.variables["cell_i"][:], dtype=np.int32),
            "rep_j": np.asarray(ds.variables["cell_j"][:], dtype=np.int32),
            "rep_lat": np.asarray(ds.variables["cell_lat"][:], dtype=np.float64),
            "rep_lon": np.asarray(ds.variables["cell_lon"][:], dtype=np.float64),
            "rep_elev_m": np.asarray(ds.variables["cell_elev_m"][:], dtype=np.float64),
            "rep_tile_id": np.asarray(ds.variables["cell_rep_tile_id"][:], dtype=np.int32),
            "rep_frac_cell": np.asarray(ds.variables["cell_rep_frac_cell"][:], dtype=np.float32),
            "nx": np.int32(getattr(ds, "grid_nx")),
            "ny": np.int32(getattr(ds, "grid_ny")),
        }

        attrs = {
            "model_threshold_scf": str(getattr(ds, "model_threshold_scf", DEFAULT_SCF_THRESHOLD)),
            "ims_var": str(getattr(ds, "ims_var", DEFAULT_IMS_VAR)),
            "ims_snow_codes": str(getattr(ds, "ims_snow_codes", ",".join(str(x) for x in sorted(DEFAULT_IMS_SNOW_CODES)))),
            "ims_no_snow_codes": str(
                getattr(ds, "ims_no_snow_codes", ",".join(str(x) for x in sorted(DEFAULT_IMS_NO_SNOW_CODES)))
            ),
            "ims_fill_values": str(
                getattr(ds, "ims_fill_values", ",".join(str(x) for x in sorted(DEFAULT_IMS_FILL_VALUES)))
            ),
        }

    return cell_counts, rep, attrs


def apply_min_ims_snow_days_filter_to_counts(
    cell_counts: dict[str, np.ndarray],
    scope_meta: pd.DataFrame,
    exp_keys: list[str],
    min_days: int,
    eligibility_experiment: str = "OL",
) -> np.ndarray:
    """
    Apply IMS snow-day eligibility to count arrays in-place.

    Eligibility uses ALL_PERIOD obs-snow days at each representative cell:
      obs_snow_days = A + C
    """
    n_cell = int(cell_counts["A"].shape[2])
    if int(min_days) <= 0:
        return np.ones(n_cell, dtype=bool)

    if "ALL_PERIOD" not in set(scope_meta["scope"]):
        raise KeyError("scope_meta does not contain ALL_PERIOD row")
    sid = int(scope_meta.loc[scope_meta["scope"] == "ALL_PERIOD", "scope_id"].iloc[0])

    if eligibility_experiment not in exp_keys:
        raise KeyError(f"eligibility_experiment={eligibility_experiment} not in exp_keys={exp_keys}")
    exp_idx = int(exp_keys.index(eligibility_experiment))

    obs_snow_days = (
        np.asarray(cell_counts["A"][exp_idx, sid, :], dtype=np.int64)
        + np.asarray(cell_counts["C"][exp_idx, sid, :], dtype=np.int64)
    )
    eligible = obs_snow_days >= int(min_days)

    ineligible = ~eligible
    if np.any(ineligible):
        for name in COUNT_ORDER:
            arr = cell_counts[name]
            arr[:, :, ineligible] = 0
            cell_counts[name] = arr

    return eligible


def build_comparison_table_from_cell_counts(
    cell_counts: dict[str, np.ndarray],
    scope_meta: pd.DataFrame,
    exp_keys: list[str],
) -> pd.DataFrame:
    """
    Build comparison table directly from per-cell scope counts.

    This path does not have per-day rows, so bootstrap CI columns are written as NaN.
    """
    if "OL" not in exp_keys or "DA" not in exp_keys:
        raise KeyError(f"exp_keys must contain OL and DA, got {exp_keys}")

    i_ol = int(exp_keys.index("OL"))
    i_da = int(exp_keys.index("DA"))

    rows = []
    for r in scope_meta.itertuples(index=False):
        sid = int(r.scope_id)
        A_ol = int(np.sum(cell_counts["A"][i_ol, sid, :], dtype=np.int64))
        B_ol = int(np.sum(cell_counts["B"][i_ol, sid, :], dtype=np.int64))
        C_ol = int(np.sum(cell_counts["C"][i_ol, sid, :], dtype=np.int64))
        D_ol = int(np.sum(cell_counts["D"][i_ol, sid, :], dtype=np.int64))
        A_da = int(np.sum(cell_counts["A"][i_da, sid, :], dtype=np.int64))
        B_da = int(np.sum(cell_counts["B"][i_da, sid, :], dtype=np.int64))
        C_da = int(np.sum(cell_counts["C"][i_da, sid, :], dtype=np.int64))
        D_da = int(np.sum(cell_counts["D"][i_da, sid, :], dtype=np.int64))

        m_ol = snow_scores_from_counts(A_ol, B_ol, C_ol, D_ol)
        m_da = snow_scores_from_counts(A_da, B_da, C_da, D_da)

        for metric in METRIC_ORDER:
            ol = m_ol[metric]
            da = m_da[metric]
            delta = (da - ol) if (np.isfinite(ol) and np.isfinite(da)) else np.nan
            rows.append(
                {
                    "scope": str(r.scope),
                    "year": int(r.year),
                    "season": str(r.season),
                    "metric": str(metric),
                    "ol": float(ol) if np.isfinite(ol) else np.nan,
                    "ol_ci_lo": np.nan,
                    "ol_ci_hi": np.nan,
                    "da": float(da) if np.isfinite(da) else np.nan,
                    "da_ci_lo": np.nan,
                    "da_ci_hi": np.nan,
                    "delta_da_minus_ol": float(delta) if np.isfinite(delta) else np.nan,
                    "delta_ci_lo": np.nan,
                    "delta_ci_hi": np.nan,
                    "n_days": np.nan,
                    "A_ol": int(A_ol),
                    "B_ol": int(B_ol),
                    "C_ol": int(C_ol),
                    "D_ol": int(D_ol),
                    "N_ol": int(A_ol + B_ol + C_ol + D_ol),
                    "A_da": int(A_da),
                    "B_da": int(B_da),
                    "C_da": int(C_da),
                    "D_da": int(D_da),
                    "N_da": int(A_da + B_da + C_da + D_da),
                }
            )

    return pd.DataFrame(rows)


def build_empty_daily_table() -> pd.DataFrame:
    """Return an empty daily table with expected columns."""
    cols = [
        "date",
        "year",
        "doy",
        "season",
        "experiment",
        "exp_name",
        "ims_file",
        "model_file",
        "model_file_found",
        "model_read_ok",
        "model_var",
        "A",
        "B",
        "C",
        "D",
        "N_valid",
        "accuracy",
        "hit_rate",
        "miss_rate",
        "false_alarm_ratio",
        "correct_rejection_rate",
        "N_land_mask",
        "N_ims_obs_valid",
        "N_ims_unknown_codes",
        "obs_scf_mean",
        "mod_scf_mean",
        "paired_common_mask_used",
        "N_common_mask",
        "mask_type",
    ]
    return pd.DataFrame(columns=cols)


def build_empty_pair_daily_table() -> pd.DataFrame:
    """Return an empty paired-daily table with expected columns."""
    cols = [
        "date",
        "year",
        "doy",
        "season",
        "A_OL",
        "B_OL",
        "C_OL",
        "D_OL",
        "N_valid_OL",
        "A_DA",
        "B_DA",
        "C_DA",
        "D_DA",
        "N_valid_DA",
    ]
    return pd.DataFrame(columns=cols)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Build IMS-vs-OL/DA map-ready per-cell metrics and summary tables."
    )

    ap.add_argument("--domain", default=DEFAULT_DOMAIN)
    ap.add_argument("--year-start", type=int, default=DEFAULT_YEAR_START)
    ap.add_argument("--year-end", type=int, default=DEFAULT_YEAR_END)

    ap.add_argument("--ims-regrid-dir", default=str(DEFAULT_IMS_REGRID_DIR))
    ap.add_argument("--ims-regrid-template", default=DEFAULT_IMS_REGRID_TEMPLATE)
    ap.add_argument("--ims-var", default=DEFAULT_IMS_VAR)

    ap.add_argument("--ol-exp-name", default=DEFAULT_EXPERIMENTS["OL"]["exp_name"])
    ap.add_argument("--ol-run-root", default=str(DEFAULT_EXPERIMENTS["OL"]["run_root"]))
    ap.add_argument("--da-exp-name", default=DEFAULT_EXPERIMENTS["DA"]["exp_name"])
    ap.add_argument("--da-run-root", default=str(DEFAULT_EXPERIMENTS["DA"]["run_root"]))

    ap.add_argument("--tilecoord", default=None, help="Optional explicit tilecoord path")

    ap.add_argument(
        "--model-var-candidates",
        default=",".join(DEFAULT_MODEL_SCF_VAR_CANDIDATES),
        help="Comma-separated candidate model SCF variable names",
    )
    ap.add_argument("--scf-threshold", type=float, default=DEFAULT_SCF_THRESHOLD)

    ap.add_argument(
        "--auto-infer-ims-codes",
        dest="auto_infer_ims_codes",
        action="store_true",
        default=DEFAULT_AUTO_INFER_IMS_CODES,
        help="Infer IMS snow/no-snow codes from sample categories",
    )
    ap.add_argument(
        "--no-auto-infer-ims-codes",
        dest="auto_infer_ims_codes",
        action="store_false",
        help="Disable auto inference and use explicit code sets",
    )
    ap.add_argument("--ims-snow-codes", default=",".join(str(x) for x in sorted(DEFAULT_IMS_SNOW_CODES)))
    ap.add_argument("--ims-no-snow-codes", default=",".join(str(x) for x in sorted(DEFAULT_IMS_NO_SNOW_CODES)))
    ap.add_argument("--ims-fill-values", default=",".join(str(x) for x in sorted(DEFAULT_IMS_FILL_VALUES)))

    ap.add_argument("--n-bootstrap", type=int, default=DEFAULT_N_BOOTSTRAP)
    ap.add_argument("--bootstrap-seed", type=int, default=DEFAULT_BOOTSTRAP_SEED)
    ap.add_argument("--ci-low", type=float, default=DEFAULT_CI_LOW)
    ap.add_argument("--ci-high", type=float, default=DEFAULT_CI_HIGH)
    ap.add_argument(
        "--min-ims-snow-days",
        type=int,
        default=DEFAULT_MIN_IMS_SNOW_DAYS,
        help="Keep only cells with at least this many IMS observed-snow days (A+C over ALL_PERIOD).",
    )
    ap.add_argument(
        "--reuse-cell-counts-nc",
        default=None,
        help="Fast path: read existing ims_ol_da_cell_counts_metrics_*.nc4 and regenerate stats without model re-read.",
    )
    ap.add_argument(
        "--reuse-scope-meta-csv",
        default=None,
        help="Optional scope metadata CSV for --reuse-cell-counts-nc mode. If omitted, read scope metadata from NetCDF.",
    )

    ap.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: repo/projects/IMS/output)",
    )
    ap.add_argument(
        "--output-tag",
        default=None,
        help="Optional custom filename tag. Default is domain/year/threshold based.",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs",
    )

    return ap.parse_args()


def main() -> None:
    args = parse_args()

    if args.year_end < args.year_start:
        raise ValueError(f"Invalid year range: {args.year_start}..{args.year_end}")

    model_var_candidates = parse_str_tuple(args.model_var_candidates)
    ims_snow_codes = parse_int_set(args.ims_snow_codes)
    ims_no_snow_codes = parse_int_set(args.ims_no_snow_codes)
    ims_fill_values = parse_int_set(args.ims_fill_values)

    repo_root = find_repo_root(Path(__file__).resolve())
    sys.path.insert(0, str(repo_root / "common/python/io"))
    from read_GEOSldas import read_tilecoord  # type: ignore

    experiments = {
        "OL": {
            "exp_name": str(args.ol_exp_name),
            "run_root": Path(args.ol_run_root),
        },
        "DA": {
            "exp_name": str(args.da_exp_name),
            "run_root": Path(args.da_run_root),
        },
    }
    exp_keys = ["OL", "DA"]
    exp_index = {k: i for i, k in enumerate(exp_keys)}

    ims_regrid_dir = Path(args.ims_regrid_dir)

    if args.output_dir is None:
        output_dir = repo_root / "projects" / "IMS" / "output"
    else:
        output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_tag = args.output_tag
    if out_tag is None:
        out_tag = f"{args.domain}_{args.year_start}_{args.year_end}_thr{args.scf_threshold:.2f}".replace(".", "p")
        if int(args.min_ims_snow_days) > 0:
            out_tag = f"{out_tag}_imsSnowDaysGe{int(args.min_ims_snow_days)}"

    daily_parquet = output_dir / f"ims_ol_da_daily_counts_{out_tag}.parquet"
    daily_csv = output_dir / f"ims_ol_da_daily_counts_{out_tag}.csv"
    pair_daily_parquet = output_dir / f"ims_ol_da_pair_daily_{out_tag}.parquet"
    pair_daily_csv = output_dir / f"ims_ol_da_pair_daily_{out_tag}.csv"
    comparison_parquet = output_dir / f"ims_ol_da_comparison_table_{out_tag}.parquet"
    comparison_csv = output_dir / f"ims_ol_da_comparison_table_{out_tag}.csv"
    scope_meta_csv = output_dir / f"ims_ol_da_scope_metadata_{out_tag}.csv"
    cell_counts_nc = output_dir / f"ims_ol_da_cell_counts_metrics_{out_tag}.nc4"

    out_targets = [daily_csv, pair_daily_csv, comparison_csv, scope_meta_csv, cell_counts_nc]
    if (not args.overwrite) and any(p.exists() for p in out_targets):
        existing = [str(p) for p in out_targets if p.exists()]
        raise FileExistsError(
            "Output exists. Use --overwrite or choose --output-tag/--output-dir. Existing: " + ", ".join(existing)
        )

    # Fast postprocess mode: start from an existing cell-counts NetCDF, apply eligibility filter,
    # then regenerate downstream files used by the notebook pipeline.
    if args.reuse_cell_counts_nc is not None:
        src_nc = Path(args.reuse_cell_counts_nc)
        if not src_nc.exists():
            raise FileNotFoundError(f"--reuse-cell-counts-nc does not exist: {src_nc}")

        print(f"Fast mode: reuse existing cell-counts NetCDF: {src_nc}")
        print(f"Min IMS snow days: {int(args.min_ims_snow_days)}")

        cell_counts, rep, src_attrs = read_counts_and_rep_from_cell_counts_nc(src_nc)
        if args.reuse_scope_meta_csv:
            scope_meta = pd.read_csv(Path(args.reuse_scope_meta_csv))
        else:
            scope_meta = read_scope_meta_from_cell_counts_nc(src_nc)

        eligible = apply_min_ims_snow_days_filter_to_counts(
            cell_counts=cell_counts,
            scope_meta=scope_meta,
            exp_keys=exp_keys,
            min_days=int(args.min_ims_snow_days),
            eligibility_experiment="OL",
        )
        n_keep = int(np.sum(eligible))
        n_total = int(eligible.size)
        print(f"Eligible cells kept: {n_keep}/{n_total}")

        # Regenerate comparison table from filtered per-scope counts.
        comparison_df = build_comparison_table_from_cell_counts(
            cell_counts=cell_counts,
            scope_meta=scope_meta,
            exp_keys=exp_keys,
        )
        maybe_write_parquet(comparison_df, comparison_parquet, "comparison table")
        comparison_df.to_csv(comparison_csv, index=False)
        print(f"Wrote comparison table csv: {comparison_csv}")

        # Keep table file contract for downstream notebook even though fast mode does not have per-day rows.
        daily_df = build_empty_daily_table()
        pair_daily = build_empty_pair_daily_table()
        maybe_write_parquet(daily_df, daily_parquet, "daily counts (empty fast-mode placeholder)")
        daily_df.to_csv(daily_csv, index=False)
        print(f"Wrote daily counts csv (empty fast-mode placeholder): {daily_csv}")
        maybe_write_parquet(pair_daily, pair_daily_parquet, "paired daily (empty fast-mode placeholder)")
        pair_daily.to_csv(pair_daily_csv, index=False)
        print(f"Wrote paired daily csv (empty fast-mode placeholder): {pair_daily_csv}")

        scope_meta.to_csv(scope_meta_csv, index=False)
        print(f"Wrote scope metadata csv: {scope_meta_csv}")

        src_threshold = float(src_attrs["model_threshold_scf"])
        src_ims_var = str(src_attrs["ims_var"])
        src_snow_codes = parse_int_set(src_attrs["ims_snow_codes"])
        src_no_snow_codes = parse_int_set(src_attrs["ims_no_snow_codes"])
        src_fill_values = parse_int_set(src_attrs["ims_fill_values"])

        write_cell_counts_netcdf(
            out_nc=cell_counts_nc,
            rep=rep,
            exp_keys=exp_keys,
            scope_meta=scope_meta,
            cell_counts=cell_counts,
            threshold=src_threshold,
            ims_var=src_ims_var,
            snow_codes=src_snow_codes,
            no_snow_codes=src_no_snow_codes,
            fill_values=src_fill_values,
            eligible_mask=eligible,
            extra_note=f"Fast postprocess mode with min_ims_snow_days={int(args.min_ims_snow_days)}.",
        )
        print(f"Wrote map-ready per-cell counts+metrics NetCDF: {cell_counts_nc}")

        print("\nDone (fast mode).")
        print("Bootstrap CI columns are NaN in fast mode because per-day paired rows are not re-read.")
        print(f"Comparison rows: {len(comparison_df)}")
        print(f"Cell dimensions: experiment=2, scope={int(scope_meta.shape[0])}, cell={int(rep['rep_idx'].size)}")
        return

    # Tilecoord from OL is sufficient because OL/DA share the same M36 grid.
    tilecoord = locate_tilecoord_file(
        run_root=Path(experiments["OL"]["run_root"]),
        exp_name=str(experiments["OL"]["exp_name"]),
        domain=str(args.domain),
        explicit=Path(args.tilecoord) if args.tilecoord else None,
    )

    print(f"Repo root: {repo_root}")
    print(f"Domain: {args.domain}")
    print(f"Years: {args.year_start}..{args.year_end}")
    print(f"IMS dir/template: {ims_regrid_dir} / {args.ims_regrid_template}")
    print(f"IMS var: {args.ims_var}")
    print(f"Model SCF var candidates: {model_var_candidates}")
    print(f"Threshold: {args.scf_threshold}")
    print(f"Tilecoord: {tilecoord}")
    for k in exp_keys:
        print(f"{k}: exp_name={experiments[k]['exp_name']} run_root={experiments[k]['run_root']}")

    tc = read_tilecoord(str(tilecoord))
    rep = choose_representative_tile_per_cell(tc)

    ny = int(rep["ny"])
    nx = int(rep["nx"])
    n_cell = int(rep["rep_idx"].size)
    rep_i = rep["rep_i"].astype(np.int64)
    rep_j = rep["rep_j"].astype(np.int64)

    print(f"Representative land cells: {n_cell}")
    print(f"Grid shape from tilecoord: ny={ny}, nx={nx}")

    scope_meta, scope_maps = build_scope_metadata(args.year_start, args.year_end)
    n_scope = int(scope_meta.shape[0])
    print(f"Scopes configured: {n_scope}")

    cell_counts = init_cell_count_arrays(n_exp=len(exp_keys), n_scope=n_scope, n_cell=n_cell)

    # Cache detected model variable names to avoid repeated probing once resolved.
    exp_model_var: dict[str, str | None] = {k: None for k in exp_keys}

    records: list[dict] = []

    for year in range(args.year_start, args.year_end + 1):
        ims_path = ims_regrid_dir / args.ims_regrid_template.format(year=year)
        if not ims_path.exists():
            warnings.warn(f"Missing IMS regridded file for year {year}: {ims_path}")
            continue

        print(f"Year {year}: reading IMS {ims_path}")

        with Dataset(ims_path, "r") as ds_ims:
            if args.ims_var not in ds_ims.variables:
                warnings.warn(f"{ims_path} missing IMS var '{args.ims_var}'; skipping year")
                continue

            ims_var_obj = ds_ims.variables[args.ims_var]
            if getattr(ims_var_obj, "ndim", None) != 3:
                warnings.warn(f"{ims_path} {args.ims_var} is not 3D; skipping year")
                continue

            n_time = int(ims_var_obj.shape[0])
            spatial_shape = tuple(int(x) for x in ims_var_obj.shape[1:])
            transpose_input_spatial = False
            if spatial_shape == (ny, nx):
                transpose_input_spatial = False
            elif spatial_shape == (nx, ny):
                transpose_input_spatial = True
                print(f"  {year}: IMS spatial shape is (x,y); transposing slices to (y,x)")
            else:
                raise ValueError(
                    f"IMS spatial shape {spatial_shape} does not match rep shape {(ny, nx)} (or transpose {(nx, ny)})"
                )

            if "land_mask" in ds_ims.variables:
                land_mask = np.asarray(ds_ims.variables["land_mask"][:], dtype=np.int8) == 1
            else:
                land_mask = np.ones((ny, nx), dtype=bool)

            if "within_max_distance" in ds_ims.variables:
                within_mask = np.asarray(ds_ims.variables["within_max_distance"][:], dtype=np.int8) == 1
            else:
                within_mask = np.ones((ny, nx), dtype=bool)

            static_mask = land_mask & within_mask
            static_mask_rep = static_mask[rep_j, rep_i]

            _, ims_dates = get_ims_time_dim_and_dates(ds_ims, args.ims_var, year)

            if args.auto_infer_ims_codes:
                seen_codes: set[int] = set()
                n_sample = min(5, n_time)
                for si in range(n_sample):
                    sample = np.asarray(ims_var_obj[si, :, :], dtype=np.float32)
                    if transpose_input_spatial:
                        sample = sample.T
                    finite = np.isfinite(sample)
                    if np.any(finite):
                        sample_codes = np.unique(np.rint(sample[finite]).astype(np.int32))
                        for c in sample_codes:
                            ci = int(c)
                            if ci not in ims_fill_values:
                                seen_codes.add(ci)

                if seen_codes:
                    snow_codes, no_snow_codes = infer_ims_code_sets(seen_codes, ims_fill_values)
                else:
                    snow_codes = set(ims_snow_codes)
                    no_snow_codes = set(ims_no_snow_codes)
            else:
                snow_codes = set(ims_snow_codes)
                no_snow_codes = set(ims_no_snow_codes)

            print(
                f"  Year {year}: snow_codes={sorted(snow_codes)} no_snow_codes={sorted(no_snow_codes)}"
            )

            for ti in range(n_time):
                day = pd.Timestamp(ims_dates[ti]).normalize()
                if day.year < args.year_start or day.year > args.year_end:
                    continue

                ims_cat = np.asarray(ims_var_obj[ti, :, :], dtype=np.float32)
                if transpose_input_spatial:
                    ims_cat = ims_cat.T

                obs_bin, obs_valid_mask, unknown_mask, _ = ims_category_to_binary_scf(
                    ims_cat,
                    snow_codes=snow_codes,
                    no_snow_codes=no_snow_codes,
                    fill_values=ims_fill_values,
                )

                base_valid_mask = static_mask & obs_valid_mask

                n_land_mask = int(np.sum(static_mask))
                n_obs_valid = int(np.sum(base_valid_mask))
                n_unknown = int(np.sum(static_mask & unknown_mask))

                if (ti + 1) % 30 == 0 or (ti + 1) == n_time:
                    print(f"  IMS day {ti + 1}/{n_time}: {day.strftime('%Y-%m-%d')}")

                # Representative-cell vectors used for fast per-day processing.
                obs_rep = obs_bin[rep_j, rep_i]
                base_valid_rep = base_valid_mask[rep_j, rep_i]

                rec_by_exp: dict[str, dict] = {}
                model_rep_by_exp: dict[str, np.ndarray] = {}

                for exp_key in exp_keys:
                    cfg = experiments[exp_key]
                    run_root = Path(cfg["run_root"])
                    exp_name = str(cfg["exp_name"])

                    rec = {
                        "date": day,
                        "year": int(day.year),
                        "doy": int(day.strftime("%j")),
                        "season": season_name(day),
                        "experiment": exp_key,
                        "exp_name": exp_name,
                        "ims_file": str(ims_path),
                        "model_file": None,
                        "model_file_found": 0,
                        "model_read_ok": 0,
                        "model_var": exp_model_var.get(exp_key),
                        "A": 0,
                        "B": 0,
                        "C": 0,
                        "D": 0,
                        "N_valid": 0,
                        "accuracy": np.nan,
                        "hit_rate": np.nan,
                        "miss_rate": np.nan,
                        "false_alarm_ratio": np.nan,
                        "correct_rejection_rate": np.nan,
                        "N_land_mask": n_land_mask,
                        "N_ims_obs_valid": n_obs_valid,
                        "N_ims_unknown_codes": n_unknown,
                        "obs_scf_mean": np.nan,
                        "mod_scf_mean": np.nan,
                        "paired_common_mask_used": 0,
                        "N_common_mask": 0,
                        "mask_type": "none",
                    }

                    model_path = locate_daily_cat_file(run_root, exp_name, str(args.domain), day)
                    if model_path is None:
                        rec_by_exp[exp_key] = rec
                        continue

                    rec["model_file"] = str(model_path)
                    rec["model_file_found"] = 1

                    try:
                        rep_vals, used_var = read_model_scf_rep_values(
                            model_path,
                            rep_idx=rep["rep_idx"],
                            candidates=model_var_candidates,
                            forced_var_name=exp_model_var.get(exp_key),
                        )
                        rec["model_var"] = used_var
                        rec["model_read_ok"] = 1
                        if exp_model_var.get(exp_key) is None:
                            exp_model_var[exp_key] = used_var
                        model_rep_by_exp[exp_key] = rep_vals
                    except Exception as exc:
                        warnings.warn(f"{exp_key} {day.strftime('%Y-%m-%d')} read failed: {exc}")

                    rec_by_exp[exp_key] = rec

                both_ok = all(rec_by_exp[e]["model_read_ok"] == 1 for e in exp_keys)

                if both_ok:
                    common_mask = base_valid_rep.copy()
                    common_mask &= np.isfinite(model_rep_by_exp["OL"])
                    common_mask &= np.isfinite(model_rep_by_exp["DA"])
                    n_common = int(np.sum(common_mask))

                    season = season_name(day)
                    scope_ids = [
                        int(scope_maps["all"]["ALL_PERIOD"]),
                        int(scope_maps["season"][season]),
                        int(scope_maps["year"][int(day.year)]),
                        int(scope_maps["year_season"][(int(day.year), season)]),
                    ]

                    for exp_key in exp_keys:
                        A, B, C, D, valid = contingency_masks_from_vectors(
                            obs_bin_rep=obs_rep,
                            mod_scf_rep=model_rep_by_exp[exp_key],
                            valid_mask=common_mask,
                            threshold=float(args.scf_threshold),
                        )
                        counts = counts_from_masks(A, B, C, D)
                        m = snow_scores_from_counts(counts["A"], counts["B"], counts["C"], counts["D"])

                        rec = rec_by_exp[exp_key]
                        rec["A"] = int(counts["A"])
                        rec["B"] = int(counts["B"])
                        rec["C"] = int(counts["C"])
                        rec["D"] = int(counts["D"])
                        rec["N_valid"] = int(counts["N"])
                        rec["accuracy"] = m["accuracy"]
                        rec["hit_rate"] = m["hit_rate"]
                        rec["miss_rate"] = m["miss_rate"]
                        rec["false_alarm_ratio"] = m["false_alarm_ratio"]
                        rec["correct_rejection_rate"] = m["correct_rejection_rate"]
                        rec["paired_common_mask_used"] = 1
                        rec["N_common_mask"] = n_common
                        rec["mask_type"] = "common_pair"

                        if np.any(valid):
                            rec["obs_scf_mean"] = float(np.nanmean(obs_rep[valid]))
                            rec["mod_scf_mean"] = float(np.nanmean(model_rep_by_exp[exp_key][valid]))

                        rec_by_exp[exp_key] = rec

                        # Update map-ready per-cell accumulators for this experiment/day.
                        update_cell_counts_for_day(
                            cell_counts=cell_counts,
                            exp_idx=exp_index[exp_key],
                            scope_ids=scope_ids,
                            A=A,
                            B=B,
                            C=C,
                            D=D,
                            N_valid=valid,
                        )
                else:
                    # Keep optional individual diagnostics in daily table for troubleshooting.
                    for exp_key in exp_keys:
                        rec = rec_by_exp[exp_key]
                        if rec["model_read_ok"] == 1:
                            indiv_mask = base_valid_rep & np.isfinite(model_rep_by_exp[exp_key])
                            A, B, C, D, valid = contingency_masks_from_vectors(
                                obs_bin_rep=obs_rep,
                                mod_scf_rep=model_rep_by_exp[exp_key],
                                valid_mask=indiv_mask,
                                threshold=float(args.scf_threshold),
                            )
                            counts = counts_from_masks(A, B, C, D)
                            m = snow_scores_from_counts(counts["A"], counts["B"], counts["C"], counts["D"])

                            rec["A"] = int(counts["A"])
                            rec["B"] = int(counts["B"])
                            rec["C"] = int(counts["C"])
                            rec["D"] = int(counts["D"])
                            rec["N_valid"] = int(counts["N"])
                            rec["accuracy"] = m["accuracy"]
                            rec["hit_rate"] = m["hit_rate"]
                            rec["miss_rate"] = m["miss_rate"]
                            rec["false_alarm_ratio"] = m["false_alarm_ratio"]
                            rec["correct_rejection_rate"] = m["correct_rejection_rate"]
                            rec["mask_type"] = "individual"

                            if np.any(valid):
                                rec["obs_scf_mean"] = float(np.nanmean(obs_rep[valid]))
                                rec["mod_scf_mean"] = float(np.nanmean(model_rep_by_exp[exp_key][valid]))

                            rec_by_exp[exp_key] = rec

                for e in exp_keys:
                    records.append(rec_by_exp[e])

    daily_df = pd.DataFrame.from_records(records)
    if daily_df.empty:
        raise RuntimeError("No daily records generated. Check paths and year range.")

    daily_df["date"] = pd.to_datetime(daily_df["date"])
    daily_df = daily_df.sort_values(["date", "experiment"]).reset_index(drop=True)

    maybe_write_parquet(daily_df, daily_parquet, "daily counts")
    daily_df.to_csv(daily_csv, index=False)
    print(f"Wrote daily counts csv: {daily_csv}")

    paired_rows = daily_df[
        (daily_df["paired_common_mask_used"] == 1)
        & (daily_df["model_read_ok"] == 1)
        & (daily_df["N_valid"] > 0)
    ].copy()

    if paired_rows.empty:
        raise RuntimeError("No paired common-mask rows available for OL-vs-DA comparison")

    pair_daily = (
        paired_rows.pivot_table(
            index=["date", "year", "doy", "season"],
            columns="experiment",
            values=[
                "A",
                "B",
                "C",
                "D",
                "N_valid",
                "obs_scf_mean",
                "mod_scf_mean",
                "model_file_found",
                "model_read_ok",
                "paired_common_mask_used",
            ],
            aggfunc="first",
        )
    )

    pair_daily.columns = [f"{v}_{e}" for v, e in pair_daily.columns]
    pair_daily = pair_daily.reset_index()

    need_cols = [
        "A_OL",
        "B_OL",
        "C_OL",
        "D_OL",
        "N_valid_OL",
        "A_DA",
        "B_DA",
        "C_DA",
        "D_DA",
        "N_valid_DA",
    ]
    for c in need_cols:
        if c not in pair_daily.columns:
            pair_daily[c] = np.nan

    valid_pair_mask = (pair_daily["N_valid_OL"].fillna(0) > 0) & (pair_daily["N_valid_DA"].fillna(0) > 0)
    pair_daily = pair_daily[valid_pair_mask].copy()
    if pair_daily.empty:
        raise RuntimeError("No paired valid days left after OL/DA pivot/filter")

    pair_daily = pair_daily.sort_values(["date"]).reset_index(drop=True)
    maybe_write_parquet(pair_daily, pair_daily_parquet, "paired daily")
    pair_daily.to_csv(pair_daily_csv, index=False)
    print(f"Wrote paired daily csv: {pair_daily_csv}")

    comparison_parts = []

    def add_scope(df_scope: pd.DataFrame, scope_name: str, year_val=np.nan, season_val=np.nan, seed_offset=0):
        if df_scope.empty:
            return
        tbl = bootstrap_compare_table_from_pair_days(
            df_scope,
            exp_ol="OL",
            exp_da="DA",
            n_boot=int(args.n_bootstrap),
            seed=int(args.bootstrap_seed + seed_offset),
            ci_low=float(args.ci_low),
            ci_high=float(args.ci_high),
        )
        if tbl.empty:
            return
        tbl["scope"] = scope_name
        tbl["year"] = year_val
        tbl["season"] = season_val
        comparison_parts.append(tbl)

    add_scope(pair_daily, "ALL_PERIOD", seed_offset=0)

    for s in SEASON_ORDER:
        add_scope(
            pair_daily[pair_daily["season"] == s],
            "SEASON_ALL_YEARS",
            season_val=s,
            seed_offset=10 + SEASON_ORDER.index(s),
        )

    for y in sorted(pair_daily["year"].unique()):
        yi = int(y)
        add_scope(pair_daily[pair_daily["year"] == yi], "YEAR", year_val=yi, seed_offset=1000 + yi)

    for y in sorted(pair_daily["year"].unique()):
        yi = int(y)
        for s in SEASON_ORDER:
            sub = pair_daily[(pair_daily["year"] == yi) & (pair_daily["season"] == s)]
            if sub.empty:
                continue
            add_scope(
                sub,
                "YEAR_SEASON",
                year_val=yi,
                season_val=s,
                seed_offset=2000 + yi * 10 + SEASON_ORDER.index(s),
            )

    if not comparison_parts:
        raise RuntimeError("No comparison rows produced")

    comparison_df = pd.concat(comparison_parts, ignore_index=True, sort=False)
    comparison_df = comparison_df[
        [
            "scope",
            "year",
            "season",
            "metric",
            "ol",
            "ol_ci_lo",
            "ol_ci_hi",
            "da",
            "da_ci_lo",
            "da_ci_hi",
            "delta_da_minus_ol",
            "delta_ci_lo",
            "delta_ci_hi",
            "n_days",
            "A_ol",
            "B_ol",
            "C_ol",
            "D_ol",
            "N_ol",
            "A_da",
            "B_da",
            "C_da",
            "D_da",
            "N_da",
        ]
    ].copy()

    maybe_write_parquet(comparison_df, comparison_parquet, "comparison table")
    comparison_df.to_csv(comparison_csv, index=False)
    print(f"Wrote comparison table csv: {comparison_csv}")

    scope_meta.to_csv(scope_meta_csv, index=False)
    print(f"Wrote scope metadata csv: {scope_meta_csv}")

    write_cell_counts_netcdf(
        out_nc=cell_counts_nc,
        rep=rep,
        exp_keys=exp_keys,
        scope_meta=scope_meta,
        cell_counts=cell_counts,
        threshold=float(args.scf_threshold),
        ims_var=str(args.ims_var),
        snow_codes=snow_codes,
        no_snow_codes=no_snow_codes,
        fill_values=ims_fill_values,
    )
    print(f"Wrote map-ready per-cell counts+metrics NetCDF: {cell_counts_nc}")

    print("\nDone.")
    print(f"Daily rows: {len(daily_df)}")
    print(f"Paired daily rows: {len(pair_daily)}")
    print(f"Comparison rows: {len(comparison_df)}")
    print(f"Cell dimensions: experiment=2, scope={n_scope}, cell={n_cell}")


if __name__ == "__main__":
    main()
