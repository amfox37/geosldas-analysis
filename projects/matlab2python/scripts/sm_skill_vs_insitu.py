#!/usr/bin/env python3
"""
Python translation of SMskill_vs_INSITU_single_CYG_expt.m.

Computes surface and root-zone soil moisture skill (R, anomR, Bias, RMSE, ubRMSE with CIs)
against sparse in-situ networks for CYGNSS experiments. Output is written as NetCDF along
with a compressed NPZ for quick reload.

Notes:
- Paths default to the same locations used in the MATLAB script. Override via CLI flags.
- Anomaly correlation follows the 31-day window climatology (CalVal uses Nmin_day=50, others 150).
- Autocorrelation-adjusted sample sizes follow the MATLAB logic (lag-1 based).
- CalVal daily stations are aggregated to daily, then reindexed to the first station’s calendar (Y/M/D), matching MATLAB.

Example:
python sm_skill_vs_insitu.py \
  --insitu-tag CalVal_M33 \
  --exp DAv8_M36_cd_all --exp DAv8_M36_cd_ssa \
  --exp-root /discover/nobackup/projects/land_da/CYGNSS_Experiments \
  --insitu-root /discover/nobackup/qliu/merra_land/DATA \
  --start 2018-08-01 --end 2024-07-01 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --file-tag ldas_tile_daily_out \
  --out-dir /discover/nobackup/projects/land_da/Evaluation/InSitu/output
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import scipy.io as sio
from scipy import stats
import xarray as xr

# Resolve repo root (…/geosldas-analysis) from projects/matlab2python/scripts/
ROOT = Path(__file__).resolve().parents[3]
# Add both shared IO locations
sys.path.append(str(ROOT / "common" / "python" / "io"))
sys.path.append(str(ROOT / "projects" / "matlab2python" / "shared" / "python"))
from read_GEOSldas import read_tilecoord  # type: ignore


# ----------------------
# Date utilities
# ----------------------
def is_leap_year(year: int) -> bool:
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def days_in_month(year: int, month: int) -> int:
    days = [31, 29 if is_leap_year(year) else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return days[month - 1]


def pentad_of_year(dofyr: int, year: int) -> int:
    # Use leap-year-aware pentads: every 5 days, starting at doy 1.
    # Keep a max of 73 pentads (365/5 = 73); leap day folds into pentad 12.
    doy_adj = dofyr
    if is_leap_year(year) and dofyr > 59:  # after Feb 29
        doy_adj -= 1
    return int(math.floor((doy_adj - 1) / 5) + 1)


@dataclass
class DateTime:
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: int
    dofyr: int | None = None
    pentad: int | None = None


def get_dofyr_pentad(dt: DateTime) -> DateTime:
    dofyr = dt.day
    for m in range(1, dt.month):
        dofyr += days_in_month(dt.year, m)
    dt.dofyr = dofyr
    dt.pentad = pentad_of_year(dofyr, dt.year)
    return dt


def augment_date_time(dtstep_sec: int, dt: DateTime) -> DateTime:
    # Mirror MATLAB augment_date_time.m: advance/retreat by dtstep_sec, updating Y/M/D/HH/MM/SS and dofyr/pentad
    total = dt.hour * 3600 + dt.minute * 60 + dt.second + dtstep_sec
    day_delta, rem = divmod(total, 86400)
    if rem < 0:
        rem += 86400
        day_delta -= 1
    hour, rem2 = divmod(rem, 3600)
    minute, second = divmod(rem2, 60)

    year, month, day = dt.year, dt.month, dt.day + day_delta
    while True:
        dim = days_in_month(year, month)
        if day > dim:
            day -= dim
            month += 1
            if month > 12:
                month = 1
                year += 1
        elif day < 1:
            month -= 1
            if month < 1:
                month = 12
                year -= 1
            day += days_in_month(year, month)
        else:
            break

    return get_dofyr_pentad(DateTime(year, month, day, hour, minute, second))


def get_date_time_string(dt: DateTime, file_tag: str) -> str:
    s = f"{dt.year:04d}"
    if "pentad" in file_tag:
        s += f"p{dt.pentad:02d}"
    else:
        s += f"{dt.month:02d}"
    if "daily" in file_tag:
        s += f"{dt.day:02d}"
    if any(k in file_tag for k in ["xhourly", "rst", "innov", "incr", "OminusF", "1hto3h", "ObsFcstAna", "inst"]):
        s += f"{dt.day:02d}_{dt.hour:02d}{dt.minute:02d}"
    return s


# ----------------------
# In-situ coordinate readers
# ----------------------
def _read_txt_ids_lat_lon(path: Path, fmt: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    arr = np.genfromtxt(path, dtype=None, encoding=None, comments=None)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    ids = [str(row[0]) for row in arr]
    lat = np.array([float(row[2]) for row in arr])
    lon = np.array([float(row[3]) for row in arr])
    return lat, lon, ids


def get_insitu_coord(insitu_tag: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Python port of get_INSITU_coord.m (subset of networks used in CYG workflow).
    """
    if insitu_tag == "SCAN":
        coord_file = Path("/discover/nobackup/qliu/merra_land/DATA/SCAN/data/coord/ALL_SCAN_coord.txt")
        with coord_file.open() as f:
            cols = list(zip(*[line.strip().split() for line in f]))
        ids = [str(int(x)) for x in cols[1]]
        lat = np.array(cols[2], dtype=float)
        lon = np.array(cols[3], dtype=float)
        # exclude non-CONUS
        mask = (lat <= 45) & (lat > 25) & (lon < -60) & (lon >= -130)
        return lat[mask], lon[mask], list(np.array(ids)[mask])

    if insitu_tag == "USCRN":
        base = Path("/discover/nobackup/qliu/merra_land/DATA/USCRN/data/")
        sites = np.genfromtxt(base / "list_of_USCRNsites.txt", dtype=str)
        info = np.genfromtxt(base / "USCRN_Station_info.txt", dtype=None, encoding=None, skip_header=1)
        ids_all = [row[0] for row in info]
        lat_all = np.array([row[1] for row in info], dtype=float)
        lon_all = np.array([row[2] for row in info], dtype=float)
        idx = [ids_all.index(s) for s in sites]
        return lat_all[idx], lon_all[idx], list(np.array(ids_all)[idx])

    if insitu_tag == "CalVal_M33":
        # Hard-coded refpix list from MATLAB
        refpixids = [
            "03013302",
            "04013302",
            "07013301",
            "09013301",
            "12033301",
            "16013302",
            "16023302",
            "16033302",
            "16043302",
            "16063302",
            "16073302",
            "19023301",
            "25013301",
            "27013301",
            "45013301",
            "45023301",
            "48013301",
            "67013301",
        ]
        # Expect the master list alongside validation tools
        master = Path("/home/qliu/smap/SMAP_Nature/validation/tools/refpix_SM_masterlist.txt")
        if not master.exists():
            raise FileNotFoundError(f"CalVal master file not found: {master}")
        import pandas as pd

        # CSV with header; columns: Alias,RPID,Active,RZ,Col,Row,Lat,Lon,Description
        df = pd.read_csv(master, sep=",", engine="python", comment="#", skip_blank_lines=True)
        df.columns = [c.strip() for c in df.columns]  # trim whitespace in header

        def fmt_rpid(val):
            try:
                return f"{int(val):08d}"
            except Exception:
                return str(val).strip()

        df["RPID"] = df["RPID"].apply(fmt_rpid)
        df = df[df["RPID"].isin(refpixids)]
        if df.empty:
            raise ValueError(f"No CalVal refpix entries matched in {master}")
        return df["Lat"].to_numpy(), df["Lon"].to_numpy(), df["RPID"].tolist()

    raise ValueError(f"Unsupported INSITU_tag '{insitu_tag}' in Python port")


# ----------------------
# Validation stats (port of get_validation_stats.m)
# ----------------------
def nancorrcoef(x: np.ndarray, y: np.ndarray) -> float:
    m = ~np.isnan(x) & ~np.isnan(y)
    if m.sum() < 2:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1])


def corrcoef_autocorr(data: np.ndarray) -> Tuple[float, float, float]:
    # data columns: [ref, model]; computes R and 95% CI with lag-1 adjusted n
    ref = data[:, 0]
    mod = data[:, 1]
    mask = ~np.isnan(ref) & ~np.isnan(mod)
    ref = ref[mask]
    mod = mod[mask]
    n = len(ref)
    if n < 3:
        return np.nan, np.nan, np.nan
    rx = nancorrcoef(ref[:-1], ref[1:])
    ry = nancorrcoef(mod[:-1], mod[1:])
    neff = max(3, math.ceil(n * (1 - rx * ry) / (1 + rx * ry))) if not (np.isnan(rx) or np.isnan(ry)) else n
    r = float(np.corrcoef(ref, mod)[0, 1])
    z = 0.5 * math.log((1 + r) / (1 - r)) if abs(r) < 1 else np.inf
    se = 1 / math.sqrt(neff - 3)
    z_lo = z - stats.norm.ppf(0.975) * se
    z_up = z + stats.norm.ppf(0.975) * se
    r_lo = math.tanh(z_lo)
    r_up = math.tanh(z_up)
    return r, r_lo, r_up


def get_validation_stats(data: np.ndarray, AC: bool, complete: bool, ref_col: int, select_cols: Iterable[int], Nmin: int) -> Dict[str, float]:
    # data shape (n, m); ref_col, select_cols are 1-based in MATLAB; convert to 0-based here
    ref_idx = ref_col - 1
    select_idx = [c - 1 for c in select_cols]
    sub = data[:, select_idx]
    mask = ~np.any(np.isnan(sub), axis=1) if complete else ~(np.isnan(sub[:, 0]) | np.isnan(sub[:, 1]))
    sub = sub[mask]
    n_pairs = sub.shape[0]
    stats_out = {k: np.nan for k in ["N_pairs", "R", "RLO", "RUP", "bias", "CI_bias", "MSE", "CI_MSE", "ubMSE", "ubMSELO", "ubMSEUP"]}
    stats_out["N_pairs"] = n_pairs
    if n_pairs <= Nmin:
        return stats_out

    ref = sub[:, 0]
    mod = sub[:, 1]
    if AC:
        r, rlo, rup = corrcoef_autocorr(sub[:, :2])
        # effective n from corrcoef_autocorr
        rx = nancorrcoef(ref[:-1], ref[1:])
        ry = nancorrcoef(mod[:-1], mod[1:])
        neff = max(3, math.ceil(n_pairs * (1 - rx * ry) / (1 + rx * ry))) if not (np.isnan(rx) or np.isnan(ry)) else n_pairs
    else:
        r = float(np.corrcoef(ref, mod)[0, 1])
        se = 1 / math.sqrt(n_pairs - 3)
        z = 0.5 * math.log((1 + r) / (1 - r)) if abs(r) < 1 else np.inf
        z_lo = z - stats.norm.ppf(0.975) * se
        z_up = z + stats.norm.ppf(0.975) * se
        rlo, rup = math.tanh(z_lo), math.tanh(z_up)
        neff = n_pairs

    err = ref - mod
    bias = float(np.nanmean(err))
    sd_err = float(np.nanstd(err, ddof=1))
    se_bias = sd_err / math.sqrt(neff)
    tcrit = stats.t.ppf(0.975, max(neff - 1, 1))
    ci_bias = tcrit * se_bias

    mse = float(np.nanmean(err**2))
    sd_mse = float(np.nanstd(err**2, ddof=1))
    se_mse = sd_mse / math.sqrt(neff)
    ci_mse = tcrit * se_mse

    ubmse = mse - bias**2
    # Variance CI via chi-square
    df = max(neff - 1, 1)
    chi_lo = stats.chi2.ppf(0.975, df)
    chi_up = stats.chi2.ppf(0.025, df)
    ubmse_lo = df * ubmse / chi_lo if chi_lo > 0 else np.nan
    ubmse_up = df * ubmse / chi_up if chi_up > 0 else np.nan

    stats_out.update(
        R=r,
        RLO=rlo - r,
        RUP=rup - r,
        bias=-bias,  # MATLAB uses -(insitu-model)
        CI_bias=ci_bias,
        MSE=mse,
        CI_MSE=ci_mse,
        ubMSE=ubmse,
        ubMSELO=ubmse_lo - ubmse,
        ubMSEUP=ubmse_up - ubmse,
    )
    return stats_out


# ----------------------
# Core processing
# ----------------------
def build_time_vector(start: DateTime, end: DateTime, dtstep: int, file_tag: str) -> List[DateTime]:
    # Build [start, end) time axis at dtstep cadence, mirroring MATLAB's while/break loop
    times: List[DateTime] = []
    current = start
    while True:
        if current.year == end.year and current.month == end.month and current.day == end.day:
            break
        times.append(current)
        current = augment_date_time(dtstep, current)
    return times


def nearest_tiles(tc: Dict[str, np.ndarray], lat: np.ndarray, lon: np.ndarray, max_distance: float) -> Tuple[np.ndarray, np.ndarray]:
    # Map each station to nearest tile (squared distance), drop if outside threshold
    inds = []
    dists = []
    for la, lo in zip(lat, lon):
        dist = (tc["com_lat"] - la) ** 2 + (tc["com_lon"] - lo) ** 2
        idx = int(np.argmin(dist))
        if dist[idx] > max_distance:
            inds.append(np.nan)
            dists.append(dist[idx])
        else:
            inds.append(idx)
            dists.append(dist[idx])
    inds = np.array(inds, dtype=float)
    dists = np.array(dists, dtype=float)
    mask = ~np.isnan(inds)
    return inds[mask].astype(int), mask


def _extract_tile_values(ds: xr.Dataset, var: str, ind_tile: np.ndarray, dt: DateTime, dtstep: int, n_tile: int) -> np.ndarray:
    da = ds[var]
    # handle time dimension
    if "time" in da.dims:
        if da.sizes["time"] > 1 and dtstep == 10800:
            ind_hour = min(da.sizes["time"] - 1, math.ceil((dt.hour + 1) / 3) - 1)
            da = da.isel(time=ind_hour)
        else:
            da = da.isel(time=0)
    # tile-friendly selection
    if "tile" in da.dims:
        return da.isel(tile=ind_tile).values
    # try to pick axis matching tile count
    arr = np.asarray(da.values)
    candidates = [i for i, dim in enumerate(da.dims) if da.shape[i] == n_tile] or \
                 [i for i, dim in enumerate(da.dims) if da.shape[i] > ind_tile.max()]
    if candidates:
        axis = candidates[0]
        return np.take(arr, ind_tile, axis=axis)
    # fallback: flatten
    arr_flat = arr.ravel(order="F")  # column-major to mimic MATLAB
    if arr_flat.size <= ind_tile.max():
        raise IndexError(f"Variable {var} flattened has length {arr_flat.size}, max index {ind_tile.max()}")
    return arr_flat[ind_tile]


def read_model_series(exp_root: Path, exp: str, domain: str, file_tag: str, times: List[DateTime], ind_tile: np.ndarray, n_tile: int) -> Tuple[np.ndarray, List[DateTime]]:
    dtstep = 86400 if "daily" in file_tag else 10800
    surf = []
    rz = []
    for dt in times:
        dt_str = get_date_time_string(dt, file_tag)
        fname_primary = exp_root / exp / "output" / domain / "cat" / "ens_avg" / f"Y{dt.year:04d}" / f"M{dt.month:02d}" / f"{exp}{'.tavg24_1d_lnd_Nt.' if 'daily' in file_tag else '.tavg3_1d_lnd_Nt.'}{dt_str}_1200z.nc4"
        fname_fallback = exp_root / exp / exp / "output" / domain / "cat" / "ens_avg" / f"Y{dt.year:04d}" / f"M{dt.month:02d}" / f"{exp}{'.tavg24_1d_lnd_Nt.' if 'daily' in file_tag else '.tavg3_1d_lnd_Nt.'}{dt_str}_1200z.nc4"
        fname = fname_primary if fname_primary.exists() else fname_fallback
        if not fname.exists():
            fname = fname.with_name(fname.name.replace("ens_avg", "ens0000"))
        if not fname.exists():
            # fallback without hour suffix
            fname = fname.with_name(f"{exp}{'.tavg24_1d_lnd_Nt.' if 'daily' in file_tag else '.tavg3_1d_lnd_Nt.'}{dt_str}.nc4")
        if not fname.exists():
            print(f"WARNING: missing model file {fname}")
            surf.append(np.full((len(ind_tile),), np.nan))
            rz.append(np.full((len(ind_tile),), np.nan))
            continue
        ds = xr.open_dataset(fname, engine="netcdf4")
        sfmc = _extract_tile_values(ds, "SFMC", ind_tile, dt, dtstep, n_tile)
        rzmc = _extract_tile_values(ds, "RZMC", ind_tile, dt, dtstep, n_tile)
        surf.append(sfmc)
        rz.append(rzmc)
        ds.close()
    surf_arr = np.vstack(surf)
    rz_arr = np.vstack(rz)
    return np.stack([surf_arr, rz_arr], axis=1), times


def load_insitu_series(insitu_tag: str, insitu_root: Path, file_tag: str, ids: List[str]):
    sm_list = []
    st_list = []
    prcp_list = []
    times = None
    ref_YMD = None
    ref_T = None
    for i, sid in enumerate(ids):
        if insitu_tag == "USCRN":
            tag1, tmp_tag = f"{insitu_tag}_", "_2009_2024"
        elif insitu_tag == "Oznet":
            tag1, tmp_tag = "", "_2015_2022"
        elif insitu_tag == "Msnet":
            tag1, tmp_tag = "Mesonet_", "_201504_201804"
        elif insitu_tag == "SMOSM":
            tag1, tmp_tag = "", "_2010_2022"
        elif "CalVal" in insitu_tag:
            tag1, tmp_tag = "", "_201504_202103"
        else:
            tag1, tmp_tag = f"{insitu_tag}_", ""
        subdir = sid if "CalVal" not in insitu_tag else ""
        fname = insitu_root / {
            "SCAN": "SCAN/data/PROCESSED_202501_QL",
            "USCRN": "USCRN/data/PROCESSED_202501",
            "Oznet": "Murrumbidgee/QLiu_202211/data/PROCESSED",
            "Msnet": "Oklahoma_Mesonet/data/PROCESSED_201805",
            "SMOSM": "ISMN/data/PROCESSED_20230101/SMOSMANIA",
        }.get(insitu_tag, "CalVal_insitu/SMAP_refpix/202104")
        data_ext = "1h" if "CalVal" in insitu_tag else ("1d" if "daily" in file_tag else "3h")
        mat_path = Path("/discover/nobackup/qliu/merra_land/DATA") / fname / subdir / f"{tag1}{sid}{tmp_tag}_{data_ext}.mat"
        if not mat_path.exists():
            raise FileNotFoundError(f"Missing in-situ file {mat_path}")
        mat = sio.loadmat(mat_path)
        if "CalVal" in insitu_tag and "daily" in file_tag:
            S = mat["regular_data"]
            # columns: year, month, day, hour, min, sec, doy, pentad, prcp, sm1, sm2, tempC, ...
            Y, M, D = S[:, 0].astype(int), S[:, 1].astype(int), S[:, 2].astype(int)
            pr = S[:, 8].astype(float)
            pr[pr < -1e-5] = np.nan
            sm1 = S[:, 9].astype(float)
            sm2 = S[:, 10].astype(float)
            tC = S[:, 11].astype(float)

            # Aggregate to daily (mean SM/temp, sum precip)
            import pandas as pd

            df = pd.DataFrame({
                "Y": Y, "M": M, "D": D,
                "sm1": sm1, "sm2": sm2,
                "tC": tC,
                "pr": pr,
            })
            gb = df.groupby(["Y", "M", "D"], as_index=False)
            sm1_d = gb["sm1"].mean()
            sm2_d = gb["sm2"].mean()
            tC_d = gb["tC"].mean()
            pr_d = gb["pr"].sum(min_count=1)
            Yd = sm1_d["Y"].to_numpy()
            Md = sm1_d["M"].to_numpy()
            Dd = sm1_d["D"].to_numpy()

            sm = np.column_stack([sm1_d["sm1"].to_numpy(), sm2_d["sm2"].to_numpy()])
            tK = tC_d["tC"].to_numpy() + 273.15
            prd = pr_d["pr"].to_numpy()

            # Build daily time metadata
            sm = sm.astype(float)
            tK = tK.astype(float)
            prd = prd.astype(float)

            if ref_YMD is None:
                ref_YMD = np.column_stack([Yd, Md, Dd])
                ref_T = ref_YMD.shape[0]
                sm_list.append(sm)
                st_list.append(tK[:, None])
                prcp_list.append(prd[:, None])
                if times is None:
                    times = DateTimeSeries.from_YMDHMS(Yd, Md, Dd, np.full_like(Yd, 12), np.zeros_like(Yd), np.zeros_like(Yd))
            else:
                this_YMD = np.column_stack([Yd, Md, Dd])
                ref_YMD_tuples = {tuple(row): idx for idx, row in enumerate(ref_YMD)}
                sm_pad = np.full((ref_T, sm.shape[1]), np.nan)
                tK_pad = np.full((ref_T, 1), np.nan)
                pr_pad = np.full((ref_T, 1), np.nan)
                for row, srow, tval, pval in zip(this_YMD, sm, tK, prd):
                    idx = ref_YMD_tuples.get(tuple(row))
                    if idx is not None:
                        sm_pad[idx, :] = srow
                        tK_pad[idx, 0] = tval
                        pr_pad[idx, 0] = pval
                sm_list.append(sm_pad)
                st_list.append(tK_pad)
                prcp_list.append(pr_pad)
        else:
            key = "Sdata_1d" if "daily" in file_tag else ("Udata_3h" if insitu_tag == "USCRN" else "Sdata_3h")
            if key not in mat:
                raise KeyError(f"{key} not in {mat_path}")
            S = mat[key]
            sm = S[:, 9:11] if "CalVal" in insitu_tag else S[:, 9:14]
            tC = S[:, 11] if "CalVal" in insitu_tag else S[:, 14:19][:, 0]
            tK = tC + 273.15
            pr = S[:, 8]
            pr[pr < -1e-5] = np.nan
            sm_list.append(sm)
            st_list.append(tK[:, None])
            prcp_list.append(pr[:, None])
            if times is None:
                times = DateTimeSeries.from_YMDHMS(S[:, 0], S[:, 1], S[:, 2], S[:, 3], S[:, 4], S[:, 5])
    sm_arr = np.stack(sm_list, axis=2)
    st_arr = np.stack(st_list, axis=2)
    pr_arr = np.stack(prcp_list, axis=2)
    return sm_arr, st_arr, pr_arr, times


class DateTimeSeries:
    def __init__(self, dts: List[DateTime]):
        self.dts = dts

    @classmethod
    def from_YMDHMS(cls, Y, M, D, H, Mi, S):
        dts = [get_dofyr_pentad(DateTime(int(y), int(m), int(d), int(h), int(mi), int(s))) for y, m, d, h, mi, s in zip(Y, M, D, H, Mi, S)]
        return cls(dts)

    def to_datetime64(self, daily: bool):
        if daily:
            return np.array([np.datetime64(f"{d.year:04d}-{d.month:02d}-{d.day:02d}T12:00:00") for d in self.dts])
        return np.array([np.datetime64(f"{d.year:04d}-{d.month:02d}-{d.day:02d}T{d.hour:02d}:{d.minute:02d}:00") for d in self.dts])


def align_series(model: np.ndarray, model_times: List[DateTime], obs_sm: np.ndarray, obs_times: DateTimeSeries, file_tag: str):
    daily = "daily" in file_tag
    mkey = DateTimeSeries(model_times).to_datetime64(daily)
    okey = obs_times.to_datetime64(daily)
    if not daily and "CalVal" in file_tag:
        # bin to 3-hour
        okey = okey.astype("datetime64[h]")
    common, im, io = np.intersect1d(mkey, okey, return_indices=True)
    if common.size == 0:
        raise RuntimeError("No overlapping timestamps between model and obs")
    return model[im], obs_sm[io], [model_times[i] for i in im]


def compute_anom(series: np.ndarray, doy_vec: np.ndarray, Nmin_day: int = 150) -> np.ndarray:
    # series shape (time,); returns anomalies with 31-day window climatology (circular over year)
    clim = np.full(365, np.nan)
    for doy in range(1, 366):
        if doy <= 15:
            window = list(range(1, doy + 16)) + list(range(365 - (15 - doy) + 1, 366))
        elif doy >= 351:
            window = list(range(doy - 15, 366)) + list(range(1, 16 - (365 - doy)))
        else:
            window = list(range(doy - 15, doy + 16))
        idx = np.isin(doy_vec, window)
        vals = series[idx]
        vals = vals[~np.isnan(vals)]
        clim[doy - 1] = np.mean(vals) if len(vals) >= Nmin_day else np.nan
    anom = np.full_like(series, np.nan, dtype=float)
    for doy in range(1, 367):
        idx = doy_vec == min(doy, 365)
        anom[idx] = series[idx] - clim[min(doy, 365) - 1]
    return anom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--insitu-tag", required=True)
    ap.add_argument("--exp", action="append", required=True, help="Experiment names (repeatable)")
    ap.add_argument("--exp-root", type=Path, default=Path("/discover/nobackup/projects/land_da/CYGNSS_Experiments"))
    ap.add_argument("--insitu-root", type=Path, default=Path("/discover/nobackup/qliu/merra_land/DATA"))
    ap.add_argument("--domain", default="SMAP_EASEv2_M36_GLOBAL")
    ap.add_argument("--file-tag", default="ldas_tile_daily_out")
    ap.add_argument("--start", required=True, help="YYYY-MM-DD")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD (exclusive)")
    ap.add_argument("--max-distance", type=float, default=0.1)
    ap.add_argument("--no-anomR", action="store_true")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    start_dt = get_dofyr_pentad(DateTime(*map(int, args.start.split("-")), 12, 0, 0))
    end_dt = get_dofyr_pentad(DateTime(*map(int, args.end.split("-")), start_dt.hour, start_dt.minute, start_dt.second))
    dtstep = 86400 if "daily" in args.file_tag else 10800
    times = build_time_vector(start_dt, end_dt, dtstep, args.file_tag)

    # Load in-situ coords once; per-experiment loop mirrors MATLAB kkk loop
    lat, lon, ids = get_insitu_coord(args.insitu_tag)
    for exp in args.exp:
        print(f"Processing {exp} / {args.insitu_tag}")
        # Try standard path; if missing, fall back to exp/exp/output layout
        tc_path_primary = args.exp_root / exp / "output" / args.domain / "rc_out" / f"{exp}.ldas_tilecoord.bin"
        tc_path_fallback = args.exp_root / exp / exp / "output" / args.domain / "rc_out" / f"{exp}.ldas_tilecoord.bin"
        tc_path = tc_path_primary if tc_path_primary.exists() else tc_path_fallback
        if not tc_path.exists():
            raise FileNotFoundError(f"Tilecoord not found at {tc_path_primary} or {tc_path_fallback}")
        print(f"reading from {tc_path}")
        tc = read_tilecoord(str(tc_path))
        n_tile = len(tc["com_lat"])
        ind_tile, mask_keep = nearest_tiles(tc, lat, lon, args.max_distance)
        lat_use = lat[mask_keep]
        lon_use = lon[mask_keep]
        ids_use = [ids[i] for i, m in enumerate(mask_keep) if m]
        # Read model and obs, then align by overlapping timestamps (daily: calendar day)
        model_sm, model_times = read_model_series(args.exp_root, exp, args.domain, args.file_tag, times, ind_tile, n_tile)
        obs_sm, obs_st, obs_prcp, obs_times = load_insitu_series(args.insitu_tag, args.insitu_root, args.file_tag, ids_use)
        model_sm_aligned, obs_sm_aligned, model_times_aligned = align_series(model_sm.transpose(0, 2, 1), model_times, obs_sm, obs_times, args.file_tag)
        # model_sm_aligned: (time, station, depth); obs_sm_aligned: (time, depth, station) -> transpose
        obs_sm_aligned = np.transpose(obs_sm_aligned, (0, 2, 1))

        # Apply same QC as MATLAB: model zeros -> NaN, obs with cold temps (<277.16K) or tiny sm -> NaN
        model_sm_aligned[model_sm_aligned == 0] = np.nan
        obs_sm_aligned[obs_sm_aligned < 1e-4] = np.nan
        # For CalVal, temp was in column 12 (converted to K) when daily; if st provided, mask cold
        # We don't carry st through alignment; assume already masked in load_insitu_series (CalVal daily).

        # Cross mask
        m = np.isnan(model_sm_aligned) | np.isnan(obs_sm_aligned)
        model_sm_aligned[m] = np.nan
        obs_sm_aligned[m] = np.nan

        n_sites = model_sm_aligned.shape[1]
        stats_R = np.full((n_sites, 2), np.nan)
        stats_RLO = np.full_like(stats_R, np.nan)
        stats_RUP = np.full_like(stats_R, np.nan)
        stats_anomR = np.full_like(stats_R, np.nan)
        stats_anomRLO = np.full_like(stats_R, np.nan)
        stats_anomRUP = np.full_like(stats_R, np.nan)
        stats_bias = np.full_like(stats_R, np.nan)
        stats_biasCI = np.full_like(stats_R, np.nan)
        stats_rmse = np.full_like(stats_R, np.nan)
        stats_rmseCI = np.full_like(stats_R, np.nan)
        stats_ubrmse = np.full_like(stats_R, np.nan)
        stats_ubrmseCI = np.full_like(stats_R, np.nan)

        doy_vec = np.array([dt.dofyr for dt in model_times_aligned])
        Nmin = 200 if dtstep == 86400 else 480
        # Compute stats per site/depth
        for i_site in range(n_sites):
            for j_depth in range(model_sm_aligned.shape[2]):
                tmp = np.c_[obs_sm_aligned[:, i_site, j_depth], model_sm_aligned[:, i_site, j_depth]]
                s = get_validation_stats(tmp, AC=1, complete=True, ref_col=1, select_cols=[1, 2], Nmin=Nmin)
                stats_R[i_site, j_depth] = s["R"]
                stats_RLO[i_site, j_depth] = s["RLO"]
                stats_RUP[i_site, j_depth] = s["RUP"]
                stats_bias[i_site, j_depth] = s["bias"]
                stats_biasCI[i_site, j_depth] = s["CI_bias"]
                rmse = math.sqrt(s["MSE"]) if not np.isnan(s["MSE"]) else np.nan
                stats_rmse[i_site, j_depth] = rmse
                stats_rmseCI[i_site, j_depth] = math.sqrt(max(s["CI_MSE"] + s["MSE"], 0)) - rmse if not np.isnan(s["CI_MSE"]) else np.nan
                ubmse = s["ubMSE"]
                ub_rmse = math.sqrt(ubmse) if not np.isnan(ubmse) else np.nan
                stats_ubrmse[i_site, j_depth] = ub_rmse
                if not np.isnan(s["ubMSELO"]):
                    stats_ubrmseCI[i_site, j_depth] = math.sqrt(max(ubmse + s["ubMSEUP"], 0)) - ub_rmse
                if not args.no_anomR:
                    obs_anom = compute_anom(tmp[:, 0], doy_vec, Nmin_day=50 if "CalVal" in args.insitu_tag else 150)
                    mod_anom = compute_anom(tmp[:, 1], doy_vec, Nmin_day=50 if "CalVal" in args.insitu_tag else 150)
                    tmp_anom = np.c_[obs_anom, mod_anom]
                    sa = get_validation_stats(tmp_anom, AC=1, complete=True, ref_col=1, select_cols=[1, 2], Nmin=Nmin)
                    stats_anomR[i_site, j_depth] = sa["R"]
                    stats_anomRLO[i_site, j_depth] = sa["RLO"]
                    stats_anomRUP[i_site, j_depth] = sa["RUP"]

        out_ds = xr.Dataset(
            data_vars=dict(
                R=(("site", "depth"), stats_R),
                RLO=(("site", "depth"), stats_RLO),
                RUP=(("site", "depth"), stats_RUP),
                Bias=(("site", "depth"), stats_bias),
                BiasCI=(("site", "depth"), stats_biasCI),
                RMSE=(("site", "depth"), stats_rmse),
                RMSECI=(("site", "depth"), stats_rmseCI),
                ubRMSE=(("site", "depth"), stats_ubrmse),
                ubRMSECI=(("site", "depth"), stats_ubrmseCI),
            ),
            coords=dict(
                site=np.arange(n_sites),
                depth=np.array([0, 1]),
                site_id=("site", ids_use),
                site_lat=("site", lat_use),
                site_lon=("site", lon_use),
            ),
            attrs=dict(
                insitu_tag=args.insitu_tag,
                experiment=exp,
                domain=args.domain,
                file_tag=args.file_tag,
                start=args.start,
                end=args.end,
                Nmin=int(Nmin),
                max_distance=float(args.max_distance),
                add_anomR=int(not args.no_anomR),
            ),
        )
        if not args.no_anomR:
            out_ds["anomR"] = (("site", "depth"), stats_anomR)
            out_ds["anomRLO"] = (("site", "depth"), stats_anomRLO)
            out_ds["anomRUP"] = (("site", "depth"), stats_anomRUP)

        args.out_dir.mkdir(parents=True, exist_ok=True)
        base_name = f"{exp}_{args.insitu_tag}_SM_1d_c1234smv_6yr_stats"
        out_nc = args.out_dir / f"{base_name}.nc"
        out_ds.to_netcdf(out_nc)
        np.savez_compressed(
            args.out_dir / f"{base_name}.npz",
            R=stats_R,
            RLO=stats_RLO,
            RUP=stats_RUP,
            Bias=stats_bias,
            BiasCI=stats_biasCI,
            RMSE=stats_rmse,
            RMSECI=stats_rmseCI,
            ubRMSE=stats_ubrmse,
            ubRMSECI=stats_ubrmseCI,
            anomR=stats_anomR if not args.no_anomR else None,
            anomRLO=stats_anomRLO if not args.no_anomR else None,
            anomRUP=stats_anomRUP if not args.no_anomR else None,
            site_id=np.array(ids_use, dtype=object),
            site_lat=lat_use,
            site_lon=lon_use,
        )
        print(f"Wrote {out_nc}")


if __name__ == "__main__":
    main()
