"""Shared helpers for analysis C (snow skill of the five VIIRS-transition arms vs IMS, SNOTEL, GHCN).

Everything arm-, date- and path-specific comes from config/arms.json. Run with the GEOSpyD py3.13 python (-I).
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path("/discover/nobackup/projects/land_da/geosldas-analysis")
ARMS_JSON = REPO / "projects/viirs_snow/config/arms.json"
VAL = Path("/discover/nobackup/projects/land_da/VIIRS_testing/validation_data")
IMS_OLD = REPO / "projects/IMS/output"
SNOTEL_OLD = REPO / "projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet"
OUT = Path("/discover/nobackup/projects/land_da/VIIRS_testing/skill_data")  # compact intermediates (land_da)
FIG = REPO / "projects/viirs_snow/docs/viirs_transition_report_figures"
sys.path.insert(0, str(REPO / "common/python/io"))

MAX_DIST_KM = 40.0
MAX_ABS_ELEV_DIFF_M = 500.0
SEASONS = ["DJF", "MAM", "JJA", "SON"]


def load_arms():
    return json.load(open(ARMS_JSON))


def arm_root(cfg, arm):
    a = cfg["arms"][arm]
    return Path(cfg["project_dir"]) / a["exp_id"] / "output" / cfg["domain"]


def daily_file(cfg, arm, day):
    exp = cfg["arms"][arm]["exp_id"]
    return arm_root(cfg, arm) / "cat/ens0000" / f"Y{day:%Y}" / f"M{day:%m}" / f"{exp}.tavg24_1d_lnd_Nt.{day:%Y%m%d}_1200z.nc4"


def arm_days(cfg, arm):
    a = cfg["arms"][arm]
    return pd.date_range(a["beg"], pd.Timestamp(cfg["end"]) - pd.Timedelta(days=1), freq="D")  # 2026-07-31 last


def snow_year_index(dates):
    """0..5 for SY1..SY6 (Aug-Jul); -1 outside."""
    d = pd.DatetimeIndex(dates)
    sy = np.where(d.month >= 8, d.year - 2020, d.year - 2021)
    return np.where((sy >= 0) & (sy <= 5), sy, -1)


def season_index(dates):
    m = pd.DatetimeIndex(dates).month
    return np.select([np.isin(m, [12, 1, 2]), np.isin(m, [3, 4, 5]), np.isin(m, [6, 7, 8])], [0, 1, 2], 3)


def read_tilecoord_quiet():
    from read_GEOSldas import read_tilecoord
    p = Path(load_arms()["project_dir"]) / "viirs_OL/output/SMAP_EASEv2_M36_GLOBAL/rc_out/viirs_OL.ldas_tilecoord.bin"
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return read_tilecoord(str(p))


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.deg2rad(lat1), np.deg2rad(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.deg2rad(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def representative_cells(tc):
    i = np.asarray(tc["i_indg"], np.int64); j = np.asarray(tc["j_indg"], np.int64)
    tid = np.asarray(tc["tile_id"], np.int64); fr = np.asarray(tc["frac_cell"], np.float64)
    fr = np.where(np.isfinite(fr), fr, -np.inf)
    nx = int(i.max()) + 1
    code = j * nx + i
    order = np.lexsort((tid, -fr, code))
    cs = code[order]
    first = np.r_[True, cs[1:] != cs[:-1]]
    rep = order[first]
    return rep, i[rep], j[rep]


def map_to_tiles(tc, lat, lon, elev_m):
    """Nearest tile (centre of mass) per station; returns DataFrame of tile index, distance, elev diff."""
    tlat = np.asarray(tc["com_lat"], float); tlon = np.asarray(tc["com_lon"], float); telev = np.asarray(tc["elev"], float)
    tix, dist = [], []
    for la, lo in zip(lat, lon):
        d = haversine_km(la, lo, tlat, tlon)
        j = int(np.argmin(d)); tix.append(j); dist.append(d[j])
    tix = np.array(tix)
    return pd.DataFrame({"tile": tix, "dist_km": dist, "tile_elev_m": telev[tix], "elev_diff_m": telev[tix] - np.asarray(elev_m, float)})


def load_snotel():
    """Return obs (date, station, swe_mm, snwd_mm) for 2020-08-01..2026-07-31 and station table; inches -> mm."""
    cols = ["date", "stationTriplet", "WTEQ", "SNWD", "latitude", "longitude", "elevation"]
    old = pd.read_parquet(SNOTEL_OLD, columns=cols)
    old = old[old["date"] >= "2020-08-01"]
    new = pd.read_parquet(VAL / "snotel/snotel_daily_wteq_snwd_20251001_20260731.parquet", columns=cols)
    df = pd.concat([old, new], ignore_index=True)
    df["date"] = pd.to_datetime(df["date"]).dt.normalize()
    df = df[(df["date"] >= "2020-08-01") & (df["date"] <= "2026-07-31")]
    df = df.drop_duplicates(["stationTriplet", "date"], keep="last")
    df["swe_mm"] = pd.to_numeric(df["WTEQ"], errors="coerce") * 25.4      # in -> mm (= kg m-2)
    df["snwd_mm"] = pd.to_numeric(df["SNWD"], errors="coerce") * 25.4     # in -> mm
    st = (df.groupby("stationTriplet").agg(lat=("latitude", "first"), lon=("longitude", "first"), elev_ft=("elevation", "first")).reset_index()
          .rename(columns={"stationTriplet": "station"}))
    for c in ("lat", "lon"): st[c] = pd.to_numeric(st[c], errors="coerce")
    st["elev_m"] = pd.to_numeric(st["elev_ft"], errors="coerce") * 0.3048
    st = st[np.isfinite(st.lat) & np.isfinite(st.lon)].reset_index(drop=True)
    df = df.rename(columns={"stationTriplet": "station"})[["date", "station", "swe_mm", "snwd_mm"]]
    return df, st


def load_ghcn():
    """Return obs (date, station, snwd_mm) (already mm) and station table for stations with any snow in window."""
    cov = pd.read_parquet(VAL / "ghcn/ghcn_snwd_station_coverage.parquet")
    cov = cov[np.isfinite(cov["lat"]) & np.isfinite(cov["lon"]) & (cov["n_days"] >= 180) & (cov["max_snwd_mm"] > 0)].copy()
    st = cov.rename(columns={"station_id": "station", "elev": "elev_m"})[["station", "lat", "lon", "elev_m", "n_days", "max_snwd_mm"]].reset_index(drop=True)
    parts = [pd.read_parquet(VAL / f"ghcn/ghcn_snwd_{y}.parquet", columns=["station_id", "date", "snwd_mm"]) for y in range(2020, 2027)]
    df = pd.concat(parts, ignore_index=True).rename(columns={"station_id": "station"})
    df["date"] = pd.to_datetime(df["date"]).dt.normalize()
    df = df[df["station"].isin(set(st["station"]))]
    return df, st
