"""Extract, for one arm, (1) per-cell IMS contingency counts by snow year x season and (2) model snow at candidate
SNOTEL / GHCN station tiles for every day. Usage: skill_extract.py <arm> [--ndays N] (N for timing tests).
Output: <OUT>/skill_<arm>.npz (compact)."""
import argparse, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc

ap = argparse.ArgumentParser(); ap.add_argument("arm"); ap.add_argument("--ndays", type=int, default=None)
ap.add_argument("--tag", default="")
args = ap.parse_args()
cfg = sc.load_arms(); arm = args.arm
days = sc.arm_days(cfg, arm)
if args.ndays: days = days[:args.ndays]
tc = sc.read_tilecoord_quiet()
rep, rep_i, rep_j = sc.representative_cells(tc)
ncell = rep.size

# station candidates + tile map (same for every arm: identical tilecoord)
_, st_sn = sc.load_snotel(); _, st_gh = sc.load_ghcn()
m_sn = sc.map_to_tiles(tc, st_sn.lat, st_sn.lon, st_sn.elev_m); m_gh = sc.map_to_tiles(tc, st_gh.lat, st_gh.lon, st_gh.elev_m)
tiles_sn = m_sn.tile.values; tiles_gh = m_gh.tile.values

nd = len(days)
shp_sn = (nd, len(tiles_sn)); shp_gh = (nd, len(tiles_gh))
sn_snod = np.full(shp_sn, np.nan, np.float32); sn_swe = np.full(shp_sn, np.nan, np.float32)
gh_snod = np.full(shp_gh, np.nan, np.float32); gh_scf = np.full(shp_gh, np.nan, np.float32)
cnt = np.zeros((4, 6, 4, ncell), np.int16)   # A,B,C,D x SY x season x cell
found = np.zeros(nd, bool)

def ims_for_year(y):
    for p in (sc.IMS_OLD / f"ims_snowcover_24km_{y}_on_m36_nearest.nc4", sc.VAL / f"ims/ims_snowcover_24km_{y}_on_m36_nearest.nc4"):
        if p.exists(): break
    else: return None
    with Dataset(p) as ds:
        cat = np.asarray(ds["ims_category"][:]).astype(np.int16)
        t = np.asarray(ds["time"][:]); base = pd.Timestamp(f"{y}-01-01")
        dates = base + pd.to_timedelta(np.rint(t).astype(int), unit="D")
        lm = np.asarray(ds["land_mask"][:]) == 1; wm = np.asarray(ds["within_max_distance"][:]) == 1
    return cat, pd.DatetimeIndex(dates), lm & wm

ims = {}
syi = sc.snow_year_index(days); sei = sc.season_index(days)
t0 = time.time()
for k, day in enumerate(days):
    f = sc.daily_file(cfg, arm, day)
    if not f.exists(): continue
    with Dataset(f) as ds:
        snod = np.asarray(ds["SNODPLAND"][0, :], np.float32); swe = np.asarray(ds["SNOMASLAND"][0, :], np.float32)
        scf = np.asarray(ds["FRLANDSNO"][0, :], np.float32)
    for a in (snod, swe, scf): a[(a > 1e14) | (a < 0)] = np.nan
    found[k] = True
    sn_snod[k] = snod[tiles_sn]; sn_swe[k] = swe[tiles_sn]
    gh_snod[k] = snod[tiles_gh]; gh_scf[k] = scf[tiles_gh]
    y = day.year
    if y not in ims: ims[y] = ims_for_year(y)
    if ims[y] is None or syi[k] < 0: continue
    cat, dates, static = ims[y]
    ti = np.where(dates == day)[0]
    if ti.size == 0: continue
    c = cat[ti[0]][rep_j, rep_i]
    valid = static[rep_j, rep_i] & np.isin(c, [0, 1, 2, 3, 4]) & np.isfinite(scf[rep])
    so = c == 4; sm = scf[rep] >= 0.5
    for q, msk in enumerate((sm & so, sm & ~so, ~sm & so, ~sm & ~so)):
        cnt[q, syi[k], sei[k]] += (valid & msk).astype(np.int16)
    if (k + 1) % 200 == 0: print(f"{arm} {day:%Y-%m-%d} {time.time()-t0:.0f}s", flush=True)

out = sc.OUT / f"skill_{arm}{args.tag}.npz"
np.savez_compressed(out, days=days.values.astype("datetime64[D]"), found=found, cnt=cnt, rep=rep, rep_i=rep_i, rep_j=rep_j,
                    sn_snod=sn_snod, sn_swe=sn_swe, gh_snod=gh_snod, gh_scf=gh_scf)
m_sn.assign(station=st_sn.station.values).to_csv(sc.OUT / "snotel_candidate_tiles.csv", index=False)
m_gh.assign(station=st_gh.station.values).to_csv(sc.OUT / "ghcn_candidate_tiles.csv", index=False)
print(f"done {arm}: {found.sum()}/{nd} days found, {time.time()-t0:.0f}s -> {out}")
