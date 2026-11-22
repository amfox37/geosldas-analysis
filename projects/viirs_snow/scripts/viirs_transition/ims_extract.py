"""IMS re-extraction with the 2024+ vertical-flip fix. For one arm: per-cell contingency counts (A hit, B false alarm,
C miss, D correct-negative) by snow year x season (same layout as skill_extract 'cnt', IMS corrected) and daily
domain counts for BOTH the original (as scored in analysis C) and the corrected IMS, for NH>=30N and all valid cells.
Why: IMS 24 km arrays for 2024-01-01 onward are stored in the file's native bottom-up row order, 2000-2023 were not;
the nearest-neighbour source indices assume the 2000-2023 order, so 2024+ must be flipped vertically (flipud) first
(see ims_diag_orient.py). Usage: ims_extract.py <arm> [--ndays N]. Output: <OUT>/imsfix/ims_<arm>.npz."""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc

ap = argparse.ArgumentParser(); ap.add_argument("arm"); ap.add_argument("--ndays", type=int, default=None)
args = ap.parse_args()
cfg = sc.load_arms(); arm = args.arm
days = sc.arm_days(cfg, arm)
if args.ndays: days = days[:args.ndays]
tc = sc.read_tilecoord_quiet()
rep, rep_i, rep_j = sc.representative_cells(tc)
ncell = rep.size
with Dataset(sc.IMS_OLD / "ims_snowcover_24km_2023_on_m36_nearest.nc4") as ds:
    sy = np.asarray(ds["source_y"][:]); sx = np.asarray(ds["source_x"][:])
    static = (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1)
    lat = np.asarray(ds["lat"][:])
st_rep = static[rep_j, rep_i]; nh_rep = lat[rep_j, rep_i] >= 30
sy_rep = sy[rep_j, rep_i]; sx_rep = sx[rep_j, rep_i]
ok_src = sy_rep >= 0

def open_year(y):
    m36 = sc.IMS_OLD / f"ims_snowcover_24km_{y}_on_m36_nearest.nc4"
    if not m36.exists(): m36 = sc.VAL / f"ims/ims_snowcover_24km_{y}_on_m36_nearest.nc4"
    raw = sc.IMS_OLD / f"ims_snowcover_24km_{y}.nc4"
    if not raw.exists(): raw = sc.VAL / f"ims/ims_snowcover_24km_{y}.nc4"
    d1 = Dataset(m36); t = np.asarray(d1["time"][:]); dates = pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(np.rint(t).astype(int), unit="D")
    d2 = Dataset(raw) if (y >= 2024 and raw.exists()) else None
    return d1, d2, pd.DatetimeIndex(dates)

yrs = {}
nd = len(days)
cnt = np.zeros((4, 6, 4, ncell), np.int16)
dst = np.zeros((nd, 2, 2, 4), np.int32)   # day, version(0 as-scored, 1 fixed), region(0 NH>=30, 1 all), A..D
nv = np.zeros((nd, 2), np.int32)          # valid cells (all region)
found = np.zeros(nd, bool)
syi = sc.snow_year_index(days); sei = sc.season_index(days)
t0 = time.time()
for k, day in enumerate(days):
    f = sc.daily_file(cfg, arm, day)
    if not f.exists(): continue
    with Dataset(f) as ds: scf = np.asarray(ds["FRLANDSNO"][0, :], np.float32)
    scf[(scf > 1e14) | (scf < 0)] = np.nan
    y = day.year
    if y not in yrs: yrs[y] = open_year(y)
    d1, d2, dates = yrs[y]
    ti = np.where(dates == day)[0]
    if ti.size == 0: continue
    c_old = np.asarray(d1["ims_category"][ti[0]]).astype(np.int16)[rep_j, rep_i]
    if d2 is not None:
        raw = np.asarray(d2["ims_snowcover"][ti[0]])[::-1]
        c_new = np.where(ok_src, raw[np.where(ok_src, sy_rep, 0), np.where(ok_src, sx_rep, 0)], -1).astype(np.int16)
    else:
        c_new = c_old
    found[k] = True
    sm = scf[rep] >= 0.5; fin = np.isfinite(scf[rep])
    for v, c in enumerate((c_old, c_new)):
        valid = st_rep & np.isin(c, [0, 1, 2, 3, 4]) & fin
        so = c == 4
        for q, msk in enumerate((sm & so, sm & ~so, ~sm & so, ~sm & ~so)):
            sel = valid & msk
            dst[k, v, 0, q] = int((sel & nh_rep).sum()); dst[k, v, 1, q] = int(sel.sum())
            if v == 1 and syi[k] >= 0: cnt[q, syi[k], sei[k]] += sel.astype(np.int16)
        nv[k, v] = int(valid.sum())
    if (k + 1) % 300 == 0: print(f"{arm} {day:%Y-%m-%d} {time.time()-t0:.0f}s", flush=True)
(sc.OUT / "imsfix").mkdir(exist_ok=True)
out = sc.OUT / "imsfix" / f"ims_{arm}.npz"
np.savez_compressed(out, days=days.values.astype("datetime64[D]"), found=found, cnt=cnt, rep=rep, rep_i=rep_i, rep_j=rep_j, dst=dst, nv=nv)
print(f"done {arm}: {found.sum()}/{nd} days, {time.time()-t0:.0f}s -> {out}")
