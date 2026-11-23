"""E2 helper: IMS contingency counts (A,B,C,D x season x cell) for 2023-08-01..2023-08-28 only, per arm, same recipe as
skill_extract.py, so skill_augmask.py can subtract them from SY4 x JJA. Usage: skill_augmask_ims.py <arm>.
Output: <OUT>/augmask_ims_cnt_<arm>.npz. Run with python -I."""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
arm = sys.argv[1]
cfg = sc.load_arms()
days = pd.date_range("2023-08-01", "2023-08-28")
tc = sc.read_tilecoord_quiet()
rep, rep_i, rep_j = sc.representative_cells(tc)
z = np.load(sc.OUT / f"skill_{arm}.npz")
assert (z["rep"] == rep).all()
p = sc.IMS_OLD / "ims_snowcover_24km_2023_on_m36_nearest.nc4"
if not p.exists(): p = sc.VAL / "ims/ims_snowcover_24km_2023_on_m36_nearest.nc4"
with Dataset(p) as ds:
    cat = np.asarray(ds["ims_category"][:]).astype(np.int16)
    t = np.asarray(ds["time"][:]); dates = pd.Timestamp("2023-01-01") + pd.to_timedelta(np.rint(t).astype(int), unit="D")
    static = (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1)
cnt = np.zeros((4, rep.size), np.int16); nfound = 0
for day in days:
    f = sc.daily_file(cfg, arm, day)
    if not f.exists(): continue
    with Dataset(f) as ds: scf = np.asarray(ds["FRLANDSNO"][0, :], np.float32)
    scf[(scf > 1e14) | (scf < 0)] = np.nan
    ti = np.where(dates == day)[0]
    if ti.size == 0: continue
    nfound += 1
    c = cat[ti[0]][rep_j, rep_i]
    valid = static[rep_j, rep_i] & np.isin(c, [0, 1, 2, 3, 4]) & np.isfinite(scf[rep])
    so = c == 4; sm = scf[rep] >= 0.5
    for q, msk in enumerate((sm & so, sm & ~so, ~sm & so, ~sm & ~so)): cnt[q] += (valid & msk).astype(np.int16)
np.savez_compressed(sc.OUT / f"augmask_ims_cnt_{arm}.npz", cnt=cnt, ndays=nfound)
print(arm, "days with IMS+model:", nfound, "total valid-cell-days:", int(cnt.sum()))
