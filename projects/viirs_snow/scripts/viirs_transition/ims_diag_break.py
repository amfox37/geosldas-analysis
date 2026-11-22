"""Where/how does the IMS category field change at 2024-01-01? Compares last/first days around the break on the M36 grid
and in the raw 24 km source files."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
P = lambda y: sc.IMS_OLD / f"ims_snowcover_24km_{y}_on_m36_nearest.nc4"
def day(y, ymd):
    with Dataset(P(y)) as ds:
        t = np.asarray(ds["time"][:]); dates = pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(np.rint(t).astype(int), unit="D")
        i = np.where(dates == pd.Timestamp(ymd))[0][0]
        return np.asarray(ds["ims_category"][i]).astype(int), np.asarray(ds["lat"][:]), np.asarray(ds["lon"][:]), (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1), np.asarray(ds["source_y"][:]), np.asarray(ds["source_x"][:])
a, lat, lon, st, sya, sxa = day(2023, "2023-12-31"); b, _, _, _, syb, sxb = day(2024, "2024-01-01")
print("source maps identical 2023 vs 2024:", np.array_equal(sya, syb) and np.array_equal(sxa, sxb))
for y in (2023, 2024):
    for d in (f"{y}-12-30" if y == 2023 else "2024-01-02",):
        pass
# daily transition around break
import itertools
for ymd in ["2023-12-29", "2023-12-30", "2023-12-31", "2024-01-01", "2024-01-02", "2024-01-03"]:
    c, *_ = day(int(ymd[:4]), ymd)
    print(ymd, {k: int((st & (c == k)).sum()) for k in range(5)})
print("transition matrix 2023-12-31 -> 2024-01-01 over static cells (rows old, cols new, cats 0-4):")
M = np.zeros((5, 5), int)
for i in range(5):
    for j in range(5): M[i, j] = int((st & (a == i) & (b == j)).sum())
print(M)
ch = st & (a != b)
for nm, sel in [("c1 new (land->'1')", st & (b == 1) & (a != 1))]:
    print(nm, int(sel.sum()), "lat quantiles", np.percentile(lat[sel], [5, 25, 50, 75, 95]).round(1), "lon quantiles", np.percentile(lon[sel], [5, 25, 50, 75, 95]).round(1))
    print("   old category counts:", {k: int((sel & (a == k)).sum()) for k in range(5)})
# raw source file
for y in (2023, 2024):
    f = sc.IMS_OLD / f"ims_snowcover_24km_{y}.nc4"
    if not f.exists(): print("no raw", f); continue
    with Dataset(f) as ds:
        v = ds["ims_snowcover"]; print(y, "raw", v.shape, {k: getattr(ds, k) for k in ds.ncattrs()})
        t = np.asarray(ds["time"][:]); 
        for k in ([int(np.argmin(abs(t - (364 if y == 2023 else 0))))] ):
            x = np.asarray(v[k]); u, n = np.unique(x[np.isfinite(x)] if x.dtype.kind == "f" else x, return_counts=True); print("  day idx", k, dict(zip(u.tolist(), n.tolist())))
