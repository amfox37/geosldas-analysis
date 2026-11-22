"""IMS-only diagnostics (no model): per year x month, share of each IMS category among land cells on the M36 grid,
mean snow fraction by latitude band, and number of valid cells. Detects product/format changes between years."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc

rows = []
for y in range(2020, 2027):
    for p in (sc.IMS_OLD / f"ims_snowcover_24km_{y}_on_m36_nearest.nc4", sc.VAL / f"ims/ims_snowcover_24km_{y}_on_m36_nearest.nc4"):
        if p.exists(): break
    else: continue
    with Dataset(p) as ds:
        cat = np.asarray(ds["ims_category"][:]).astype(np.int16)
        t = np.asarray(ds["time"][:]); dates = pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(np.rint(t).astype(int), unit="D")
        st = (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1)
        lat = np.asarray(ds["lat"][:])
    nh = st & (lat >= 30)
    for k, d in enumerate(dates):
        c = cat[k]
        v = st & np.isin(c, [0, 1, 2, 3, 4])
        r = dict(year=y, month=d.month, day=d, n_static=int(st.sum()), n_valid=int(v.sum()),
                 n_fill=int((st & ~np.isin(c, [0, 1, 2, 3, 4])).sum()),
                 snow=int((v & (c == 4)).sum()), snow_nh=int((nh & (c == 4)).sum()), valid_nh=int((nh & v).sum()),
                 file=p.name)
        for q in range(5): r[f"c{q}"] = int((st & (c == q)).sum())
        rows.append(r)
df = pd.DataFrame(rows)
df.to_parquet(sc.OUT / "ims_diag_product_daily.parquet")
m = df.groupby(["year", "month"]).agg(n_valid=("n_valid", "mean"), n_fill=("n_fill", "mean"), snow=("snow", "mean"),
        snow_nh=("snow_nh", "mean"), valid_nh=("valid_nh", "mean"), c0=("c0","mean"), c1=("c1","mean"), c2=("c2","mean"), c3=("c3","mean"), c4=("c4","mean"), days=("day","count")).reset_index()
m["snow_frac_nh"] = m.snow_nh / m.valid_nh
m.to_csv(sc.FIG / "IMS_product_monthly.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(m.round(3).to_string())
