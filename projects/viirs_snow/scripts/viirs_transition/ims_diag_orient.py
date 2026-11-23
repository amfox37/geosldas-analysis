"""Test orientation of the raw 24 km IMS arrays vs the lat/lon file: for each candidate transform of the raw 2024 array
(and of the 2025 array), map to M36 land cells via the stored nearest source indices and report the fraction of M36 land
cells labelled 'sea' (code 1) and total snow cells. 2023 (pre-break) value ~6% is the reference."""
import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
with Dataset(sc.IMS_OLD / "ims_snowcover_24km_2023_on_m36_nearest.nc4") as ds:
    sy = np.asarray(ds["source_y"][:]); sx = np.asarray(ds["source_x"][:])
    st = (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1) & (sy >= 0)
    ref = np.asarray(ds["ims_category"][0 if False else 364])
    print("2023-12-31 sea frac of land cells:", round(float((ref[st] == 1).mean()), 3), "snow", int((ref[st] == 4).sum()))
def stats(raw, name):
    T = {"as-is": raw, "flipud": raw[::-1], "fliplr": raw[:, ::-1], "rot180": raw[::-1, ::-1], "T": raw.T, "T+flipud": raw.T[::-1], "T+fliplr": raw.T[:, ::-1], "T+rot180": raw.T[::-1, ::-1]}
    for k, a in T.items():
        c = a[sy[st], sx[st]]
        print(f"  {name:12s} {k:9s} sea-frac {float((c==1).mean()):.3f}  snow {int((c==4).sum()):6d}  land {int((c==2).sum()):6d}")
for y, idx in ((2024, 0), (2024, 200), (2025, 0)):
    p = sc.IMS_OLD / f"ims_snowcover_24km_{y}.nc4" if y == 2024 else sc.VAL / f"ims/ims_snowcover_24km_{y}.nc4"
    with Dataset(p) as ds: raw = np.asarray(ds["ims_snowcover"][idx]).astype(int)
    stats(raw, f"{y} idx{idx}")
