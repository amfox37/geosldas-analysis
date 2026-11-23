import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
with Dataset(sc.IMS_OLD / "ims_snowcover_24km_2024.nc4") as ds:
    r1 = np.asarray(ds["ims_snowcover"][0]); r7 = np.asarray(ds["ims_snowcover"][200])
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
cm = matplotlib.colors.ListedColormap(["#999999", "#4a90d9", "#8fbf6a", "#c9a6e0", "#ffffff"])
for a, r, t in zip(ax, (r1, r7), ("raw 24 km IMS 2024-01-01", "raw 24 km IMS 2024-07-19")):
    a.imshow(r, cmap=cm, vmin=-0.5, vmax=4.5, origin="upper"); a.set_title(t); a.set_xticks([]); a.set_yticks([])
fig.tight_layout(); fig.savefig(sc.FIG / "IMS_2_raw_2024_orientation_check.png", dpi=90)
