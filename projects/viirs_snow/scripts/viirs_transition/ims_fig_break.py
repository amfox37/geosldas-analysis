import sys
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
def day(y, ymd):
    with Dataset(sc.IMS_OLD / f"ims_snowcover_24km_{y}_on_m36_nearest.nc4") as ds:
        t = np.asarray(ds["time"][:]); dates = pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(np.rint(t).astype(int), unit="D")
        i = np.where(dates == pd.Timestamp(ymd))[0][0]
        return np.asarray(ds["ims_category"][i]).astype(int), np.asarray(ds["lat"][:]), np.asarray(ds["lon"][:]), (np.asarray(ds["land_mask"][:]) == 1) & (np.asarray(ds["within_max_distance"][:]) == 1)
a, lat, lon, st = day(2023, "2023-12-31"); b, *_ = day(2024, "2024-01-01")
fig, ax = plt.subplots(1, 3, figsize=(18, 4.2))
cm = matplotlib.colors.ListedColormap(["#999999", "#4a90d9", "#8fbf6a", "#c9a6e0", "#ffffff"])
for k, (c, t) in enumerate([(a, "IMS on M36 land cells, 2023-12-31"), (b, "2024-01-01")]):
    x = np.where(st, c, np.nan); im = ax[k].imshow(x, cmap=cm, vmin=-0.5, vmax=4.5, origin="upper", interpolation="nearest"); ax[k].set_title(t); ax[k].set_xticks([]); ax[k].set_yticks([])
d = np.where(st & (b == 1) & (a != 1), 1.0, np.where(st, 0.0, np.nan))
ax[2].imshow(d, cmap="Reds", origin="upper", interpolation="nearest"); ax[2].set_title("land cells that become IMS code 1 ('sea') on 2024-01-01"); ax[2].set_xticks([]); ax[2].set_yticks([])
fig.suptitle("IMS categories (0 outside, 1 sea=blue, 2 land=green, 3 sea ice=purple, 4 snow=white) on the M36 grid"); fig.tight_layout()
fig.savefig(sc.FIG / "IMS_1_break_2024-01-01_maps.png", dpi=110)
