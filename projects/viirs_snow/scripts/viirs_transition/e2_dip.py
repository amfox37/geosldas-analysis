"""E2 item 1: the 2023 week-8 J1 obs dip. Daily J1/J2 obs counts (obs_availability, TA_to_J1J2 arm) and the input file sizes
(/discover/nobackup/projects/gmao/polarm/JPSS/VJ{1,2}10C1/<yyyy>) for 15 Sep - 10 Oct of 2023-25, plus a scan of missing /
small files over 2023-08..2026-07. Writes E2_dip_obs_vs_files.png and E2_dip_files.csv. Run with python -I."""
import os, re, sys
from pathlib import Path
import numpy as np, pandas as pd, netCDF4
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent)); import skill_common as sc
REPO = sc.REPO / "projects/viirs_snow"; JPSS = Path("/discover/nobackup/projects/gmao/polarm/JPSS")
ds = netCDF4.Dataset(REPO / "output/obs_availability/obs_availability_TA_to_J1J2.nc")
d = pd.to_datetime(ds["date"][:].astype(str), format="%Y%m%d"); n = ds["n_obs"][:].sum(axis=2)
df = pd.DataFrame(n, index=d, columns=[str(s) for s in ds["species_name"][:]])
def fsize(prod, day):
    p = JPSS / prod / f"{day.year}" / f"{prod}.A{day.year}{day.dayofyear:03d}.002.h5"
    return p.stat().st_size / 1e6 if p.exists() else 0.0
rows = []
for prod in ("VJ110C1", "VJ210C1"):
    s = pd.Series({pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=int(re.match(rf"{prod}\.A\d{{4}}(\d{{3}})", f)[1]) - 1): (JPSS / prod / str(y) / f).stat().st_size / 1e6
                   for y in range(2023, 2027) for f in os.listdir(JPSS / prod / str(y)) if re.match(rf"{prod}\.A\d{{4}}\d{{3}}\.002\.h5$", f)}).sort_index()["2023-08-01":"2026-07-31"]
    med = s.rolling(15, center=True, min_periods=5).median()
    for t in pd.date_range("2023-08-01", "2026-07-31").difference(s.index): rows.append(dict(product=prod, date=t.date(), kind="missing", MB=0.0))
    for t, v in s[s < 0.7 * med].items(): rows.append(dict(product=prod, date=t.date(), kind="small(<70% of 15-d median)", MB=round(v, 2)))
pd.DataFrame(rows).to_csv(sc.FIG / "E2_dip_files.csv", index=False)
fig, ax = plt.subplots(2, 3, figsize=(13, 5.5), sharey="row")
for j, y in enumerate((2023, 2024, 2025)):
    x = df.loc[f"{y}-09-15":f"{y}-10-10"]
    ax[0, j].plot(x.index, x.J1 / 1e3, "o-", ms=3, label="J1"); ax[0, j].plot(x.index, x.J2 / 1e3, "s-", ms=3, label="J2")
    ax[0, j].set_title(f"{y}: obs per day (ObsFcstAna, all steps)"); ax[0, j].set_ylim(0, None)
    t = x.index
    ax[1, j].bar(t - pd.Timedelta(hours=5), [fsize("VJ110C1", a) for a in t], width=0.4, label="VJ110C1 file MB")
    ax[1, j].bar(t + pd.Timedelta(hours=5), [fsize("VJ210C1", a) for a in t], width=0.4, label="VJ210C1 file MB")
    ax[1, j].set_title(f"{y}: input file size")
    for a in ax[:, j]: a.tick_params(axis="x", rotation=45)
ax[0, 0].set_ylabel("obs / 1000"); ax[1, 0].set_ylabel("MB"); ax[0, 0].legend(); ax[1, 0].legend()
fig.suptitle("E2: J1 input file missing on 2023-09-30 and half-size on 2023-10-01 -> J1 obs dip (J2 unaffected)")
fig.tight_layout(); fig.savefig(sc.FIG / "E2_dip_obs_vs_files.png", dpi=130)
print(pd.DataFrame(rows).to_string())
