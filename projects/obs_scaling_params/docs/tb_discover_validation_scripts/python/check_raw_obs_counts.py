import sys
sys.path.insert(0, "/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params")
from pathlib import Path
from collections import Counter
import numpy as np
from obs_scaling.io import read_obs_fcst_ana

EXP = Path("/discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL")
counts = Counter()
n_files = 0
n_obs_total = 0
for month_dir in sorted((EXP / "ana" / "ens_avg" / "Y2020").glob("M0[1-2]")):
    for f in sorted(month_dir.glob("*.ldas_ObsFcstAna.*.nc4")):
        rec = read_obs_fcst_ana(f)
        if rec is None:
            continue
        n_files += 1
        mask = rec.obs_species == 6
        n_obs_total += int(mask.sum())
        for t in rec.obs_tilenum[mask]:
            counts[int(t)] += 1

print(f"files: {n_files}, total species-6 obs: {n_obs_total}")
if counts:
    vals = np.array(list(counts.values()))
    print(f"tiles with >=1 obs: {len(counts)}")
    print(f"max count per tile: {vals.max()}, mean: {vals.mean():.2f}, median: {np.median(vals)}")
    print(f"tiles with >=20 obs: {(vals>=20).sum()}")
    print(f"tiles with >=10 obs: {(vals>=10).sum()}")
