import sys
sys.path.insert(0, "/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params")
sys.path.insert(0, "/gpfsm/dnb34/tdirs/login/discover34.252777.amfox/claude-236601614/-gpfsm-dnb06-projects-p284-geosldas-analysis/bad106a1-1ee2-4c2d-96cc-d1c5e084a6fd/scratchpad")
from pathlib import Path
import numpy as np
from obs_scaling.io import read_obs_fcst_ana
from obs_scaling.tile_io import read_tilecoord, read_tilegrids
from obs_scaling.seqbin import read_tb_scaling_file

EXP = Path("/discover/nobackup/projects/land_da/obs_scaling_test/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL")
tile_coord = read_tilecoord(EXP / "rc_out" / "LS_OLv8_M36.ldas_tilecoord.bin")

SPECIES = {6, 8}  # SMAP_L1C_Tbh_D, SMAP_L1C_Tbv_D

observed_tilenums = set()
files_read = 0
for month_dir in sorted((EXP / "ana" / "ens_avg" / "Y2016").glob("M0[1-3]")):
    for f in sorted(month_dir.glob("*.ldas_ObsFcstAna.*.bin")):
        rec = read_obs_fcst_ana(f)
        if rec is None:
            continue
        files_read += 1
        mask = np.isin(rec.obs_species, list(SPECIES))
        observed_tilenums.update(rec.obs_tilenum[mask].astype(np.int64).tolist())

print(f"Files read: {files_read}")
print(f"Total model tiles (tilecoord): {tile_coord.n_tile}")
print(f"Unique observed obs_tilenum values (species {SPECIES}): {len(observed_tilenums)}")

# obs_tilenum is 1-based model tile index -> corresponding tile_id
observed_tilenums_arr = np.array(sorted(observed_tilenums), dtype=np.int64)
model_idx = observed_tilenums_arr - 1
out_of_range = (model_idx < 0) | (model_idx >= tile_coord.n_tile)
print(f"obs_tilenum out of model tile range: {int(out_of_range.sum())}")
observed_tile_ids = set(tile_coord.tile_id[model_idx[~out_of_range]].astype(np.int64).tolist())

scale_file = read_tb_scaling_file(
    str(EXP / "stats" / "z_score_clim" /
        "PYTHON_TEST_M21C_SMAP_Tb_2016_p1_2016_p18_hscale_0.00_W_15p_Nmin_20_D_p09.bin")
)
scaling_tile_ids = scale_file.tile_id.astype(np.int64)
scaling_tile_id_set = set(scaling_tile_ids.tolist())

print(f"Total M36 administering tiles (scaling file): {len(scaling_tile_id_set)}")
absent = observed_tile_ids - scaling_tile_id_set
print(f"Observed tiles absent from scaling file: {len(absent)}")
if absent:
    print("  sample absent:", list(absent)[:10])

dup_count = scaling_tile_ids.size - len(scaling_tile_id_set)
print(f"Duplicate scaling tile IDs: {dup_count}")
