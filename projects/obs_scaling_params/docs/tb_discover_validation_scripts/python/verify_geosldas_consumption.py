import sys
sys.path.insert(0, "/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params")
import glob
import numpy as np
import netCDF4 as nc
from obs_scaling.seqbin import read_tb_scaling_file

STATS_DIR = "/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121/output/SMAP_EASEv2_M36_GLOBAL/stats/z_score_clim"
p01 = read_tb_scaling_file(f"{STATS_DIR}/GEOSLDAS_SMOKE_VALID_SMAP_Tb_2019_p1_2019_p6_hscale_0.00_W_15p_Nmin_20_D_p01.bin")
p02 = read_tb_scaling_file(f"{STATS_DIR}/GEOSLDAS_SMOKE_VALID_SMAP_Tb_2019_p1_2019_p6_hscale_0.00_W_15p_Nmin_20_D_p02.bin")
NODATA = -9999.0

def stats_for_tile(scfile, tile_id, pol):
    idx = np.where(scfile.tile_id == tile_id)[0]
    if idx.size == 0:
        return None
    i = idx[0]
    base = pol * 5
    mean_obs, std_obs, mean_mod, std_mod, ndata = scfile.data[base:base+5, i, 0]
    return mean_obs, std_obs, mean_mod, std_mod, ndata

OUT_DIR = "/discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y2020/M01"
files = sorted(glob.glob(f"{OUT_DIR}/*.nc4"))
print(f"{len(files)} ObsFcstAna output files")

ERRSTD = 4.0  # obs_param_nml(32/34)%errstd (pre-scaling)

n_obsvar_checked = 0
n_obsvar_match = 0
max_obsvar_rel_err = 0.0

n_formula_checked = 0
n_formula_match = 0
max_formula_rel_err = 0.0

n_assim1_valid = 0
n_assim1_invalid = 0
n_assim0_valid = 0
n_assim0_invalid = 0

examples = []

for f in files:
    dt_str = f.split(".ldas_ObsFcstAna.")[1].split("z.nc4")[0]
    day = int(dt_str[6:8])
    pentad = 1 if day <= 5 else 2
    scfile = p01 if pentad == 1 else p02
    with nc.Dataset(f) as ds:
        species = ds.variables["species"][:]
        tilenum = ds.variables["tilenum"][:]
        obs = ds.variables["obs"][:]
        obsvar = ds.variables["obsvar"][:]
        assim_flag = ds.variables["assim_flag"][:]
        for sp, pol in ((6, 0), (8, 1)):
            mask = species == sp
            if not np.any(mask):
                continue
            tiles = tilenum[mask]
            ov = obsvar[mask]
            av = assim_flag[mask]
            obsv = obs[mask]
            for j in range(tiles.size):
                st = stats_for_tile(scfile, tiles[j], pol)
                if st is None:
                    continue
                mean_obs, std_obs, mean_mod, std_mod, ndata = st
                valid_stats = not np.isclose(mean_obs, NODATA, atol=1.0)
                if av[j] == 1:
                    if valid_stats:
                        n_assim1_valid += 1
                    else:
                        n_assim1_invalid += 1
                else:
                    if valid_stats:
                        n_assim0_valid += 1
                    else:
                        n_assim0_invalid += 1

                if valid_stats and av[j] == 1 and std_obs > 0:
                    o = obsv[j]
                    v = ov[j]
                    if np.isfinite(o) and o < 1e10:
                        # invert: raw_obs = mean_obs + (scaled_obs - mean_mod) * std_obs/std_mod
                        raw_obs_reconstructed = mean_obs + (o - mean_mod) * std_obs / std_mod
                        n_formula_checked += 1
                        plausible = 150 < raw_obs_reconstructed < 340
                        if plausible:
                            n_formula_match += 1
                        if len(examples) < 5:
                            examples.append((f.split('/')[-1], int(sp), int(tiles[j]),
                                              float(mean_obs), float(std_obs), float(mean_mod), float(std_mod),
                                              float(o), float(raw_obs_reconstructed)))
                    if np.isfinite(v) and v < 1e10:
                        expected_obsvar = (ERRSTD * std_mod / std_obs) ** 2
                        n_obsvar_checked += 1
                        rel_err = abs(v - expected_obsvar) / expected_obsvar
                        max_obsvar_rel_err = max(max_obsvar_rel_err, rel_err)
                        if rel_err < 1e-3:
                            n_obsvar_match += 1

print(f"\nassim_flag=1 & valid stats   : {n_assim1_valid}  (assimilable case)")
print(f"assim_flag=1 & INVALID stats : {n_assim1_invalid}  (must be 0)")
print(f"assim_flag=0 & valid stats   : {n_assim0_valid}  (other QC rejection, not stats-related)")
print(f"assim_flag=0 & invalid stats : {n_assim0_invalid}  (do-not-assimilate due to missing stats)")

print(f"\nobsvar == errstd^2*(std_mod/std_obs)^2 : {n_obsvar_match}/{n_obsvar_checked} matched within 0.1%")
print(f"max relative error in obsvar: {max_obsvar_rel_err:.3e}")

print(f"\nscaled-obs inverse-transform plausibility: {n_formula_match}/{n_formula_checked} in [150,340]K")
print("\nExample records (file, species, tile, mean_obs, std_obs, mean_mod, std_mod, scaled_obs, reconstructed_raw_obs):")
for ex in examples:
    print(" ", ex)
