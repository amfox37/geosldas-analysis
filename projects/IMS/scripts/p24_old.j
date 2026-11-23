#!/bin/bash
#SBATCH --job-name=ims_p24_old
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/p24_old.%j.out
bash /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/run2024.sh /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/regrid_old /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out_old2024 SMAP_EASEv2_M36_GLOBAL_2024_2024_thr0p50
