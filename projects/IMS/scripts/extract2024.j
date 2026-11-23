#!/bin/bash
#SBATCH --job-name=ims_fix_extract
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:30:00
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/extract2024.%j.out
set -u
W=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix; R=/discover/nobackup/projects/land_da/M21C_land_sweeper
for spec in "LS_DAv8_M36:$R/LS_DAv8_M36_v3/LS_DAv8_M36" "LS_OLv8_M36:$R/LS_OLv8_M36_v2/LS_OLv8_M36"; do
  exp=${spec%%:*}; root=${spec#*:}
  d=$W/model/$exp/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg; mkdir -p $d
  ( time tar -xf $root/output/SMAP_EASEv2_M36_GLOBAL/SMAP_EASEv2_M36_GLOBAL_cat_tars/Y2024.tar -C $d --wildcards '*tavg24_1d_lnd_Nt.2024*_1200z.nc4' ; echo $exp done rc=$? ) &
done
wait
for e in LS_DAv8_M36 LS_OLv8_M36; do echo $e $(find $W/model/$e -name '*.nc4' | wc -l) files; done
