#!/bin/bash
#SBATCH --job-name=ims_regen
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:20:00
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/regen.%j.out
P=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
S=/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py
for spec in "0:SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50:SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_regen" "10:SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50:SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10"; do
  m=${spec%%:*}; rest=${spec#*:}; src=${rest%%:*}; tag=${rest#*:}
  $P -I $S --domain SMAP_EASEv2_M36_GLOBAL --year-start 2000 --year-end 2024 \
    --reuse-cell-counts-nc /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out/ims_ol_da_cell_counts_metrics_$src.nc4 --reuse-scope-meta-csv /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out/ims_ol_da_scope_metadata_$src.csv \
    --min-ims-snow-days $m --output-dir /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out --output-tag $tag --overwrite
done
