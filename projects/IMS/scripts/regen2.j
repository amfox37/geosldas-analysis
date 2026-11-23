#!/bin/bash
#SBATCH --job-name=ims_regen2
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:25:00
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/regen2.%j.out
/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3 -I /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py --domain SMAP_EASEv2_M36_GLOBAL --year-start 2000 --year-end 2024 \
  --reuse-cell-counts-nc /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50.nc4 --reuse-scope-meta-csv /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out/ims_ol_da_scope_metadata_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50.csv \
  --min-ims-snow-days 10 --output-dir /discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/out --output-tag SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10 --overwrite
