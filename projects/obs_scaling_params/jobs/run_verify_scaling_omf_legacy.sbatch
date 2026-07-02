#!/bin/bash
#SBATCH --job-name=verify_omf_legacy
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=verify_scaling_omf_legacy.%j.out
#SBATCH --error=verify_scaling_omf_legacy.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
STATS_DIR="/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121/output/SMAP_EASEv2_M36_GLOBAL/stats"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/verify_scaling_omf.py \
  --stats-file "${STATS_DIR}/legacy_ASCAT_python_zscore_0p25deg/OLv7_M36_legacy_ASCAT_zscore_stats_2015_doy91_2021_doy120_W_75d_Nmin_20_sp_ALL_all_pentads.nc4" \
  --species 5,6,7 \
  --species-names "ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM" \
  --start-year 2017 --end-year 2020
