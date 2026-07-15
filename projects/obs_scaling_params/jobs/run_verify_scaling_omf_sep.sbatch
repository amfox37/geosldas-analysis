#!/bin/bash
#SBATCH --job-name=verify_omf_sep
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=verify_scaling_omf_sep.%j.out
#SBATCH --error=verify_scaling_omf_sep.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
STATS_DIR="/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121/output/SMAP_EASEv2_M36_GLOBAL/stats/H121_ASCAT_python_zscore_0p25deg"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD

for entry in "8:ASCAT_HSAF_META_SM" "9:ASCAT_HSAF_METB_SM" "10:ASCAT_HSAF_METC_SM"; do
  sp="${entry%%:*}"
  name="${entry#*:}"
  echo "=== ${name} (species ${sp}) ==="
  python -u scripts/verify_scaling_omf.py \
    --stats-file "${STATS_DIR}/OLv7_M36_H121_ASCAT_zscore_stats_2015_doy91_2021_doy120_W_75d_Nmin_20_sp_${name}_all_pentads.nc4" \
    --species "${sp}" \
    --species-names "${name}" \
    --start-year 2017 --end-year 2020
done
