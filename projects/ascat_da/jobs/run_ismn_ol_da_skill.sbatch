#!/bin/bash
#SBATCH --job-name=ismn_ol_da_skill
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=ismn_ol_da_skill.%j.out
#SBATCH --error=ismn_ol_da_skill.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_DIR="${OUTPUT_DIR:-/discover/nobackup/projects/land_da/Evaluation/ISMN/ascat_da_ol_da_skill}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# Every ISMN network in the local archive (no predetermined network list),
# one window spanning the full experiment period. All four runs are read from
# SMAP_L4_SM_gph sm_surface/sm_rootzone, daily-averaged over each file's eight
# 3-hourly instants -- the analysis state, and the only collection the
# DA_SMAP_comb_fp_scaled run carries.
python -u projects/ascat_da/scripts/run_ismn_ol_da_skill.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --run "OL=/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121" \
  --run "DA_H121=/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_H121" \
  --run "DA_legacy=/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_legacy" \
  --run "DA_SMAP_comb_fp_scaled=/discover/nobackup/amfox/Experiments/DAv7_M36_SMAP_type_13_comb_fp_scaled/DAv7_M36_SMAP_type_13_comb_fp_scaled" \
  --reference-run OL \
  --output-dir "${OUTPUT_DIR}"
