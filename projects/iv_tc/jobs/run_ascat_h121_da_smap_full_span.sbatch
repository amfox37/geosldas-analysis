#!/bin/bash
#SBATCH --job-name=ascat_h121_da_smap_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=ascat_h121_da_smap_full_span.%j.out
#SBATCH --error=ascat_h121_da_smap_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"
DA_RUN="${DA_RUN:-/discover/nobackup/amfox/Experiments/DAv7_M36_SMAP_type_13_comb_fp_scaled/DAv7_M36_SMAP_type_13_comb_fp_scaled}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

python -u projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor ascat_h121 \
  --run "DA_SMAP_comb_fp_scaled=${DA_RUN}" \
  --collection ".SMAP_L4_SM_gph." \
  --model-variable sm_surface \
  --model-time-reduce mean \
  --output-root "${OUTPUT_ROOT}"
