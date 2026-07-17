#!/bin/bash
#SBATCH --job-name=ascat_h121_ol_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=ascat_h121_ol_full_span.%j.out
#SBATCH --error=ascat_h121_ol_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"
OL_RUN="${OL_RUN:-/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

python -u projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor ascat_h121 \
  --run "OL=${OL_RUN}" \
  --output-root "${OUTPUT_ROOT}"
