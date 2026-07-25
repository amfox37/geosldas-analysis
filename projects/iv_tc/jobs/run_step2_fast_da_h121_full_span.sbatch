#!/bin/bash
#SBATCH --job-name=step2_fast_da_h121_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=step2_fast_da_h121_full_span.%j.out
#SBATCH --error=step2_fast_da_h121_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
DA_RUN="${DA_RUN:-/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_H121}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# smosic + smap_l3 over the full OL span; cygnss_l3 only exists from
# 2018-08-01 onward.
python -u projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date 2015-04-01 \
  --end-date 2021-03-31 \
  --sensor smosic \
  --sensor smap_l3 \
  --run "DA_H121=${DA_RUN}" \
  --output-root "${OUTPUT_ROOT}"

python -u projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date 2018-08-01 \
  --end-date 2021-03-31 \
  --sensor cygnss_l3 \
  --run "DA_H121=${DA_RUN}" \
  --output-root "${OUTPUT_ROOT}"
