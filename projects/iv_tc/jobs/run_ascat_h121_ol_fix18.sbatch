#!/bin/bash
#SBATCH --job-name=ascat_h121_ol_fix18
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=ascat_h121_ol_fix18.%j.out
#SBATCH --error=ascat_h121_ol_fix18.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
OL_RUN="${OL_RUN:-/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121}"

# The 18 days previously blocked by corrupt H121 source files (masked behind
# the first 35 found by job 57189168, fixed via redownload_h121_corrupt_batch2.sh)
DATES=(
  2015-07-05
  2015-07-06
  2016-07-10
  2016-07-11
  2017-06-25
  2017-06-26
  2018-06-27
  2018-10-21
  2018-10-22
  2018-10-23
  2019-06-25
  2019-06-26
  2019-09-01
  2019-10-22
  2019-10-24
  2019-10-25
  2020-09-04
  2020-09-05
)

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

for d in "${DATES[@]}"; do
  python -u projects/iv_tc/scripts/generate_daily_pairs.py \
    --start-date "${d}" \
    --end-date "${d}" \
    --sensor ascat_h121 \
    --run "OL=${OL_RUN}" \
    --output-root "${OUTPUT_ROOT}"
done
