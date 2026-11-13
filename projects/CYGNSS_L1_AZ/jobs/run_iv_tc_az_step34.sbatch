#!/bin/bash
#SBATCH --job-name=iv_tc_az_step34
#SBATCH --account=s3208
#SBATCH --partition=compute
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --output=iv_tc_az_step34.%A_%a.out
#SBATCH --error=iv_tc_az_step34.%A_%a.err

# Steps 3-4 for one arm per array task: pentad climatology, IV (smap_l3, ascat_h121),
# and TC (primary smap_l3, secondary ascat_h121, annual = full period).
# Per-tile results are kept as written; the 552-tile L1-domain mask and any
# domain means are applied afterwards in post-processing, not here.
# Run after run_iv_tc_az_step2.sbatch finishes (--dependency=afterok:<jobid>).

set -euo pipefail

# sbatch runs a spooled copy of this file, so source the settings by absolute path.
source "${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}/projects/CYGNSS_L1_AZ/jobs/iv_tc_az_arms.sh"

cd "${REPO_ROOT}"
module purge
module load python/GEOSpyD

echo "arm=${ARM_NAME} ${START_DATE}..${END_DATE} root=${OUT_ROOT}"

echo "=== step3 climatology ==="
python -u projects/iv_tc/scripts/compute_climatology.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor smap_l3 \
  --sensor ascat_h121 \
  --run "${ARM_NAME}" \
  --pair-root "${OUT_ROOT}" \
  --output-root "${OUT_ROOT}"

echo "=== step4 IV ==="
python -u projects/iv_tc/scripts/compute_iv.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor smap_l3 \
  --sensor ascat_h121 \
  --run "${ARM_NAME}" \
  --pair-root "${OUT_ROOT}" \
  --climatology-root "${OUT_ROOT}" \
  --output-root "${OUT_ROOT}"

echo "=== step4 TC: smap_l3 + model + ascat_h121 (min-count 20, annual) ==="
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor smap_l3 \
  --secondary-sensor ascat_h121 \
  --run "${ARM_NAME}" \
  --pair-root "${OUT_ROOT}" \
  --climatology-root "${OUT_ROOT}" \
  --output-root "${OUT_ROOT}" \
  --allow-missing
