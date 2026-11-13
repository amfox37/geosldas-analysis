#!/bin/bash
#SBATCH --job-name=iv_tc_az_step2
#SBATCH --account=s3208
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-2
#SBATCH --output=iv_tc_az_step2.%A_%a.out
#SBATCH --error=iv_tc_az_step2.%A_%a.err

# Step 2: daily obs/model pairs for smap_l3 and ascat_h121, one array task per arm
# (0=OL, 1=L1_filt275, 2=L3; see iv_tc_az_arms.sh).
# Submit from this jobs/ directory so the logs land here:
#   sbatch run_iv_tc_az_step2.sbatch
#   sbatch --dependency=afterok:<jobid> run_iv_tc_az_step34.sbatch
# Resumable: existing valid pair files are skipped (no --overwrite).

set -euo pipefail

# sbatch runs a spooled copy of this file, so source the settings by absolute path.
source "${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}/projects/CYGNSS_L1_AZ/jobs/iv_tc_az_arms.sh"

cd "${REPO_ROOT}"
module purge
module load python/GEOSpyD

echo "arm=${ARM_NAME} run_dir=${ARM_DIR} ${START_DATE}..${END_DATE} -> ${OUT_ROOT}"
mkdir -p "${OUT_ROOT}"

python -u projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor smap_l3 \
  --sensor ascat_h121 \
  --run "${ARM_NAME}=${ARM_DIR}" \
  --output-root "${OUT_ROOT}"
