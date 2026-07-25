#!/bin/bash
#SBATCH --job-name=tc_smap_smosic_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_smap_smosic_full_span.%j.out
#SBATCH --error=tc_smap_smosic_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# smap_l3 (primary) + smosic (secondary) + model, annual, full OL/DA_H121
# span. Same triplet already spot-checked live against MATLAB
# (2018-08-01..2019-07-31): N_sm 0 mismatches, R2/correlation/covariance
# clean at 1e-4; sigma2's near-zero-denominator noise is float32 precision,
# not a coverage or engine defect. --allow-missing covers smap_l3's known
# 37-day safe-mode gap (2019-06-20..07-22) and smosic's 1-day gap.
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor smap_l3 \
  --secondary-sensor smosic \
  --run OL \
  --run DA_H121 \
  --pair-root "${ROOT}" \
  --allow-missing
