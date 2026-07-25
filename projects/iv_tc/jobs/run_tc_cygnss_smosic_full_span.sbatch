#!/bin/bash
#SBATCH --job-name=tc_cygnss_smosic_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_cygnss_smosic_full_span.%j.out
#SBATCH --error=tc_cygnss_smosic_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2018-08-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# smosic's existing step3 climatology only covers the full 2015-04-01
# span; cygnss_l3 (primary sensor below) only exists from 2018-08-01
# onward, so its climatology filename won't match. Build a smosic
# climatology on cygnss_l3's own window so both source climatology
# files share the same start/end date.
echo "=== step3 climatology: smosic on cygnss_l3 window (${START_DATE}..${END_DATE}) ==="
python -u projects/iv_tc/scripts/compute_climatology.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor smosic \
  --run OL \
  --run DA_H121 \
  --pair-root "${ROOT}" \
  --output-root "${ROOT}"

# cygnss_l3 (primary) + smosic (secondary) + model, annual.
echo "=== step4 TC: cygnss_l3 + model + smosic (${START_DATE}..${END_DATE}) ==="
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor cygnss_l3 \
  --secondary-sensor smosic \
  --run OL \
  --run DA_H121 \
  --pair-root "${ROOT}" \
  --allow-missing
