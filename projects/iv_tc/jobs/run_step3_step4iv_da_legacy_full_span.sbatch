#!/bin/bash
#SBATCH --job-name=iv_tc_step34_da_legacy_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=iv_tc_step34_da_legacy_full_span.%j.out
#SBATCH --error=iv_tc_step34_da_legacy_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# Full DA_legacy span (2015-04-01..2021-03-31), the 4 sensors with
# (near-)complete coverage over that whole window.
echo "=== step3 climatology: smosic, smap_l3, ascat_h121, ascat_h119_h120 (2015-04-01..2021-03-31) ==="
python -u projects/iv_tc/scripts/compute_climatology.py \
  --start-date 2015-04-01 \
  --end-date 2021-03-31 \
  --sensor smosic \
  --sensor smap_l3 \
  --sensor ascat_h121 \
  --sensor ascat_h119_h120 \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --output-root "${ROOT}"

# CYGNSS only exists from 2018-08-01 onward -- its own window/Nday_min.
echo "=== step3 climatology: cygnss_l3 (2018-08-01..2021-03-31) ==="
python -u projects/iv_tc/scripts/compute_climatology.py \
  --start-date 2018-08-01 \
  --end-date 2021-03-31 \
  --sensor cygnss_l3 \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --output-root "${ROOT}"

echo "=== step4 IV: smosic, smap_l3, ascat_h121, ascat_h119_h120 (2015-04-01..2021-03-31) ==="
python -u projects/iv_tc/scripts/compute_iv.py \
  --start-date 2015-04-01 \
  --end-date 2021-03-31 \
  --sensor smosic \
  --sensor smap_l3 \
  --sensor ascat_h121 \
  --sensor ascat_h119_h120 \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --climatology-root "${ROOT}" \
  --output-root "${ROOT}"

echo "=== step4 IV: cygnss_l3 (2018-08-01..2021-03-31) ==="
python -u projects/iv_tc/scripts/compute_iv.py \
  --start-date 2018-08-01 \
  --end-date 2021-03-31 \
  --sensor cygnss_l3 \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --climatology-root "${ROOT}" \
  --output-root "${ROOT}"
