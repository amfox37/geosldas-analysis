#!/bin/bash
#SBATCH --job-name=tc_cygnss_smosic_da_legacy_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_cygnss_smosic_da_legacy_full_span.%j.out
#SBATCH --error=tc_cygnss_smosic_da_legacy_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2018-08-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# DA_legacy completes the cygnss_l3+smosic combo already run for OL/DA_H121
# (run_tc_cygnss_smosic_full_span.sbatch). DA_SMAP is intentionally excluded:
# it assimilates SMAP, and smap_l3/cygnss_l3 are treated as circular
# verification obs for that run elsewhere in this project. DA_legacy
# assimilates the legacy EUMETSAT ASCAT BUFR product, which is independent
# of both cygnss_l3 and smosic.
#
# smosic's existing step3 climatology only covers the full 2015-04-01 span;
# cygnss_l3 (primary sensor below) only exists from 2018-08-01 onward, so
# build a matching-window smosic climatology for DA_legacy (mirrors the
# OL/DA_H121 climatology already built by run_tc_cygnss_smosic_full_span.sbatch).
echo "=== step3 climatology: smosic on cygnss_l3 window (${START_DATE}..${END_DATE}) ==="
python -u projects/iv_tc/scripts/compute_climatology.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --sensor smosic \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --output-root "${ROOT}"

echo "=== step4 TC: cygnss_l3 + model + smosic, DA_legacy (${START_DATE}..${END_DATE}) ==="
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor cygnss_l3 \
  --secondary-sensor smosic \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --allow-missing
