#!/bin/bash
#SBATCH --job-name=tc_smap_smosic_da_legacy_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_smap_smosic_da_legacy_full_span.%j.out
#SBATCH --error=tc_smap_smosic_da_legacy_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# DA_legacy completes the smap_l3+smosic combo already run for OL/DA_H121
# (run_tc_smap_smosic_full_span.sbatch). DA_SMAP is intentionally excluded:
# it assimilates SMAP, and smap_l3 is a circular verification obs for that
# run. DA_legacy assimilates the legacy EUMETSAT ASCAT BUFR product, which
# is independent of both smap_l3 and smosic. Both sensors already have
# full-span climatology for DA_legacy (built by
# run_step3_step4iv_da_legacy_full_span.sbatch), so no extra climatology
# step is needed here.
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor smap_l3 \
  --secondary-sensor smosic \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --allow-missing
