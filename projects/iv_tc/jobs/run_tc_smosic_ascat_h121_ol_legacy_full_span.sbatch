#!/bin/bash
#SBATCH --job-name=tc_smosic_ascat_h121_ol_legacy_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_smosic_ascat_h121_ol_legacy_full_span.%j.out
#SBATCH --error=tc_smosic_ascat_h121_ol_legacy_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# OL and DA_legacy complete the smosic+ascat_h121 combo already run for
# DA_SMAP (run_tc_smosic_ascat_h121_da_smap_full_span.sbatch). OL has
# nothing assimilated, so it's a valid baseline for any pair. DA_legacy
# assimilates the legacy EUMETSAT ASCAT BUFR product (h119/h120), a
# different sensor than ascat_h121 (H121 CDR), so it's independent of
# both smosic and ascat_h121 here. DA_H121 is excluded: it assimilates
# ascat_h121 directly, which would make this pair circular for that run.
# Both sensors already have full-span climatology for OL and DA_legacy, so
# no extra climatology step is needed here.
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor smosic \
  --secondary-sensor ascat_h121 \
  --run OL \
  --run DA_legacy \
  --pair-root "${ROOT}" \
  --allow-missing
