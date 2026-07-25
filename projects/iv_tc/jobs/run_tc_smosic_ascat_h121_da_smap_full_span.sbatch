#!/bin/bash
#SBATCH --job-name=tc_smosic_ascat_h121_da_smap_full_span
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=tc_smosic_ascat_h121_da_smap_full_span.%j.out
#SBATCH --error=tc_smosic_ascat_h121_da_smap_full_span.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
ROOT="${ROOT:-/discover/nobackup/projects/land_da/Evaluation/IVs/output_python}"
START_DATE="${START_DATE:-2015-04-01}"
END_DATE="${END_DATE:-2021-03-31}"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# smosic (primary) + ascat_h121 (secondary) + model, annual,
# DA_SMAP_comb_fp_scaled. This run assimilates SMAP, so smap_l3/cygnss_l3
# are excluded as verification obs (circularity); smosic + ascat_h121 are
# both independent of the assimilated species. Both sensors already get
# full-span climatology from run_step3_step4iv_da_smap_full_span.sbatch.
python -u projects/iv_tc/scripts/compute_tc.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --primary-sensor smosic \
  --secondary-sensor ascat_h121 \
  --run DA_SMAP_comb_fp_scaled \
  --pair-root "${ROOT}" \
  --allow-missing
