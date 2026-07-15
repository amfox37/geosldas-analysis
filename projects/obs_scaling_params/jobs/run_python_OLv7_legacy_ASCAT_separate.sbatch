#!/bin/bash
#SBATCH --job-name=obs_scaling_legacy_ASCAT_sep
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=obs_scaling_legacy_ASCAT_sep.%j.out
#SBATCH --error=obs_scaling_legacy_ASCAT_sep.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/run_scaling_params.py \
  --exp-path /discover/nobackup/projects/land_da/hsaf_cdr_test \
  --exp-run OLv7_M36_MULTI_type_13_H121 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2015-04 \
  --end 2021-04 \
  --species ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM \
  --no-combine-species \
  --grid-resolution 0.25 \
  --window-days 75 \
  --ndata-min 20 \
  --obsfcstana-format nc4 \
  --prefix OLv7_M36_legacy_ASCAT_zscore_stats_ \
  --out-dir legacy_ASCAT_python_zscore_0p25deg
