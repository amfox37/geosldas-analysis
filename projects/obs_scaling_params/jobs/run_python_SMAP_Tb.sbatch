#!/bin/bash
#SBATCH --job-name=obs_scaling_SMAP_Tb
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=obs_scaling_SMAP_Tb.%j.out
#SBATCH --error=obs_scaling_SMAP_Tb.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
EXP_PATH="${EXP_PATH:?Set EXP_PATH to the parent directory containing the experiment}"
EXP_RUN="${EXP_RUN:?Set EXP_RUN to the GEOS-LDAS experiment id}"
START_MONTH="${START_MONTH:?Set START_MONTH as YYYY-MM}"
END_MONTH="${END_MONTH:?Set END_MONTH as YYYY-MM}"
PREFIX="${PREFIX:?Set PREFIX to the scaling filename prefix}"
DOMAIN="${DOMAIN:-SMAP_EASEv2_M09_GLOBAL}"
ORBIT="${ORBIT:-D}"
OFA_FORMAT="${OFA_FORMAT:-bin}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/run_tb_scaling_params.py \
  --exp-path "${EXP_PATH}" \
  --exp-run "${EXP_RUN}" \
  --domain "${DOMAIN}" \
  --start "${START_MONTH}" \
  --end "${END_MONTH}" \
  --orbit "${ORBIT}" \
  --obsfcstana-format "${OFA_FORMAT}" \
  --prefix "${PREFIX}"
