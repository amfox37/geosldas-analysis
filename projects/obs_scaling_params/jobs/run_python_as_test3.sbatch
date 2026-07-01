#!/bin/bash
#SBATCH --job-name=obs_scaling_python_as_test3
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=obs_scaling_python_as_test3.%j.out
#SBATCH --error=obs_scaling_python_as_test3.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/discover/nobackup/${USER}/obs_scaling_as_test3}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/run_scaling_params.py \
  --exp-path "${SCRATCH_ROOT}" \
  --exp-run LS_DAv8_M36_as_test3 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2020-05 \
  --end 2020-05 \
  --run-months 5 \
  --species SMOS_fit_Tbh_A,SMOS_fit_Tbh_D,SMOS_fit_Tbv_A,SMOS_fit_Tbv_D,SMAP_L1C_Tbh_A,SMAP_L1C_Tbh_D,SMAP_L1C_Tbv_A,SMAP_L1C_Tbv_D,ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM,MYD10C1,MOD10C1,CYGNSS_SM_6hr \
  --grid-resolution 1.0 \
  --window-days 3 \
  --ndata-min 1 \
  --dt-assim 10800 \
  --t0-assim 0 \
  --prefix as_test3_python_zscore_stats_ \
  --out-dir python_as_test3_zscore_1deg \
  --obsfcstana-format nc4 \
  --no-all-pentads
