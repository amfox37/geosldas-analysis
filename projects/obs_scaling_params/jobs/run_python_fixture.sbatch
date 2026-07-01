#!/bin/bash
#SBATCH --job-name=obs_scaling_python_fixture
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=obs_scaling_python_fixture.%j.out
#SBATCH --error=obs_scaling_python_fixture.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
FIXTURE_ROOT="${FIXTURE_ROOT:-/discover/nobackup/${USER}/obs_scaling_fixture}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/prepare_fixture_case.py --fixture-root "${FIXTURE_ROOT}"
python -u scripts/run_scaling_params.py \
  --exp-path "${FIXTURE_ROOT}" \
  --exp-run LS_DAv8_M36_as_test2 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2020-05 \
  --end 2020-05 \
  --run-months 5 \
  --species ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM \
  --grid-resolution 1.0 \
  --window-days 1 \
  --ndata-min 1 \
  --dt-assim 86400 \
  --t0-assim 32400 \
  --prefix fixture_python_zscore_stats_ \
  --out-dir python_fixture_zscore_1deg \
  --no-all-pentads
