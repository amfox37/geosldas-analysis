#!/bin/bash
#SBATCH --job-name=obs_scaling_matlab_fixture
#SBATCH --account=s2519
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=obs_scaling_matlab_fixture.%j.out
#SBATCH --error=obs_scaling_matlab_fixture.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
FIXTURE_ROOT="${FIXTURE_ROOT:-/discover/nobackup/${USER}/obs_scaling_fixture}"
MATLAB_MODULE="${MATLAB_MODULE:-matlab}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/prepare_fixture_case.py --fixture-root "${FIXTURE_ROOT}"

module load "${MATLAB_MODULE}"
export OBS_SCALING_REPO_ROOT="${REPO_ROOT}"
export OBS_SCALING_FIXTURE_ROOT="${FIXTURE_ROOT}"

matlab -batch "run('${REPO_ROOT}/projects/obs_scaling_params/matlab_reference/run_fixture_case.m')"
