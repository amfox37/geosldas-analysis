#!/bin/bash
#SBATCH --job-name=obs_scaling_matlab_as_test3
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=obs_scaling_matlab_as_test3.%j.out
#SBATCH --error=obs_scaling_matlab_as_test3.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/discover/nobackup/${USER}/obs_scaling_as_test3}"
MATLAB_MODULE="${MATLAB_MODULE:-matlab}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load "${MATLAB_MODULE}"
export OBS_SCALING_REPO_ROOT="${REPO_ROOT}"
export OBS_SCALING_SCRATCH_ROOT="${SCRATCH_ROOT}"

matlab -batch "run('${REPO_ROOT}/projects/obs_scaling_params/matlab_reference/run_as_test3_case.m')"
