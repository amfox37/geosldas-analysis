#!/bin/bash
#SBATCH --job-name=obs_scaling_params
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --account=g0610
#SBATCH --time=07:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --constraint=mil

# Run run_obs_scaling_params.py on a compute node; all arguments are passed through, e.g.
#
#   sbatch jobs/run_obs_scaling_params.sbatch --method latlon --exp-path ... --exp-run ... \
#       --domain ... --start 2019-02 --end 2020-08 --species ... --prefix ...

set -euo pipefail

export OMP_NUM_THREADS=1

module purge
module load python/GEOSpyD

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
cd "${REPO_ROOT}/projects/obs_scaling_params"

python -u run_obs_scaling_params.py "$@" | grep -v '^Processing'
