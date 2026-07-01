#!/bin/bash
#SBATCH --job-name=scaling_params
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --account=g0610
#SBATCH --time=07:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --constraint=mil

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load python/GEOSpyD

# Submit from the repo root, or set REPO_ROOT explicitly in the sbatch command.
REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
cd "${REPO_ROOT}/projects/obs_scaling_params"

srun -n 1 python -u scripts/run_scaling_params.py "$@"
