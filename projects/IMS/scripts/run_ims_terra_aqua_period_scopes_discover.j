#!/bin/bash
#SBATCH --job-name=ims_terra_aqua_scopes
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --account=g0610
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module purge
module load python/GEOSpyD

REPO_ROOT=/discover/nobackup/projects/land_da/geosldas-analysis

srun -n 1 bash "${REPO_ROOT}/projects/IMS/scripts/run_ims_terra_aqua_period_scopes_discover.sh"

exit 0
