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

# --- good HPC hygiene: keep BLAS/OpenMP from oversubscribing cores ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load the Discover Python stack
module purge
module load python/GEOSpyD

# Go to the script directory
cd /discover/nobackup/projects/land_da/GEOSldas_CYGNSS/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSldas_App/util/inputs/obs_scaling_params/python

# (Optional) activate a conda env if you have one inside GEOSpyD
# source activate geospyd

# Run the script (use absolute paths for any inputs/outputs if the script expects them)
srun -n 1 python -u run_get_model_and_obs_clim_stats_latlon_grid.py
# Example with args:
# srun -n 1 python -u run_get_model_and_obs_clim_stats_latlon_grid.py --years 2000-2024 --grid M36

exit 0
