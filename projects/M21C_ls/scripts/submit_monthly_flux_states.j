#!/usr/bin/env bash
#SBATCH --job-name=monthly_flux_states
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=monthly_flux_states.%j.out
#SBATCH --error=monthly_flux_states.%j.err

module purge
module load python/GEOSpyD/25.3.1-0/3.13

REPO="/discover/nobackup/projects/land_da/geosldas-analysis"

REPO="${REPO}" \
PYTHON="python3" \
OUTPUT_DIR="${REPO}/projects/M21C_ls/output/monthly_flux_states" \
bash "${REPO}/projects/M21C_ls/scripts/build_monthly_flux_states_discover.sh"
