#!/usr/bin/env bash
#SBATCH --job-name=tr_extract
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --time=00:08:00
#SBATCH --array=2023-2025
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/tr_logs/tr_extract.%A_%a.out
#SBATCH --error=/discover/nobackup/projects/land_da/VIIRS_testing/tr_logs/tr_extract.%A_%a.err
# Analysis E step 1: one array task per year (SLURM_ARRAY_TASK_ID = year). Test: 10 days x 3 arms = 3.8 s incl. startup; est. 1-3 min per year.
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/tr_extract.py" --year "${SLURM_ARRAY_TASK_ID}"
