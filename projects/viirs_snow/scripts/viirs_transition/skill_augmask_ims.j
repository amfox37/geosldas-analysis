#!/usr/bin/env bash
#SBATCH --job-name=augims
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:06:00
#SBATCH --array=0-3
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/skill_data/augmask_ims.%A_%a.out
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS=(T_to_J2 TA_to_J1J2 T_control TA_control)
"${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/skill_augmask_ims.py" "${ARMS[${SLURM_ARRAY_TASK_ID}]}"
