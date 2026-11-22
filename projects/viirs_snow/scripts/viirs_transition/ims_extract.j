#!/usr/bin/env bash
#SBATCH --job-name=ims_ext
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:05:00
#SBATCH --array=0-3
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/skill_data/imsfix/ims_ext.%A_%a.out
# One task per arm. 200-day test = 103 s, but 1096 days took 1518 s under load (1.4 s/day) -> ~51 min for a 6-year arm; first run with 25 min timed out.
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS=($("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms']))"))
"${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/ims_extract.py" "${ARMS[${SLURM_ARRAY_TASK_ID}]}"
