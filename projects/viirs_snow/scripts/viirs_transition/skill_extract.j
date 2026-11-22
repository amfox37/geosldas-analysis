#!/usr/bin/env bash
#SBATCH --job-name=skill_ext
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:12:00
#SBATCH --array=0-4
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/skill_data/skill_ext.%A_%a.out
# One task per arm (arms.json order). Timing test: 30 days = 6 s -> ~7.5 min for a 6-year arm + ~100 s station->tile mapping.
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS=($("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms']))"))
"${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/skill_extract.py" "${ARMS[${SLURM_ARRAY_TASK_ID}]}"
