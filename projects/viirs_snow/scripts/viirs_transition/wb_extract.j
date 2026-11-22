#!/usr/bin/env bash
#SBATCH --job-name=wb_extract
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --array=0-4
#SBATCH --output=wb_extract.%A_%a.out
#SBATCH --error=wb_extract.%A_%a.err
# Analysis D step 1: one array task per arm (order of arms.json). Measured ~0.5 s per day (10-day test), so ~18 min
# for the 6-year arms. Submit: cd <repo>/projects/viirs_snow/scripts/viirs_transition && sbatch wb_extract.j
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS_ALL=($("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms']))"))
ARM="${ARMS_ALL[${SLURM_ARRAY_TASK_ID}]}"
PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/wb_extract.py" \
  --arm "${ARM}" --out "${REPO}/projects/viirs_snow/output/water_balance/wb_${ARM}.nc"
