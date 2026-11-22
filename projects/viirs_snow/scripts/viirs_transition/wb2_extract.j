#!/usr/bin/env bash
#SBATCH --job-name=wb2_extract
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --time=00:08:00
#SBATCH --array=0-4
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/wb2_logs/wb2_extract.%A_%a.out
#SBATCH --error=/discover/nobackup/projects/land_da/VIIRS_testing/wb2_logs/wb2_extract.%A_%a.err
# 34 days x ~1 s per arm. Submit: sbatch wb2_extract.j
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS_ALL=($("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms']))"))
ARM="${ARMS_ALL[${SLURM_ARRAY_TASK_ID}]}"
PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/wb2_extract.py" \
  --arm "${ARM}" --out "${REPO}/projects/viirs_snow/output/water_balance/wb2_gap_${ARM}.nc"
