#!/usr/bin/env bash
#SBATCH --job-name=obs_avail
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --array=0-4
#SBATCH --output=obs_avail.%A_%a.out
#SBATCH --error=obs_avail.%A_%a.err

# One array task per arm (order of arms.json). 6-year arms took ~10 min in a one-month timing test
# (8 s per arm-month x 72), the two controls (3 years) about half. The first full run (job 59362332) took 11:11-11:24
# for the 6-year arms against a 12 min limit, hence 15 min now. Submit with:
#   cd <repo>/projects/viirs_snow/scripts/viirs_transition && sbatch obs_availability.j
# (the .out/.err files land in the submit directory and are not versioned; delete them when done)

REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
ARMS_ALL=($("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms']))"))

REPO="${REPO}" PYTHON="${PYTHON}" ARMS="${ARMS_ALL[${SLURM_ARRAY_TASK_ID}]}" \
  bash "${REPO}/projects/viirs_snow/scripts/viirs_transition/obs_availability_discover.sh"
