#!/usr/bin/env bash
#SBATCH --job-name=incr_stats
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --array=0-4
#SBATCH --output=incr_stats.%A_%a.out
#SBATCH --error=incr_stats.%A_%a.err

# Analysis B: one array task per arm (order of arms.json); runs incr_stats.py (catch_progn_incr, ~0.3 s/day) then
# incr_innov.py (ObsFcstAna, ~0.4-0.5 s/day in a 3-day test). Limit is a first guess: 2191 days x ~0.8 s = ~30 min for
# the 6-year arms before measuring; refine from the first run's elapsed time. Submit with:
#   cd <repo>/projects/viirs_snow/scripts/viirs_transition && sbatch incr_stats.j
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
HERE="${REPO}/projects/viirs_snow/scripts/viirs_transition"
ARM=$("${PYTHON}" -I -c "import json; print(list(json.load(open('${REPO}/projects/viirs_snow/config/arms.json'))['arms'])[${SLURM_ARRAY_TASK_ID}])")
export PYTHONDONTWRITEBYTECODE=1
"${PYTHON}" -I "${HERE}/incr_stats.py" --arm "${ARM}"
"${PYTHON}" -I "${HERE}/incr_innov.py" --arm "${ARM}"
