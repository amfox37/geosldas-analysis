#!/usr/bin/env bash
#SBATCH --job-name=wb2_analyze
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/wb2_logs/wb2_analyze.%j.out
#SBATCH --error=/discover/nobackup/projects/land_da/VIIRS_testing/wb2_logs/wb2_analyze.%j.err
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -I "${REPO}/projects/viirs_snow/scripts/viirs_transition/wb2_analyze.py"
