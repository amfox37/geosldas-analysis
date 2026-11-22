#!/usr/bin/env bash
#SBATCH --job-name=ims_score
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:06:00
#SBATCH --output=/discover/nobackup/projects/land_da/VIIRS_testing/skill_data/imsfix/ims_score.%j.out
REPO="/discover/nobackup/projects/land_da/geosldas-analysis"
PYTHON="/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3"
cd "${REPO}/projects/viirs_snow/scripts/viirs_transition" && "${PYTHON}" -I ims_score.py
