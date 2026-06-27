#!/usr/bin/env bash
#SBATCH --job-name=raw_catch_incr
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=raw_catch_incr.%j.out
#SBATCH --error=raw_catch_incr.%j.err

module purge
module load python/GEOSpyD/25.3.1-0/3.13

REPO="/discover/nobackup/projects/land_da/geosldas-analysis"

REPO="${REPO}" \
PYTHON="python3" \
OUTPUT_DIR="${REPO}/projects/M21C_ls/output/raw_increment_summaries" \
bash "${REPO}/projects/M21C_ls/scripts/build_raw_catch_progn_increment_summaries_discover.sh"
