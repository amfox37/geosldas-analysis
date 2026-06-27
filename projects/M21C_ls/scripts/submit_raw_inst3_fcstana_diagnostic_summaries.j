#!/usr/bin/env bash
#SBATCH --job-name=raw_inst3_diag
#SBATCH --account=g0610
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=raw_inst3_diag.%j.out
#SBATCH --error=raw_inst3_diag.%j.err

module purge
module load python/GEOSpyD/25.3.1-0/3.13

REPO="/discover/nobackup/projects/land_da/geosldas-analysis"

REPO="${REPO}" \
PYTHON="python3" \
OUTPUT_DIR="${REPO}/projects/M21C_ls/output/raw_increment_summaries" \
bash "${REPO}/projects/M21C_ls/scripts/build_raw_inst3_fcstana_diagnostic_summaries_discover.sh"
