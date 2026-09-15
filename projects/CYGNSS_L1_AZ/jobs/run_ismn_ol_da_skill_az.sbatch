#!/bin/bash
#SBATCH --job-name=ismn_ol_da_skill_az
#SBATCH --account=s3208
#SBATCH --partition=compute
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=ismn_ol_da_skill_az.%j.out
#SBATCH --error=ismn_ol_da_skill_az.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/projects/CYGNSS_L1_AZ/output/ismn_dense075_coh05_ol_da_skill}"
START_DATE="${START_DATE:-2020-01-01}"
END_DATE="${END_DATE:-2021-12-31}"

CYGL1_ROOT="/discover/nobackup/projects/land_da/cygl1_operator_test"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# Reuses the ascat_da ISMN skill driver (network-agnostic ISMN reader, any
# GEOSldas domain/collection/variable pair via CLI flags -- no project-specific
# code needed). This run points it at the AZ-domain OL and dense075_coh05
# CYGNSS L1 DA experiments and reads the tavg24_1d_lnd_Nt daily collection
# (single file/day, SFMC/RZMC) instead of ascat_da's SMAP_L4_SM_gph. --nmin is
# lowered from the 1000-day default (tuned for multi-year hsaf_cdr_test runs)
# since this is a 24-month window.
python -u projects/ascat_da/scripts/run_ismn_ol_da_skill.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --collection ".tavg24_1d_lnd_Nt." \
  --model-surface-variable SFMC \
  --model-rz-variable RZMC \
  --run "OL=${CYGL1_ROOT}/OLv8_M36_all_sensors_AZ" \
  --run "DA_dense075_coh05=${CYGL1_ROOT}/DAv8_M36_AZ_paired_cygl1_dense075_coh05" \
  --reference-run OL \
  --nmin 90 \
  --output-dir "${OUTPUT_DIR}"
