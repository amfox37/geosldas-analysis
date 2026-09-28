#!/bin/bash
#SBATCH --job-name=ismn_ol_da_skill_az
#SBATCH --account=s3208
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=ismn_ol_da_skill_az.%j.out
#SBATCH --error=ismn_ol_da_skill_az.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
START_DATE="${START_DATE:-2020-01-01}"
END_DATE="${END_DATE:-2022-12-31}"
# Space-separated NAME=EXP_ID pairs (dirs under CYGL1_ROOT); the one named OL is the reference.
# Default = post-fix fixed-operator OL + the four arms run through 2022.
RUNS="${RUNS:-OL=OLv8_M36_AZ_fixedop L1coh=DA_L1_full_xc015_coh040216_fixedop L1coh_err39=DA_L1_full_xc015_coh040216_err39_fixedop L3=DA_L3_fixedop SMAP=DA_SMAP_fixedop}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/projects/CYGNSS_L1_AZ/output/ismn_fixedop_${START_DATE//-/}_${END_DATE//-/}}"

CYGL1_ROOT="/discover/nobackup/projects/land_da/cygl1_operator_test"

cd "${REPO_ROOT}"

module purge
module load python/GEOSpyD

# Reuses the ascat_da ISMN skill driver (network-agnostic ISMN reader, any
# GEOSldas domain/collection/variable pair via CLI flags -- no project-specific
# code needed). It points it at the AZ-domain runs in RUNS and reads the tavg24_1d_lnd_Nt daily collection
# (single file/day, SFMC/RZMC) instead of ascat_da's SMAP_L4_SM_gph. --nmin is
# lowered from the 1000-day default (tuned for multi-year hsaf_cdr_test runs).
# 2026-09-09 pre-fix run (dense075_coh05, 2020-2021, job 58298184) used
# RUNS="OL=OLv8_M36_all_sensors_AZ DA_dense075_coh05=DAv8_M36_AZ_paired_cygl1_dense075_coh05" END_DATE=2021-12-31.
RUN_ARGS=()
for r in ${RUNS}; do RUN_ARGS+=(--run "${r%%=*}=${CYGL1_ROOT}/${r#*=}"); done
python -u projects/ascat_da/scripts/run_ismn_ol_da_skill.py \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --collection ".tavg24_1d_lnd_Nt." \
  --model-surface-variable SFMC \
  --model-rz-variable RZMC \
  "${RUN_ARGS[@]}" \
  --reference-run OL \
  --nmin 90 \
  --output-dir "${OUTPUT_DIR}"
