#!/usr/bin/env bash
set -euo pipefail

# Run on Discover after pulling geosldas-analysis.
# Builds monthly diagnostic ANA-FCST burden metrics from raw daily inst3 files.

REPO="${REPO:-/discover/nobackup/projects/land_da/geosldas-analysis}"
PYTHON="${PYTHON:-python3}"
DOMAIN="${DOMAIN:-SMAP_EASEv2_M36_GLOBAL}"
RUN_BASE="${RUN_BASE:-/discover/nobackup/projects/land_da/M21C_land_sweeper}"

DA_ROOT="${DA_ROOT:-${RUN_BASE}/LS_DAv8_M36_v3/LS_DAv8_M36/output/${DOMAIN}/cat/ens_avg}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO}/projects/M21C_ls/output/raw_increment_summaries}"
OUTPUT_FILE="${OUTPUT_FILE:-${OUTPUT_DIR}/inst3_fcstana_raw_monthly_diagnostic_200006_202405.nc}"

SCRIPT="${REPO}/projects/M21C_ls/scripts/build_raw_inst3_fcstana_diagnostic_summaries.py"

mkdir -p "${OUTPUT_DIR}"

"${PYTHON}" "${SCRIPT}" \
  --root "${DA_ROOT}" \
  --file-prefix LS_DAv8_M36 \
  --start 2000-06-01 \
  --end 2024-05-31 \
  --output "${OUTPUT_FILE}" \
  --state-vars SFMC RZMC PRMC TSURF TSOIL1

echo "Wrote raw inst3 ANA-FCST diagnostic summaries to ${OUTPUT_FILE}"
