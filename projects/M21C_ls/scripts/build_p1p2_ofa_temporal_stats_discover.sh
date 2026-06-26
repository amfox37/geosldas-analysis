#!/usr/bin/env bash
set -euo pipefail

# Run this on Discover from any directory. Override REPO, PYTHON, MONTHLY_BASE,
# OL_MONTHLY_ROOT, DA_MONTHLY_ROOT, or OUTPUT_DIR if needed.

REPO="${REPO:-/discover/nobackup/projects/land_da/geosldas-analysis}"
PYTHON="${PYTHON:-python}"
DOMAIN="${DOMAIN:-SMAP_EASEv2_M36_GLOBAL}"
MONTHLY_BASE="${MONTHLY_BASE:-/discover/nobackup/qliu/SMAP_diag}"

OL_MONTHLY_ROOT="${OL_MONTHLY_ROOT:-${MONTHLY_BASE}/LS_OLv8_M36/output/${DOMAIN}/ana/ens_avg}"
DA_MONTHLY_ROOT="${DA_MONTHLY_ROOT:-${MONTHLY_BASE}/LS_DAv8_M36/output/${DOMAIN}/ana/ens_avg}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO}/projects/M21C_ls/output/ofa_temporal_stats}"

SCRIPT="${REPO}/projects/M21C_ls/scripts/build_ofa_temporal_stats.py"

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${OL_MONTHLY_ROOT}" \
  --expid LS_OLv8_M36 \
  --exp-tag OL \
  --start 2000-06-01 \
  --end 2002-06-30 \
  --output-dir "${OUTPUT_DIR}" \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${DA_MONTHLY_ROOT}" \
  --expid LS_DAv8_M36 \
  --exp-tag DA \
  --start 2000-06-01 \
  --end 2002-06-30 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-suffix _stats.nc4 \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${OL_MONTHLY_ROOT}" \
  --expid LS_OLv8_M36 \
  --exp-tag OL \
  --start 2002-07-01 \
  --end 2007-05-31 \
  --output-dir "${OUTPUT_DIR}" \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${DA_MONTHLY_ROOT}" \
  --expid LS_DAv8_M36 \
  --exp-tag DA \
  --start 2002-07-01 \
  --end 2007-05-31 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-suffix _stats.nc4 \
  --overwrite

echo "Wrote P1/P2 OFA temporal stats to ${OUTPUT_DIR}"
