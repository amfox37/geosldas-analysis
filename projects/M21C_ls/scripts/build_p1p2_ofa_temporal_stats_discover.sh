#!/usr/bin/env bash
set -euo pipefail

# Run this on Discover from any directory. Override REPO, PYTHON, RUN_BASE,
# OL_MONTHLY_ROOT, DA_MONTHLY_ROOT, or OUTPUT_DIR if needed.

REPO="${REPO:-/discover/nobackup/projects/land_da/geosldas-analysis}"
PYTHON="${PYTHON:-python}"
DOMAIN="${DOMAIN:-SMAP_EASEv2_M36_GLOBAL}"
RUN_BASE="${RUN_BASE:-/discover/nobackup/projects/land_da/M21C_land_sweeper}"

OL_MONTHLY_ROOT="${OL_MONTHLY_ROOT:-${RUN_BASE}/LS_OLv8_M36_v2/LS_OLv8_M36/output/${DOMAIN}/ana/ens_avg}"
DA_MONTHLY_ROOT="${DA_MONTHLY_ROOT:-${RUN_BASE}/LS_DAv8_M36_v3/LS_DAv8_M36/output/${DOMAIN}/ana/ens_avg}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO}/projects/M21C_ls/output/ofa_temporal_stats}"

SCRIPT="${REPO}/projects/M21C_ls/scripts/build_ofa_temporal_stats.py"
MONTHLY_SUMS_PATTERN="{tag}.ens_avg.ldas_ObsFcstAna_sums.{ym}.nc4"

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${OL_MONTHLY_ROOT}" \
  --expid LS_OLv8_M36 \
  --exp-tag OL \
  --start 2000-06-01 \
  --end 2002-06-30 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-pattern "${MONTHLY_SUMS_PATTERN}" \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${DA_MONTHLY_ROOT}" \
  --expid LS_DAv8_M36 \
  --exp-tag DA \
  --start 2000-06-01 \
  --end 2002-06-30 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-pattern "${MONTHLY_SUMS_PATTERN}" \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${OL_MONTHLY_ROOT}" \
  --expid LS_OLv8_M36 \
  --exp-tag OL \
  --start 2002-07-01 \
  --end 2007-05-31 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-pattern "${MONTHLY_SUMS_PATTERN}" \
  --overwrite

"${PYTHON}" "${SCRIPT}" \
  --monthly-root "${DA_MONTHLY_ROOT}" \
  --expid LS_DAv8_M36 \
  --exp-tag DA \
  --start 2002-07-01 \
  --end 2007-05-31 \
  --output-dir "${OUTPUT_DIR}" \
  --monthly-pattern "${MONTHLY_SUMS_PATTERN}" \
  --overwrite

echo "Wrote P1/P2 OFA temporal stats to ${OUTPUT_DIR}"
