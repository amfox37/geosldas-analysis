#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}
OUT_DIR=${OUT_DIR:-${REPO_ROOT}/projects/IMS/output}

RAW_TAG=SMAP_EASEv2_M36_GLOBAL_2000_2007_thr0p50_terraAquaScopes_raw
FINAL_TAG=SMAP_EASEv2_M36_GLOBAL_2000_2007_thr0p50_imsSnowDaysGe10_terraAquaScopes

python "${REPO_ROOT}/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py" \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start 2000 \
  --year-end 2007 \
  --ims-regrid-dir "${OUT_DIR}" \
  --ims-regrid-template 'ims_snowcover_24km_{year}_on_m36_nearest.nc4' \
  --ims-var ims_category \
  --ol-exp-name LS_OLv8_M36 \
  --ol-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36 \
  --da-exp-name LS_DAv8_M36 \
  --da-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36 \
  --scf-threshold 0.5 \
  --min-ims-snow-days 0 \
  --custom-period-scope P1_MODIS_Terra_SCF:2000-06-01:2002-06-30 \
  --custom-period-scope P2_MODIS_Terra_Aqua_SCF:2002-07-01:2007-05-31 \
  --output-dir "${OUT_DIR}" \
  --output-tag "${RAW_TAG}" \
  --overwrite

python "${REPO_ROOT}/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py" \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start 2000 \
  --year-end 2007 \
  --reuse-cell-counts-nc "${OUT_DIR}/ims_ol_da_cell_counts_metrics_${RAW_TAG}.nc4" \
  --reuse-scope-meta-csv "${OUT_DIR}/ims_ol_da_scope_metadata_${RAW_TAG}.csv" \
  --min-ims-snow-days 10 \
  --output-dir "${OUT_DIR}" \
  --output-tag "${FINAL_TAG}" \
  --overwrite
