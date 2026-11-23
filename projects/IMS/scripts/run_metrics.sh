#!/usr/bin/env bash
# usage: run_metrics.sh <year_start> <year_end> <outdir> <tagbase>
set -euo pipefail
W=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix
REPO=/discover/nobackup/projects/land_da/geosldas-analysis
Y0=$1; Y1=$2; OUT=$3; TAG=$4
P=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
$P -I $REPO/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start $Y0 --year-end $Y1 --ims-regrid-dir $W/regrid --ims-var ims_category \
  --scf-threshold 0.5 --min-ims-snow-days 0 --output-dir $OUT --output-tag $TAG --overwrite
$P -I $REPO/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start $Y0 --year-end $Y1 --reuse-cell-counts-nc $OUT/ims_ol_da_cell_counts_metrics_${TAG}.nc4 \
  --reuse-scope-meta-csv $OUT/ims_ol_da_scope_metadata_${TAG}.csv --min-ims-snow-days 10 \
  --output-dir $OUT --output-tag ${TAG}_imsSnowDaysGe10 --overwrite
