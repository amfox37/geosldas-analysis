#!/usr/bin/env bash
# usage: run2024.sh <regrid_dir> <outdir> <tag>   (pass 1 only, 2024 only, model dailies extracted under model/)
set -euo pipefail
W=/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix
P=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
TC=/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_OLv8_M36.ldas_tilecoord.bin
$P -I /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/scripts/run_ims_ol_da_cell_metrics.py --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start 2024 --year-end 2024 --ims-regrid-dir $1 --ims-var ims_category --scf-threshold 0.5 --min-ims-snow-days 0 \
  --ol-run-root $W/model/LS_OLv8_M36 --da-run-root $W/model/LS_DAv8_M36 --tilecoord $TC \
  --output-dir $2 --output-tag $3 --overwrite
