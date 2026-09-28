#!/bin/bash
#SBATCH --job-name=obs_scaling_H121_ASCAT_peatqc
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=obs_scaling_H121_ASCAT_peatqc.%j.out
#SBATCH --error=obs_scaling_H121_ASCAT_peatqc.%j.err

set -euo pipefail

# Regenerates the H121 ASCAT HSAF z-score scaling stats from the SAME
# OLv7_M36_MULTI_type_13_H121 OFA archive used by
# run_python_OLv7_H121_ASCAT.sbatch, but excluding observations at tiles
# flagged by the retroactive peat-footprint QC mask (built from static BCS
# porosity + tile geometry, FOV=25.0 km matching this OL's own obs_param
# config -- see /gpfsm/dnb06/projects/p284/hsaf_cdr_test/
# build_hsaf_peat_obs_tilenum_mask.py for the derivation).
#
# Writes to a SEPARATE --out-dir/--prefix so the original (un-QC'd) stats
# file is left untouched for comparison.

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"

cd "${REPO_ROOT}/projects/obs_scaling_params"

module load python/GEOSpyD
python -u scripts/run_scaling_params.py \
  --exp-path /discover/nobackup/projects/land_da/hsaf_cdr_test \
  --exp-run OLv7_M36_MULTI_type_13_H121 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2015-04 \
  --end 2021-04 \
  --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM \
  --grid-resolution 0.25 \
  --window-days 75 \
  --ndata-min 20 \
  --obsfcstana-format nc4 \
  --exclude-tiles-csv /discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_peat_obs_tilenum_mask.csv \
  --exclude-tiles-column reject \
  --prefix OLv7_M36_H121_ASCAT_zscore_stats_peatqc_ \
  --out-dir H121_ASCAT_python_zscore_0p25deg_peatqc_masked
