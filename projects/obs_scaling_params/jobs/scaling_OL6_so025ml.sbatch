#!/bin/bash
#SBATCH --job-name=scaling_OL6_so025ml
#SBATCH --output=/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params/jobs/%x.%j.out
#SBATCH --account=g0610
#SBATCH --qos=allnccs
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --constraint=mil
set -euo pipefail
PY=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
cd /discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params
export OMP_NUM_THREADS=8
echo "== scaling $(date)"
$PY -u scripts/run_scaling_params.py --exp-path /discover/nobackup/projects/land_da/M21C_testing --exp-run M21C_test_CF0360_OL6_superob025ml_clim --domain CF0360x6C_GLOBAL --start 2019-05 --end 2021-04 --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM --grid-resolution 0.25 --window-days 75 --ndata-min 20 --obsfcstana-format nc4 --prefix CF0360_OL6_HSAF_so025ml_zscore_stats_ --out-dir zclim_hsaf_so025ml
echo "== done $(date)"; ls -la /discover/nobackup/projects/land_da/M21C_testing/M21C_test_CF0360_OL6_superob025ml_clim/output/CF0360x6C_GLOBAL/stats/zclim_hsaf_so025ml
