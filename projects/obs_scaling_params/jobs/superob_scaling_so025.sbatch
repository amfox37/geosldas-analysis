#!/bin/bash
#SBATCH --job-name=superob_scaling_so025
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
echo "== aggregation $(date)"
$PY -u scripts/make_superob_obsfcstana.py --src-ana /gpfsm/dnb10/projects/p15/iau/merra_land/M21C_land/m21cland_excl_ens6_OL1/output/CF0360x6C_GLOBAL/ana/ens_avg --dst-ana /discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params/test_runs_so/m21cland_excl_ens6_OL1/output/CF0360x6C_GLOBAL/ana/ens_avg --species 12,13,14 --grid-res 0.25 --start 2019-02 --end 2020-08 --workers 32
echo "== scaling $(date)"
export OMP_NUM_THREADS=8
$PY -u scripts/run_scaling_params.py --exp-path /discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params/test_runs_so --exp-run m21cland_excl_ens6_OL1 --domain CF0360x6C_GLOBAL --start 2019-02 --end 2020-08 --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM --grid-resolution 0.25 --window-days 75 --ndata-min 20 --obsfcstana-format nc4 --prefix CF0360_OL1_HSAF_zscore_stats_ --out-dir zclim_so025
echo "== done $(date)"; ls -la /discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params/test_runs_so/m21cland_excl_ens6_OL1/output/CF0360x6C_GLOBAL/stats/zclim_so025
