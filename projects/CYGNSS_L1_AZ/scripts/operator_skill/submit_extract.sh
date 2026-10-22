#!/bin/bash
# Build the per-obs L1 operator table, one SLURM job per year.
#   submit_extract.sh 2020 2021 2022          # full years
#   MONTHS=6-6 submit_extract.sh 2020         # quick test month
# Extra args after the years go to extract_obs_table.py via EXTRA, e.g. EXTRA=--overwrite
P=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
HERE=/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ
LOG=$HERE/output/operator_skill/logs
mkdir -p $LOG
for y in "$@"; do
  sbatch --account=s3208 --time=${WALL:-2:00:00} --ntasks=1 --cpus-per-task=4 --mem=40G \
         ${QOS:+--qos=$QOS} --job-name=l1op_$y${MONTHS:+_$MONTHS} --output=$LOG/extract_${y}${MONTHS:+_$MONTHS}_%j.out \
         --wrap="$P $HERE/scripts/operator_skill/extract_obs_table.py --year $y --months ${MONTHS:-1-12} $EXTRA"
done
