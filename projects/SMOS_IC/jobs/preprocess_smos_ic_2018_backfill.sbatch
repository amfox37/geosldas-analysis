#!/bin/bash
#SBATCH --job-name=smos_ic_2018_backfill
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=smos_ic_2018_backfill.%j.out
#SBATCH --error=smos_ic_2018_backfill.%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"
SMOS_ROOT="${SMOS_ROOT:-/discover/nobackup/projects/land_da/SMOS_IC}"
OUT_DIR="${OUT_DIR:-/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily}"
TILECOORD="${TILECOORD:-/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121/output/SMAP_EASEv2_M36_GLOBAL/rc_out/OLv7_M36_MULTI_type_13_H121.ldas_tilecoord.bin}"
START_DATE="${START_DATE:-2018-01-01}"
END_DATE="${END_DATE:-2018-07-31}"
FORCE_EXTRACT="${FORCE_EXTRACT:-0}"

cd "${REPO_ROOT}"

module load python/GEOSpyD

extract_if_needed() {
  local archive="$1"
  local dirname="$2"
  local archive_path="${SMOS_ROOT}/${archive}"
  local dir_path="${SMOS_ROOT}/${dirname}"

  if [[ ! -f "${archive_path}" ]]; then
    echo "Missing archive: ${archive_path}" >&2
    exit 1
  fi

  if [[ "${FORCE_EXTRACT}" != "1" && -d "${dir_path}" ]]; then
    if find "${dir_path}" -maxdepth 1 -name '*.nc' -print -quit | grep -q .; then
      echo "Using existing extracted directory: ${dir_path}"
      return
    fi
  fi

  echo "Extracting ${archive_path} into ${SMOS_ROOT}"
  tar -xzf "${archive_path}" -C "${SMOS_ROOT}"
}

extract_if_needed "SMOS_IC_V2_ASC_2018_2021.tar.gz" "SMOS_IC_V2_ASC_2018_2021"
extract_if_needed "SMOS_IC_V2_DES_2018_2021.tar.gz" "SMOS_IC_V2_DES_2018_2021"
extract_if_needed "SMOS_IC_V2_ASC.tar.gz" "SMOS_IC_V2_ASC"
extract_if_needed "SMOS_IC_OPER_DES.tar.gz" "SMOS_IC_OPER_DES"

python -u projects/SMOS_IC/scripts/preprocess_smos_ic_daily_to_m36.py \
  --smos-root "${SMOS_ROOT}" \
  --tilecoord "${TILECOORD}" \
  --out-dir "${OUT_DIR}" \
  --start-date "${START_DATE}" \
  --end-date "${END_DATE}" \
  --scene-flag-max 1 \
  --tb-rmse-max 8.0
