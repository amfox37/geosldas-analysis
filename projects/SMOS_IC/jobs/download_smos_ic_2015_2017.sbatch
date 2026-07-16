#!/bin/bash
#SBATCH --job-name=smos_ic_download_2015_2017
#SBATCH --account=g0610
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=smos_ic_download_2015_2017.%j.out
#SBATCH --error=smos_ic_download_2015_2017.%j.err

set -euo pipefail

SMOS_ROOT="${SMOS_ROOT:-/discover/nobackup/projects/land_da/SMOS_IC}"

if ! command -v lftp >/dev/null 2>&1; then
  echo "lftp is not available on PATH. Load the needed module or install lftp before running." >&2
  exit 1
fi

mkdir -p "${SMOS_ROOT}/SMOS_IC_V2_ASC_2015_2017"
mkdir -p "${SMOS_ROOT}/SMOS_IC_V2_DES_2015_2017"

cd "${SMOS_ROOT}"

lftp <<'EOF'
set ssl:verify-certificate false
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
open -u anonymous,anonymous ftp://ib.remote-sensing.inrae.fr
mirror -c \
  --include-glob 'SM_RE07_MIR_CDF3SA_2015*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SA_2016*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SA_2017*.nc' \
  /SMOS_IC_V2__SM_2010-2024/ASC ./SMOS_IC_V2_ASC_2015_2017
mirror -c \
  --include-glob 'SM_RE07_MIR_CDF3SD_2015*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SD_2016*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SD_2017*.nc' \
  /SMOS_IC_V2__SM_2010-2024/DES ./SMOS_IC_V2_DES_2015_2017
bye
EOF
