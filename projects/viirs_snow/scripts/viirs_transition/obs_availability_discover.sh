#!/usr/bin/env bash
set -euo pipefail

# Analysis A (docs/viirs_transition_plan.md): SCF observation availability per arm.
# Run on Discover from any directory. Reads config/arms.json and the 3-hourly ObsFcstAna files of each arm;
# writes one NetCDF per arm to output/obs_availability/ (gitignored). Override with environment variables:
#   ARMS="TA_to_J1J2 TA_control" ./obs_availability_discover.sh     # subset of arms (default: all in arms.json)
#   BEG=2024-01-01 END=2024-01-31 ARMS=TA_to_J1J2 OUT_DIR=... ./obs_availability_discover.sh   # test range

REPO="${REPO:-/discover/nobackup/projects/land_da/geosldas-analysis}"
PYTHON="${PYTHON:-/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3}"
PROJECT="${REPO}/projects/viirs_snow"
OUT_DIR="${OUT_DIR:-${PROJECT}/output/obs_availability}"
SCRIPT="${PROJECT}/scripts/viirs_transition/obs_availability.py"

if [[ -z "${ARMS:-}" ]]; then
  ARMS="$("${PYTHON}" -I -c "import json; print(' '.join(json.load(open('${PROJECT}/config/arms.json'))['arms']))")"
fi

mkdir -p "${OUT_DIR}"
for arm in ${ARMS}; do
  echo "=== ${arm} ==="
  extra=()
  [[ -n "${BEG:-}" ]] && extra+=(--beg "${BEG}")
  [[ -n "${END:-}" ]] && extra+=(--end "${END}")
  PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -I "${SCRIPT}" --arm "${arm}" \
    --out "${OUT_DIR}/obs_availability_${arm}.nc" "${extra[@]}"
done
echo "Wrote ${OUT_DIR}"
