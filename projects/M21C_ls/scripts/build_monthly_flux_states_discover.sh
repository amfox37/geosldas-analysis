#!/usr/bin/env bash
set -euo pipefail

# Run this on Discover from any directory after pulling geosldas-analysis.
# It builds themed monthly auxiliary OL/DA files so Analysis 2 can be enabled
# without remaking or copying the larger monthly state bundles. Override GROUPS
# with a space-separated subset, for example:
#   GROUPS="flux_core water_budget" ./build_monthly_flux_states_discover.sh

REPO="${REPO:-/discover/nobackup/projects/land_da/geosldas-analysis}"
PYTHON="${PYTHON:-python}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO}/projects/M21C_ls/output/monthly_flux_states}"

SCRIPT="${REPO}/projects/discover_JH/read_monthly_states.py"

mkdir -p "${OUTPUT_DIR}"

GROUPS="${GROUPS:-flux_core latent_components water_budget energy_context}"

build_group() {
  local group="$1"
  shift
  local vars=("$@")

  echo "=== Building ${group}: ${vars[*]} ==="
  "${PYTHON}" "${SCRIPT}" \
    --vars "${vars[@]}" \
    --out-ol "${OUTPUT_DIR}/OLv8_${group}_2000_2024_compressed.nc" \
    --out-da "${OUTPUT_DIR}/DAv8_${group}_2000_2024_compressed.nc"
}

for group in ${GROUPS}; do
  case "${group}" in
    flux_core)
      build_group "${group}" EVLAND LHLAND SHLAND
      ;;
    latent_components)
      build_group "${group}" LHLANDTRNS LHLANDSOIL LHLANDINTR LHLANDSBLN
      ;;
    water_budget)
      build_group "${group}" RUNSURFLAND BASEFLOWLAND QINFILLAND SMLAND TWLAND WCHANGELAND
      ;;
    energy_context)
      build_group "${group}" SWLAND LWLAND GHLAND
      ;;
    *)
      echo "Unknown group: ${group}" >&2
      exit 2
      ;;
  esac
done

echo "Wrote monthly auxiliary state files to ${OUTPUT_DIR}"
