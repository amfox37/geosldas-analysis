#!/bin/bash
# Shared settings for the CYGNSS L1 AZ IV/TC jobs (sourced, not run).
#
# Triplet: model + SMAP L3 + ASCAT H121.  Neither product is assimilated in
# the arms below, and SMAP L3 and ASCAT are independent of CYGNSS.
# DA_SMAP_fixedop is deliberately NOT listed: it assimilates SMAP Tb, so
# SMAP L3 is not an independent reference for it.

PROJ="${PROJ:-/discover/nobackup/projects/land_da/cygl1_operator_test}"
REPO_ROOT="${REPO_ROOT:-/discover/nobackup/projects/land_da/geosldas-analysis}"

# Everything the jobs write goes in the land_da project space (never nobackup).
OUT_ROOT="${OUT_ROOT:-${PROJ}/iv_tc_output}"

START_DATE="${START_DATE:-2020-01-01}"
END_DATE="${END_DATE:-2022-12-31}"

# Index order = SLURM_ARRAY_TASK_ID.  NAME is the run label used in file names.
ARM_NAMES=(OL L1_filt275 L3)
ARM_DIRS=(
  "${PROJ}/OLv8_M36_AZ_fixedop"
  "${PROJ}/DA_L1_full_xc015_coh040216_fixedop"
  "${PROJ}/DA_L3_fixedop"
)

ARM_IDX="${SLURM_ARRAY_TASK_ID:-${ARM_IDX:-}}"
if [[ -z "${ARM_IDX}" || "${ARM_IDX}" -ge "${#ARM_NAMES[@]}" ]]; then
  echo "ERROR: need SLURM_ARRAY_TASK_ID (or ARM_IDX) in 0..$(( ${#ARM_NAMES[@]} - 1 ))" >&2
  exit 2
fi
ARM_NAME="${ARM_NAMES[${ARM_IDX}]}"
ARM_DIR="${ARM_DIRS[${ARM_IDX}]}"
