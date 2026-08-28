#!/usr/bin/env bash

set -euo pipefail

HOST="ftphsaf.meteoam.it"
PRODUCT="h121"
FORMAT="netcdf"
SATELLITE="all"
YEAR=""
MONTH=""
PARALLEL=4
DRY_RUN=0

usage() {
    cat <<'EOF'
Download H SAF ASCAT soil-moisture files with lftp.

Usage:
  download_hsaf_ascat.sh --user USER --dest DIR [options]

Required:
  -u, --user USER          H SAF FTP username
  -d, --dest DIR           Local product root (for example, /data/H121)

Options:
  -p, --product PRODUCT    h121, h139, or h29 (default: h121)
  -s, --satellite SAT      metop_a, metop_b, metop_c, or all (default: all)
  -y, --year YYYY          Limit the download to one year
  -m, --month MM           Limit the download to one month (requires --year)
  -j, --parallel N         Number of parallel transfers (default: 4)
      --dry-run            Show what lftp would transfer without downloading
  -h, --help               Show this help

Examples:
  # Test one month from Metop-A; lftp prompts for the password.
  download_hsaf_ascat.sh -u USER -d /data/H121 -s metop_a -y 2020 -m 01 --dry-run

  # Download that month.
  download_hsaf_ascat.sh -u USER -d /data/H121 -s metop_a -y 2020 -m 01

  # Resume/synchronize the complete H121 NetCDF archive.
  download_hsaf_ascat.sh -u USER -d /data/H121
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 2
}

USER_NAME=""
DEST=""

while (($#)); do
    case "$1" in
        -u|--user) USER_NAME=${2:-}; shift 2 ;;
        -d|--dest) DEST=${2:-}; shift 2 ;;
        -p|--product) PRODUCT=${2:-}; shift 2 ;;
        -s|--satellite) SATELLITE=${2:-}; shift 2 ;;
        -y|--year) YEAR=${2:-}; shift 2 ;;
        -m|--month) MONTH=${2:-}; shift 2 ;;
        -j|--parallel) PARALLEL=${2:-}; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$USER_NAME" ]] || die "--user is required"
[[ -n "$DEST" ]] || die "--dest is required"
command -v lftp >/dev/null 2>&1 || die "lftp is not installed (macOS: brew install lftp)"

PRODUCT=$(printf '%s' "$PRODUCT" | tr '[:upper:]' '[:lower:]')
FORMAT=$(printf '%s' "$FORMAT" | tr '[:upper:]' '[:lower:]')
SATELLITE=$(printf '%s' "$SATELLITE" | tr '[:upper:]' '[:lower:]')

case "$PRODUCT" in
    h121|h139|h29) ;;
    *) die "--product must be h121, h139, or h29" ;;
esac

case "$SATELLITE" in
    all|metop_a|metop_b|metop_c) ;;
    *) die "--satellite must be metop_a, metop_b, metop_c, or all" ;;
esac

[[ -z "$YEAR" || "$YEAR" =~ ^[0-9]{4}$ ]] || die "--year must be YYYY"
[[ -z "$MONTH" || "$MONTH" =~ ^(0[1-9]|1[0-2])$ ]] || die "--month must be MM (01-12)"
[[ -z "$MONTH" || -n "$YEAR" ]] || die "--month requires --year"
[[ "$PARALLEL" =~ ^[1-9][0-9]*$ ]] || die "--parallel must be a positive integer"
[[ "$USER_NAME" != *$'\n'* && "$USER_NAME" != *$'\r'* ]] || die "invalid username"
[[ "$DEST" != *$'\n'* && "$DEST" != *$'\r'* ]] || die "invalid destination"

mkdir -p "$DEST"
DEST=$(cd "$DEST" && pwd -P)

# Escape the only characters significant inside an lftp double-quoted string.
LFTP_DEST=${DEST//\\/\\\\}
LFTP_DEST=${LFTP_DEST//\"/\\\"}

if [[ "$SATELLITE" == "all" ]]; then
    SATELLITES=(metop_a metop_b metop_c)
else
    SATELLITES=("$SATELLITE")
fi

LFTP_COMMANDS="set cmd:fail-exit yes; set ftp:passive-mode yes; lcd \"$LFTP_DEST\"; cd /$PRODUCT/$FORMAT;"
for sat in "${SATELLITES[@]}"; do
    remote_path=$sat
    [[ -z "$YEAR" ]] || remote_path+="/$YEAR"
    [[ -z "$MONTH" ]] || remote_path+="/$MONTH"

    mirror_opts="--continue --only-newer --parallel=$PARALLEL"
    ((DRY_RUN == 0)) || mirror_opts+=" --dry-run"
    LFTP_COMMANDS+=" mirror $mirror_opts $remote_path $remote_path;"
done
LFTP_COMMANDS+=" bye"

printf 'Server:      ftp://%s\n' "$HOST"
printf 'Remote root: /%s/%s\n' "$PRODUCT" "$FORMAT"
printf 'Local root:  %s\n' "$DEST"
printf 'Selection:   %s%s%s\n' "$SATELLITE" "${YEAR:+/$YEAR}" "${MONTH:+/$MONTH}"
((DRY_RUN == 0)) || printf 'Mode:        dry run\n'
printf '\nlftp will prompt for the H SAF password.\n\n'

lftp -u "$USER_NAME" "ftp://$HOST" -e "$LFTP_COMMANDS"
