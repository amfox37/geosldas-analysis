#!/usr/bin/env bash

set -euo pipefail

COLLECTION="CYGNSS_L3_SOIL_MOISTURE_V3.2"
EXTENSIONS="36km.a32.d33.nc"
YEAR=""
MONTH=""
LIMIT=""
DRY_RUN=0

usage() {
    cat <<'EOF'
Download CYGNSS L3 soil-moisture granules with podaac-data-downloader.

Usage:
  download_cygnss_l3.sh --dest DIR --year YYYY [options]

Required:
  -d, --dest DIR           Local product root (for example, /data/CYGNSS_L3)
  -y, --year YYYY          Year to download

Options:
  -m, --month MM            Limit to one month (default: all 12 months)
  -c, --collection NAME     PO.DAAC collection short name (default: CYGNSS_L3_SOIL_MOISTURE_V3.2)
  -e, --extensions REGEX    Filename regex to download (default: 36km.a32.d33.nc;
                             the collection also serves a 9km grid, use ".nc" for both)
  -n, --limit N              Max granules per month (useful for testing)
      --dry-run             Search and list granules without downloading
  -h, --help                Show this help

Files land in DEST/Y<yyyy>/M<mm>, matching this repo's obs directory convention.

Requires a ~/.netrc entry for machine urs.earthdata.nasa.gov (chmod 600).

Examples:
  # Test one month without downloading.
  download_cygnss_l3.sh -d /data/CYGNSS_L3 -y 2018 -m 08 --dry-run

  # Download that month for real.
  download_cygnss_l3.sh -d /data/CYGNSS_L3 -y 2018 -m 08

  # Download the whole year, one month at a time.
  download_cygnss_l3.sh -d /data/CYGNSS_L3 -y 2018
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 2
}

DEST=""

while (($#)); do
    case "$1" in
        -d|--dest) DEST=${2:-}; shift 2 ;;
        -y|--year) YEAR=${2:-}; shift 2 ;;
        -m|--month) MONTH=${2:-}; shift 2 ;;
        -c|--collection) COLLECTION=${2:-}; shift 2 ;;
        -e|--extensions) EXTENSIONS=${2:-}; shift 2 ;;
        -n|--limit) LIMIT=${2:-}; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$DEST" ]] || die "--dest is required"
[[ -n "$YEAR" ]] || die "--year is required"
command -v podaac-data-downloader >/dev/null 2>&1 || die "podaac-data-downloader is not installed (pip install podaac-data-subscriber)"
grep -q 'machine urs.earthdata.nasa.gov' ~/.netrc 2>/dev/null || die "no urs.earthdata.nasa.gov entry in ~/.netrc"

[[ "$YEAR" =~ ^[0-9]{4}$ ]] || die "--year must be YYYY"
[[ -z "$MONTH" || "$MONTH" =~ ^(0[1-9]|1[0-2])$ ]] || die "--month must be MM (01-12)"
[[ -z "$LIMIT" || "$LIMIT" =~ ^[1-9][0-9]*$ ]] || die "--limit must be a positive integer"

mkdir -p "$DEST"
DEST=$(cd "$DEST" && pwd -P)

if [[ -n "$MONTH" ]]; then
    MONTHS=("$MONTH")
else
    MONTHS=(01 02 03 04 05 06 07 08 09 10 11 12)
fi

for mm in "${MONTHS[@]}"; do
    start="${YEAR}-${mm}-01T00:00:00Z"
    end=$(date -u -d "${YEAR}-${mm}-01 +1 month -1 second" +%Y-%m-%dT%H:%M:%SZ)
    out_dir="$DEST/Y$YEAR/M$mm"
    mkdir -p "$out_dir"

    args=(-c "$COLLECTION" -d "$out_dir" -sd "$start" -ed "$end" -e "$EXTENSIONS")
    [[ -z "$LIMIT" ]] || args+=(--limit "$LIMIT")
    ((DRY_RUN == 0)) || args+=(--dry-run)

    printf 'Collection: %s\n' "$COLLECTION"
    printf 'Window:     %s to %s\n' "$start" "$end"
    printf 'Local dir:  %s\n' "$out_dir"
    ((DRY_RUN == 0)) || printf 'Mode:       dry run\n'
    printf '\n'

    podaac-data-downloader "${args[@]}"
done
