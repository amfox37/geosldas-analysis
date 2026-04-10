# SNOTEL Downloader + Validation

Download daily SNOTEL snow observations from USDA NRCS AWDB:

- `WTEQ` (SWE)
- `SNWD` (snow depth)

The downloader writes per-station CSV files, station metadata, and a combined parquet file used by downstream validation notebooks.

## Scripts

- `scripts/download_snotel_swe_snwd.py`
- `scripts/read_parquet_example.py`

## Notebook

- `notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb`
  - Compares OL/DA daily model snow variables against SNOTEL SWE + snow depth.
  - Uses station lat/lon from the parquet output for model-to-site extraction.

## Requirements

- Python 3.10+
- `pandas`
- `requests`
- `tqdm`
- `zeep`
- `pyarrow` (for parquet output)

Conda/mamba example (`regrid` environment):

```bash
mamba install -n regrid pandas requests tqdm pyarrow
mamba run -n regrid pip install zeep
```

## Quick Start

From repo root (`/Users/amfox/Desktop/geosldas-analysis`):

Download all stations with defaults:

```bash
python projects/SNOTEL/scripts/download_snotel_swe_snwd.py \
  --out projects/SNOTEL/output
```

Download only 10 stations:

```bash
python projects/SNOTEL/scripts/download_snotel_swe_snwd.py \
  --out projects/SNOTEL/output \
  --limit 10
```

Overwrite existing station CSVs:

```bash
python projects/SNOTEL/scripts/download_snotel_swe_snwd.py \
  --out projects/SNOTEL/output \
  --overwrite
```

Custom date range:

```bash
python projects/SNOTEL/scripts/download_snotel_swe_snwd.py \
  --out projects/SNOTEL/output \
  --begin 2000-06-01 \
  --end 2025-09-30
```

Read parquet quickly:

```bash
python projects/SNOTEL/scripts/read_parquet_example.py \
  --parquet projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet
```

## CLI Options

- `--begin` default `1999-10-01`
- `--end` default `2025-09-30`
- `--out` output directory (required)
- `--limit` number of stations (default: all)
- `--sleep` delay between REST calls, default `0.2` seconds
- `--overwrite` re-download even if per-station CSV exists

## Output Layout

```text
<out>/
  snotel_station_triplets.txt
  stations_metadata.csv
  stations/
    <stationTriplet>.csv
    <stationTriplet>.ERROR.txt   # only for failed stations
  all_stations_daily_wteq_snwd.parquet
```

Downstream usage:

- `all_stations_daily_wteq_snwd.parquet` is the standard input for model validation notebooks.
- Station metadata (`latitude`, `longitude`, `elevation`) are merged into the parquet for extraction/mapping.

## What Data Is Downloaded

### Data sources

- Station discovery + metadata: AWDB SOAP/WSDL service
- Daily time series: AWDB REST v1 `/data` endpoint

The script requests network code `SNTL` (SNOTEL), `duration=DAILY`, and elements:

- `WTEQ` (snow water equivalent)
- `SNWD` (snow depth)

### Station identifiers

- Station IDs are AWDB triplets, e.g. `1234:ID:SNTL`
- Per-station CSV filenames replace `:` with `_`, e.g. `1234_ID_SNTL.csv`

### File contents

- `snotel_station_triplets.txt`
  - One station triplet per line for this run.

- `stations_metadata.csv`
  - Metadata returned by SOAP for requested stations:
  - `stationTriplet`, `name`, `latitude`, `longitude`, `elevation`, `state`, `countyName`, `huc`

- `stations/<station>.csv`
  - Daily series for one station.
  - Columns: `date`, `WTEQ`, `SNWD`, `stationTriplet`
  - One row per day where data exists in AWDB response (outer-merged across elements).

- `all_stations_daily_wteq_snwd.parquet`
  - Combined daily data across stations.
  - Includes per-row station metadata columns merged in:
  - `name`, `latitude`, `longitude`, `elevation`

### Date coverage

- Controlled by `--begin` and `--end` (inclusive).
- Defaults:
  - `--begin 1999-10-01`
  - `--end 2025-09-30`

### Missing/empty data behavior

- Missing element/day values are written as empty/NaN in tabular outputs.
- If a station returns no daily data in the requested window, an empty CSV is still written.
- Request failures for individual stations are logged as:
  - `stations/<station>.ERROR.txt`

### Units

- Values are stored exactly as returned by AWDB (no unit conversion in downloader).
- For SNOTEL, `WTEQ` and `SNWD` are commonly interpreted in inches in downstream workflows; convert as needed for model comparison.

Common conversion used in notebook workflows:

- inches -> meters: `m = inches * 0.0254`

## Resume Behavior

If a per-station CSV already exists and `--overwrite` is not used:

- the station is skipped for download
- the existing CSV is still included in the combined parquet build

## Read the Combined Parquet

```bash
python projects/SNOTEL/scripts/read_parquet_example.py \
  --parquet projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet
```

With optional filters:

```bash
python projects/SNOTEL/scripts/read_parquet_example.py \
  --parquet projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet \
  --station "1234:ID:SNTL" \
  --start 2015-10-01 \
  --end 2016-09-30
```

## Typical Workflow

1. Run `download_snotel_swe_snwd.py` to build station CSVs + combined parquet.
2. Sanity-check parquet with `read_parquet_example.py`.
3. Run `notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb` for OL/DA validation.

## Latest OL/DA validation workflow (quick reference)

1. Build or refresh SNOTEL source parquet:
   - Script: `projects/SNOTEL/scripts/download_snotel_swe_snwd.py`
   - Required output:
     - `projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet`

2. Run validation notebook:
   - `projects/SNOTEL/notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb`
   - Reads:
     - SNOTEL parquet above (`stationTriplet`, `WTEQ`, `SNWD`, station metadata)
     - GEOSldas OL/DA daily files (`SNODPLAND`, `SNOMASLAND`) + tilecoord

3. Review notebook outputs:
   - `projects/SNOTEL/outputs_snotel_ol_da_validation/`
   - Includes:
     - `snotel_raw_timeseries_*.nc`
     - `snotel_station_metrics_*.csv/.parquet`
     - `snotel_domain_metrics_*.csv/.parquet`
     - figures/tables subfolders

## Notes

- The script uses retry/backoff for SOAP and REST calls.
- AWDB values are preserved as returned by NRCS; unit conversion is handled in downstream analysis notebooks.
