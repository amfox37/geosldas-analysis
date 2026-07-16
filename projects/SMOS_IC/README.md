# SMOS-IC Workflow

This project stages SMOS-IC daily soil-moisture files for GEOSldas OL/DA validation workflows.

## Goal

Build a reusable, quality-controlled SMOS-IC daily cache on the GEOSldas M36 grid so we do **not** repeat expensive read/QC/regridding for every experiment run.

## Recommended workflow

1. Preprocess SMOS-IC observations once (obs-only):
   - Script: `scripts/preprocess_smos_ic_daily_to_m36.py`
   - Input: SMOS-IC ASC/DES daily NetCDF files
   - QC used:
     - `Scene_Flags <= 1`
     - `RMSE <= 8 K` (SM default)
     - `0 <= Soil_Moisture <= 1`
   - Output: daily sparse M36 NetCDF files (`.nc`) + mapping cache (`.npz`) +
     manifest CSV

2. Pair with model daily data (OL or DA):
   - Read preprocessed daily SMOS files by date.
   - Read model daily surface SM for same date.
   - Build paired outputs (`sm_obs`, `sm_mod`, `idx_EASEv2_lonxlat`) in the same style as other step2 workflows.

3. Continue with existing downstream steps:
   - Climatology
   - IVs / TC stats
   - Rdiff / figure notebooks

## Remapping Method

The current preprocessing script remaps SMOS-IC from its native 25 km EASE grid
to the GEOSldas M36 EASE grid with cached conservative overlap weights.

- Build one representative M36 land tile per EASE cell from the GEOSldas
  tilecoord file: choose the tile with maximum `frac_cell`, tie-break by
  smallest `tile_id`.
- Build conservative overlap weights once and cache them:
  - Source SMOS-IC 25 km cell edges come from the NetCDF `crs.GeoTransform`.
  - Target M36 cell edges come from fixed EASE2 M36 geometry
    (`dx=36032.220840584 m`, `964x406`).
  - For each target cell, compute geometric overlap area with source cells in
    projected EASE2 x/y space.
  - Save sparse mapping arrays: `map_tgt`, `map_src`, and `map_area`.
- For each day:
  - QC ASC and DES separately with `Scene_Flags <= 1`, `RMSE <= 8 K`, and
    `0 <= Soil_Moisture <= 1`.
  - Merge ASC+DES on the source grid using the mean where both passes are
    valid, otherwise the available pass.
  - Apply conservative remapping with area-weighted means:
    `num_t = sum(A_ts * SM_s)`, `den_t = sum(A_ts)`,
    `SM_t = num_t / den_t`.
  - Keep only target cells with coverage above `--min-coverage-frac`.
  - Save sparse daily outputs with `idx_EASEv2_lonxlat`, `sm_obs`, and
    `coverage_frac`.

This avoids per-day interpolation/triangulation and keeps the source-to-target
geometry fixed across the full run.

## Typical command

```bash
source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate /Users/amfox/mamba/envs/regrid

python projects/SMOS_IC/scripts/preprocess_smos_ic_daily_to_m36.py \
  --smos-root /Users/amfox/Desktop/SMOS_IC \
  --tilecoord /Users/amfox/Desktop/geosldas-analysis/outputs_daily_snow_seasonal_ol_da/LS_OLv8_M36.ldas_tilecoord.bin \
  --out-dir projects/SMOS_IC/output/preprocessed_m36_daily \
  --scene-flag-max 1 \
  --tb-rmse-max 8.0
```

## Discover Backfill

The existing Discover cache at
`/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily` was
originally produced from local Mac paths and synced to Discover. The raw SMOS-IC
inputs currently available on Discover are four tar archives under
`/discover/nobackup/projects/land_da/SMOS_IC/`; after extraction, the
preprocessing script expects these directories directly under that same root:

- `SMOS_IC_V2_ASC_2018_2021`
- `SMOS_IC_V2_DES_2018_2021`
- `SMOS_IC_V2_ASC`
- `SMOS_IC_OPER_DES`

The useful immediate backfill from already-present archives is 2018-01-01
through 2018-07-31. For the full HSAF comparison window, the missing 2015-2017
raw files must be downloaded separately.

Because Discover compute nodes may not be able to reach the INRAE FTP server,
the most reliable path is to download raw files locally, transfer either raw
files or preprocessed output to Discover, and run the remaining preprocessing
there if needed.

Local download for 2015-2017 ASC/DES:

```bash
cd /Users/amfox/Desktop/geosldas-analysis
SMOS_ROOT=/Users/amfox/Desktop/SMOS_IC \
  bash projects/SMOS_IC/jobs/download_smos_ic_2015_2017.sbatch
```

Transfer the raw 2015-2017 folders to Discover:

```bash
scp -r /Users/amfox/Desktop/SMOS_IC/SMOS_IC_V2_ASC_2015_2017 \
  amfox@discover:/discover/nobackup/projects/land_da/SMOS_IC/

scp -r /Users/amfox/Desktop/SMOS_IC/SMOS_IC_V2_DES_2015_2017 \
  amfox@discover:/discover/nobackup/projects/land_da/SMOS_IC/
```

Then run the full missing-window backfill on Discover:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
sbatch projects/SMOS_IC/jobs/preprocess_smos_ic_2015_2018_backfill.sbatch
```

The 2015-2017 download job is also available for environments where outbound
FTP works:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
sbatch projects/SMOS_IC/jobs/download_smos_ic_2015_2017.sbatch
```

For only the 2018-01-01 through 2018-07-31 backfill from already-present
archives:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
sbatch projects/SMOS_IC/jobs/preprocess_smos_ic_2018_backfill.sbatch
```

The job extracts archives only when the expected extracted directory is missing
or empty. Set `FORCE_EXTRACT=1` to re-extract. It writes daily sparse NetCDF
files to:

```text
/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily
```

## Data download provenance (SMOS-IC v2)

SMOS-IC v2 product documentation/source:

- https://ib.remote-sensing.inrae.fr/index.php/smos-ic-v2-product-documentation/

The original raw cache was downloaded locally under `/Users/amfox/Desktop/SMOS_IC`
on `gs6101-madrid` with `lftp`. The following commands document that local
workflow.

### Operational ASC files (`SMOS_IC_V2_ASC`)

```bash
lftp <<'EOF'
set ssl:verify-certificate false
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
open -u anonymous,anonymous ftp://ib.remote-sensing.inrae.fr
mirror -c \
  --include-glob 'SM_OPER_MIR_CDF3SA_*.nc' \
  /SMOS_IC_V2__SM_2010-2024/ASC ./SMOS_IC_V2_ASC
bye
EOF
```

### Operational DES files (`SMOS_IC_OPER_DES`)

```bash
lftp <<'EOF'
set ssl:verify-certificate false
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
open -u anonymous,anonymous ftp://ib.remote-sensing.inrae.fr
mirror -c \
  --include-glob 'SM_OPER_MIR_CDF3SD_*.nc' \
  /SMOS_IC_V2__SM_2010-2024/DES ./SMOS_IC_OPER_DES
bye
EOF
```

### Reprocessed ASC subset for 2018-2021 (`SMOS_IC_V2_ASC_2018_2021`)

```bash
lftp <<'EOF'
set ssl:verify-certificate false
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
open -u anonymous,anonymous ftp://ib.remote-sensing.inrae.fr
mirror -c \
  --include-glob 'SM_RE07_MIR_CDF3SA_2018*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SA_2019*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SA_2020*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SA_2021*.nc' \
  /SMOS_IC_V2__SM_2010-2024/ASC ./SMOS_IC_V2_ASC_2018_2021
bye
EOF
```

### Reprocessed DES subset for 2018-2021 (`SMOS_IC_V2_DES_2018_2021`)

```bash
lftp <<'EOF'
set ssl:verify-certificate false
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
open -u anonymous,anonymous ftp://ib.remote-sensing.inrae.fr
mirror -c \
  --include-glob 'SM_RE07_MIR_CDF3SD_2018*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SD_2019*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SD_2020*.nc' \
  --include-glob 'SM_RE07_MIR_CDF3SD_2021*.nc' \
  /SMOS_IC_V2__SM_2010-2024/DES ./SMOS_IC_V2_DES_2018_2021
bye
EOF
```

### Reprocessed ASC/DES subset for 2015-2017

Use the script wrapper so ASC and DES stay symmetric:

```bash
cd /Users/amfox/Desktop/geosldas-analysis
SMOS_ROOT=/Users/amfox/Desktop/SMOS_IC \
  bash projects/SMOS_IC/jobs/download_smos_ic_2015_2017.sbatch
```

## Local Cache Provenance

Local files currently present in this workspace:

- Raw SMOS-IC root: `/Users/amfox/Desktop/SMOS_IC`
- Preprocessed cache: `projects/SMOS_IC/output/preprocessed_m36_daily`
- Existing preprocessed date range: `2018-08-01` through `2024-06-30`
- Existing cache tarball: `projects/SMOS_IC/output/preprocessed_m36_daily.tar.gz`

The existing Discover preprocessed cache was produced from this local
preprocessed directory and synced to:

```text
/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily
```

To package and transfer a locally regenerated preprocessed cache:

```bash
cd /Users/amfox/Desktop/geosldas-analysis/projects/SMOS_IC/output
tar -czf preprocessed_m36_daily.tar.gz preprocessed_m36_daily
scp preprocessed_m36_daily.tar.gz \
  amfox@discover:/discover/nobackup/projects/land_da/SMOS_IC/
```

Then on Discover:

```bash
cd /discover/nobackup/projects/land_da/SMOS_IC
tar -xzf preprocessed_m36_daily.tar.gz
```

## Notes

- Source folders expected under `--smos-root`:
  - `SMOS_IC_V2_ASC_2015_2017`
  - `SMOS_IC_V2_DES_2015_2017`
  - `SMOS_IC_V2_ASC_2018_2021`
  - `SMOS_IC_V2_DES_2018_2021`
  - `SMOS_IC_V2_ASC`
  - `SMOS_IC_OPER_DES`
- Missing days are allowed; the script writes only dates with at least one valid pass.
- For VOD-oriented workflows, use `--tb-rmse-max 6.0`.
