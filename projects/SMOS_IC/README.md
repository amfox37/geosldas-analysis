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
   - Output: daily sparse M36 files (`.npz`) + manifest CSV

2. Pair with model daily data (OL or DA):
   - Read preprocessed daily SMOS files by date.
   - Read model daily surface SM for same date.
   - Build paired outputs (`sm_obs`, `sm_mod`, `idx_EASEv2_lonxlat`) in the same style as other step2 workflows.

3. Continue with existing downstream steps:
   - Climatology
   - IVs / TC stats
   - Rdiff / figure notebooks

## Regridding proposal (efficient + consistent)

Use a one-time nearest-neighbor lookup from fixed SMOS-IC source cells to fixed M36 representative land cells.

- Build one representative tile per M36 cell (max `frac_cell`, tie-break min `tile_id`).
- Build one KDTree on representative M36 cell centers (`com_lat`, `com_lon`).
- Query all SMOS source cell centers once and cache source->target mapping.
- For each day:
  - apply QC on ASC and DES separately,
  - merge passes,
  - aggregate source values to target cells via `np.bincount` mean.

This avoids per-day `griddata`/triangulation and is much faster for long runs.

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

To download 2015-2017 ASC/DES raw files on a Discover compute node:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
sbatch projects/SMOS_IC/jobs/download_smos_ic_2015_2017.sbatch
```

That job writes to:

```text
/discover/nobackup/projects/land_da/SMOS_IC/SMOS_IC_V2_ASC_2015_2017
/discover/nobackup/projects/land_da/SMOS_IC/SMOS_IC_V2_DES_2015_2017
```

The preprocessing script currently indexes fixed raw directory names. After the
download, it will also index these two 2015-2017 directories directly.

After the download finishes, run the full missing-window backfill:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
sbatch projects/SMOS_IC/jobs/preprocess_smos_ic_2015_2018_backfill.sbatch
```

On Discover:

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

The following are the exact `lftp` commands used to download DES files used in this project.

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

## Notes

- Source folders expected under `--smos-root`:
  - `SMOS_IC_V2_ASC_2018_2021`
  - `SMOS_IC_V2_DES_2018_2021`
  - `SMOS_IC_V2_ASC`
  - `SMOS_IC_OPER_DES`
- Missing days are allowed; the script writes only dates with at least one valid pass.
- For VOD-oriented workflows, use `--tb-rmse-max 6.0`.
