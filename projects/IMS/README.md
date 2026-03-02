# IMS Workflow (Download -> Regrid -> Analysis)

This project has three main steps:

1. Download IMS daily ASCII files and build annual NetCDF files.
2. Regrid annual IMS files to the GEOSldas M36 EASE grid (nearest-neighbor categorical mapping).
3. Compare IMS vs OL/DA in the notebook using contingency metrics (hit/miss/FAR/etc.) with bootstrap CIs.

The analysis notebook expects IMS already regridded to M36 EASE.

## Project Layout

- `scripts/download_ims_ascii_to_nc.py`: download + convert `.asc.gz` to annual NetCDF.
- `scripts/regrid_ims_to_m36_nearest.py`: regrid annual IMS files to M36 EASE.
- `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`: OL vs DA comparison against regridded IMS.
- `output/`: local output/cache staging area.
- `notebooks/outputs_ims_ol_da_validation/`: notebook output tables.

## Environment

Recommended Python packages:

- `numpy`
- `pandas`
- `xarray`
- `netCDF4`
- `scipy`
- `requests`
- `jupyter`
- `pyarrow` (for parquet output)

Example with mamba into `regrid`:

```bash
mamba install -n regrid numpy pandas xarray netcdf4 scipy requests jupyter pyarrow
```

## 1) Download IMS Annual NetCDF

Run from `projects/IMS/scripts`:

```bash
python download_ims_ascii_to_nc.py \
  --year 2024 \
  --latlon-file ims_24km_latlon.nc4 \
  --out-file ../output/ims_snowcover_24km_2024.nc4 \
  --cache-dir ../output/gz \
  --allow-missing-days
```

Notes:

- Output variable is `ims_snowcover` (raw IMS category code values).
- If a day is missing/unreadable and `--allow-missing-days` is set, that day is written as `NaN`.
- Filename convention expected upstream: `imsYYYYDDD_00UTC_24km_v1.3.asc.gz`.

## 2) Regrid Annual IMS to M36 EASE

Run from `projects/IMS/scripts`:

```bash
python regrid_ims_to_m36_nearest.py \
  --ims-dir /gpfsm/dnb06/projects/p163/IMS \
  --out-dir ../output \
  --year-start 2000 \
  --year-end 2023 \
  --ims-template 'ims_snowcover_24km_{year}.nc4' \
  --out-template 'ims_snowcover_24km_{year}_on_m36_nearest.nc4' \
  --tilecoord /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_OLv8_M36.ldas_tilecoord.bin \
  --copy-time-values \
  --skip-missing-input
```

Important behavior:

- Regridding is nearest-neighbor only (categorical-safe, no averaging).
- One representative tile per EASE cell is used (`max frac_cell`, tie-break `min tile_id`).
- Output categorical variable is `ims_category`.
- Output includes mapping diagnostics (`source_y`, `source_x`, `nearest_distance_km`, `within_max_distance`) and `rep_elev_m`.

## 3) Run OL vs DA Analysis Notebook

Notebook:

- `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`

It compares regridded IMS (`ims_category`) with OL/DA daily cat files and writes:

- daily extraction cache
- paired OL-vs-DA daily table (same collocated samples)
- final comparison table with bootstrap CIs

Configure in the notebook:

- `YEAR_START`, `YEAR_END`
- `IMS_REGRID_DIR`, `IMS_REGRID_TEMPLATE`
- `EXPERIMENTS` run roots and `DOMAIN`
- `SCF_THRESHOLD` (default `0.5`)
- IMS snow/no-snow code mapping (`AUTO_INFER_IMS_CODES`, `IMS_SNOW_CODES`, `IMS_NO_SNOW_CODES`)
- bootstrap settings (`N_BOOTSTRAP`, seed, CI bounds)

Fair comparison note:

- OL and DA are scored using the same per-day collocated sample mask.
- Final table reports contingency metrics only:
  - `accuracy`
  - `hit_rate`
  - `miss_rate`
  - `false_alarm_ratio`
  - `correct_rejection_rate`
- Includes `OL`, `DA`, and `delta_da_minus_ol` with bootstrap CIs.

## Typical End-to-End Order

1. Build annual raw files with `download_ims_ascii_to_nc.py`.
2. Regrid annual files with `regrid_ims_to_m36_nearest.py`.
3. Run notebook `ims_daily_seasonal_ol_da_scf.ipynb` to generate OL-vs-DA verification tables.
