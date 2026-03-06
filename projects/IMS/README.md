# IMS Workflow (Download -> Regrid -> Validation)

IMS processing in this project has three core stages:

1. Download daily IMS ASCII snow-cover and write annual NetCDF.
2. Regrid annual IMS to M36 EASE with categorical-safe nearest neighbor.
3. Compare IMS with OL/DA using contingency metrics (accuracy/hit/miss/FAR/CRR).

## Project Layout

- `scripts/download_ims_ascii_to_nc.py`: download + convert `.asc.gz` to annual NetCDF (`ims_snowcover`).
- `scripts/regrid_ims_to_m36_nearest.py`: regrid annual IMS files to M36 EASE (`ims_category`).
- `scripts/run_ims_ol_da_cell_metrics.py`: batch OL/DA vs IMS metrics and cache outputs.
- `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`: direct end-to-end validation notebook.
- `notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`: plotting/tables from precomputed script outputs.
- `notebooks/ims_regrid_qc_colorado.ipynb`: map QC for raw IMS vs regridded IMS (plus model snapshots).
- `output/`: local staging/output area.

## Environment

Recommended packages:

- `numpy`
- `pandas`
- `xarray`
- `netCDF4`
- `scipy`
- `requests`
- `jupyter`
- `pyarrow`

Example:

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

- Expected source filename pattern: `imsYYYYDDD_00UTC_24km_v1.3.asc.gz`.
- If `--allow-missing-days` is set, missing days are written as `NaN`.

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

Behavior:

- Nearest-neighbor only (no averaging of categorical IMS values).
- One representative tile per EASE cell (`max frac_cell`, tie-break `min tile_id`).
- Writes diagnostics like `source_y/source_x`, distance, and representative elevation.

## 3A) Direct Notebook Validation

- `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`

Uses regridded IMS and OL/DA daily files directly. Produces paired daily metrics and bootstrap comparison tables.

## 3B) Precomputed Fast Pipeline (Recommended for reruns/maps)

Script:

- `scripts/run_ims_ol_da_cell_metrics.py`

Example:

```bash
python run_ims_ol_da_cell_metrics.py \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --year-start 2000 \
  --year-end 2024 \
  --ims-regrid-dir /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output \
  --ims-regrid-template 'ims_snowcover_24km_{year}_on_m36_nearest.nc4' \
  --ol-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36 \
  --da-run-root /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36 \
  --output-dir /discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output
```

Important options:

- `--min-ims-snow-days` (default `10`): keep only cells where IMS reports snow on at least N days.
- `--fast-mode`: skips model re-read and recomputes stats from cached cell-count files.

Fast-mode note:

- Comparison CIs are recomputed with cell bootstrap.
- Day-level paired table is not reconstructed in fast mode (placeholder output is written for file compatibility).

Main outputs:

- `ims_ol_da_cell_counts_metrics_*.nc4`
- `ims_ol_da_scope_metadata_*.csv`
- `ims_ol_da_comparison_table_*.parquet/csv`
- `ims_ol_da_pair_daily_*.parquet/csv`
- `ims_ol_da_daily_counts_*.parquet/csv`

Plotting notebook for precomputed outputs:

- `notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`

## Typical End-to-End Order

1. `download_ims_ascii_to_nc.py`
2. `regrid_ims_to_m36_nearest.py`
3. Either:
   - run `ims_daily_seasonal_ol_da_scf.ipynb` directly, or
   - run `run_ims_ol_da_cell_metrics.py` then `ims_maps_and_tables_from_precomputed_outputs.ipynb`.
