# IMS Validation Workflow

This project validates GEOS-LDAS OL/DA snow cover against IMS daily snow-cover products on the M36 grid.

## Latest workflow run order

1. Build annual IMS NetCDF from daily ASCII:
   - `scripts/download_ims_ascii_to_nc.py`

2. Regrid IMS annual files to M36 EASE (nearest-neighbor categorical mapping):
   - `scripts/regrid_ims_to_m36_nearest.py`

3. Build OL/DA vs IMS validation metrics (recommended scripted pipeline):
   - `scripts/run_ims_ol_da_cell_metrics.py`

4. Make maps/tables from precomputed outputs:
   - `notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`

5. Optional direct notebook workflow (no precompute step):
   - `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`

6. Optional regridding QC map checks:
   - `notebooks/ims_regrid_qc_colorado.ipynb`

## Inputs

- IMS daily source data (`imsYYYYDDD_00UTC_24km_v1.3.asc.gz`) and IMS 24-km lat/lon file.
- M36 tilecoord (used by regridding script).
- GEOS-LDAS daily model outputs for:
  - `LS_OLv8_M36`
  - `LS_DAv8_M36`

## Outputs

- IMS annual files (example):
  - `projects/IMS/output/ims_snowcover_24km_YYYY.nc4`
- Regridded IMS annual files:
  - `projects/IMS/output/ims_snowcover_24km_YYYY_on_m36_nearest.nc4`
- Validation outputs:
  - `projects/IMS/output/outputs_ims_ol_da_validation/`
  - includes metric caches/tables such as:
    - `ims_ol_da_cell_counts_metrics_*.nc4`
    - `ims_ol_da_scope_metadata_*.csv`
    - `ims_ol_da_comparison_table_*.parquet/.csv`
    - `ims_ol_da_pair_daily_*.parquet/.csv`
    - `ims_ol_da_daily_counts_*.parquet/.csv`

## Notes

- Current default rerun path is script-first (`run_ims_ol_da_cell_metrics.py`) plus plotting notebook.
- `run_ims_ol_da_cell_metrics.py --fast-mode` reuses cached counts and skips model re-read.
- `--min-ims-snow-days` controls analysis domain filtering (default `10`).
