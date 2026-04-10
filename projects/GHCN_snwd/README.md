# GHCN SNWD Validation Workflow

This project validates GEOS-LDAS OL/DA snow depth against GHCN-Daily SNWD observations.

## Latest notebook run order

1. Build GHCN SNWD intermediate products:
   - `notebooks/ghcn_snwd_global_1998_present_build.ipynb`
   - Produces build outputs (manifest, station metadata/coverage parquet, per-year parquet partitions).

2. Run baseline OL/DA validation:
   - `notebooks/ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb`
   - Current default notebook for summary bars/maps and DA-OL metric deltas.

## Inputs used by baseline notebook

- GHCN build outputs (from step 1), typically:
  - `ghcn_snwd_manifest.json`
  - `ghcn_snwd_station_metadata.parquet`
  - `ghcn_snwd_station_coverage.parquet`
  - `ghcn_snwd_parquet/`
- GEOSldas daily model files + tilecoord for:
  - `LS_OLv8_M36`
  - `LS_DAv8_M36`
- Model snow depth is computed as:
  - `SNODPLAND * SCF` (SCF fallback candidates are handled in-notebook).

## Outputs

- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/`
- Typical files:
  - `ghcn_snwd_raw_timeseries_*.nc`
  - `ghcn_station_metrics_baseline_core_*.csv`
  - `ghcn_domain_metrics_baseline_core_*.csv`
  - figures/maps tables written by notebook cells

## Notes

- The notebook is configured to run both on local and Discover-style paths via top-cell overrides.
- The older notebook is archived at:
  - `notebooks/legacy/ghcn_snwd_daily_seasonal_ol_da_snwd.ipynb`
- `ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb` is the maintained baseline workflow.
