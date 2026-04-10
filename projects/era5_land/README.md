# ERA5 / ERA5-Land Comparison Workflow

This project compares GEOS-LDAS OL/DA monthly outputs against ERA5 and ERA5-Land references, then builds periodized metrics and figures.

## Latest notebook run order

1. Prepare reference monthly files (as needed):
   - `notebooks/download_ERA5_monthly.ipynb`
   - `notebooks/download_ERA5_land.ipynb`

2. Build strict regridded model-vs-reference summaries:
   - `notebooks/compare_with_reanalysis_strict.ipynb`
   - Set:
     - `MODEL_KIND` = `OL` or `DA`
     - `REF_KIND` = `ERA5` or `ERA5-Land`
   - Produces:
     - `projects/era5_land/notebooks/ERA5_vs_OLv8_M36_strict_summary.nc`
     - `projects/era5_land/notebooks/ERA5_vs_DAv8_M36_strict_summary.nc`
     - `projects/era5_land/notebooks/ERA5L_vs_OLv8_M36_strict_summary.nc`
     - `projects/era5_land/notebooks/ERA5L_vs_DAv8_M36_strict_summary.nc`

3. Build periodized tables and figures:
   - `notebooks/plot_ERA5_comparison.ipynb`
   - Outputs include:
     - `ERA5_periodized_metrics_summary.csv`
     - `ERA5_Land_periodized_metrics_summary.csv`
     - map/time-series comparison figures

## Inputs used by strict workflow

- GEOS-LDAS monthly model files for OL and DA (configured in notebook).
- ERA5 merged monthly NetCDF.
- ERA5-Land merged monthly NetCDF.
- Regridding weights (reused or created):
  - `weights_era5_to_m36_consnormed.nc`
  - `weights_era5l_to_m36_consnormed.nc`

## Notes

- `plot_ERA5_comparison.ipynb` currently includes logic to treat ERA5-Land no-snow NaNs as zeros for SWE/snow depth before statistics.
- Legacy notebooks are archived under `notebooks/legacy/` (including `compare_with_ERA5*.ipynb` and `compare_with_ERA5_Land_*` variants).
- `compare_with_reanalysis_strict.ipynb` is the current default workflow.
