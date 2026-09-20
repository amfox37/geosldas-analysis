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

4. Build the manuscript figures (Figures 11, 12 and 13):
   - `projects/M21C_ls/notebooks/paper_figures_unified.ipynb`, cells 33-39
     (cell 33 = ERA5-Land helpers, 35 = Fig. 11, 37 = Fig. 12, 39 = Fig. 13).
     That notebook is the single generator for all manuscript figures; it writes
     PNG + PDF + a per-figure stats CSV into `projects/M21C_ls/output/paper_figures/`
     and records provenance in `paper_figures_manifest.csv`. The copies under
     `projects/M21C_ls/docs/paper_figures/` (PNG only) are what the manuscript
     references, so sync them after regenerating.
   - It reads the cached periodized metrics written by
     `notebooks/plot_ERA5L_comparison_bars.ipynb` (`cache/era5l_periodized_metrics_bars/`),
     so run that notebook first if the cache is absent.

   **Common support.** Figures 11 and 13 average OL and DA over the cells at which
   *both* experiments yield a finite metric (`mean_se_shared` in cell 33), applied per
   period, per layer and per metric. This is necessary because `mask_both` in the strict
   summaries is built from each experiment's own soil temperature and snow-cover fraction,
   so DA and OL do not share a comparison mask; averaging each run over its own finite set
   mixes the assimilation signal with a difference in support. No land-fraction or area
   weighting is applied anywhere in this pathway.

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
