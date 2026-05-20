# ERA5-Land comparison methods

This note documents the ERA5/ERA5-Land reanalysis-comparison workflow behind
the manuscript ERA5-Land figures:

- `projects/era5_land/notebooks/compare_with_reanalysis_strict.ipynb`
- `projects/era5_land/notebooks/plot_ERA5L_comparison_bars.ipynb`
- `projects/era5_land/notebooks/plot_ERA5L_postage_stamp_maps.ipynb`
- `projects/era5_land/notebooks/figures_era5l_bars/bars_surface_rz_sm_combined_3x3_era5l_only.png`
- `projects/era5_land/notebooks/figures_era5l_postage_stamps/postage_stamp_da_minus_ol_surface_rz_sm_6x3.png`
- `projects/era5_land/notebooks/figures_era5l_bars/bars_snow_combined_3x3_era5l_only.png`

The strict comparison notebook can build both ERA5 and ERA5-Land summaries.
However, the paper figures shown here use ERA5-Land only. The final plotting
notebooks read:

- `projects/era5_land/notebooks/ERA5L_vs_OLv8_M36_strict_summary.nc`
- `projects/era5_land/notebooks/ERA5L_vs_DAv8_M36_strict_summary.nc`

The ERA5 summaries are documented below because they were part of the workflow
and are useful provenance, but they are not the reference used for the paper
bar and postage-stamp figures described in this note.

## Reference data provenance

ERA5 and ERA5-Land monthly reference fields were downloaded from the Copernicus
Climate Data Store using:

- `projects/era5_land/notebooks/download_ERA5_monthly.ipynb`
- `projects/era5_land/notebooks/download_ERA5_land.ipynb`

The ERA5 notebook uses CDS dataset:

- `reanalysis-era5-single-levels-monthly-means`

and requests monthly averaged fields for:

- `soil_temperature_level_1`
- `snow_depth`
- `snow_density`
- `volumetric_soil_water_layer_1`
- `volumetric_soil_water_layer_2`
- `volumetric_soil_water_layer_3`

The ERA5-Land notebook uses CDS dataset:

- `reanalysis-era5-land-monthly-means`

and the strict comparison file requires the merged ERA5-Land fields represented
with short names:

- `swvl1`, `swvl2`, `swvl3`: volumetric soil water layers
- `stl1`: soil temperature layer 1
- `snowc`: snow cover, stored as percent in the source file used here
- `sd`: snow depth water equivalent, used as SWE in meters of water equivalent
- `sde`: snow depth in meters

The ERA5-Land download/merge workflow includes yearly files for 2000-2025 and
the strict comparison used:

- `/Users/amfox/Desktop/geosldas-analysis/ERA5L_monthly_merged_rebuilt.nc`

The ERA5 comparison used:

- `/Users/amfox/Desktop/geosldas-analysis/era5_monthly_nc/ERA5_monthly_merged.nc`

ERA5 and ERA5-Land are living CDS archives in this workflow rather than locally
versioned products. For manuscript wording, report the data as CDS ERA5-Land
monthly averaged reanalysis fields from the local merged snapshot used in the
April 2026 analysis. Add the exact CDS access date if needed by the journal.

## Model data

The model inputs are monthly GEOS-LDAS OL and DA files:

- OL: `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/OLv8_land_variables_2000_2024_compressed.nc`
- DA: `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/DAv8_land_variables_2000_2024_compressed.nc`

The model grid is `SMAP_EASEv2_M36_GLOBAL`. The strict notebook reads the
GEOS-LDAS tile-coordinate file and maps tile data to the M36 grid. The code
aggregates finite tiles within each M36 cell using an unweighted mean, then
masks cells with no finite tile values.

Model variables used in the ERA5-Land figures are:

- `SFMC`: surface soil moisture
- `RZMC`: root-zone soil moisture
- `FRLANDSNO`: snow-cover fraction
- `TSOIL1`: upper-layer soil temperature, used for warm-soil masking
- `SNOMASLAND`: SWE in `kg m-2`, converted to meters water equivalent by
  dividing by 1000
- `SNODPLAND`: snow depth in meters; the plotted grid-cell snow depth is
  `SNODPLAND * FRLANDSNO`

## Regridding and collocation

Reference fields are regridded to the M36 grid using xESMF
`conservative_normed` regridding with periodic longitude handling. The strict
workflow creates or reuses source-specific weights:

- ERA5: `weights_era5_to_m36_consnormed.nc`
- ERA5-Land: `weights_era5l_to_m36_consnormed.nc`

ERA5-Land regridding is chunked over time in 6-month blocks to keep memory use
manageable. Source latitude order is normalized, longitude is wrapped to
`[0, 360)`, and the model and reference records are intersected by monthly
period. The aligned strict outputs contain:

- `time=288`, from `2000-06-30` to `2024-05-31`
- `y=379`, `x=964`

Temporal collocation is exact monthly collocation after intersecting common
model/reference months. No interpolation in time is applied.

## ERA5 versus ERA5-Land variable handling

The strict notebook has explicit branches for ERA5 and ERA5-Land:

- ERA5 snow cover and snow depth are derived from `sd` and `rsn`.
- ERA5-Land snow cover is read directly from `snowc` and divided by 100.
- ERA5-Land SWE is `sd`.
- ERA5-Land snow depth is read directly from `sde`.

For both products, root-zone soil moisture is computed from the first three
soil layers over 0-100 cm:

- `RZ = 0.07 * swvl1 + 0.21 * swvl2 + 0.72 * swvl3`

The final manuscript figures use only the ERA5-Land branch.

## Masks and missing data

Surface and root-zone soil moisture metrics use a joint warm, snow-free mask:

- model soil temperature `TSOIL1 > 275.15 K`
- reference soil temperature `stl1 > 275.15 K`
- model snow-cover fraction `< 0.01`
- reference snow-cover fraction `< 0.01`

The mask used for soil-moisture metrics is the intersection of the model and
reference masks. This removes frozen or snow-covered conditions from the soil
moisture comparison.

For the final ERA5-Land plotting notebooks, snow metrics use an ever-snow
domain mask over the full record:

- include a grid cell if either model or ERA5-Land snow-cover fraction exceeds
  `1e-6` in any month

The ERA5-Land plotting notebooks also fill ERA5-Land no-snow missing values in
`SWE_era` and `SNWD_era` with zero before computing SWE and snow-depth metrics.
All metrics use finite paired model-reference samples only. A grid-cell metric
is retained only if at least 24 valid monthly pairs are available.

## Periods

The ERA5-Land plotting notebooks compute periodized metrics for:

- `2000-06-01` to `2007-05-31`
- `2007-06-01` to `2015-03-31`
- `2015-04-01` to `2024-05-31`
- `2000-06-01` to `2024-05-31`

The soil-moisture bar and postage-stamp figures use the first three periods.
The snow bar figure uses the full period, `2000-06-01` to `2024-05-31`.

## Metrics

Metrics are computed per M36 grid cell from monthly paired samples.

For surface and root-zone soil moisture, the paper figures use:

- `R`: Pearson correlation of raw monthly values
- `anomR`: anomaly correlation after removing month-of-year climatology within
  each period
- `ubRMSE`: root-mean-square difference of the monthly anomalies

For snow-cover fraction, SWE, and snow depth, the paper snow bar figure uses:

- `RMSE = sqrt(mean((model - ERA5-Land)^2))`
- `ubRMSE`: root-mean-square difference of monthly anomalies
- `Bias = mean(model - ERA5-Land)`

The anomaly definition is month-of-year climatology removal within the selected
period. Metrics are computed independently for OL and DA against ERA5-Land.

## Figure-specific summaries

The surface/root-zone soil-moisture bar figure shows spatial means with
standard errors over finite M36 grid cells. Surface bars are solid; root-zone
bars are hatched. The valid finite-cell counts printed in the current figure
range from about 69,700 to 75,200 depending on metric and period.

The postage-stamp map figure shows `DA - OL` for surface and root-zone
soil-moisture metrics across the three periods. Positive values indicate higher
DA correlation or anomaly correlation; negative values indicate lower DA
ubRMSE. Panel annotations are land-fraction-weighted spatial means, computed
from `frac_cell` in the GEOS-LDAS tile-coordinate file. The maps are shown from
60 S to 90 N using the canonical M36 EASE latitude/longitude arrays.

The snow bar figure shows spatial means with standard errors for the full
period. Rows are snow-cover fraction, SWE, and snow depth; columns are RMSE,
ubRMSE, and bias. The current figure reports `n=52,499` cells for snow-cover
fraction and `n=62,294` cells for SWE and snow depth.

## Current output files

Strict ERA5-Land summary files:

- `projects/era5_land/notebooks/ERA5L_vs_OLv8_M36_strict_summary.nc`,
  last modified locally `2026-04-03 15:10:51`
- `projects/era5_land/notebooks/ERA5L_vs_DAv8_M36_strict_summary.nc`,
  last modified locally `2026-04-03 15:27:56`

Strict ERA5 summary files, documented but not used for the final paper figures:

- `projects/era5_land/notebooks/ERA5_vs_OLv8_M36_strict_summary.nc`,
  last modified locally `2026-04-03 14:53:44`
- `projects/era5_land/notebooks/ERA5_vs_DAv8_M36_strict_summary.nc`,
  last modified locally `2026-04-03 14:46:34`

Figure files:

- `projects/era5_land/notebooks/figures_era5l_bars/bars_surface_rz_sm_combined_3x3_era5l_only.png`,
  last modified locally `2026-04-12 16:27:03`
- `projects/era5_land/notebooks/figures_era5l_postage_stamps/postage_stamp_da_minus_ol_surface_rz_sm_6x3.png`,
  last modified locally `2026-04-13 14:27:22`
- `projects/era5_land/notebooks/figures_era5l_bars/bars_snow_combined_3x3_era5l_only.png`,
  last modified locally `2026-04-12 16:31:55`

## Draft methods text

Monthly GEOS-LDAS OL and DA fields were evaluated against ERA5-Land monthly
averaged reanalysis fields from the Copernicus Climate Data Store. Although the
workflow can also generate ERA5 comparisons, the manuscript figures shown here
use ERA5-Land only. ERA5-Land variables were regridded to the GEOS-LDAS
`SMAP_EASEv2_M36_GLOBAL` grid using conservative-normalized xESMF regridding.
Model tile outputs were aggregated to the same M36 cells, and model and
reference fields were matched by common monthly period from June 2000 through
May 2024.

Surface soil moisture was compared using GEOS-LDAS `SFMC` and ERA5-Land
`swvl1`. Root-zone soil moisture was compared using GEOS-LDAS `RZMC` and an
ERA5-Land 0-100 cm estimate computed as `0.07*swvl1 + 0.21*swvl2 +
0.72*swvl3`. Soil-moisture statistics were computed only where both model and
ERA5-Land indicated warm, snow-free conditions, using a soil-temperature
threshold of 275.15 K and snow-cover-fraction threshold of 0.01.

Snow-cover fraction, SWE, and snow depth were evaluated against ERA5-Land
`snowc`, `sd`, and `sde`, respectively. Model SWE was `SNOMASLAND / 1000`, and
model grid-cell snow depth was `SNODPLAND * FRLANDSNO`. ERA5-Land no-snow
missing values in SWE and snow depth were treated as zero in the plotting
workflow. Snow metrics were restricted to grid cells with snow in either the
model or ERA5-Land during at least one month of the full record.

For each grid cell and period, metrics were computed from finite paired monthly
samples, requiring at least 24 valid months. Correlations were computed from raw
monthly values. Anomaly correlations and unbiased RMSE were computed after
removing the month-of-year climatology within each period. Bar plots show
spatial means and standard errors over finite grid cells. Map panels show
periodized `DA - OL` differences, so positive values indicate higher DA
correlation metrics and negative values indicate lower DA ubRMSE.

## Items still worth verifying

- Add the formal CDS ERA5-Land citation and exact CDS access date required by
  the manuscript or journal.
- Confirm final figure numbering, since the manuscript draft has shifted figure
  numbers across the snow and reanalysis sections.
