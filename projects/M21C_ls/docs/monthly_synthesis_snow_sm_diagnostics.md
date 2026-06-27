# Monthly Synthesis Snow/Soil-Moisture Diagnostics

Companion notes for:

- `projects/M21C_ls/notebooks/monthly_synthesis_snow_sm_diagnostics.ipynb`
- output directory:
  `projects/M21C_ls/output/monthly_synthesis_diagnostics/`

## Purpose

This exploratory workflow addresses coauthor questions about whether MODIS
snow-cover assimilation connects to the broader observing-system evolution
story. The notebook tests three monthly/seasonal diagnostics:

1. whether snow DA affects seasonal soil moisture before microwave
   soil-moisture DA begins;
2. whether snow DA and soil-moisture DA project onto ET or evaporative fraction
   when monthly flux variables are available;
3. whether snow DA activity changes the later monthly soil-moisture DA work
   required during microwave assimilation periods.

These diagnostics are exploratory. They are intended to identify possible
main-text or supplemental figure candidates, not to provide final validation.

## Monthly-Only Constraint

The local products used here are monthly. The notebook must not be interpreted
as a melt-out timing or event-scale diagnostic.

Use language such as:

- seasonal propagation;
- spring-to-summer carryover;
- monthly/seasonal association;
- spring snow-season DA activity;
- later warm-season soil-moisture response.

Avoid language such as:

- 0-30 days after melt-out;
- post-melt daily response;
- precise melt timing.

## Local Inputs

Current local monthly files under
`/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`:

- `OLv8_land_variables_2000_2024_compressed.nc`
- `DAv8_land_variables_2000_2024_compressed.nc`
- `LS_monthly_increments_2000_2024.nc`
- `spatial_stats_OL_200006_202405.pkl`
- `spatial_stats_DA_200006_202405.pkl`

The compressed OL/DA land-variable files currently include:

- `SFMC`: surface soil moisture, `m3 m-3`;
- `RZMC`: root-zone soil moisture, `m3 m-3`;
- `PRECTOTCORRLAND`: total precipitation over land, `kg m-2 s-1`;
- `FRLANDSNO`: fractional snow-covered land area, unitless;
- `TSOIL1`: soil temperature layer 1, `K`;
- `SNOMASLAND`: total snow storage over land, `kg m-2`;
- `SNODPLAND`: snow depth within snow-covered land fraction, `m`;
- coordinates: `time`, `lat`, `lon`.

The monthly increment file currently includes:

- `SFMC_INC`: surface soil-moisture increment, `m3 m-3`;
- `RZMC_INC`: root-zone soil-moisture increment, `m3 m-3`;
- `SNOWMASS_INCR`: snow-mass increment, `kg m-2`.

The current local compressed monthly land-variable files do not include ET,
latent heat, or sensible heat variables. Therefore Analysis 2 is explicitly
skipped until monthly flux files are available.

## Periods and Windows

Main periods:

| ID | Dates | Meaning |
| --- | --- | --- |
| `modis_only` | 2000-06-01 to 2007-05-31 | MODIS snow-cover DA only; no microwave soil-moisture DA |
| `pre_smap_mw` | 2007-06-01 to 2015-03-31 | ASCAT active; SMOS begins January 2010 |
| `smap_era` | 2015-04-01 to 2024-05-31 | SMAP-era microwave period in the available monthly files |

Seasonal windows:

- predictor snow window: `MAM`;
- soil-moisture response windows: `AMJ`, `MJJ`, and exploratory `JJA`;
- soil-moisture DA work windows: `AMJ`, `MJJ`, and `JJA`.

For Analysis 1, the notebook uses years 2001-2006. This avoids using AMJ/MJJ
response months after ASCAT begins in June 2007.

For Analysis 3, the first pass uses:

- `pre_smap_mw`: years 2008-2014;
- `smap_era`: years 2016-2023.

## Masks

The first-pass Northern Hemisphere seasonal-snow mask is:

- latitude > 20 N;
- maximum monthly OL/DA `FRLANDSNO` over the full record > 0.05;
- mean JJA OL/DA `FRLANDSNO` over the full record < 0.20.

This intentionally removes most snow-free and permanent-snow/glacier-like
cells. Thresholds are simple first-pass values and should be revisited if a
figure moves toward publication.

## Sign Conventions

- Signed state response is `DA - OL`.
- Positive `DA - OL SFMC` or `DA - OL RZMC` means DA is wetter than OL.
- Positive `DA - OL FRLANDSNO` means DA has more snow-covered area than OL.
- Positive `DA - OL SNOMASLAND` means DA has more snow storage than OL.
- Absolute snow DA activity uses `abs(DA - OL FRLANDSNO)` or
  `abs(DA - OL SNOMASLAND)`.
- Soil-moisture DA work in Analysis 3 uses monthly increments:
  `SFMC_INC` and `RZMC_INC`. Absolute work uses the absolute value of the
  monthly increment before seasonal averaging, so opposite-signed monthly
  corrections do not cancel.

`DA - OL` snow variables are proxies for snow assimilation impact, not actual
analysis increments. The increment file does include `SNOWMASS_INCR`, but the
first notebook pass keeps the predictor tied to DA-minus-OL snow-state
differences so the diagnostic remains comparable with Analysis 1 maps.

## Outputs

Inventory and notes:

- `monthly_synthesis_input_inventory.csv`
- `monthly_synthesis_availability_summary.csv`
- `monthly_synthesis_mask_summary.csv`

Analysis 1:

- `analysis1_gridcell_year_snow_to_sm_table.csv`
- `analysis1_binned_snow_activity_vs_sm_response.csv`
- `analysis1_seasonal_snow_to_sm_maps.nc`
- `figures/analysis1_maps_snow_to_sm_response.{png,pdf}`
- `figures/analysis1_binned_abs_scf_vs_sm_response.{png,pdf}`

Analysis 2:

- `analysis2_flux_response_status.csv`

Analysis 3:

- `analysis3_gridcell_year_snow_activity_vs_sm_work_table.csv`
- `analysis3_binned_snow_activity_vs_sm_da_work.csv`
- `analysis3_smapera_snow_activity_and_sm_work_maps.nc`
- `figures/analysis3_binned_abs_scf_vs_rzmc_increment_work.{png,pdf}`
- `figures/analysis3_smapera_maps_snow_activity_and_rzmc_work.{png,pdf}`

Recommendation scaffold:

- `monthly_synthesis_recommendation_quicklook.csv`

## Interpretation Guide

Analysis 1 is the cleanest causal seasonal diagnostic. If MAM snow DA activity
corresponds to AMJ/MJJ SFMC or RZMC differences during 2001-2006, that supports
seasonal propagation from MODIS snow-cover DA into soil moisture.

Analysis 3 directly addresses whether snow DA activity changes later
soil-moisture DA work. If high snow-activity bins have smaller later increments,
snow DA may reduce later microwave correction burden. If high snow-activity
bins have larger later increments, snow DA may identify places where snow
process errors remain hydrologically important. If the relationship is flat,
snow DA and microwave SM DA are largely decoupled at monthly scale.

Analysis 2 should remain skipped or clearly labeled incomplete until monthly
flux files with ET/LE/H are copied in and their units/sign conventions are
confirmed.

Do not overclaim DA-minus-OL differences as improvement without independent
validation.
