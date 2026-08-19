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

## Analysis A Robustness Workflow

The publication-oriented Analysis A now has a separate falsification workflow:

- contract: `config/analysis_a_robustness.json`;
- implementation: `scripts/analysis_a_robustness.py`;
- derived output: `output/monthly_synthesis_diagnostics/analysisA_robustness/`;
- report section: `docs/monthly_synthesis_report_out.md`.

This workflow does not replace or alter the original notebook products. It
first rebuilds the original 2001-2006 tile-year sample and asserts that its
eight-bin Figure 2 statistics match the saved product. It also fails if a date
escapes the MODIS-only period or if MAM `snow_abs_netpack` is smaller than the
absolute MAM `snow_net` value for any finite tile-year.

For each of six responses, the workflow uses one common complete-case sample
for its magnitude and signed formulations and requires at least four of the six
years per tile. In the current inputs all selected values are finite, so the
restriction retains all 48,067 seasonal-snow tiles and 288,402 tile-years.

The fitted sequence is:

- `M0`: pooled tile-year relationship;
- `M1`: tile-demeaned relationship;
- `M2`: `M1` plus year fixed effects;
- `M3`: `M2` plus the within-tile anomaly of OL MAM `SNOMASLAND`.

The between-tile coefficient uses tile means weighted by valid-year count.
Magnitude fits pair MAM `snow_abs_netpack` with the original absolute `DA - OL`
response definition. Signed fits pair MAM `snow_net` with signed `DA - OL`
responses. M4 is intentionally omitted because the local monthly inputs do not
contain a clean pre-assimilation, tile-level MODIS availability count.

The practical identifying-variation threshold was fixed before fitting the
responses at a within-tile predictor SD of `0.1 kg m-2`. Identification is
declared adequate only when at least half the eligible tiles exceed that value
and within-tile variance is at least 5% of total variance.

Uncertainty uses 1,000 spatial-block bootstrap replicates that retain all years
and tiles in each sampled block. The primary blocks are approximately 5 by 5
degrees; 10-degree blocks are a sensitivity. Additional sensitivities remove
predictor anomalies below the 1st or above the 99th percentile and compare
M1, M2, and M3 directly. Individual six-year tile correlations are descriptive
and are never assigned tile-level significance.

The predeclared classification is applied only to MJJ ET, total runoff, total
water, and RZMC. The current production result is A: all four have positive M3
coefficients, 5-degree block intervals above zero, retained fractions above
40%, and directionally consistent positive signed M3 coefficients. This is
evidence for a physically coherent modeled propagation pathway in six
MODIS-only seasons, not causal proof or independent validation.

## Water-Year Differential Budget

The follow-on workflow is fixed by
`config/water_year_snow_da_budget.json` and implemented in
`scripts/water_year_snow_da_budget.py`. It uses conventional October-September
WY2001-WY2006 and the same 48,067-tile Northern Hemisphere seasonal-snow mask.
Every tile-year retains the native signed `snow_net` input, absolute snow
activity, snowmelt, infiltration, ET, surface runoff, baseflow, total runoff,
storage change, residual, and monthly SFMC/RZMC trajectories.

The closing equation is `I_snow = dET + dRunoff + dStorage + dPeatFreeWater +
residual`. `dStorage` is the change in `DA - OL TWLAND` between instantaneous
00Z 1 October restart endpoints. `dPeatFreeWater` is the change in PEATCLSM
free-standing surface water, a store that `TWLAND` excludes by construction. Integrated `DA - OL WCHANGELAND` is retained separately because the
source audit showed that it closes the model-process tendency balance but
omits the discontinuous analysis injection. Snowmelt and infiltration are
pathway diagnostics rather than additional closing terms. SFMC and RZMC are
state diagnostics and are also excluded because total land-water storage
already contains soil water.

The compressed OL and DA precipitation fields are conceptually common forcing
but are not bit-identical. The production gate allows at most `0.2 kg m-2` for
an annual tile discrepancy and reports the observed maximum (`0.102 kg m-2`)
and maximum annual area-weighted domain discrepancy (`0.000127 kg m-2`). This
retains the known float32 compression artifact without silently accepting a
material forcing mismatch.

For snow-addition tile-years, the six-year direct partition is 64.3% runoff
(43.1% surface runoff and 21.2% baseflow), 35.9% ET, 4.2% storage change,
-2.7% peatland free-standing water, and -1.7% residual. The 5-degree spatial-block 95% interval for total
runoff is 61.1-67.4%; 10-degree blocks are retained as a sensitivity. Annual
domain runoff fractions span 55.6-70.1% across the six water years.

RZMC and SFMC describe the intermediate state response. Snow-addition
tile-years have an area-weighted peak RZMC difference of `0.0189 m3 m-3`, MJJ
mean difference of `0.0108 m3 m-3`, and September difference of `0.0082 m3
m-3`. Persistence is right-censored in 67.2% of cases at the predeclared
`0.001 m3 m-3` threshold. Because the mean RZMC difference is already positive
in October, persistence includes inherited DA-OL state differences and is not
the residence time of only the current water-year snow increment.

The derived tile archive and full CSV outputs live under
`output/monthly_synthesis_diagnostics/water_year_snow_da_budget/`. The archive
is reused by default; pass `--rebuild-tile-budget` only when the source files or
mask definition change. Permanent figures and the interpreted result are in
`docs/monthly_synthesis_report_out.md` and its PDF rendering.
