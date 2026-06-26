# P1 vs P2 Comparison Diagnostics

This note documents the dedicated P1/P2 diagnostics notebook:

- `projects/M21C_ls/notebooks/p1_p2_comparison_diagnostics.ipynb`
- outputs: `projects/M21C_ls/output/p1_p2_comparison/`

The notebook is an exploratory comparison layer for the MODIS Terra-only to
Terra+Aqua transition. It does not replace the manuscript figure notebook. It
pulls together the products currently available for the two early MODIS snow
periods and applies common-support comparisons wherever possible.

## Periods

| Period | Dates | Meaning | Days | Months |
| --- | --- | --- | ---: | ---: |
| P1 | 2000-06-01 to 2002-06-30 | MODIS Terra SCF | 760 | 25 |
| P2 | 2002-07-01 to 2007-05-31 | MODIS Terra+Aqua SCF | 1796 | 59 |

These labels are intentionally aligned with the plain-language period naming in
the manuscript workflow. P1 and P2 together correspond to the broader
`SCF-only period` used in the main ISMN and ERA5-Land figures.

## Inputs

The notebook writes an input inventory table:

- `output/p1_p2_comparison/input_availability.csv`

Current input products:

- OFA period-split temporal stats:
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/temporal_stats_{OL,DA}_20000601_20020630.nc4`
  and
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/temporal_stats_{OL,DA}_20020701_20070531.nc4`
- M36 tile coordinates:
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/LS_OLv8_M36.ldas_tilecoord.bin`
- IMS Terra/Terra+Aqua scope product:
  `projects/IMS/output/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2007_thr0p50_imsSnowDaysGe10_terraAquaScopes.nc4`
- SNOTEL raw time-series cache:
  `projects/SNOTEL/outputs_snotel_ol_da_validation/snotel_raw_timeseries_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.nc`
- GHCN snow-depth raw time-series cache:
  `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_snwd_raw_timeseries_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.nc`
- ERA5-Land strict monthly aligned summaries:
  `projects/era5_land/notebooks/ERA5L_vs_OLv8_M36_strict_summary.nc`
  and
  `projects/era5_land/notebooks/ERA5L_vs_DAv8_M36_strict_summary.nc`

## Common-Support Rules

The intent is to avoid comparing P1 and P2 over different spatial or station
samples unless that is unavoidable.

- OFA: common tile support across P1, P2, OL, and DA, with MODIS weighted
  observation support and finite OmF standard deviation. The current common
  support has 55,918 tiles.
- IMS: metric-specific valid cells in each period from the Terra/Terra+Aqua
  scope NetCDF, with the snow-domain mask inherited from the upstream IMS
  product.
- SNOTEL and GHCN: stations must meet the paired-sample threshold in both
  periods and both experiments. The current threshold is 100 paired days.
- ERA5-Land: grid cells must be valid for P1, P2, OL, and DA for each variable
  and metric. The current minimum paired monthly support is 18 months for P1
  and 36 months for P2.

Positive values mean DA is better throughout the notebook:

- `DA - OL` for higher-is-better metrics such as `R`, IMS accuracy, IMS hit
  rate, and IMS correct rejection rate.
- `OL - DA` for lower-is-better metrics such as RMSE, ubRMSE, absolute bias,
  IMS miss rate, and IMS false alarm ratio.
- OFA OmF standard-deviation improvement is `(OL - DA) / OL * 100`.

## Output Tables

Primary tables:

- `p1_p2_cross_dataset_rollup.csv`: compact cross-dataset summary.
- `p1_p2_ofa_summary.csv`: MODIS observation counts and OmF standard-deviation
  diagnostics on common tile support.
- `p1_p2_ims_summary.csv`: IMS categorical snow-cover skill improvements.
- `p1_p2_station_metrics_common_support.csv`: paired OL/DA station metrics for
  SNOTEL and GHCN on common support.
- `p1_p2_station_delta_summary.csv`: station-level DA-vs-OL delta summaries.
- `p1_p2_era5land_delta_summary.csv`: ERA5-Land gridcell DA-vs-OL delta
  summaries.

Figures are written as PNG and PDF under:

- `output/p1_p2_comparison/figures/`

Current figure products:

- `p1_p2_ofa_modis_maps.*`
- `p1_p2_ims_skill_maps.*`
- `p1_p2_snotel_swe_ubrmse_station_maps.*`
- `p1_p2_snotel_snow_depth_ubrmse_station_maps.*`
- `p1_p2_ghcn_snow_depth_ubrmse_station_maps.*`
- `p1_p2_era5land_sm_ubrmse_maps.*`
- `p1_p2_era5land_rz_ubrmse_maps.*`
- `p1_p2_era5land_scf_rmse_maps.*`
- `p1_p2_era5land_swe_rmse_maps.*`
- `p1_p2_era5land_snwd_rmse_maps.*`

## Main Results

### MODIS OFA

| Metric | P1 | P2 |
| --- | ---: | ---: |
| DA observations per day on common tiles | 20,477/day | 40,030/day |
| Spatial mean OmF StdDev, OL | 0.2329 | 0.2475 |
| Spatial mean OmF StdDev, DA | 0.2100 | 0.2087 |
| Spatial mean OmF StdDev improvement | 8.28% | 13.93% |

Interpretation: P2 has almost twice the MODIS observation density per day on
the same common tiles, and DA reduces MODIS OmF standard deviation more strongly
in P2. This is the cleanest direct signal that the Terra+Aqua period provides
more useful SCF information to the assimilation.

### IMS Snow-Cover Skill

| Metric | P1 mean improvement | P1 cells improved | P2 mean improvement | P2 cells improved |
| --- | ---: | ---: | ---: | ---: |
| Accuracy | +0.0082 | 51.3% | +0.0128 | 62.3% |
| Hit rate | +0.0433 | 63.0% | +0.0655 | 81.5% |
| False alarm ratio | +0.0061 | 32.2% | +0.0126 | 30.9% |
| Correct rejection rate | -0.0061 | 18.8% | -0.0101 | 19.9% |

Interpretation: DA improves snow detection more in P2, especially hit rate.
The tradeoff is that correct rejection rate is lower in both periods and false
alarm ratio is spatially mixed. The IMS story is therefore not simply "better
everywhere"; it is more specifically better detection/hit behavior with some
cost in no-snow classification.

### SNOTEL and GHCN

| Dataset | Variable | Metric | P1 mean improvement | P2 mean improvement |
| --- | --- | --- | ---: | ---: |
| SNOTEL | SWE | R | +0.0223 | +0.0271 |
| SNOTEL | SWE | RMSE | +0.6409 kg m-2 | +0.7870 kg m-2 |
| SNOTEL | SWE | ubRMSE | +0.2830 kg m-2 | +0.2880 kg m-2 |
| SNOTEL | SWE | abs bias | +0.6892 kg m-2 | +0.9239 kg m-2 |
| SNOTEL | Snow depth | R | +0.0045 | -0.0050 |
| SNOTEL | Snow depth | RMSE | +0.0053 m | +0.0054 m |
| SNOTEL | Snow depth | ubRMSE | +0.0001 m | -0.0005 m |
| SNOTEL | Snow depth | abs bias | +0.0084 m | +0.0100 m |
| GHCN | Snow depth | R | +0.0132 | +0.0113 |
| GHCN | Snow depth | RMSE | +0.7621 mm | +1.6661 mm |
| GHCN | Snow depth | ubRMSE | +0.8071 mm | +1.5894 mm |
| GHCN | Snow depth | abs bias | +0.1872 mm | +0.6698 mm |

Interpretation: SNOTEL SWE is consistently improved in both periods and is
slightly stronger in P2. GHCN snow-depth improvements are also larger in P2,
especially RMSE and ubRMSE. SNOTEL snow depth is small and mixed: bulk RMSE and
absolute bias improve, but P2 correlation and ubRMSE are slightly worse.

### ERA5-Land

| Variable | Metric | P1 mean improvement | P2 mean improvement |
| --- | --- | ---: | ---: |
| Surface soil moisture | R | +0.0034 | +0.0072 |
| Surface soil moisture | RMSE | +0.00169 m3 m-3 | +0.00342 m3 m-3 |
| Surface soil moisture | ubRMSE | +0.00015 m3 m-3 | +0.00010 m3 m-3 |
| Root-zone soil moisture | R | -0.0034 | +0.0002 |
| Root-zone soil moisture | RMSE | +0.00120 m3 m-3 | +0.00236 m3 m-3 |
| Root-zone soil moisture | ubRMSE | -0.00013 m3 m-3 | -0.00010 m3 m-3 |
| Snow-cover fraction | R | +0.0177 | +0.0259 |
| Snow-cover fraction | RMSE | +0.00332 | +0.00457 |
| SWE | R | +0.0029 | +0.0114 |
| SWE | RMSE | +0.00014 m | +0.00003 m |
| Snow depth | R | +0.0047 | +0.0117 |
| Snow depth | RMSE | +0.00273 m | +0.00358 m |

Interpretation: ERA5-Land generally agrees that P2 is at least as favorable as
P1, especially for surface soil moisture, snow-cover fraction, and snow-depth
RMSE/bias. The caveat is that ubRMSE is flat or slightly negative for some
root-zone and snow variables. ERA5-Land therefore supports a cautious story:
DA improves bulk agreement and bias/RMSE structure more clearly than it reduces
all random-error components.

## Overall Interpretation

The overall picture is that the Terra+Aqua period is stronger than the
Terra-only period:

1. P2 has about twice as many MODIS observations per day on common tiles.
2. The direct MODIS OFA diagnostic shows stronger OmF standard-deviation
   reduction in P2.
3. IMS snow-cover skill improves more in P2, particularly hit rate.
4. Independent snow validation mostly supports the same direction, especially
   GHCN snow depth and SNOTEL SWE.
5. ERA5-Land provides a more nuanced but broadly consistent check: several
   RMSE, bias, and correlation metrics improve more in P2, while some ubRMSE
   metrics are flat or slightly worse.

This is not a "DA improves every metric everywhere" result. The more defensible
manuscript interpretation is: adding Aqua MODIS increases useful snow-cover
constraint density, strengthens the DA reduction in MODIS OmF variance, and
improves snow detection and several independent validation metrics, but with
metric-specific tradeoffs in no-snow classification and some variance-style
diagnostics.

## Notebook Behavior

The notebook is intended to be readable interactively:

- tables render inline with `display(...)`;
- figures render inline in Jupyter via `show_figure(fig)`;
- all tables and figures are also written to disk for provenance.

For batch/headless execution, the notebook detects non-notebook execution and
uses the `Agg` backend so saved-figure runs do not hang on GUI windows.

The ERA5-Land map cell uses finite-coordinate scatter plotting rather than
`pcolormesh`, because the aligned EASE-grid coordinate arrays can include
masked/non-finite cells that Cartopy refuses as pcolormesh coordinates.

## Reproducibility Notes

The full notebook was executed locally after the current inputs were copied
into place. The final run produced:

- 14 OFA summary rows,
- 10 IMS summary rows,
- 24 station delta-summary rows,
- 31,044 station metric rows,
- 40 ERA5-Land summary rows,
- 76 cross-dataset rollup rows.

Representative OFA, GHCN, and ERA5-Land maps were visually checked after the
final successful run.
