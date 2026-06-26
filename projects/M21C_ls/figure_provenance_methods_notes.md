# M21C LS manuscript figure provenance and methods notes

Working note for tracing figures in:

- `/Users/amfox/Documents/Publications/SM and SCF paper/draft_033026.docx`

The Word document embeds resized PNGs, so exact file hashes do not match the
local figure files. The mappings below are based on captions, aspect ratios,
notebook save names, and project READMEs.

## Figure provenance

| Draft figure | Draft content | Likely source workflow | Likely output or notes |
| --- | --- | --- | --- |
| Figure 2 | Mean number of assimilated observations per day for MODIS SCF, ASCAT SM, SMOS Tb, SMAP Tb, and CYGNSS SM — five global postage-stamp maps. | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` | Confirmed: notebook contains exact day-count values (MODIS 8700, ASCAT 6121, SMOS 5047, SMAP 3207, CYGNSS 2131) and title string "Obs/day (mean over N days)". |
| Figure 3 | Temporal evolution of number of assimilated observations: total (top panel) and by observing system (bottom panel), monthly, 2000–2024. | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` | Confirmed: notebook contains title "Number of Observations (DA): LS_DAv8_M36_200006–202405" and "Total Obs / Month". |
| Figure 4 | Spatial maps of O-F StdDev improvement, `(OL−DA)/OL` [%], by sensor (MODIS, ASCAT, SMOS, SMAP, CYGNSS). | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Updated unified workflow makes positive/red values mean lower DA O-F StdDev than OL and masks tiny OL denominators. |
| Figure 5 | Monthly time series of O-F StdDev improvement, `(OL−DA)/OL` [%], for MODIS, ASCAT, SMOS, SMAP, and CYGNSS, 2000–2024. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Updated unified workflow uses P1-P9 shading, positive-improvement convention, 5-month centered smoothing over monthly values, and a period-mean summary strip. |
| Figure 6 | O-F StdDev improvement, `(OL−DA)/OL` [%], by P1-P9 observing-system period and sensor. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Updated unified workflow uses a 9-row x 5-column availability-aware grid for MODIS, ASCAT, SMOS, SMAP, and CYGNSS; positive/red means lower DA O-F StdDev than OL. |
| Figure 7 | In situ soil-moisture skill changes between DA and OL for surface and root-zone soil moisture, grouped by broader validation periods mapped to P1-P9. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Updated unified workflow reads the hybrid ISMN table and labels validation periods as V1 SCF-only period (P1-P2), V2 microwave pre-SMAP period (P3-P5), and V3 SMAP-era microwave period (P6-P9); positive values mean DA better for all metrics. |
| Figure 8 | IMS categorical snow-cover skill differences between DA and OL. | `projects/IMS/scripts/run_ims_ol_da_cell_metrics.py` plus `projects/IMS/notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb` | Plotting notebook saves an `ims_all_period_delta_metrics_*_nh_robinson_2x3.png` style figure from precomputed IMS metric outputs. |
| Figure 9 | Seasonal SWE performance at SNOTEL stations: top-row domain bars and bottom-row station maps of DA improvement. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Unified workflow reads cached SNOTEL station metrics/top-row CI tables and plots map values as `OL - DA` for RMSE, ubRMSE, and absolute bias, so positive/red means DA improved. |
| Figure 10 | GHCN snow depth validation: seasonal station bars and NH station improvement maps. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb` | Rebuilt from cached GHCN baseline-core station metrics. Map colors use `OL - DA`, so positive/red means DA improved lower-is-better snow-depth metrics. Legacy source figure used DA minus OL. |
| Figures 11-13 | ERA5-Land comparison figures for soil moisture and snow variables. | `projects/M21C_ls/notebooks/paper_figures_unified.ipynb`, using ERA5-Land metric caches derived from `projects/era5_land/notebooks/compare_with_reanalysis_strict.ipynb`. | Unified outputs are `fig11_era5land_soil_moisture_bars.*`, `fig12_era5land_soil_moisture_improvement_maps.*`, and `fig13_era5land_snow_comparison_bars.*`. |

## Project README inventory

These project areas have local README files with current workflow notes:

- `projects/M21C_ls/README.md`
- `projects/era5_land/README.md`
- `projects/GHCN_snwd/README.md`
- `projects/IMS/README.md`
- `projects/ISMN/README.md`
- `projects/SNOTEL/README.md`

## Methods details to carry forward

### P1 vs P2 Terra-only / Terra+Aqua diagnostics

- Dedicated notebook:
  `projects/M21C_ls/notebooks/p1_p2_comparison_diagnostics.ipynb`.
- Detailed written summary:
  `projects/M21C_ls/docs/p1_p2_comparison_summary.md`.
- Output directory:
  `projects/M21C_ls/output/p1_p2_comparison/`.
- This notebook isolates the early MODIS snow-cover transition:
  P1 `2000-06-01` to `2002-06-30` (`MODIS Terra SCF`) and P2
  `2002-07-01` to `2007-05-31` (`MODIS Terra+Aqua SCF`).
- It is a diagnostic companion to the manuscript figure workflow, not a
  replacement for `paper_figures_unified.ipynb`. P1+P2 together correspond to
  the broader `SCF-only period` used in the main ISMN and ERA5-Land period
  figures.
- The notebook compares P1 and P2 on common support where possible:
  common MODIS OFA tiles, common station support for SNOTEL/GHCN, and common
  ERA5-Land grid cells for each variable/metric.
- It writes a compact cross-dataset rollup:
  `projects/M21C_ls/output/p1_p2_comparison/p1_p2_cross_dataset_rollup.csv`,
  plus OFA, IMS, station, and ERA5-Land component tables.
- Current summary: P2 has about twice the MODIS DA observation density per day
  on common tiles, stronger MODIS OmF StdDev reduction (13.93% vs 8.28%),
  stronger IMS snow hit-rate improvement, and generally stronger independent
  validation support from GHCN snow depth and SNOTEL SWE. ERA5-Land is
  supportive but more nuanced: RMSE/bias/correlation often improve more in P2,
  while some ubRMSE metrics are flat or slightly negative.
- Interpretive caveat: the P1/P2 result is not "DA improves every metric
  everywhere." It is more specifically that Aqua MODIS increases useful SCF
  constraint density and strengthens several snow-detection and bulk-error
  diagnostics, with tradeoffs in no-snow classification and some variance-style
  metrics.

### M21C_ls / ISMN soil moisture skill

- Main manuscript-relevant in situ workflow appears to be
  `projects/M21C_ls/notebooks/insitu_skill_cached_batch_figures.ipynb`.
- The notebook uses cached raw time-series files to compare OL and DA against
  in situ soil moisture networks.
- Figure-style outputs include paired DA minus OL skill changes for surface and
  root-zone soil moisture.
- Metrics used in the combined figure are correlation (`R`), anomaly
  correlation (`anomR`), and unbiased RMSE (`ubRMSE`).
- Current batch outputs are written under
  `projects/M21C_ls/outputs_ismn_network_skill/batch_figures/`.

### IMS snow-cover validation

- Current default rerun path is script-first:
  `projects/IMS/scripts/run_ims_ol_da_cell_metrics.py`, followed by
  `projects/IMS/notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`.
- Detailed manuscript-methods notes for the IMS table and Figure 8 maps are in
  `projects/IMS/IMS_methods_readme.md`.
- IMS daily snow-cover data are regridded to the M36 EASE grid using
  nearest-neighbor categorical mapping.
- OL and DA are compared against IMS binary snow/no-snow fields.
- The categorical metrics used in the manuscript-style maps include accuracy,
  hit rate, miss rate, false alarm ratio, and correct rejection rate.
- The validation script writes gridded count/metric NetCDF output and comparison
  tables; the plotting notebook builds the final map/table figures from those
  precomputed outputs.
- Unified Figure 8 uses a positive-improvement convention for every metric:
  `DA - OL` for accuracy, hit rate, and correct rejection rate; `OL - DA` for
  miss rate and false alarm ratio.
- The Terra-only/Terra+Aqua supplemental comparison uses the Discover rerun
  product with custom scopes `P1_MODIS_Terra_SCF` and
  `P2_MODIS_Terra_Aqua_SCF`.

### SNOTEL SWE and snow-depth validation

- Main workflow:
  `projects/SNOTEL/notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb`.
- Detailed manuscript-methods notes for the SNOTEL SWE panel are in
  `projects/SNOTEL/SNOTEL_methods_readme.md`.
- Source observations are daily NRCS AWDB SNOTEL `WTEQ` and `SNWD`.
- The downloader preserves AWDB units; downstream validation converts SNOTEL
  `SNWD` from inches to meters and `WTEQ` from inches to kg m-2.
- Model variables used include `SNODPLAND` for snow depth and `SNOMASLAND` for
  SWE.
- Station-to-model matching uses station latitude/longitude and M36 tile
  coordinates, with a maximum horizontal distance of 40 km.
- The manuscript Figure 9 candidate applies an elevation-difference filter of
  `|dz| < 500 m`.
- Seasonal summaries include all non-JJA periods plus SON, DJF, and MAM.
- Plotted skill quantities include RMSE, unbiased RMSE, and absolute bias, shown
  as DA minus OL changes.

### GHCN snow-depth validation

- Main workflow:
  `projects/GHCN_snwd/notebooks/ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb`.
- Detailed manuscript-methods notes for the GHCN SNWD Figure 10 panel are in
  `projects/GHCN_snwd/GHCN_methods_readme.md`.
- The build notebook creates station metadata, coverage tables, a manifest, and
  yearly parquet partitions from GHCN-Daily SNWD observations.
- Baseline validation compares GHCN snow depth against model snow depth computed
  as `SNODPLAND * SCF`.
- SCF fallback variables are handled in the notebook.
- Station filters include minimum valid days, minimum average snow days per
  year, maximum model-tile distance of 40 km, maximum elevation difference of
  500 m, and a minimum paired-day requirement.
- Metrics include `R`, `anomR`, RMSE, `ubRMSE`, bias, and absolute bias.
- The unified Figure 10 output covers the Northern Hemisphere and reports
  station-level `OL - DA` improvements for lower-is-better metrics, so
  positive/red values indicate lower DA error.

### ERA5 / ERA5-Land comparison

- Current default workflow:
  `projects/era5_land/notebooks/compare_with_reanalysis_strict.ipynb`, followed
  by `projects/era5_land/notebooks/plot_ERA5_comparison.ipynb` and companion
  bar/postage-stamp plotting notebooks.
- Detailed manuscript-methods notes for the ERA5-Land paper figures are in
  `projects/era5_land/ERA5L_methods_readme.md`.
- Monthly GEOS-LDAS OL and DA outputs are compared against ERA5 and ERA5-Land
  reference products on the M36 grid.
- The strict workflow can generate both ERA5 and ERA5-Land outputs, but the
  paper figures currently use ERA5-Land only.
- ERA5/ERA5-Land fields are regridded to M36 using conservative-normalized
  xESMF weights.
- Variables include surface soil moisture, root-zone soil moisture, snow-cover
  fraction, SWE, snow depth, and soil temperature.
- ERA5-Land no-snow NaNs are treated as zeros for SWE and snow-depth statistics
  in the plotting workflow.
- Periodized soil-moisture metrics use the manuscript validation periods:
  V1 `SCF-only period` (P1-P2), V2 `microwave pre-SMAP period` (P3-P5),
  and V3 `SMAP-era microwave period` (P6-P9). Metrics include `R`,
  `anomR`, `ubRMSE`, RMSE, and bias.

### ISMN project area

- The standalone `projects/ISMN` area is mainly a support/inventory workflow for
  local ISMN data access and network coverage checks.
- It documents use of `ismn.interface.ISMN_Interface`, metadata cache creation,
  network/station listing, and network date-range summaries.
- Manuscript soil-moisture skill figures appear to be produced from the
  `projects/M21C_ls` in situ skill notebooks rather than directly from the
  standalone ISMN notebooks.

### Figure 7 / ISMN in-situ soil-moisture trace

- The current figure notebook is
  `projects/M21C_ls/notebooks/insitu_skill_cached_batch_figures.ipynb`.
  It is cache-only: it reads `outputs_ismn_network_skill/*_raw_timeseries.nc`
  and does not reread ISMN archives or model daily files.
- The raw cache files were built by
  `projects/M21C_ls/notebooks/insitu_skill_ismn_network_ol_da.ipynb`, one
  network/window-mode run at a time. That notebook can either rebuild from ISMN
  and model files or, once the cache exists, reload the cache.
- The ISMN source data path documented in the standalone ISMN project is a local
  downloaded/extracted ISMN archive, with metadata generated by
  `ismn.interface.ISMN_Interface` as `python_metadata/ISMN_data.csv`. The local
  laptop path in `projects/ISMN/README.md` is `/Users/amfox/Desktop/ISMN_data`;
  the M21C build notebook was configured for the Discover path
  `/discover/nobackup/projects/land_da/ISMN_data`.
- The standalone `projects/ISMN` notebooks provide the inventory/metadata
  support:
  `interface.ipynb` initializes the ISMN reader and metadata cache,
  `ismn_network_listing.ipynb` lists available networks/stations, and
  `ismn_network_date_ranges_from_metadata.ipynb` writes
  `projects/ISMN/notebooks/ismn_network_timerange_summary.csv`.
- Networks represented in the current hybrid Figure-7-style table are SNOTEL,
  SCAN, USCRN, SMOSMANIA, OZNET, and ARM.
- The hybrid figure uses three windows for SNOTEL, SCAN, and ARM
  (`pre-ASCAT`: 2000-06-01 to 2007-05-31; `pre-SMAP`: 2007-06-01 to
  2015-03-31; `SMAP-era`: 2015-04-01 to 2024-06-01), and two windows for
  USCRN, SMOSMANIA, and OZNET (`pre-SMAP` and `SMAP-era`).
- ISMN observation processing in the M21C build notebook:
  all soil-moisture sensors at a station are inspected, only records with ISMN
  soil-moisture quality flag `G` are retained, duplicate timestamps are averaged,
  surface soil moisture is the available sensor depth nearest 0.05 m, and daily
  means are shifted by 12 hours.
- Root-zone observations use the `matlab_strict` setting in the current caches.
  The cache filename records the active network-specific strict rule, e.g.
  `snotel_n3_c123smv`, `scan_n4_c1234smv`, `uscrn_n4_c1234smv`,
  `smosm_n3_c3smv`, and `oznet_n3_c3smv`.
- Model values are daily GEOS-LDAS `SFMC` and `RZMC` extracted from the OL and
  DA runs on `SMAP_EASEv2_M36_GLOBAL`. Stations are mapped to nearest M36 tiles
  using squared latitude/longitude degree distance with threshold
  `MAX_DISTANCE_DEG2=0.1`; the cache also stores the corresponding haversine
  distance.
- Validation requires at least 1000 observation days per window. Per-station
  statistics are computed for `R`, `anomR`, and `ubRMSE` using
  `projects/matlab2python/scripts/sm_skill_vs_insitu.py`; anomaly correlation
  uses day-of-year climatology with `NMIN_DAY=30`.
- Network summaries enforce common station membership across windows for each
  network/domain/metric. The plotted deltas are paired site-level mean changes:
  `DA-OL` for `R` and `anomR`, and `OL-DA` for `ubRMSE` so positive values mean
  improved lower error.
- Source ISMN output products:
  `outputs_ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE.png`
  and the machine-readable table
  `outputs_ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE_table.csv`.
- Unified manuscript outputs:
  `output/paper_figures/fig07_ismn_skill_by_validation_period.png`,
  `output/paper_figures/fig07_ismn_skill_by_validation_period.pdf`, and
  `output/paper_figures/fig07_ismn_skill_by_validation_period_table.csv`.

## Unified paper-figure notebook plan

Recommended new notebook:

- `projects/M21C_ls/notebooks/paper_figures_unified.ipynb`

The notebook should be a paper-figure assembly layer, not a full raw-data
processing workflow. By default it should read existing analysis products and
regenerate final manuscript figures with shared style, period definitions, sign
conventions, and output naming. Heavy upstream rebuilds should remain in the
existing source notebooks/scripts and be referenced in a figure registry.

### Notebook structure

1. Shared configuration:
   repository root, local/Discover path resolution, output directory, manuscript
   figure list, and a `RUN_FIGURES` switch for rerunning a subset.
2. Canonical observing-system periods:
   one table with period IDs, dates, labels, sensor availability, and display
   colors. Use this table for Fig. 1 and all later shaded periods/labels.
3. Shared plotting conventions:
   global map projection helpers, period shading, panel labels, color maps, and
   an explicit improvement convention. The coauthor-feedback target is red =
   improvement wherever possible.
4. Dataset registry:
   a machine-readable dictionary mapping each figure to required files, upstream
   workflow, and missing-file message.
5. Figure sections:
   one section per manuscript figure, each reading only the files it needs and
   writing to a single paper-figure output directory such as
   `projects/M21C_ls/output/paper_figures/`.
6. Manifest:
   write a CSV or JSON manifest with figure filename, source files, source-file
   mtimes, periods, key settings, and upstream workflow.

### Dataset registry by figure

| Figure | Main inputs for unified notebook | Upstream workflow if products need rebuilding |
| --- | --- | --- |
| Fig. 1 observing-system timeline | Manual/canonical period table in the notebook or a small CSV/YAML. Should include MODIS Terra start, Terra+Aqua transition, ASCAT/SMOS/SMAP/CYGNSS starts and major ASCAT platform transitions. | None, unless timeline dates are derived from OFA counts. |
| Fig. 2 obs/day maps | `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/temporal_stats_DA_20000601_20240531.nc4`, `spatial_stats_DA_200006_202405.pkl`, and `LS_OLv8_M36.ldas_tilecoord.bin`. The NetCDF has `tile=112573`, `species=14`, and variables including `N_data`, `OmF_stdv`, and `OmF_norm_stdv`; the pickle has monthly arrays such as `N_data`, `OmF_stdv`, and `date_vec` with shape `(288, 14)`. | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` or the underlying ObsFcstAna stats generation. |
| Fig. 3 obs/month time series | Same full-period DA temporal/spatial OFA stats as Fig. 2, especially monthly `N_data` grouped by species. | Same as Fig. 2. |
| Fig. 4 full-period OmF maps | Full-period OL and DA temporal stats: `temporal_stats_OL_20000601_20240531.nc4`, `temporal_stats_DA_20000601_20240531.nc4`, plus tilecoord. Uses variables such as `OmF_stdv` and `N_data`. | `LS_ofa_figures_refactor_20260327.ipynb`. |
| Fig. 5 monthly OmF evolution | Full-period OL/DA spatial-stat pickles: `spatial_stats_OL_200006_202405.pkl` and `spatial_stats_DA_200006_202405.pkl`. Uses monthly `OmF_stdv`, `N_data`, and `date_vec`. | `paper_figures_unified.ipynb`; outputs `fig05_monthly_omf_stddev_improvement_evolution.*`, monthly time-series CSV, and P1-P9 period-mean CSV. |
| Fig. 6 OmF maps by period and sensor | P1-P9 period-split OL/DA temporal stats under `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`, named `temporal_stats_{OL,DA}_YYYYMMDD_YYYYMMDD.nc4`, plus tilecoord. | `paper_figures_unified.ipynb`; outputs `fig06_period_sensor_omf_stddev_improvement_maps.*` and `fig06_period_sensor_omf_stddev_improvement_summary.csv`. |
| Fig. 7 ISMN skill | Assembled local table `projects/M21C_ls/output/ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE_table.csv`; raw caches in `projects/M21C_ls/output/ismn_network_skill/*_raw_timeseries.nc` remain the rebuild source. The table has network/domain/metric/window/delta columns and 108 rows. | `paper_figures_unified.ipynb`; outputs `fig07_ismn_skill_by_validation_period.*` and a relabeled plot table. Raw caches from `insitu_skill_ismn_network_ol_da.ipynb`; source table from `insitu_skill_cached_batch_figures.ipynb`. |
| Fig. 8 IMS snow-cover skill | Local precomputed IMS products under `projects/IMS/output/`: full-record count NetCDF/table/metadata `ims_ol_da_*_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.*`, plus Terra/Aqua custom-scope products `ims_ol_da_*_SMAP_EASEv2_M36_GLOBAL_2000_2007_thr0p50_imsSnowDaysGe10_terraAquaScopes.*`. The Terra/Aqua files were generated on Discover by `projects/IMS/scripts/run_ims_terra_aqua_period_scopes_discover.sh`. | `paper_figures_unified.ipynb`; outputs `fig08_ims_snow_cover_skill_maps.*`, `fig08_ims_snow_cover_skill_map_summary.csv`, `fig08_supp_ims_terra_aqua_scope_maps.*`, `fig08_supp_ims_terra_aqua_scope_map_summary.csv`, `fig08_supp_ims_terra_aqua_scope_bars.*`, and `fig08_supp_ims_terra_aqua_scope_bar_values.csv`. Red/positive means DA improved for every metric; miss rate and false alarm ratio are sign-flipped to `OL - DA`. |
| Fig. 9 SNOTEL SWE | `projects/SNOTEL/outputs_snotel_ol_da_validation/snotel_station_metrics_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.csv` plus `tables/snotel_swe_toprow_bar_values_ci_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.csv`. | `paper_figures_unified.ipynb`; outputs `fig09_snotel_swe_skill.*`, `fig09_snotel_swe_bar_values_ci.csv`, `fig09_snotel_swe_station_improvement_table.csv`, and `fig09_snotel_swe_station_improvement_summary.csv`. Map colors use `OL - DA`, so positive/red means DA improved lower-is-better SWE metrics. |
| Fig. 10 GHCN snow depth | `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_metrics_baseline_core_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`. | `paper_figures_unified.ipynb`; outputs `fig10_ghcn_snow_depth_skill.*`, `fig10_ghcn_snow_depth_bar_values_ci.csv`, `fig10_ghcn_snow_depth_station_improvement_table.csv`, and `fig10_ghcn_snow_depth_station_improvement_summary.csv`. Map colors use `OL - DA`, so positive/red means DA improved lower-is-better snow-depth metrics. Raw GHCN build inputs come from `ghcn_snwd_global_1998_present_build.ipynb`; legacy plotting was in `ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb`. |
| Fig. 11 ERA5-Land soil-moisture bars | ERA5-Land periodized metric caches in `projects/era5_land/cache/era5l_periodized_metrics_bars/`, currently `ERA5L_vs_OLv8_M36_strict_summary__e5c443a4ae1b.nc` and `ERA5L_vs_DAv8_M36_strict_summary__27bfb9f0eb98.nc`, derived from strict OL/DA summaries. | `paper_figures_unified.ipynb`; outputs `fig11_era5land_soil_moisture_bars.*` and `fig11_era5land_soil_moisture_bar_stats.csv`. Uses V1/V2/V3 period labels matching Fig. 7; surface/root-zone are distinguished with blue/orange shade pairs. |
| Fig. 12 ERA5-Land soil-moisture maps | Same ERA5-Land periodized metric caches as Fig. 11, plus canonical finite M36 EASE grid coordinates for map plotting. | `paper_figures_unified.ipynb`; outputs `fig12_era5land_soil_moisture_improvement_maps.*` and `fig12_era5land_soil_moisture_map_summary.csv`. Red/positive means DA improved: `DA - OL` for `R`/`anomR`, `OL - DA` for `ubRMSE`. |
| Fig. 13 ERA5-Land snow comparison | Same ERA5-Land periodized metric caches as Fig. 11, using full-period snow metrics over the snow-possible domain. | `paper_figures_unified.ipynb`; outputs `fig13_era5land_snow_comparison_bars.*` and `fig13_era5land_snow_comparison_stats.csv`. Labels use OL/DA consistently; RMSE and ubRMSE improvement are `OL - DA`, while bias annotation reports absolute-bias improvement. |

### Implementation notes

- Keep raw-data rebuilding out of the unified notebook unless a small figure
  truly needs it. For example, ISMN should use the hybrid table first; SNOTEL
  and GHCN should use station-metric CSVs; ERA5-Land should use strict summary
  NetCDFs; IMS should use the precomputed cell-count NetCDF/table products.
- Add a `resolve_path()` helper that checks local laptop paths first and then
  Discover paths. This is essential for IMS, whose manuscript-ready precomputed
  outputs are documented on Discover but are not present in the local checkout.
- Add a `require_files(fig_id, paths)` helper that prints exactly which upstream
  notebook/script must be run when a product is missing.
- Prefer regenerating figures from tables/NetCDFs over embedding existing PNGs.
  Existing PNGs remain useful for visual regression checks.

### Period decisions

- Use `2024-05-31` as the common manuscript end date for figure labeling.
- The OFA monthly count file
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/spatial_stats_DA_200006_202405.pkl`
  stores `date_vec` as monthly `YYYYMM` values from `200006` through `202405`.
- In that monthly count product, MODIS species `12` is present from
  `2000-06-01` through `2024-05-01`, while MODIS species `11` first appears in
  `2002-07-01` and continues through `2024-05-01`.
- The companion obs-species notebook labels the MODIS products as `MYD10C1`
  and `MOD10C1`; the timing therefore supports using `2002-07-01` as the
  Terra-only to Terra+Aqua transition for assimilation/count figures.
- Use two linked period definitions:
  a fine-grain observing-system timeline for Fig. 1/Fig. 3/Fig. 5/Fig. 6, and
  a three-period validation summary for Fig. 7/Fig. 11/Fig. 12.

Fine-grain observing-system periods:

| ID | Dates | Label |
| --- | --- | --- |
| P1 | 2000-06-01 to 2002-06-30 | MODIS Terra SCF |
| P2 | 2002-07-01 to 2007-05-31 | MODIS Terra+Aqua SCF |
| P3 | 2007-06-01 to 2010-04-30 | SCF + ASCAT-A |
| P4 | 2010-05-01 to 2013-03-31 | SCF + ASCAT-A + SMOS |
| P5 | 2013-04-01 to 2015-03-31 | SCF + ASCAT-A/B + SMOS |
| P6 | 2015-04-01 to 2018-07-31 | SCF + ASCAT-A/B + SMOS + SMAP |
| P7 | 2018-08-01 to 2019-10-31 | SCF + ASCAT-A/B + SMOS + SMAP + CYGNSS |
| P8 | 2019-11-01 to 2021-11-30 | SCF + ASCAT-A/B/C + SMOS + SMAP + CYGNSS |
| P9 | 2021-12-01 to 2024-05-31 | SCF + ASCAT-B/C + SMOS + SMAP + CYGNSS |

P1/P2 OFA temporal-stat regeneration:

- The unified notebook expects Fig. 6 period files named
  `temporal_stats_{OL,DA}_YYYYMMDD_YYYYMMDD.nc4` in
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`.
- Local laptop inputs now include the P1/P2 split temporal-stat files generated
  on Discover from the monthly ObsFcstAna sum files in the M21C run output
  directories using
  `projects/M21C_ls/scripts/build_p1p2_ofa_temporal_stats_discover.sh`.
- Default Discover input roots in that driver are
  `/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg`
  and
  `/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v3/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg`;
  override `OL_MONTHLY_ROOT` or `DA_MONTHLY_ROOT` if the monthly stats are
  stored elsewhere.
- The M21C run-local monthly files used for this build are the non-deduplicated
  sum products, e.g. `OL.ens_avg.ldas_ObsFcstAna_sums.200006.nc4` and
  `DA.ens_avg.ldas_ObsFcstAna_sums.200006.nc4`. The P1/P2 driver passes this
  filename pattern explicitly and does not use the `_dedup` products.
- The driver writes to
  `/discover/nobackup/projects/land_da/geosldas-analysis/projects/M21C_ls/output/ofa_temporal_stats/`
  by default. Copy the four regenerated files back to the local diagnostic
  directory before rerunning the unified notebook availability checks.

Three-period validation summaries:

| ID | Dates | Label | Fine-period mapping |
| --- | --- | --- | --- |
| V1 | 2000-06-01 to 2007-05-31 | SCF-only period | P1 + P2 |
| V2 | 2007-06-01 to 2015-03-31 | microwave pre-SMAP period | P3 + P4 + P5 |
| V3 | 2015-04-01 to 2024-05-31 | SMAP-era microwave period | P6 + P7 + P8 + P9 |

## Items to verify before final manuscript text

- Confirm whether Figure 7 uses the `custom`, `hybrid`, or `two_period`
  all-network output.
- Confirm whether the final text should call ERA5-Land a reanalysis reference
  or a land-surface analysis reference in Figures 11-13 captions.
