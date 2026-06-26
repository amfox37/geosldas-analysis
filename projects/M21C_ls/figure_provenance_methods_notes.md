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
| Figure 4 | Spatial maps of (DA−OL)/OL O-F StdDev [%] by sensor (MODIS, ASCAT, SMOS, SMAP, CYGNSS). | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` | Confirmed: notebook contains "Normalized Percent Difference (OL - DA)" spatial figure logic for all five sensors. |
| Figure 5 | Monthly time series of normalized percent difference in O-F StdDev (OL−DA) for ASCAT, SMOS, SMAP, and CYGNSS, 2000–2024. | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` | Confirmed: title string "Normalized Percent Difference (OL - DA): LS_DAv8_M36_200006_202405" with per-sensor lines. |
| Figure 6 | (DA−OL)/OL O-F StdDev (%) by Period and Sensor — 7-row × 4-column map grid (periods I–VII, ASCAT/SMOS/SMAP/CYGNSS). | `projects/M21C_ls/notebooks/LS_ofa_figures_refactor_20260327.ipynb` | Confirmed: only notebook containing the title string "by Period and Sensor". Not in the older LS_ofa_figures.ipynb. |
| Figure 7 | In situ soil moisture skill differences between DA and OL for surface and root-zone soil moisture, grouped by observing-system period. | `projects/M21C_ls/notebooks/insitu_skill_cached_batch_figures.ipynb` | Likely `projects/M21C_ls/outputs_ismn_network_skill/batch_figures/all_networks_custom_OL_DA_delta_surface_rz_R_anomR_ubRMSE.png` or the corresponding `all_networks_hybrid_...` version. |
| Figure 8 | IMS categorical snow-cover skill differences between DA and OL. | `projects/IMS/scripts/run_ims_ol_da_cell_metrics.py` plus `projects/IMS/notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb` | Plotting notebook saves an `ims_all_period_delta_metrics_*_nh_robinson_2x3.png` style figure from precomputed IMS metric outputs. |
| Figure 9 | Seasonal SWE performance at SNOTEL stations: top-row domain bars and bottom-row station maps for DA minus OL. | `projects/SNOTEL/notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb` | Strong match to `projects/SNOTEL/outputs_snotel_ol_da_validation/figures/snotel_swe_2x3_bars_allsites_maps_da_minus_ol_elevfilt500_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.png`. |
| Figure 10 | GHCN snow depth validation: seasonal/domain bars and DA minus OL maps. | `projects/GHCN_snwd/notebooks/ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb` | Strong match to `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/figures/ghcn_baseline_core_snodpland_2x3_bars_maps_nh_da_minus_ol_ALL_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.png`. |
| Figures 11-14 | ERA5/ERA5-Land comparison figures for soil moisture and snow variables. | `projects/era5_land/notebooks/compare_with_reanalysis_strict.ipynb` plus plotting notebooks in `projects/era5_land/notebooks/` | Candidate outputs include `figures_era5l_postage_stamps/postage_stamp_da_minus_ol_surface_rz_sm_6x3.png` and bar figures under `figures_era5l_bars/`. |

## Project README inventory

These project areas have local README files with current workflow notes:

- `projects/M21C_ls/README.md`
- `projects/era5_land/README.md`
- `projects/GHCN_snwd/README.md`
- `projects/IMS/README.md`
- `projects/ISMN/README.md`
- `projects/SNOTEL/README.md`

## Methods details to carry forward

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
- The likely Figure 10 output covers the Northern Hemisphere and reports DA
  minus OL metric differences for the full analysis period and seasons.

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
- Periodized metrics are computed for pre-SMAP, SMAP-era, and full-period
  windows. Metrics include `R`, `anomR`, `ubRMSE`, RMSE, and bias.

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
- Current output products:
  `outputs_ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE.png`
  and the machine-readable table
  `outputs_ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE_table.csv`.

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
| Fig. 5 monthly OmF evolution | Full-period OL/DA spatial-stat pickles: `spatial_stats_OL_200006_202405.pkl` and `spatial_stats_DA_200006_202405.pkl`. Uses monthly `OmF_stdv`, `N_data`, and `date_vec`. | `LS_ofa_figures_refactor_20260327.ipynb`. |
| Fig. 6 OmF maps by period and sensor | Period-split OL/DA temporal stats under `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`, e.g. `temporal_stats_{OL,DA}_20070601_20100430.nc4`, `20100501_20130331.nc4`, `20130401_20150331.nc4`, `20150401_20180731.nc4`, `20180801_20191031.nc4`, `20191101_20211130.nc4`, and `20211201_20240531.nc4`, plus tilecoord. | `LS_ofa_figures_refactor_20260327.ipynb`. |
| Fig. 7 ISMN skill | Prefer the already assembled local table `projects/M21C_ls/output/ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE_table.csv` and, if needed, the raw caches in `projects/M21C_ls/output/ismn_network_skill/*_raw_timeseries.nc`. The table has network/domain/metric/window/delta columns and 108 rows. | `projects/M21C_ls/notebooks/insitu_skill_cached_batch_figures.ipynb`; raw caches from `insitu_skill_ismn_network_ol_da.ipynb`. |
| Fig. 8 IMS snow-cover skill | Precomputed IMS products on Discover: `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.nc4` and `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_comparison_table_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.csv`. The existing Fig. 8-style PNG is `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/figures_ims_maps_and_tables/ims_all_period_delta_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10_nh_robinson_2x3.png`. Local checkout only has partial 2024 IMS raw/conversion files. | `projects/IMS/scripts/run_ims_ol_da_cell_metrics.py`, then `projects/IMS/notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`. |
| Fig. 9 SNOTEL SWE | `projects/SNOTEL/outputs_snotel_ol_da_validation/snotel_station_metrics_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.csv` or parquet, plus `tables/snotel_swe_toprow_bar_values_ci_...csv` if reusing existing bootstrapped bar values. | `projects/SNOTEL/notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb`. |
| Fig. 10 GHCN snow depth | `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_metrics_baseline_core_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`, plus station selection/tile-map CSVs if rebuilding map annotations. | `projects/GHCN_snwd/notebooks/ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb`; build inputs from `ghcn_snwd_global_1998_present_build.ipynb`. |
| Fig. 11 ERA5-Land soil-moisture bars | `projects/era5_land/notebooks/ERA5L_vs_OLv8_M36_strict_summary.nc` and `ERA5L_vs_DAv8_M36_strict_summary.nc`. Each strict summary has `time=288`, `y=379`, `x=964`, and variables `SM_model`, `SM_era`, `RZ_model`, `RZ_era`, `SCF_model`, `SCF_era`, `SWE_model`, `SWE_era`, `SNWD_model`, `SNWD_era`, and `mask_both`. | `projects/era5_land/notebooks/compare_with_reanalysis_strict.ipynb`, then `plot_ERA5L_comparison_bars.ipynb`. |
| Fig. 12 ERA5-Land soil-moisture maps | Same strict ERA5-Land summaries as Fig. 11, plus `LS_OLv8_M36.ldas_tilecoord.bin` or another resolved M36 tilecoord file for weighted spatial means/map masking. | `projects/era5_land/notebooks/plot_ERA5L_postage_stamp_maps.ipynb`. |
| Fig. 13 ERA5-Land snow comparison | Same strict ERA5-Land summaries as Fig. 11, using `SCF_*`, `SWE_*`, and `SNWD_*`. | `projects/era5_land/notebooks/plot_ERA5L_comparison_bars.ipynb`. |

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
- Local laptop inputs currently include the P3-P9 temporal-stat files, but not
  the P1/P2 split temporal-stat files. Regenerate P1/P2 on Discover from the
  monthly ObsFcstAna sum files in the M21C run output directories using
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
| V1 | 2000-06-01 to 2007-05-31 | SCF-only era | P1 + P2 |
| V2 | 2007-06-01 to 2015-03-31 | pre-SMAP microwave era | P3 + P4 + P5 |
| V3 | 2015-04-01 to 2024-05-31 | SMAP-era multi-sensor era | P6 + P7 + P8 + P9 |

## Items to verify before final manuscript text

- Confirm whether Figure 7 uses the `custom`, `hybrid`, or `two_period`
  all-network output.
- Confirm the exact ERA5/ERA5-Land files used for Figures 11-14 by matching the
  draft panels against the `figures_era5l_bars/` and
  `figures_era5l_postage_stamps/` outputs.
