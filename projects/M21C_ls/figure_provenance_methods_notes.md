# M21C LS manuscript figure provenance and methods notes

Working note for tracing figures in:

- `/Users/amfox/Documents/Publications/SM and SCF paper/draft_033026.docx`

The Word document embeds resized PNGs, so exact file hashes do not match the
local figure files. The mappings below are based on captions, aspect ratios,
notebook save names, and project READMEs.

## Figure provenance

| Draft figure | Draft content | Likely source workflow | Likely output or notes |
| --- | --- | --- | --- |
| Figure 2 | Mean number of assimilated observations per day for MODIS SCF, ASCAT SM, SMOS Tb, SMAP Tb, and CYGNSS SM. | `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb` | CYGNSS/combined assimilation figure notebook; output figures are written to `paper_figures/`. |
| Figure 3 | Temporal evolution of number of assimilated observations, total and by observing system. | `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb` | Same combined assimilation figure workflow. |
| Figure 4 | Relative DA vs OL difference in time-series standard deviation of OmF residuals, by sensor. | `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb` | Same combined assimilation figure workflow. |
| Figure 5 | Monthly time series of normalized percent difference in OmF residual standard deviation for ASCAT, SMOS, SMAP, and CYGNSS. | `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb` | Same combined assimilation figure workflow. |
| Figure 6 | Spatial distribution of relative difference in OmF residual standard deviation, stratified by observing-system period and sensor. | `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb` or nearby CYGNSS/ObsFcstAna workflow | Mapping is less certain than Figures 7-10; confirm visually before final manuscript citation. |
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

## Items to verify before final manuscript text

- Visually confirm Figures 2-6 against the CYGNSS/combined assimilation figure
  outputs because the draft figure numbering appears to have shifted relative to
  some saved figure names.
- Confirm whether Figure 7 uses the `custom`, `hybrid`, or `two_period`
  all-network output.
- Confirm the exact ERA5/ERA5-Land files used for Figures 11-14 by matching the
  draft panels against the `figures_era5l_bars/` and
  `figures_era5l_postage_stamps/` outputs.
