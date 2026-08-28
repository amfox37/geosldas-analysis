# CYGNSS_L1_AZ

Workspace for CYGNSS Level-1 AZ analyses in the GEOSldas analysis repository.

## Goal

Use this project area for scripts, notebooks, jobs, and generated outputs
specific to the CYGNSS L1 AZ workflow. The initial focus is the
`OLv8_M36_all_sensors_AZ` GEOSldas monitoring run, which exercises the new
CYGNSS L1 DDM scalar and H SAF ASCAT readers alongside existing SMOS/SMAP
brightness-temperature readers.

## Layout

- `scripts/`: reusable Python or shell helpers.
- `notebooks/`: exploratory analysis and figure notebooks.
- `jobs/`: batch scripts for Discover or other remote runs.
- `runs/`: durable run notes and status records.
- `output/`: generated local products and figures. This directory is ignored by
  Git through the repository-wide output rules.

## OFA Summary Figures

The local ObsFcstAna summary cache in `data/omf_compare_sums/OL_AZ_all_sensors`
can be plotted with:

```bash
python projects/CYGNSS_L1_AZ/scripts/plot_az_omf_summary.py
```

The script follows the CYGNSS/M21C figure convention: species are summarized
into observing-system groups with `N_data`-weighted observation and OmF
statistics, and low-count diagnostic values are excluded from those weighted
means (`--nmin 20` by default). Observation counts are plotted as raw monitored
counts because this run is open-loop-style monitoring, not assimilation. O-F
mean maps use a zero-centered diverging color scale because the sign is
physically meaningful for the unscaled observations in this OL run.

Primary outputs are written under
`projects/CYGNSS_L1_AZ/output/omf_summary_figures/`:

- `az_fig01_mean_monitored_observations_per_day_by_system.png` (includes a
  SMAP horizontal-only count panel)
- `az_fig02_monthly_monitored_observation_counts.png`
- `az_fig03_full_period_omf_stddev_by_system.png`
- `az_fig04_monthly_omf_stddev_by_system.png`
- `az_fig05_full_period_omf_mean_by_system.png`
- `az_fig06_monthly_omf_mean_by_system.png`
- `az_fig07_full_period_o_mean_by_system.png`
- `az_fig08_monthly_o_mean_by_system.png` (observed and forecast means)
- `az_fig09_full_period_f_mean_by_system.png`
- `az_fig10_monthly_cygnss_l3_l1_forecast_mean.png`
- `az_fig11_tile_elevation.png`

The script also writes per-species QC maps/time series and
`az_monthly_species_summary.csv`. Map-panel annotations and
`az_map_support_summary.csv` include both native valid-tile statistics and
statistics recomputed over the CYGNSS L1 valid-tile mask
(`CYGNSS_L1_DDM3X5_CROP_SCALAR N_data >= --nmin`) so observing systems can be
compared on common spatial support. Fig07 observed-mean and Fig09 forecast-mean
maps use matched colorbar ranges for corresponding panels with matching units;
ASCAT is left unmatched because its observations are percent and its forecast is
soil moisture. CYGNSS L1 and CYGNSS L3 are kept separate in the grouped OmF
figures because their native units differ.

## CYGNSS L1 R-Sweep Figures

A separate four-experiment set assimilating CYGNSS L1 with full, half and
quarter observation-error variance against a scaled open-loop baseline, over
calendar 2020. Stats products live in `output/stats_output/`. Plot with:

```bash
python projects/CYGNSS_L1_AZ/scripts/plot_cygl1_R_sweep.py
```

Outputs to `projects/CYGNSS_L1_AZ/output/R_sweep_figures/`:

- `R_sweep_fig01_omf_stdv_maps.png` (O-F stdv, open loop and each DA run minus OL)
- `R_sweep_fig02_omf_stdv_monthly.png`
- `R_sweep_fig03_observation_counts.png`
- `R_sweep_fig04_mean_obs_and_forecast.png`
- `R_sweep_fig05_oma_vs_omf.png` (analysis-impact check)

The same grouping and weighting conventions as `plot_az_omf_summary.py` apply.
Maps use `pcolormesh` on the native EASEv2 M36 grid rather than scatter, so
panels show true tile footprints and the unsampled mask.

Headline result and its caveats are in `runs/cygl1_assim_R_sweep.md`:
assimilating CYGNSS L1 degrades the fit to every independent observing system,
increasingly so as R is reduced, while the analysis itself behaves correctly.
The observation is well scaled but correlates with the forecast at only r ~ 0.35
against 0.64-0.93 for the other sensors, so its increments are correctly-sized
noise. The runs are otherwise undocumented — that note lists the provenance gaps.

## Runs

- `runs/OLv8_M36_all_sensors_AZ.md`: current long-run validation/monitoring
  experiment notes.
- `runs/cygl1_assim_R_sweep.md`: CYGNSS L1 assimilation observation-error
  sweep, 2020, with findings and open provenance questions.

## Notes

Keep raw data and heavy generated products out of Git. Document durable external
paths, processing commands, and product versions here as the workflow takes
shape.
