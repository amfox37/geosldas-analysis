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

## Runs

- `runs/OLv8_M36_all_sensors_AZ.md`: current long-run validation/monitoring
  experiment notes.

## Notes

Keep raw data and heavy generated products out of Git. Document durable external
paths, processing commands, and product versions here as the workflow takes
shape.
