# M21C_ls

Land Sweeper (LS) diagnostics for the M21C campaign. Notebooks track LS assimilation performance, ObsFcstAna statistics, and validation at core sites; scripts support time-series generation.

## Notebook Guide

### In-situ skill workflows

- `notebooks/insitu_skill_ismn_network_ol_da.ipynb`
  - ISMN-network OL/DA skill workflow with window switches and raw-timeseries cache support.
  - Includes strict root-zone logic options for MATLAB-aligned comparisons.

- `notebooks/insitu_skill_preSMAP_SMAPera.ipynb`
  - Pre-SMAP vs SMAP-era in-situ OL/DA skill from raw MAT timeseries.

- `notebooks/insitu_skill_preSMAP_SMAPera_mat_vs_ismn.ipynb`
  - Side-by-side MAT-obs vs ISMN-obs skill comparison for pre-SMAP and SMAP-era periods.

- `notebooks/insitu_skill_from_raw_timeseries.ipynb`
  - Full-period in-situ skill analysis from raw MAT timeseries.

- `notebooks/insitu_skill_from_raw_timeseries_date_ranges.ipynb`
  - Date-range segmented version of raw-timeseries in-situ skill analysis.

### Snow workflows

- `notebooks/snow_daily_seasonal_ol_da_from_discover.ipynb`
  - Daily/seasonal snow diagnostics (OL vs DA) from Discover model outputs.

- `notebooks/snow_da_impact_OLv8_vs_DAv8_M21C_land_sweeper.ipynb`
  - Snow DA impact analysis across OLv8 and DAv8 experiments.

### Legacy / figure notebooks

- `notebooks/LS_insitu_plotter_102625.ipynb` - broad in-situ plotting/inspection notebook.
- `notebooks/LS_ofa_figures.ipynb` - OL vs DA figure package generation.
- `notebooks/obsfcstana_ts_stats_LS_021725.ipynb` - ObsFcstAna time-series statistics.
- `notebooks/LS_compare_OL_DA.ipynb` - compact OL vs DA comparison notebook.

## Scripts

- `scripts/LS_full_TS_insitu_plotter.py`
  - Python driver to regenerate large in-situ plot sets outside notebook execution.
