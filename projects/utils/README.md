# utils

Shared notebooks, scripts, and Fortran helpers used across projects. This is the main hub for ObsFcstAna post-processing, OmF/OMA diagnostics, and misc utility plotting.

## Key notebooks
- `notebooks/obsfcstana_ts_stats_111324.ipynb` – current ObsFcstAna tile-stat workflow (multi-sensor analysis, cross-masking).
- `notebooks/ofa_figures_071725.ipynb` – produces the latest OFA figure suite for the ObsFcstAna analysis pipeline.
- `notebooks/process_experiment_111124.ipynb` – template for aggregating experiment outputs into monthly stats/graphics.

## Supporting utilities
- `scripts/mapper_functions.py` – shared plotting helpers imported by many notebooks.
- `scripts/read_bufr_files.py` – Python BUFR reader used both here and in project-specific notebooks.
- `fortran/read_obs_scaling.F90` – reference Fortran routines for scaling parameters.
