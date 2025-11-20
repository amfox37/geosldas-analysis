# discover_JH

General diagnostics sandbox for the Discover/JH environment. Contains ASCAT tools, snow QC notebooks, ObsFcstAna readers, and helper scripts that feed into multiple projects. Use this directory for ad hoc analyses before migrating stable workflows into project-specific folders.

## Representative notebooks
- `analyze_snowDA_041025.ipynb` – archived snow DA analysis workflow (kept for reference).
- `process_snow_incr_monthly_061125.ipynb` – summarizes monthly snow increment files.
- `UpT13_reader_plotter.ipynb` – quick visualization of UpT13 observation subsets.

## Key helpers
- `my_functions.py`, `helper.py`, `mapper_functions.py` – shared IO/plotting routines used by many notebooks.
- `ObsFcstAna_reader.py` – standalone ObsFcstAna reader useful for scripting without notebooks.
