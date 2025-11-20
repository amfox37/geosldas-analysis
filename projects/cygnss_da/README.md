# cygnss_da

Workspace for CYGNSS land data assimilation diagnostics. Contains notebooks for CYGNSS-only OFA figures, in situ comparisons, and proposal graphics, plus scripts for staging the raw CYGNSS observation files.

## Key notebooks
- `notebooks/CYG_insitu_plotter_100325.ipynb` – plots CYGNSS skill against in situ references.
- `notebooks/cygnss_cd_ofa_figures_081825.ipynb` – generates current control/diagnostic figure sets for CYGNSS experiments.
- `notebooks/cygnss_cd_TC_figures_101025.ipynb` – tropical cyclone-focused diagnostics.

## Supporting scripts
- `scripts/mv_cygnss_obs_files.py` – helper to organize CYGNSS observations into experiment directories.
- `scripts/pymove_files_v3.py` – generic mover used when collecting inputs on Discover.
