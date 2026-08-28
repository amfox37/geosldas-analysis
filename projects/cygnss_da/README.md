# cygnss_da

Workspace for CYGNSS land data assimilation diagnostics. Contains notebooks for CYGNSS-only OFA figures, in situ comparisons, and proposal graphics, plus scripts for staging the raw CYGNSS observation files.

## Key notebooks
- `notebooks/CYG_insitu_plotter_100325.ipynb` – plots CYGNSS skill against in situ references.
- `notebooks/cygnss_cd_ofa_figures_081825.ipynb` – generates current control/diagnostic figure sets for CYGNSS experiments.
- `notebooks/cygnss_cd_TC_figures_101025.ipynb` – tropical cyclone-focused diagnostics.

## Supporting scripts
- `scripts/download_cygnss_l3.sh` – `podaac-data-downloader` wrapper for the
  CYGNSS L3 Soil Moisture V3.2 product (PO.DAAC collection
  `CYGNSS_L3_SOIL_MOISTURE_V3.2`, https://doi.org/10.5067/CYGNU-L3S32).
  Defaults to the 36 km grid (the collection also serves a 9 km grid) and
  writes into `Y<yyyy>/M<mm>`, one month at a time. Requires a
  `~/.netrc` entry for `urs.earthdata.nasa.gov` (`chmod 600`) and
  `pip install podaac-data-subscriber`. For example, test one month with
  `scripts/download_cygnss_l3.sh -d /data/CYGNSS_L3 -y 2018 -m 08 --dry-run`;
  omit `--dry-run` to download it, or omit `--month` to fetch the whole year.
- `scripts/mv_cygnss_obs_files.py` – helper to organize CYGNSS observations into experiment directories.
- `scripts/pymove_files_v3.py` – generic mover used when collecting inputs on Discover.
