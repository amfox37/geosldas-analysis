# ascat_da

Diagnostics and utilities for ASCAT-based land data assimilation experiments. The notebooks here generate diagnostic figures and statistics for different ASCAT filtering/footprint configurations, while the scripts automate data preparation.

## Key notebooks
- `notebooks/compare_comb_fp_043024.ipynb` – latest comparison of combined FP experiments with OL, including map plots and PDF summaries.
- `notebooks/ASCAT_masking_021324.ipynb` – builds and inspects ASCAT observation masks for assimilation tests.
- `notebooks/regrid_ASCAT_mask_022724.ipynb` – regrids ASCAT masks to the analysis grid for quick QA.

## Supporting scripts
- `scripts/check_ascat_duplicates.py` – helper to spot duplicate ASCAT files before ingest.
- `scripts/filename_lister_v2.py` – utility for listing/renaming raw observation files.
