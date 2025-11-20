# era5_land

Notebooks/scripts for comparing GEOS-LDAS output with ERA5-Land. Used to download ERA5-Land data, compute aridity indices, and build comparison figures.

## Key notebooks
- `notebooks/anomaly_calculator_081525.ipynb` – computes anomaly statistics from ERA5-Land and LDAS products.
- `notebooks/compare_with_ERA5_Land.ipynb` and `..._masked.ipynb` – side-by-side comparisons of GEOS-LDAS vs ERA5-Land (global and masked subsets).
- `notebooks/download_ERA5_land.ipynb` – outlines the ERA5-Land data download workflow.

## Supporting scripts
- `scripts/compute_aridity_indices.py` – CLI helper for aridity index computation.
- `scripts/era5l_monthly_nc/unzip_and_merge_era5l.py` – merges monthly ERA5-Land files.
