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

## Notebook Runtime/Smoke-Test Runbook

Use this when developing new M21C notebooks from Codex/shell or when notebook runs fail with OpenMP shared-memory errors.

### 1) Prevent OpenMP SHM failures in shell runs

Symptom:
- `OMP: Error #179: Function Can't open SHM2 failed`

Working fix:
- Set `MKL_THREADING_LAYER=SEQUENTIAL` before Python starts.
- Also cap thread counts to reduce backend/thread contention.

Example:

```bash
MKL_THREADING_LAYER=SEQUENTIAL \
OMP_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 \
~/mamba/envs/regrid/bin/python your_script.py
```

### 2) Matplotlib inline vs pop-out windows on macOS

Issue encountered:
- `ismn.interface` can force `TkAgg` on Darwin at import time.

Recommended notebook pattern:
- Avoid top-level `from ismn.interface import ISMN_Interface`.
- Lazy-import ISMN only in code paths that require it.
- Re-enforce inline backend after lazy import.

### 3) Smoke-test command pattern (cache-only notebooks)

For quick one-file smoke tests:
- Run with the environment variables above.
- Use `~/mamba/envs/regrid/bin/python`.
- Ensure xarray can read NetCDF via `netCDF4` (present in `regrid` env).

### 4) Current batch cache notebook behavior

`notebooks/insitu_skill_cached_batch_figures.ipynb` now includes:
- Startup env guard for the OpenMP issue.
- Cache-only processing over `*_raw_timeseries.nc`.
- OL/DA bar figures and paired delta figures (Figure-4 style).
- Delta figure key/legend text with paired `n` counts (surface/rootzone) by window.
