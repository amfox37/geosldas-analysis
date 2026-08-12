# M21C_ls

Land Sweeper (LS) diagnostics for the M21C campaign. Notebooks track LS assimilation performance, ObsFcstAna statistics, and validation at core sites; scripts support time-series generation.

## Shared Observing-System Registry

`config/observing_system_registry.json` is the machine-readable source of truth
for the P1-P9 observing-system periods, V1-V3 validation periods, sensor
availability, and segment-length constraints. Both the unified paper notebook
and new trend/breakpoint work load it through `scripts/m21c_periods.py`.

See `docs/observing_system_period_registry.md` for the period rationale. P7 is
only 15 months long and is explicitly excluded from period-specific slope
claims and changepoint detection-agreement scoring under a 24-month minimum
segment.

## Notebook Guide

### Current notebooks

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

### Figure / utility notebooks

- `notebooks/LS_insitu_plotter_102625.ipynb` - broad in-situ plotting/inspection notebook.
- `notebooks/LS_ofa_figures.ipynb` - OL vs DA figure package generation.
- `notebooks/LS_compare_OL_DA.ipynb` - compact OL vs DA comparison notebook.
- `notebooks/tsoil1_climatology.ipynb` - tsoil climatology diagnostics.
- `notebooks/paper_figures_unified.ipynb` - unified manuscript figures using the shared period registry.

### Trends and breakpoints

- `docs/trends_breakpoints_methods.md` - input contract and statistical guardrails for the replacement trend workflow.
- `config/trend_breakpoint_inputs.json` - monthly file, source, dimension, date, and tile-area contract.
- `config/trend_breakpoint_variable_selection.json` - phase-1 and phase-2 selections from existing monthly-synthesis products.
- `config/trend_statistics.json` - support, seasonal-adjustment, modified-MK, and FDR settings.
- `scripts/audit_trend_breakpoint_inputs.py` - audit gate for current OL v2/DA v3 monthly inputs.
- `scripts/trend_breakpoint_series.py` - paired OL/DA tile loader, synthesis masks, and area-weighted monthly series.
- `scripts/trend_statistics.py` - exact Theil-Sen, conservative modified Mann-Kendall, and BH-FDR engine.
- `scripts/build_trend_statistics.py` - CLI for provenance-complete tile-level trend NetCDF files.
- `scripts/validate_trend_statistics.py` - fixed-seed white-noise and AR(1) false-positive/power validation.
- `tests/test_trend_breakpoint_series.py` - focused tests for pairing, units, masks, weights, signs, dates, and missing data.
- `tests/test_trend_statistics.py` - synthetic trend, autocorrelation, support, FDR, and parallel-consistency tests.
- `legacy/trends_2025/` - provenance archive of the superseded 2025 trend exploration.

### Archived notebooks

Older exploratory notebooks are kept in:

- `notebooks/legacy/Core_site_stats_plotter_042124.ipynb`
- `notebooks/legacy/Plot_validation_core_sites_010324.ipynb`
- `notebooks/legacy/obsfcstana_ts_stats_LS_021725.ipynb`
- `notebooks/legacy/pentad_scatter_compare.ipynb`
- `notebooks/legacy/process_experiment_LS_021725.ipynb`

## Scripts

- `scripts/LS_full_TS_insitu_plotter.py`
  - Python driver to regenerate large in-situ plot sets outside notebook execution.
- `scripts/m21c_periods.py`
  - Loads and validates the shared observing-system registry.
- `scripts/audit_trend_breakpoint_inputs.py`
  - Validates the existing monthly-synthesis coverage, inventory, variables,
    source experiments, and period constraints before trend calculations.
- `scripts/trend_breakpoint_series.py`
  - Loads audited monthly products with strict OL/DA pairing, M36 tile-area
    weighting, explicit monthly-rate conversion, and shared snow/warm masks.

Run the trend input gate and focused loader tests with:

```bash
python projects/M21C_ls/scripts/audit_trend_breakpoint_inputs.py --no-write
python -m pytest projects/M21C_ls/tests/test_trend_breakpoint_series.py projects/M21C_ls/tests/test_trend_statistics.py -q
python projects/M21C_ls/scripts/validate_trend_statistics.py --n-series 100 --n-jobs 2
```

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
