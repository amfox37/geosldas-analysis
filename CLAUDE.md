# CLAUDE.md — geosldas-analysis

Working notes and context for AI-assisted development in this repo.

---

## Python environments

Default environment for almost all work:

```
/Users/amfox/mamba/envs/regrid/bin/python
```

Has: `eccodes`, `cartopy`, `xarray`, `netCDF4`, `numpy`, `matplotlib`, `scipy`, `pyproj`.

Other environments (use only if specifically needed):
- `xr` — Anaconda, xarray/netCDF4 only, no eccodes
- `cygnss` — CYGNSS-specific dependencies

**OpenMP / MKL errors on macOS** (e.g. `OMP: Error #15`):
```bash
MKL_THREADING_LAYER=SEQUENTIAL OMP_NUM_THREADS=1 python script.py
```

---

## Repo structure

```
common/
  python/{io,plotting,stats}/   # shared Python helpers
  matlab/                       # shared MATLAB utilities
  fortran/                      # shared Fortran utilities
data/                           # docs + .gitignore only; large data via symlinks
  external/test_data -> /Users/amfox/Desktop/GEOSldas_diagnostics/test_data
projects/
  ascat_da/         # H121 CDR vs Legacy BUFR ASCAT integration & validation
  cygnss_da/        # CYGNSS soil moisture DA (M21C manuscript figures)
  era5_land/        # ERA5 / ERA5-Land reanalysis comparison
  snow_da/          # Snow data assimilation diagnostics
  SNOTEL/           # SNOTEL SWE/snow-depth in-situ validation
  GHCN_snwd/        # GHCN snow-depth in-situ validation
  IMS/              # IMS snow-cover categorical skill (DA vs OL)
  ISMN/             # ISMN soil-moisture in-situ skill
  SMOS_IC/          # SMOS-IC comparison
  M21C_ls/          # Combined manuscript (SM + SCF): figure provenance, ISMN batch
  utils/            # General utilities
  matlab2python/    # Migrating legacy MATLAB postprocessing to Python
  discover_JH/      # Sandbox for NASA Discover / JupyterHub work
```

Each project has its own `README.md` with current workflow notes.

---

## GEOSldas naming conventions

| Token | Meaning |
|---|---|
| `OL` | Open Loop (no data assimilation) |
| `DA` | Data Assimilation |
| `LS` | Land Sweeper (coupled snow+SM DA experiment type) |
| `v8` | GEOSldas version 8 |
| `M36` | SMAP EASE2 36 km global grid |
| `M33` | 33 km variant |
| `M21C` | Campaign/manuscript name (≠ grid size) |
| `EASEv2` | EASE-Grid 2.0 projection |

Typical experiment directory names: `LS_OLv8_M36`, `LS_DAv8_M36`, `DAv8_M36`, `hsaf_cdr_test_DAv8_M36`.

---

## ObsFcstAna (OFA) files

Located under: `<expt>/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y<YYYY>/M<MM>/`

One `.nc4` file per 3-hour assimilation window (8 per day). Variables:
- `obs` — scaled observation value (model-space units after z-score scaling)
- `fcst` — model background (O−F denominator)
- `ana` — analysis value
- `innov` — O−F (obs minus forecast); equivalently `obs - fcst`
- `incr` — A−F (analysis increment); equivalently `ana - fcst`; **only meaningful where `assim_flag == 1`**
- `assim_flag` — 1 = assimilated, 0 = monitored/rejected
- `species` — integer ID identifying the observing system

ASCAT species IDs:
```python
# Legacy EUMETSAT BUFR (ASCSMR02, ~25 km)
9:  Metop-A  (ASCAT_META_SM)
10: Metop-B  (ASCAT_METB_SM)
11: Metop-C  (ASCAT_METC_SM)
# H121 CDR (12.5 km Fibonacci DGG, NetCDF)
15: Metop-A  (ASCAT_HSAF_META_SM)
16: Metop-B  (ASCAT_HSAF_METB_SM)
17: Metop-C  (ASCAT_HSAF_METC_SM)
```

Obs units before scaling: degree-of-saturation % (`sfds`).
Model units: volumetric soil moisture m³/m³ (`sfmc`).
Scaling (z-score) handles the unit conversion — don't compare raw obs to model directly.

---

## Local data paths

| Dataset | Path |
|---|---|
| H121 CDR test OFA runs | `/Users/amfox/Desktop/geosldas-analysis/data/hsaf_cdr_test/` |
| Legacy ASCAT BUFR | `/Users/amfox/Desktop/ASCAT_SSM_CDR/legacy_bufr/` |
| H121 NetCDF | `/Users/amfox/Desktop/ASCAT_SSM_CDR/H121/` |
| SMOS-IC | `/Users/amfox/Desktop/SMOS_IC` |
| ISMN | `/Users/amfox/Desktop/ISMN_data` |
| GEOSldas test data | `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data` |

Long-term / full experiment data lives on NASA Discover (HPC). The `discover_JH` project has notes on that environment.

---

## Notebook hygiene

**Always clear outputs before committing.** The pre-commit hook at `.githooks/pre-commit` does this automatically on staged notebooks, but verify it fires:

```bash
git config core.hooksPath .githooks
```

To clear manually:
```bash
/Users/amfox/mamba/envs/regrid/bin/jupyter nbconvert \
  --ClearOutputPreprocessor.enabled=True \
  --to notebook --inplace <notebook.ipynb>
```

nbformat 4.4 notebooks have **no cell IDs** — the `NotebookEdit` tool will fail on them. Workaround: manipulate notebook JSON directly with a Python script using `json.load` / `json.dump`.

**f-string `\n` in notebook cell source**: a literal newline inside an f-string stored in notebook JSON will cause a `SyntaxError`. Build multi-line strings as a list and `'\n'.join(...)` instead.

**ISMN notebooks**: lazy-import `ISMN_Interface` — importing at the top of a notebook triggers a TkAgg backend switch that breaks matplotlib.

---

## Active projects (brief)

**ascat_da** — Replacing Legacy BUFR ASCAT with H121 CDR in GEOSldas. Infrastructure verified (2020-01-01 single-day test). Next step: multi-year monitor-mode run to derive H121-specific z-score scaling factors.

**M21C_ls / cygnss_da** — Combined SM+SCF manuscript. Main figures tracked in `projects/M21C_ls/figure_provenance_methods_notes.md`. CYGNSS-period figures come from `projects/cygnss_da/notebooks/cygnss_da_combined_revised_paper_figures.ipynb`.

**era5_land** — ERA5 / ERA5-Land comparison against GEOSldas OL and DA on M36 grid. Uses `compare_with_reanalysis_strict.ipynb` and postage-stamp figure notebooks.

**snow_da / SNOTEL / GHCN_snwd / IMS** — Snow validation for the manuscript. IMS uses a script-first workflow (`run_ims_ol_da_cell_metrics.py`) then a plotting notebook.

---

## Git / commit checklist

1. Clear notebook outputs (hook should do it; double-check for new notebooks)
2. Don't commit large data files — `data/` directories have `.gitignore` entries
3. Don't commit `.nc4`, `.bin`, `.mat` output files
4. Commit message style: imperative, short subject line (see `git log`)
