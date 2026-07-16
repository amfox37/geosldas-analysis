# IV / TC Python port - working plan

Status as of 2026-07-15. Written for whoever, human or AI, picks this up next,
possibly in a different coding session on a different machine. This file is the
durable record; do not rely on chat history.

## Why

We need robust Python code for independent validation (IV), triple collocation
(TC), and R-diff analysis across multiple GEOSldas experiments and multiple
observation products. The immediate production driver is the HSAF ASCAT CDR
comparison on Discover:

| Role | EXP_ID | Status as of 2026-07-15 |
|---|---|---|
| OL baseline | `OLv7_M36_MULTI_type_13_H121` | finished; also source run for scaling stats |
| DA, H121 CDR ASCAT assimilated | `DAv7_M36_ASCAT_type_13_H121` | running on Discover when this was first written |
| DA, legacy EUMETSAT ASCAT assimilated | `DAv7_M36_ASCAT_type_13_legacy` | running on Discover when this was first written |

Those experiments use the per-species z-score ASCAT scaling stats validated in
`projects/obs_scaling_params/`. Once the runs finish, we need IV/TC skill
metrics for a 3-way comparison: OL vs DA-H121 vs DA-legacy.

That immediate HSAF comparison should be the first use case, but the Python port
should not be a one-off HSAF script. The real target is a maintainable IV/TC
toolkit that can handle SMAP, SMOS-IC, ASCAT, and potentially CYGNSS through the
same interfaces.

## Inventory: what exists today

### MATLAB: current canonical capability

The MATLAB pipeline under `common/matlab/` is still the broadest and most
production-like implementation.

- `common/matlab/IVs/scripts/step2/Save_SMPL3_LDAS_tavg24_nc4_daily.m`
  extracts daily SMAP L3/model pairs.
- `common/matlab/IVs/scripts/step2_asc/Save_ASCAT_LDAS_tavg24_nc4_daily.m`
  extracts daily ASCAT/model pairs.
- `common/matlab/IVs/scripts/step2_smosic/Save_SMOS_IC_LDAS_tavg24_nc4_daily.m`
  extracts daily SMOS-IC/model pairs.
- `common/matlab/IVs/scripts/step2_cyg/Save_CYGNSS_LDAS_tavg24_nc4_daily.m`
  extracts daily CYGNSS/model pairs.
- `common/matlab/IVs/scripts/step3/Compute_clim*.m` computes pentad
  climatologies.
- `common/matlab/IVs/scripts/step4/Compute_IVd_IVs*.m` computes IVd/IVS
  lag-2-day anomaly correlations.
- `common/matlab/IVs/scripts/step5/Rdiff_script_*.m` computes R-diff between
  experiments.
- `common/matlab/TC/Compute_TC_two_input_files.m` and
  `common/matlab/TC/Compute_TC_three_input_files.m` compute TC metrics from
  daily matched datasets.

SMAP L3 status: the live/current MATLAB path uses **SPL3SMP v009 / R19240**.
Confirmed code evidence:

- `common/matlab/IVs/scripts/step2/Save_SMPL3_LDAS_tavg24_nc4_daily.m` reads
  `SPL3SMP_v009`, matches `SMAP_L3_SM_P_*_R19240_*.h5`, and writes
  `SMPL3_R19...`.
- `common/matlab/IVs/scripts/step3/Compute_clim_SMPL3.m` consumes
  `SMPL3_R19...`.
- `common/matlab/IVs/scripts/step4/Compute_IVd_IVs_SMPL3.m` consumes
  `SMPL3_R19...`.
- `common/matlab/TC/Compute_TC_two_input_files.m` and
  `Compute_TC_three_input_files.m` consume `SMPL3_R19...`.

Older SMAP v008/R18290 paths still exist, for example
`common/matlab/IVs/scripts/step2/Save_SMPL3_LDAS_gph_nc4_daily.m` and old
monolithic TC code, but they are disconnected from the live step2 -> step3 ->
step4 -> TC chain. Treat them as reference only, not as parity targets.

### Python: current partial capability

Python currently has useful pieces, but not a full production port.

- `projects/SMOS_IC/scripts/preprocess_smos_ic_daily_to_m36.py` is a real
  script for preprocessing SMOS-IC daily files onto the M36 grid.
- `projects/SMOS_IC/notebooks/ivs_tc_ascat_smosic_python_workflow.ipynb` is a
  working notebook-only port for ASCAT + SMOS-IC + model. It includes step2-like
  daily matching, climatology, IV, TC, and R-diff logic.
- That notebook is hardcoded to the CYGNSS-era `OLv8_M36_cd` / `DAv8_M36_cd`
  experiments and to one ASCAT + SMOS-IC + model triplet.
- The notebook TC solver is inlined and 3-input only.
- There is no Python SMAP L3 reader yet.
- `common/python/stats/` currently has no reusable IV/TC code.
- `projects/matlab2python/matlab_postprocess/TC|IVs/` is reference MATLAB copied
  for translation, not Python production code.

## Design target

Create a dedicated Python project/package for IV and TC, rather than continuing
to grow a notebook or burying the implementation inside `matlab2python`.

Recommended layout:

```text
projects/iv_tc/
  README.md
  pyproject.toml or project-local environment notes
  iv_tc/
    __init__.py
    config.py
    dates.py
    grid.py
    pairs.py
    climatology.py
    iv.py
    tc.py
    rdiff.py
    workflow.py
    io/
      __init__.py
      model.py
      smap_l3.py
      smos_ic.py
      ascat.py
      cygnss.py
      pair_store.py
  scripts/
    build_daily_pairs.py
    compute_climatology.py
    compute_iv.py
    compute_tc.py
    compute_rdiff.py
    run_workflow.py
  tests/
```

This can later be promoted into `common/python/` if it becomes a shared library,
but the first production work should live in its own project directory where we
can document assumptions, keep fixtures, and test against MATLAB outputs without
disturbing unrelated common code.

## Core abstractions

The important design decision is to separate data adapters from statistical
engines. Sensor-specific logic should end at the daily matched-pair boundary.
After that, IV, TC, climatology, and R-diff should operate on generic arrays.

Suggested records:

```python
@dataclass
class SparseObservation:
    date: datetime.date
    sensor: str
    idx: np.ndarray          # M36/EASE/global cell index
    values: np.ndarray       # observation values, documented units
    units: str
    qc_summary: dict

@dataclass
class DailyPair:
    date: datetime.date
    sensor: str
    run: str
    idx: np.ndarray          # matched grid cells
    obs: np.ndarray          # observed values
    model: np.ndarray        # model values sampled at idx
    obs_units: str
    model_units: str
```

The contract should be explicit:

- `idx`, `obs`, and `model` are 1-D arrays after masking invalid obs/model cells.
- Missing values are NaN internally, not sentinel literals.
- Units are carried with the data and checked before computing metrics.
- All daily pair files must identify sensor, experiment/run, date, grid, units,
  and source product/version.

## Sensor adapters

Each sensor adapter should expose the same high-level method:

```python
read_daily_observation(date, config) -> SparseObservation
```

The model adapter then samples the requested GEOSldas run at the observation
cells to produce a `DailyPair`.

Initial adapters:

- **ASCAT**: first production target for HSAF CDR comparison. Support the H121
  CDR and legacy EUMETSAT paths needed by the current experiments. Preserve
  enough metadata to distinguish H121 vs legacy in filenames and outputs. H121
  must be resolved through the day-indexed `flists` tree, not by trying to infer
  dates from raw filenames.
- **SMOS-IC**: wrap the existing preprocessed M36 daily product from
  `projects/SMOS_IC/scripts/preprocess_smos_ic_daily_to_m36.py`. Exclude
  macOS AppleDouble `._*.nc` files if globbing the flat preprocessed directory.
- **SMAP L3**: port the live MATLAB SPL3SMP v009/R19240 reader. Do not port the
  disconnected R18290/v008 path except as reference while checking logic. The
  real NSIDC-downloaded path includes the extra
  `n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.009/YYYY.MM.DD/` nesting.
- **CYGNSS**: include the adapter slot and interface from the start. It does not
  need to be the first production adapter, but the design must not require a TC
  rewrite when CYGNSS is added.

Adapter-specific responsibilities:

- Locate files for a date.
- Read raw variables and QC flags.
- Apply documented QC and fill-value handling.
- Convert or validate units using product metadata or MATLAB-documented rules.
- Map observations to the analysis grid.
- Return sparse daily observations with no model or statistics logic embedded.

## Confirmed Discover data locations

These paths were confirmed on Discover on 2026-07-15 and should be used as the
starting defaults for a fixture-copy script and for production config templates.

### ASCAT H121 CDR

The H121 raw product root is defined in the GEOSldas namelist driving the
current HSAF DA runs:

```text
/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/H121/{metop_a,metop_b,metop_c}/Y{YYYY}/M{MM}/
```

Raw filenames look like:

```text
W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOP{A,B,C}-12.5km-H121_C_LIIB_<creation_ts>_<sensing_start>_<sensing_end>____.nc
```

Do **not** glob raw H121 files by date. All granules for a month are flat under
`Y{YYYY}/M{MM}/`, and the date is embedded in the sensing-time fields rather
than in a stable date-specific filename prefix. The authoritative date lookup is
the manifest tree used by the Fortran reader
`clsm_ensupd_read_obs.F90:read_obs_sm_ASCAT_HSAF`:

```text
/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/flists/Y{YYYY}/M{MM}/D{DD}/H121_METOPA.txt
/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/flists/Y{YYYY}/M{MM}/D{DD}/H121_H139_METOPB.txt
/discover/nobackup/projects/land_da/ASCAT_SSM_CDR/flists/Y{YYYY}/M{MM}/D{DD}/H121_H139_METOPC.txt
```

Each manifest is plain text with one bare filename per line. The Python adapter
should join each listed filename to the platform-specific monthly raw root.

### Legacy ASCAT H119/H120

The existing MATLAB/Python IV workflow reads the processed daily legacy ASCAT
`.mat` files here:

```text
/discover/nobackup/qliu/merra_land/DATA/ASCAT_HSAF/H119_H120_processed/Y{YYYY}/M{MM}/ASCAT_HSAF_H119_SM_{YYYYMMDD}_AD.mat
```

### SMOS-IC

The preprocessed M36 daily files are a flat directory:

```text
/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily/smos_ic_sm_m36_{YYYYMMDD}.nc
```

Ignore `._smos_ic_sm_m36_*.nc` AppleDouble sidecar files if present.

### SMAP L3

The live SMAP product is SPL3SMP v009 / R19240 under the NSIDC downloader
directory structure:

```text
/discover/nobackup/projects/land_da/Evaluation/IVs/data/SPL3SMP_v009/n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.009/{YYYY.MM.DD}/SMAP_L3_SM_P_{YYYYMMDD}_R19240_*.h5
```

Do not use the older simplified `SPL3SMP_v009/Y{YYYY}/...` pattern; that was a
MATLAB-script assumption and not the confirmed on-disk structure.

### GEOSldas model output

For the current HSAF experiments, model output and tile coordinates follow:

```text
<run_root>/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/Y{YYYY}/M{MM}/
<run_root>/output/SMAP_EASEv2_M36_GLOBAL/rc_out/<run>.ldas_tilecoord.bin
```

The current DA-H121 and DA-legacy runs already have partial output and tilecoord
files on disk while they are mid-flight.

## Model adapter

The model adapter should handle GEOSldas output layout, variable naming, and
daily averaging/sampling.

Responsibilities:

- Support standard output roots like
  `<exp>/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/...`.
- Read the model soil moisture variable used by the MATLAB IV/TC chain.
- Preserve documented units and fail loudly on unexpected units.
- Sample the model field at the sparse observation cell indices.
- Hide differences between daily/tavg24 files and any fallback 3-hourly files
  behind one interface.

The IV/TC engines should not know how model files are organized.

## Pair store

Daily matched pairs are the natural checkpoint between slow I/O and fast stats.
The first implementation can keep `.npz` support because the existing notebook
does that, but production should prefer a self-describing format:

- NetCDF4 is the safest first production store because the existing repo already
  uses NetCDF heavily.
- Zarr is worth considering if many years, sensors, and experiments make
  repeated partial reads a bottleneck.

Minimum daily pair metadata:

- `date`
- `sensor`
- `run` / `experiment`
- `grid`
- `obs_product`
- `obs_version`
- `model_variable`
- `obs_units`
- `model_units`
- QC/masking summary counts
- code version or git SHA when available

## Statistical engines

### Climatology

Port MATLAB step3 into reusable Python:

```python
compute_pentad_climatology(pairs, *, sensor, run, years, min_count, grid)
```

This should be independent of the sensor and model reader. It consumes daily
pairs and writes climatology products with enough metadata to support IV/TC
auditing.

### IV

Port MATLAB step4 into reusable Python:

```python
compute_iv(pairs, climatology, *, lag_days=2, metric="anom_corr")
```

Keep IVd and IVS names if they match the MATLAB outputs, but make the function
interfaces clear about what is being correlated and how anomalies are formed.

### R-diff

Generalize MATLAB step5 to arbitrary experiment pairs:

```python
compute_rdiff(iv_a, iv_b, *, label_a, label_b)
```

For the HSAF comparison, expected pairwise outputs are:

- DA-H121 minus OL
- DA-legacy minus OL
- DA-H121 minus DA-legacy

### TC

TC should accept named inputs instead of hardcoded ASCAT/SMOS/model variable
names. Classic TC is still a 3-way covariance calculation, but the implementation
should be generic over which three data sources are supplied:

```python
compute_tc(inputs=[
    NamedSeries("ASCAT_H121", ...),
    NamedSeries("SMOS_IC", ...),
    NamedSeries("MODEL_DA_H121", ...),
])
```

This supports:

- ASCAT + SMOS-IC + model
- ASCAT + SMAP + model
- SMOS-IC + SMAP + model
- CYGNSS + SMAP + model
- other valid triplets, provided the assumptions are documented

Do not hardcode the current notebook triplet into the TC engine. Let workflow
configuration decide which triplets to run.

## Workflow for the HSAF comparison

First production workflow:

1. Build daily pairs for each experiment and sensor needed for validation.
2. Compute pentad climatologies for each `(run, sensor)` pair.
3. Compute IV for each `(run, sensor)` pair.
4. Compute R-diff for OL vs DA-H121, OL vs DA-legacy, and DA-H121 vs DA-legacy.
5. Compute TC triplets for each experiment using the selected independent
   observation pair plus model.
6. Generate summary tables/maps only after the numeric outputs are validated.

The workflow should take a config file or command-line options listing:

- experiment names and roots
- date range
- sensors to read
- TC triplets to compute
- output root
- grid definition
- min-count thresholds

## Validation and parity plan

This needs testing at three levels.

1. Formula unit tests:
   - climatology binning
   - anomaly calculation
   - IV correlation
   - TC covariance/error-variance formulas
   - R-diff sign convention

2. Adapter tests with tiny fixtures:
   - fill values and NaNs
   - units
   - QC masks
   - grid index mapping
   - file-not-found behavior

3. MATLAB parity tests:
   - For one short date range and one or two sensors, compare Python daily pair
     files against MATLAB step2 outputs.
   - Compare climatology, IV, R-diff, and TC fields against MATLAB outputs where
     available.
   - Treat SMAP R19, ASCAT, and SMOS-IC as the highest-value parity targets.

For TC, exact bitwise parity may not be realistic if masking/order differs, but
the valid-cell mask, covariance inputs, and final metrics should agree within a
documented tolerance.

## Staged implementation

Recommended order:

1. Create `projects/iv_tc/` package skeleton and move the working notebook
   functions into modules with minimal behavior changes.
2. Add tests around the extracted climatology, IV, TC, and R-diff functions using
   synthetic arrays.
3. Implement the daily pair abstraction and pair-store writer/reader.
4. Wrap the existing ASCAT and SMOS-IC notebook logic as adapters.
5. Run Python against the existing OLv8/DAv8 notebook-era case and compare with
   the notebook outputs.
6. Add config-driven N-experiment workflows.
7. Run the HSAF OL/DA-H121/DA-legacy workflow.
8. Add the SMAP L3 v009/R19240 adapter and parity-check against MATLAB step2.
9. Add CYGNSS adapter support when its product contract and target experiment
   are clear.

## Explicit non-goals

- Do not preserve parity with dead SMAP v008/R18290 scripts unless a downstream
  consumer is found.
- Do not make CYGNSS the first production target. Do keep CYGNSS compatible with
  the adapter and TC-input design from day one.
- Do not let notebooks become the production implementation. Notebooks can drive
  or inspect outputs, but core logic belongs in tested modules.
- Do not guess units or conversions from plot appearance. Units must come from
  product metadata, GEOSldas metadata, or documented MATLAB logic.
- Follow the repo's notebook hygiene rules when notebooks are touched: clear
  outputs before commit, and watch for known f-string newline / no-cell-id
  issues described in the root project guidance.
