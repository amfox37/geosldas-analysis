# IV / TC Python Tools

Early Python port of the independent validation (IV), triple collocation (TC),
and R-diff workflow used for GEOSldas soil-moisture evaluation.

This project is intentionally small and plain at the start. The first useful
piece is a Discover fixture collector that can dry-run the known data locations
and, once checked, copy a few representative files into `test_data/inputs/`.

## Layout

- `iv_tc/`: importable Python package.
- `scripts/collect_discover_fixtures.py`: dry-run/copy helper for small Discover
  fixtures.
- `scripts/generate_daily_pairs.py`: date/run/sensor loop that writes step-2
  daily obs/model pair files.
- `scripts/compute_climatology.py`: sensor-independent step-3 climatology over
  the compact daily pair files.
- `scripts/compute_iv.py`: sensor-independent step-4 IVd/IVs statistics using
  the compact daily pairs and pentad climatology.
- `tests/`: local tests using fake files, so path logic can be checked without
  Discover access.

## Reader Status

- SMOS-IC uses the preprocessed sparse M36 NetCDF files
  (`idx_EASEv2_lonxlat`, `sm_obs`).
- SMAP L3 uses the existing MATLAB-parity SPL3SMP v009/R19240 QC from
  `Save_SMPL3_LDAS_tavg24_nc4_daily.m`: AM/PM are QC'd separately and averaged
  onto the M36 sparse index. `soil_moisture`'s `valid_min`/`valid_max` CF
  attributes **are** enforced by default (netCDF4's default `set_auto_mask`
  behavior) -- a deliberate choice to treat SMAP's own published valid range
  as real QC, not MATLAB parity: the legacy `h5read`-based reader applies no
  such attribute masking and only ever clips `sm < 0`, so it lets through
  ~1.7% more points (spot-checked), all above `valid_max`. Pass
  `--smap-l3-matlab-valid-range-parity` to reproduce that MATLAB behavior
  instead. The underlying reader/AM-PM-merge logic was verified byte-for-byte
  against a live run of the MATLAB script for 2018-10-15/`OL` in that parity
  mode (53,586/53,586 points, identical `idx`/`obs`/`mod`).
- CYGNSS L3 reads the product-supplied daily `SM_daily` field (not a local
  average of `SM_subdaily`), applies the field's valid range and optional
  finite `SIGMA_daily` mask, then uses the MATLAB step2 spatial mapping onto
  model-land M36 cells. As with MATLAB `griddata(..., 'natural')`, no
  nearest-neighbor extrapolation is performed outside the finite-data hull.
- ASCAT H121/H139 NetCDF swaths reuse the existing `projects/ascat_da` H121 QC
  and GEOS M36 super-ob logic (`QC_DEFAULT_H121`, `form_super_obs`) instead of
  duplicating bit masks here.
- H121/H139 daily pair generation resolves exact raw files from the
  `ASCAT_SSM_CDR/flists/Y*/M*/D*/*.txt` manifests. It does not use broad
  date-globs over the monthly raw-file directories.
- The ASCAT input used by the existing MATLAB IV/TC step2 scripts is the
  processed HSAF H119/H120 daily `.mat` product
  (`H119_H120_processed/.../ASCAT_HSAF_H119_SM_YYYYMMDD_AD.mat`). That is
  distinct from the new H121/H139 NetCDF swaths.
- H119/H120 `.mat` files are read with the same core QC as the MATLAB step2
  scripts (`sm_tile <= 100`, `conf_flag_tile == 0`) and interpolated onto M36
  model-land cells. SciPy `linear` interpolation is used as the local Python
  approximation to MATLAB `griddata(..., 'natural')`. `fill_nearest` defaults
  to `False`: cells outside the interpolation convex hull stay `NaN`, matching
  the MATLAB scripts (`Save_ASCAT_LDAS_tavg24_nc4_daily.m`,
  `Save_ASCAT_LDAS_gph_daily.m`, `Save_ASCAT_SMAPL3_LDAS_daily.m`), which never
  fill those gaps — they only ever *mask further* (intersecting with the model
  land-tile mask and, for the three-way TC script, other products' coverage
  too). Pass `--h119-h120-fill-nearest` to opt into the extrapolated behavior.
- The H119/H120 reader also needs
  `ASCAT_HSAF/Auxiliary/TUW_WARP5_grid_info_2_2.nc`, which supplies the land
  GPI lat/lon lookup used by the MATLAB scripts.
- The shared interpolation helper (used by both H119/H120 and CYGNSS L3)
  triangulates on **every** coordinate-finite GPI, QC-failed ones included,
  matching MATLAB's `griddata`: its Delaunay mesh is built from all `(lon,
  lat)` pairs regardless of `z`, so a query point inside a triangle with a
  NaN vertex comes out NaN too. An earlier version filtered QC-failed points
  out *before* triangulation, which builds a sparser, all-valid mesh with no
  such NaN-poisoned triangles -- for a spot-checked H119/H120 day this more
  than doubled the point count (108,634 vs MATLAB's 48,486) by manufacturing
  coverage MATLAB's own mesh never had. Verified against a live run of
  `Save_ASCAT_LDAS_tavg24_nc4_daily.m` for 2018-10-15/`OL` after the fix:
  `sm_mod` identical on the shared cells (max abs diff 0.0), confirming the
  model-side tile matching was never the issue. `idx` and `sm_obs` are
  **not** identical, though, and re-confirmed as such on 2026-07-18 against
  a second live MATLAB run of the same date: MATLAB's 48,486 cells are
  fully contained in Python's 50,118 (a persistent ~3.4% coverage superset,
  not day-to-day noise -- see the Triple Collocation section below for how
  this compounds at the climatology/TC level), and `sm_obs` on the 48,486
  shared cells has mean abs diff 0.12, max abs diff 6.86 percentage points.
  This residual is attributed to the remaining `linear`-vs-`natural`
  `griddata` method difference, not a new bug.

## Fixture Dry Run

From the repository root:

```bash
python projects/iv_tc/scripts/collect_discover_fixtures.py \
  --date 2020-05-15 \
  --run OLv7_M36_MULTI_type_13_H121=/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121
```

The default mode is dry-run: it prints found/missing source files and the
portable destination path it would use. Add `--copy` only after reviewing the
dry-run output.

## Daily Pair Generation

Step 2 pair generation resolves the observation product, model daily file, and
tilecoord for each requested date/run/sensor, reads a `DailyPair`, then writes:

```text
output_root/step2_pairs/{sensor}/{run_name}/YYYYMMDD.npz
```

Each `.npz` contains the compact downstream arrays used by the Python IV/TC/Rdiff
workflow:

- `idx0`: zero-based M36 sparse cell index (`i + j * 964`)
- `sm_obs`: observation values
- `sm_mod`: matched GEOSldas model values

Example:

```bash
python projects/iv_tc/scripts/generate_daily_pairs.py \
  --start-date 2018-10-15 \
  --end-date 2018-10-15 \
  --sensor smosic \
  --sensor ascat_h121 \
  --run OL=/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121 \
  --output-root /discover/nobackup/projects/land_da/Evaluation/IVs/output_python
```

Use explicit ASCAT product names (`ascat_h121` or `ascat_h119_h120`) so the
current CDR NetCDF and processed HSAF `.mat` paths cannot be mixed up silently.

## Pentad Climatology

Step 3 consumes the generic daily pair files, so the same implementation works
for every observation product. It follows the MATLAB `Compute_clim*.m` logic:

- each sample contributes to a circular 31-day calendar window (day +/- 15);
- February 29 contributes to February 28 on the 365-day climatology calendar;
- means are retained only where the cell/pentad count reaches `Nday_min`;
- the 365 daily windows are represented by their 73 pentad centers
  (day-of-year 3, 8, ..., 363).

Unlike MATLAB's full 964 x 406 arrays, the Python file stores only M36 cells
seen in at least one input pair. `idx0` maps each row back to the global M36
index. The remaining array names preserve the MATLAB terminology:
`obs_sm_clim`, `mod_sm_clim`, and `N_sm_clim`.

```bash
python projects/iv_tc/scripts/compute_climatology.py \
  --start-date 2018-08-01 \
  --end-date 2024-06-30 \
  --sensor smosic \
  --sensor ascat_h121 \
  --run OL \
  --pair-root /discover/nobackup/projects/land_da/Evaluation/IVs/output_python
```

The default count threshold matches MATLAB's `4 * (end_year - start_year)`,
with the command's inclusive end date converted to MATLAB's end-exclusive
boundary. Use `--min-count` to make a short test range meaningful.

**Verified bit-exact against a live MATLAB run** (`Compute_clim_SMOSIC.m`):
`smosic`/`OLv8_M36_cd`, 2018-08-01 to 2024-06-29 (2,160 daily pairs), against
the archived `SMOSIC_OLv8_M36_cd_clim_pentad_201808_202405_w31.mat`. `N_sm_clim`
matched with 0/4,580,500 mismatches; `mod_sm_clim`/`obs_sm_clim` matched with
max abs diff 0.0 across 4,409,554 valid cell-pentads.

**Also verified against `Compute_clim_ASCL4.m`** (ASCAT HSAF H119, same
2018-08-01..2024-06-29 span): 0/3,598,733 `N_sm_clim` mismatches, 0 mask
mismatches on 3,552,521 valid cell-pentads, max abs diff ~7.6e-6 (this run
used MATLAB's real step-2 `ASCL4_H119_SMSF_L4_*_QC_1_*.mat` files converted
directly to the `idx0`/`sm_obs`/`sm_mod` npz schema -- `idx_EASEv2_lonxlat`
is MATLAB's 1-based index, so subtract 1). **Do not** use the archived
`python_output/step2_pairs.tar.gz` `ascat/` pairs for this: that's a
separate, older Python port that over-covers relative to MATLAB by ~50%
(70,633 vs 46,688 cells on a spot-checked day) -- almost certainly the same
griddata/QC-filtering-before-triangulation bug already found and fixed in
this package's own reader (see the H119/H120 triangulation note above), but
in different, unfixed code. It is not a valid MATLAB-parity reference.

```bash
python projects/iv_tc/scripts/compare_climatology_matlab.py \
  --python /path/to/step3_climatology/smosic/OL/20180801_20240629_w31.npz \
  --matlab /path/to/SMOSIC_OLv8_M36_cd_clim_pentad_201808_202405_w31.mat
```

**Verified against a live MATLAB run on the current OL run**
(`OLv7_M36_MULTI_type_13_H121`, not any archived data): `smap_l3`/`OL`,
2015-04-01 to 2017-03-31 (727/731 daily pairs, matching this project's known
isolated SMAP gaps). `Save_SMPL3_LDAS_tavg24_nc4_daily.m` /
`Compute_clim_SMPL3.m` / `Compute_IVd_IVs_SMPL3.m` were copied, repointed at
`mod_path`/`mod_run` for the real current OL run and a scratch `out_path`/
`data_path`, and run live end-to-end via a SLURM job (`matlab/R2023b`,
~9.5 minutes for all three stages). `N_sm_clim` matched with 0/6,273,320
mismatches, `mod_sm_clim`/`obs_sm_clim` matched with 0 mask mismatches on
5,938,990 valid cell-pentads (max abs diff ~2.4e-7).

This required generating the Python step-2 pairs with
`--smap-l3-matlab-valid-range-parity`: the default reader QC (CF
`valid_min`/`valid_max` enforcement) is a deliberate, documented deviation
from MATLAB's own `h5read`-based QC (see Reader Status above), and using the
default pairs against a live MATLAB run reproduces that same ~0.3%
coverage gap here, not a step-3 bug -- confirmed by rerunning with the
parity flag and getting an exact match.

## Independent Validation

Step 4 follows `Compute_IVd_IVs*.m` directly. For each current day `d`, it
requires the exact calendar day `d - Nlag`, intersects their sparse M36
indices, and subtracts the observation and model climatologies for each day's
own pentad. Missing dates are skipped; they are not compressed into adjacent
samples. Defaults match MATLAB: `Nlag=2` days and `Nmin=100` samples per cell.

The population-moment equations are intentionally preserved as written in
MATLAB, including its use of the current anomaly means in the lagged
covariances. The finite-value behavior is also preserved: MATLAB selects
samples using the current model anomaly only, so a NaN in another anomaly can
propagate into that cell's sums even though `N_sm` is incremented. Undefined
`R_mod_obs` values retain MATLAB's `-9999` sentinel; undefined IVd/IVs fields
remain NaN, and valid squared correlations are capped at 1.

```bash
python projects/iv_tc/scripts/compute_iv.py \
  --start-date 2018-08-01 \
  --end-date 2024-06-29 \
  --climatology-end-date 2024-06-30 \
  --sensor smosic \
  --sensor ascat_h121 \
  --run OL \
  --pair-root /discover/nobackup/projects/land_da/Evaluation/IVs/output_python
```

Outputs are written to:

```text
output_root/step4_iv/{sensor}/{run_name}/START_END_lag2_n100.npz
```

The sparse file stores `idx0` plus the MATLAB fields `N_sm`, `R2_ivd_mod`,
`R2_ivd_obs`, `R2_ivs_mod`, `R2_ivs_obs`, `R_mod_obs`, `Nmin`, and `Nlag`,
along with run, units, date-range, climatology, and day-accounting metadata.

The Python CLI treats both start and end dates as inclusive. The MATLAB scripts
stop before `end_time`; therefore a MATLAB configuration with
`end_time = 2024-06-30` processes through 2024-06-29 and should be compared
with Python `--end-date 2024-06-29`.

For a direct Discover parity check, expand the sparse Python fields and compare
them with the full-grid MATLAB output using:

```bash
python projects/iv_tc/scripts/compare_iv_matlab.py \
  --python /path/to/step4_iv/smosic/OL/20180801_20240629_lag2_n100.npz \
  --matlab /path/to/SMOSIC_OLv8_M36_cd_IVD_IVS_stats_lag2day_201808_202405.mat
```

**Verified against the same live MATLAB run** (`Compute_IVd_IVs_SMOSIC.m`),
using the step-3 climatology above as input: `N_sm` and every field's
valid/undefined mask matched MATLAB exactly (0 mismatches across ~50-62k
cells, 5 fields). Numeric values agreed to ~1e-7 mean / ~1e-5 max absolute
difference (`--rtol 1e-4 --atol 1e-4`) rather than bit-exact, traced to this
pipeline's `float32` on-disk pair storage (`sm_obs`/`sm_mod` in step-2 NPZ
files) versus MATLAB's double-precision accumulation throughout -- not a
logic difference, since step 3 (same inputs) matched bit-exact.

This run also confirmed the `r_model_obs[~np.isfinite(...)]` zero-variance
guard (vs. MATLAB's `isnan`-only check, which leaves stray `Inf`/`-Inf` in
place) never actually diverges on real data: recomputing with an
`isnan`-only fill produced byte-identical `R_mod_obs` across all 62,276
cells, confirming the edge case is floating-point-noise-only and does not
occur with real soil-moisture variance.

**Also verified against `Compute_IVd_IVs_ASCL4.m`**, using the same
MATLAB-native step-2 conversion described above: `N_sm` and every field's
valid/undefined mask matched exactly (0 mismatches, 42,116-50,925 cells per
field), numeric agreement ~1e-7 mean / ~1e-6 max absolute difference --
tighter than the `smosic` run since this path has only one float64-to-
float32 cast (MATLAB double straight to npz) instead of two.

**Verified against a live MATLAB run on the current OL run**, same
`smap_l3`/`OL`/2015-04-01..2017-03-31 setup as the step-3 entry above (live
`matlab/R2023b` SLURM job, no archived data): `N_sm` and every field's
valid/undefined mask matched exactly (0 mismatches, 71,101-97,310 cells per
field, using `--smap-l3-matlab-valid-range-parity` pairs so the comparison
isolates step-3/step-4 math from the deliberate reader-QC deviation).
Numeric agreement ~1.5e-7 mean / ~3e-5 max absolute difference.

The command checks `N_sm` exactly and reports valid-mask and numeric
differences for every saved IV field. It exits nonzero on any mismatch beyond
the requested tolerances.

## Triple Collocation

The TC engine follows `common/matlab/TC/Compute_TC_ASCAT_SMOSIC_MOD.m` for
the two-daily-pair workflow. The three source columns are always ordered as
primary observation, model, and secondary observation. For the reference
SMOS-IC/model/ASCAT calculation this means `(SMOS-IC, model, ASCAT)`.

The implementation preserves the details of that MATLAB script:

- both daily sparse files must exist, and only their exact common M36 cells
  contribute;
- model values come from the primary pair, while the model climatology comes
  from the secondary climatology file, matching the script as written;
- ASCAT observations and their climatology are multiplied by `0.005`
  (MATLAB `/200`) before anomaly calculation;
- all three anomalies must be finite for `N_sm` to increment;
- annual, May-September summer, and October-April winter modes are available;
- `Nmin=20`, population covariance (`sum/N`), the `R2 < 0.001` mask, and the
  `R2 > 1` cap match MATLAB;
- negative variances become NaN, but negative cross-covariances are retained,
  as in the current script where the old covariance masks are commented out.

For the direct SMOS-IC/ASCAT comparison on Discover:

```bash
python projects/iv_tc/scripts/compute_tc.py \
  --start-date 2018-08-01 \
  --end-date 2024-06-29 \
  --primary-sensor smosic \
  --secondary-sensor ascat_h119_h120 \
  --primary-name SMOS \
  --secondary-name ASC \
  --run OL \
  --pair-root /discover/nobackup/projects/land_da/Evaluation/IVs/output_python
```

`--secondary-scale` defaults to `0.005` for an ASCAT sensor and `1` for other
products. Other defaults reproduce MATLAB (`--model-values-source primary`,
`--model-climatology-source secondary`, `--season annual`, `--min-count 20`).
Missing daily pairs are fatal unless `--allow-missing` is supplied.

Outputs are written to:

```text
output_root/step4_tc/{primary}__{model}__{secondary}/{run}/START_END_SEASON_n20.npz
```

The file stores sparse `idx0` and `N_sm` vectors plus `(cell, 3)` arrays:
`variance`, `R2_TC`, and `sigma2` use source-column order `(0, 1, 2)`, while
`covariance` and `correlation` use pair-column order `(0,1), (0,2), (1,2)`.
Names, units, scale factors, date range, seasonal mode, climatology/model
source choices, and day-accounting diagnostics are included as metadata.

The Python CLI uses an inclusive end date. To match the MATLAB configuration
with `end_time = 2024-06-30`, use Python `--end-date 2024-06-29`. If step 3
was built through the MATLAB boundary date, as in the reference workflow, add
`--climatology-end-date 2024-06-30` so the CLI resolves that climatology file.
Compare the result directly with the reference MAT file using:

```bash
python projects/iv_tc/scripts/compare_tc_matlab.py \
  --python /path/to/step4_tc/SMOS__model__ASC/OL/20180801_20240629_annual_n20.npz \
  --matlab /path/to/ASCL4_SMOSIC_OLv8_M36_cd_TC_stats_201808_202405.mat
```

The comparator checks `N_sm` exactly, then reports finite-mask and numeric
differences for all three R2, error-variance, pair-correlation, and
cross-covariance fields. Its default MATLAB variable mapping matches
`Compute_TC_ASCAT_SMOSIC_MOD.m`; custom comma-separated mappings support the
other three-input TC scripts.

**Known limitation, quantified against a live MATLAB run (2026-07-18):**
running the full step2→step3→step4 chain live for `smosic`/`ascat_h119_h120`
against the current OL (`OLv7_M36_MULTI_type_13_H121`, 2018-08-01..2019-07-31)
does **not** reach the tight parity the `smap_l3` and MATLAB-native-conversion
runs above achieve. This is *not* a step-3/step-4 logic bug -- it traces
entirely to the pre-existing `ascat_h119_h120` step-2 interpolation
difference (`griddata` `linear` vs. MATLAB's `natural`; see Reader Status),
which turns out to be **persistent per-cell** rather than random day-to-day
noise: specific boundary/coastal cells are covered by Python's interpolation
on every day of the year and never by MATLAB's (or vice versa), so a ~3%/day
coverage difference compounds linearly instead of averaging out.

Measured impact:
- **Step-2** (one spot-checked day, 2018-08-01): Python covers 62,557 cells
  vs. MATLAB's 60,530 (ratio 1.034, MATLAB fully contained in Python's set,
  consistent with the known residual). `sm_mod` bit-identical; `sm_obs` mean
  abs diff 0.15 percentage points, max 10.8 on the shared cells.
- **Step-3 climatology** (5 pentads spot-checked across the year): coverage
  ratio ~1.024; `N_sm_clim` (obs count per cell/pentad) differs on ~16% of
  shared cells, by up to 27 out of the 31-day window -- the compounding
  effect described above. Obs value agreement on shared cells stays small
  (mean abs diff ~0.2 percentage points).
- **Step-4 TC**: `N_sm` mismatched on ~35% of cells (28,025/79,097).
  `sigma2`'s division by a near-zero `C_SMOS_mod`/`C_SMOS_ASC`/`C_mod_ASC`
  cross-covariance (cells where the two source anomalies are essentially
  uncorrelated -- weak model skill / complex terrain) amplifies the small
  `N_sm` divergence into large-looking errors: of 75,446 common valid
  cells, the 19,090 (~25%) with `|C_SMOS_mod| < 1e-4` carry essentially all
  of it (mean abs diff 0.013, max 118.9 on `sigma2_ASC`); the other 56,356
  well-conditioned cells agree tightly (mean abs diff 0.00012, max 2.27).
  This near-zero-denominator sensitivity is an inherent property of the TC
  error-variance estimator at weakly-correlated triplets, not a code defect.

Not currently planned: closing the `ascat_h119_h120` `linear`-vs-`natural`
gap further, or masking near-zero-denominator cells out of the parity check.
If tighter `ascat_h119_h120`-driven TC/climatology parity is needed later,
start there.
