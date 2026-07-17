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
  identical `idx`/`obs`/`mod` (max abs diff 0.0).

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
