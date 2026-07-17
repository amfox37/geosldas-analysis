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
- `tests/`: local tests using fake files, so path logic can be checked without
  Discover access.

## Reader Status

- SMOS-IC uses the preprocessed sparse M36 NetCDF files
  (`idx_EASEv2_lonxlat`, `sm_obs`).
- SMAP L3 uses the existing MATLAB-parity SPL3SMP v009/R19240 QC from
  `Save_SMPL3_LDAS_tavg24_nc4_daily.m`: AM/PM are QC'd separately and averaged
  onto the M36 sparse index.
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
