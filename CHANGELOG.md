# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `projects/viirs_snow`: VIIRS snow cover (SCF) work moved out of `projects/snow_da` (reader test, obs study, `viirs_snow_lib/`, runs); the Python package `lib` was renamed `viirs_snow_lib`. Use the commit prefix `viirs_snow:` from now on.
- `projects/viirs_snow (formerly snow_da)`: VIIRS snow cover (SCF) reader test report, scripts, configs and results (`docs/viirs_scf_reader_test_report.md`).
- Imported project directories: `projects/utils`, `ascat_da`, `cygnss_da`, `era5_land`, `snow_da`, `discover_JH`, `M21C_ls`, and `matlab2python` (calc/plot notebooks, helpers, scripts).
- Added shared assets under `common/` (Matlab, Python, Fortran utilities) and lightweight `data/` scaffolding with gitignores.

### Changed
- Updated top-level `README.md` to describe the reorganized GEOS LDAS Analysis layout and notebook hygiene.
- Commented external `addpath` entries in `projects/matlab2python/matlab_postprocess/in_situ/SMskill_vs_INSITU_single_CYG_expt.m` to reduce hard-coded dependencies.

### Fixed

### Removed

### Deprecated
