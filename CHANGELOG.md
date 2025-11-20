# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Imported project directories: `projects/utils`, `ascat_da`, `cygnss_da`, `era5_land`, `snow_da`, `discover_JH`, `M21C_ls`, and `matlab2python` (calc/plot notebooks, helpers, scripts).
- Added shared assets under `common/` (Matlab, Python, Fortran utilities) and lightweight `data/` scaffolding with gitignores.

### Changed
- Updated top-level `README.md` to describe the reorganized GEOS LDAS Analysis layout and notebook hygiene.
- Commented external `addpath` entries in `projects/matlab2python/matlab_postprocess/in_situ/SMskill_vs_INSITU_single_CYG_expt.m` to reduce hard-coded dependencies.

### Fixed

### Removed

### Deprecated
