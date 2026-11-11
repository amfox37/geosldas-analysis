# viirs_obs_study

Obs-space study of the VIIRS snow cover fraction (SCF) obs for the GEOSldas reader: keep or drop `Basic_QA` class 1,
and JPSS-1 vs JPSS-2 (agreement, added coverage). Plan, periods, decision rules and log:
[`../../docs/viirs_obs_study_plan.md`](../../docs/viirs_obs_study_plan.md).

Python: `/Users/amfox/mamba/envs/regrid/bin/python` (numpy, netCDF4 with HDF4 and HDF5 support).

## Files

| File | What |
|---|---|
| `../../lib/scf_qc.py` | Python mirror of the QC in the GEOSldas readers `read_JPSS_SCF_hdf` (VIIRS) and `read_MODIS_SCF_hdf` (MODIS): raw readers, QC dicts, SCF |
| `../../lib/test_scf_qc.py` | unit tests of `scf_qc.py` on synthetic cells |
| `check_qc_vs_fortran.py` | checks `scf_qc.py` against the Fortran readers on real files |
| `fortran_check/` | builds a small driver around the two Fortran reader routines, copied verbatim from a GridComp checkout |
| `periods.py` | study periods P1–P7, products, file paths (shared by all scripts) |
| `get_obs.py` | downloads the VIIRS + MODIS files of all periods; writes `manifest.txt` |
| `../../lib/obs2tile.py` | CMG cell → GEOSldas tile (M36, CF0360) as `get_tile_num_for_obs`; cached maps in `data/viirs_obs_study/tilemap/` |
| `check_tile_map.py` | builds and checks the tile maps (M36 end to end against local ObsFcstAna; CF0360 counts) |
| `coverage.py` | Part A1: cells and tiles with / without QA class 1, per day; writes `data/viirs_obs_study/out/coverage/` |
| `coverage_report.py` | Part A1 tables (stdout) and figures A1-1, A1-2 for the report |
| `agreement.py` | Part A2: VIIRS vs MODIS Aqua and Terra, same day and cell, per QA class, band and MODIS SCF class; writes `data/viirs_obs_study/out/agreement/` |
| `agreement_report.py` | Part A2 tables (stdout) and figures A2-1..3 for the report (needs `coverage.py` output) |

Results go into [`../../docs/viirs_obs_study_report.md`](../../docs/viirs_obs_study_report.md); figures into
`../../docs/viirs_obs_study_report_figures/`. Outputs under `data/` are git-ignored.

## Part A1, coverage

```bash
cd /Users/amfox/Desktop/geosldas-analysis/projects/snow_da
PY=/Users/amfox/mamba/envs/regrid/bin/python
$PY scripts/viirs_obs_study/check_tile_map.py         # tile maps: build once (~10 s) and check
$PY scripts/viirs_obs_study/coverage.py --workers 4    # ~1 min for all 130 VIIRS files
$PY scripts/viirs_obs_study/coverage_report.py         # tables + figures
```

`check_tile_map.py` reproduces the SCF obs of GEOSldas ObsFcstAna files from the raw CMG files. It needs, under
`data/viirs_obs_study/`:

- `validation_m36/MYD10C1/2020/`: MYD10C1 2020-05-05..07
  (`../viirs_scf_reader_test/get_snow_cmg_obs.py <dir> 2020-05-05 2020-05-07 MYD10C1`); the M36 run
  `LS_DAv8_M36_as_test3` is local on the Desktop.
- `validation_cf0360/`: ObsFcstAna and rc_out of the `M21C_test_CF0360_jpss_*` runs (copied from Discover,
  `CF0360_scf_ObsFcstAna_202001.tar.gz`); `tiles/` with the runs' tilecoord/tilegrids (from
  `CF0360_tilecoord_tilegrids.tar.gz`, MD5s in `MD5SUMS_tiles.txt`); `cmg/` with VJ110C1, VNP10C1, MYD10C1,
  MOD10C1 2020-01-08..11.

Last run 2026-11-04: M36 17827/17827 MYD10C1 obs; CF0360 (run `5d4d636_OL`, 2020-01-09..10) 307578/307579 MYD10C1,
392970/392971 MOD10C1, 455427/455428 VJ110C1, 450102/450103 VNP10C1. Prints ALL OK at >= 99.99%.

Things that matter for the exact match (details in `lib/obs2tile.py` and the report): use the run's own tilecoord
(the older local `M21C_test_CF0360.ldas_tilecoord.bin` has other tile boxes), only its land tiles (typ 100), and
real*4 CMG cell centers with a fused multiply-add, as the Discover build computes them.

## Part A2, agreement with MODIS

```bash
cd /Users/amfox/Desktop/geosldas-analysis/projects/snow_da
PY=/Users/amfox/mamba/envs/regrid/bin/python
$PY scripts/viirs_obs_study/agreement.py --workers 4    # ~1.5 min (49 days, 130 VIIRS + 98 MODIS files)
$PY scripts/viirs_obs_study/agreement_report.py         # tables + figures; run coverage.py first
```

## Checking the QC mirror against the Fortran readers

Run this whenever the reader QC or `scf_qc.py` changes.

```bash
cd /Users/amfox/Desktop/geosldas-analysis/projects/snow_da

# 1. build the driver from the GridComp code to be checked (records the commit in built_from.txt)
bash scripts/viirs_obs_study/fortran_check/build.sh \
     /Users/amfox/Desktop/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp \
     /tmp/scf_check

# 2. compare on any VIIRS (V*10C1*.h5) and MODIS (M*10C1*.hdf) daily files
/Users/amfox/mamba/envs/regrid/bin/python scripts/viirs_obs_study/check_qc_vs_fortran.py \
     /tmp/scf_check/test_driver <files ...>

# 3. unit tests
/Users/amfox/mamba/envs/regrid/bin/python -m pytest -q lib/test_scf_qc.py
```

`build.sh` extracts `read_JPSS_SCF_hdf` and `read_MODIS_SCF_hdf` from `clsm_ensupd_read_obs.F90`, compiles them with
`io_hdf5.F90` and the HDF4 bridge of the same checkout and with stubs for the GEOSldas logging/abort routines
(`test_mod_header.F90`), and links against the Mac Baselibs (override with `BASELIBS=...`). The driver reads one
daily file for the whole globe and writes the kept obs as raw real*4. The check passes if both sides keep the same
cells in the same order with bit-identical SCF; lon/lat may differ by ~1e-5 deg because the Fortran computes the cell
centers in real*4.

Last run (2026-11-04, GridComp `6d0ed32`): all 9 files OK.

| File | Kept obs |
|---|---|
| VJ110C1.A2024015, .A2024016 | 3134159, 3257111 |
| VJ210C1.A2024015, .A2024016 | 3064167, 3189283 |
| VNP10C1.A2024015, .A2024016 | 3126945, 3246656 |
| VJ110C1.A2020015, VNP10C1.A2020015 | 3187407, 3182079 |
| MYD10C1.A2020015 | 2386249 |
