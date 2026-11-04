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

Analysis scripts for Parts A and B will be added here as they are written.

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
