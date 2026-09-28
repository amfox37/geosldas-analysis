# cygl1_operator_test

Development and validation of a GEOSldas observation operator for
**preprocessed CYGNSS L1 DDM 3x5 "crop scalar" coefficients**
(`obs_param_nml(56)`, `CYGNSS_L1_DDM3X5_CROP_SCALAR`, `varname = 'cygl1scal'`),
run over an Arizona box (`MINLON -118 MAXLON -106 MINLAT 29 MAXLAT 40`).
Separate from the CYGNSS L3 soil-moisture species (`obs_param_nml(54)`).

*A copy of this file is tracked in git as
`geosldas-analysis/projects/CYGNSS_L1_AZ/docs/cygl1_operator_test_project_README.md`, and all
run configs are in `runs/configs/` there. Keep the two copies in sync.*

## Status (2026-09-28)

The L1 operator bug was fixed on 2026-09-25 (details below). Everything since then was rerun from scratch on
the fixed build (the "fixedop" experiments). All pre-fix results are invalid and archived in `archive/`.

**Fixed-operator experiments** (all ended clean, with 0 `LDAS ERROR`/`forrtl`):

| EXP_ID | assimilated | period | templates/ (NML_INPUT_PATH) |
|---|---|---|---|
| `OLv8_M36_AZ_fixedop` | none (unscaled OL, all 13 species monitored) | 2020–2022 | `OL_monitor_unscaled` |
| `DA_L1_full_fixedop` | full L1, xcorr 0.625 (ill-conditioned-R test) | Jan 2020 | `DA_L1_full` |
| `DA_L1_full_xc015_fixedop` | full L1, xcorr 0.15, errstd 2.75 (benchmark) | 2020 | `DA_L1_full_xc015` |
| `DA_L1_full_xc015_err39_fixedop` | benchmark with errstd 3.9 | 2020 | `DA_L1_full_xc015_err39` |
| `DA_L1_full_xc015_coh040216_fixedop` | benchmark + coherency filter 0.40–2.16 + rebuilt L1 clim | **2020–2022** | `DA_L1_full_xc015_coh040216` |
| `DA_L1_full_xc015_coh040216_err39_fixedop` | filter + clim + errstd 3.9 (**best L1 arm by O-F**) | **2020–2022** | `DA_L1_full_xc015_coh040216_err39` |
| `DA_L1_dense075_coh05_fixedop` | thinned L1 (pre-fix "best" config) | 2020 | `DA_L1_dense075_coh05` |
| `DA_L3_fixedop` | CYGNSS L3 (`CYGNSS_SM_6hr`), full L1 monitored | **2020–2022** | `DA_L3_fullL1mon` |
| `DA_SMAP_fixedop` | SMAP L1C Tb, full L1 monitored | **2020–2022** | `DA_SMAP_fullL1mon` |

The four 2020–2022 arms were extended by editing `run/CAP.rc` END_DATE (backups `run/CAP.rc.bak_end20210101`).
Their exeinp files in the project root still show the original END_DATE.

**Key results and where they are** (in `geosldas-analysis/projects/CYGNSS_L1_AZ/`):

- `docs/cygl1_fixedop_2020_arm_comparison_report.md`: O-F arm comparison for 2020 plus the 2021–2022 extension, the noise-vs-gain
  split and next steps.
  - Filter + errstd 3.9 is the best L1 arm by O-F. Over 2020–22: SMAP −1.20%, L3 −2.08%, ASCAT +0.55% vs OL.
  - The spring degradation recurs every year, with a moving window and two mechanisms.
- `runs/cygl1_fixedop_ismn_validation_2020_2022.md`: ISMN in-situ validation, **scored over the L1 area only**. 78% of the ISMN
  stations in the box sit on tiles with no L1 obs.
  - L1 improves root-zone R by +0.03 to +0.04 and beats L3 in situ. L3 is neutral.
  - Over the L1 area, O-F shows L1 and L3 tied on Tb, and L3's O-F lead is ASCAT plus its own L3 fit.
  - It also has the SOILSCAPE per-tile section.
- `notebooks/cygl1_ismn_skill_figures.ipynb`: CYGNSS-paper-style in-situ figures.
- `docs/cygl1_obs_quality_vs_qc_report.md`: L1 obs quality vs coherency_ratio. This motivated the 0.40–2.16 filter.
- `docs/cygl1_obs_error_correlation_report.md`: L1 obs-error correlation and why the full stream needs a short xcorr.

**Transfer bundles** (project root; each has a matching `*_README.md`, also tracked in the repo's `docs/`):

- `cygl1_fixedop_2020_2022_stats.tar.gz`: pkl/nc4 O-F stats for the whole of 2020–2022 (4 arms), per-year nc4, the 2020-only arms and
  summary tables. It supersedes `cygl1_fixedop_2020_7arm_stats.tar.gz`.
- `cygl1_ismn_fixedop_2020_2022_inputs.tar.gz`: inputs for the ISMN figures notebook when running it locally.

**CYGNSS L1 coverage limit:** the preprocessed L1 stream only contains specular points within 200 km of Arizona (the
CYGNSS_operator preprocessing region, `REGION=arizona`, `REGION_BUFFER_KM=200`) and south of about 37.9° N (the CYGNSS orbit).
It stops at about 106.8° W, so 238 of the 909 domain tiles (28%) and about 90 tiles north of 37.4° N get no L1 obs. That includes
the Jornada and CO-Z1 SOILSCAPE sites. Re-preprocessing with a larger region was considered on 2026-09-28 and **not done**.

### The operator bug (fixed 2026-09-25)

All AZ experiments up to 2026-09-25 were computed with a **defective L1 operator**
and have been archived (configs + stats only, see `archive/`).

**The bug:** `cygnss_preproc_find_obs` (`cygnss_preprocessed_obs.F90`) took each
tile's support coefficients and incidence angle from the obs closest to the
tile among *all* obs in the day file(s), with no time-window or nodata filter.
The reader (`read_obs_cygnss_l1_scalar`) took the observed value from the
closest obs *inside* the 3-h window. So O and F could come from different obs
(different time, footprint and incidence angle): ~30% of tile-cycles on the
full `CYGNSS_L1` stream in a one-week test, fewer on thinned streams. It affected
every L1 O-F statistic, every run that assimilated L1, the paired runs' L1
"% vs OL" numbers (OL read the full stream, the DA arms read the thinned one),
and the L1 z-score scaling climatology.

**The fix:** `GEOSldas_GridComp` commits `cdfaaf4` + `11dfdb1` on
`feature/amfox/cygnss-ascat-hsaf-v8`. The operator now looks up exactly the
obs the reader kept (owner tile + J2000 time stamp from a shared
`cygnss_l1_obs_J2000()`, with an exact sp lon/lat match), and aborts on any mismatch.

**What is still good:** the obs (`CYGNSS_L1/`, `CYGNSS_L1_thinned_dense075_coh05/`)
and the non-L1 scaling parameters (`scaling_params/`).

## Build

Personal checkout `/gpfsm/dnb34/amfox/GEOSldas_cygnss_operator/GEOSldas/`
(same as `/discover/nobackup/amfox/GEOSldas_cygnss_operator/GEOSldas/`).
The top level is on `feature/amfox/cygnss-ascat-hsaf-v8-regtest`, with `origin/develop` (v21.0.0+) merged in,
GMAO_Shared v3.0.2, GEOS_Util v3.0.3, MAPL v2.70.0, and GEOSgcm_GridComp at
`develop`. It was a clean full rebuild on 2026-09-25 (`parallel_build.csh -account s3208 -q debug`,
~12 min), and `install/bin/GEOSldas.x` has md5 `6a70b4968b2b376b601636f780a0da2f`.
The previous build/install are kept as `build_pre_20260925/`, `install_pre_20260925/`.

**Must** use this checkout's `install/bin/ldas_setup`. The shared
`land_da/GEOSldas_develop` install doesn't have the CYGNSS L1 operator. After
`ldas_setup`, check that `<EXP_ID>/build` points at
`/gpfsm/dnb34/amfox/GEOSldas_cygnss_operator/GEOSldas/install`. The old
`cygl1_pregate_build` relink is obsolete: the hard gate is retired and the default
install is ungated.

## Layout

| Path | What |
|---|---|
| `CYGNSS_L1/` | Full CYGNSS L1 coefficient-product obs (schema 0.5), 2020-01-01 to 2022-12-31 (`old_Y2019`/`Y2019`: older vintages). |
| `CYGNSS_L1_coh040_216/` | Full stream screened to coherency_ratio 0.40–2.16, 2020–2022. The obs for the two `coh040216` arms. |
| `CYGNSS_L1_thinned_dense075_coh05/` | Thinned subset (min separation 0.75 deg on top of the intermediate tier, coherency_ratio >= 0.5), 2020-01-01 to 2021-12-31. The obs used by the paired experiments. |
| `scaling_params/` | z-score scaling params, still valid: `z_score_clim/` (SMOS/SMAP Tb, species 21-24, 31-34), `python_z_score_clim_quarter_degree/` (H SAF ASCAT 49-51, CYGNSS L3 54). Copied from the archived `OLv8_M36_all_sensors_AZ` (2020-2022 climatology). Post-fix L1 clims: `cygnss_l1_z_score_clim/` (full stream) and `cygnss_l1_z_score_clim_coh040216/` (screened stream), both built from `OLv8_M36_AZ_fixedop` for 2020–2022 (W=75 d, Nmin=20, short `AZ_CYGNSS_L1_zscore_all_pentads.nc4` symlink). |
| `templates/` | exeinp + `LDASsa_SPECIAL_inputs_ens{upd,prop}.nml` to start the new experiments from (see below). |
| `<EXP_ID>/`, `<EXP_ID>.txt` | The fixed-operator experiments (table above) and their exeinp files. |
| `archive/` | Tarballs of everything retired (see below). |
| `bat_inp_*.txt` | Batch inputs (account s3208, 24 tasks): `debug_limdom` (debug qos), `limdom_1mo_monitor` (1:30 h, monitor-only months), `limdom_1mo_cygl1assim` (2 h, L1-assim months), `limdom_long` (8 h). |

The L1 scaling clims were regenerated post-fix on 2026-09-26 (see `scaling_params/` above); the pre-fix clim is only in
`archive/OLv8_M36_all_sensors_AZ_configs_stats.tar.gz`.

## Templates

`NML_INPUT_PATH` points at the template dir itself and `scalepath` at
`scaling_params/`. `EXP_ID`/`BEG_DATE`/`END_DATE` still carry the values of the
archived run each one came from, so copy the template and set those before `ldas_setup`
(absolute exeinp path; every exeinp needs `DO_ISSM: 0`).

| Template | From | Purpose |
|---|---|---|
| `OL_monitor_unscaled/` | `OLv8_M36_all_sensors_AZ` | Monitor-only, unscaled, all 13 species on the full `CYGNSS_L1` stream: the run to regenerate the L1 z-score clim from. |
| `OL_monitor_paired/` | `OLv8_M36_AZ_paired_monitor_dense075obs` (never run) | Optional, reference only: a scaled OL on the thinned obs (clone of the L1 arm with species 56 `assim=.false.`). Not needed for scoring (see below). |
| `DA_L1_dense075_coh05/` | `DAv8_M36_AZ_paired_cygl1_dense075_coh05` | CYGNSS L1 assim (thinned dense075_coh05 obs, errstd 2.75 dB, xcorr/ycorr 0.625 deg). |
| `DA_L3/` | `DAv8_M36_AZ_paired_cygl1_dense075_coh05_L3assim` | CYGNSS L3 (`CYGNSS_SM_6hr`) assim, L1 monitor-only. |
| `DA_SMAP/` | `DAv8_M36_all_sensors_AZ_scaled_smapassim_6mo` | SMAP L1C Tb (31-34) assim, L1 monitor-only. |
| `DA_L1_full*/` | new, 2026-09-26 | Full-stream L1 arms used by the fixed-operator experiments (xcorr / errstd / coherency variants; see the table at the top). |
| `DA_L3_fullL1mon/`, `DA_SMAP_fullL1mon/` | `DA_L3/`, `DA_SMAP/` | The same as the originals but monitoring the full L1 stream, so L1 O-F is comparable across arms. These were used by `DA_L3_fixedop`/`DA_SMAP_fixedop`. |

**Scoring rule:** compare DA against the single unscaled OL (`OLv8_M36_AZ_fixedop`)
using the DA run's own (scaled) obs, cross-masked: match the two runs' ObsFcstAna on
species/tile/time/sp lon-lat, take O from the DA run, and compute O-F for both runs on that
common set. With the fixed operator, F for each obs depends only on that obs and the model state,
so no separate scaled OL or "OL on the DA's obs files" is needed.

## Archive (`archive/`)

- `*_configs_stats.tar.gz` (2026-09-25): `run/`, `scratch/`, `input/`, `rc_out/`
  (minus regenerable BCS copies), `stats/`, `da_performance_metrics/`, plus the
  experiment's exeinp and nml dir, for `DAv8_M36_AZ_paired_cygl1_dense075_coh05`,
  `..._L3assim`, `DAv8_M36_all_sensors_AZ_scaled_smapassim_6mo`,
  `OLv8_M36_all_sensors_AZ_scaled` and `OLv8_M36_all_sensors_AZ` (the latter includes the invalid
  `cygnss_l1_z_score_clim`). Model output (`ana/`, `cat/`, `rs/`) deleted.
- `old_configs_20260925.tar.gz`: exeinp files and nml dirs for all earlier
  experiments (LS_* tests, R sweeps, single-obs, thinning/coherency families),
  the old shared top-level nml, early coefficient test files, and the unused batch inputs.
- `root_analysis_files_20260925.tar.gz`: loose analysis scripts, npz, maps
  and notes from the project root (tile quality/SNR maps, coherence, R
  diagnosis, coherency spec, L3/SMAP 24-month README).
- Older tarballs (2026-09-23 and earlier): previously archived experiments and
  stats bundles. All of them are pre-fix, so none of their L1 results is valid.

Analysis scripts and write-ups live in the separate git repo
`/gpfsm/dnb06/projects/p284/geosldas-analysis` (`projects/CYGNSS_L1_AZ/`). Documents dated before 2026-09-25 there
(e.g. `runs/cygl1_dense075_coh05_ismn_validation.md`, `runs/cygl1_coh05_density_spectrum_and_24mo_result.md`) are pre-fix, so their
L1 results are invalid. The post-fix documents are the ones listed under Status.

## Known gotchas

- `ldas_setup` positional args must be absolute paths, or it fails with a misleading
  coupled-ADAS assertion error.
- `ldas_setup` bug in this checkout's `ldas.py` (~l.355-395): a restart source that is
  itself a domain-cropped experiment trips `.til.domain` handling. Restart from
  the global `LS_OLv8_M36_v2` instead.
- MetOp-A (`ASCAT_HSAF_META_SM`) ends 2021-11-15 (clean early return, not an abort).
- `END_DATE=YYYYMMDD` needs obs through that date. `CYGNSS_L1` has daily files
  through 2022-12-31, and the thinned set through 2021-12-31.
- To extend a finished run, edit `run/CAP.rc` END_DATE (exactly `END_DATE: 20230101 000000`) and resubmit `lenkf.j`. No new
  `ldas_setup` is needed.
- land_da has an 800k HARD inode limit, and each arm-month adds about 1.4k inodes. Check `showquota` before long extensions. The
  verified tar packer is `/discover/nobackup/projects/land_da/util/inode_pack/` (`pack_verify.sh`).
- Scoring jobs (`score_cygl1_arm.py`) take about 15–20 s per month per arm, and the ISMN job about 22 min. Request about 1 h of
  walltime, not 8 h: oversized requests sit in the Priority queue for hours.
- HISTORY.rc increment collections (`catch_progn_incr`, `inst3_1d_lndfcstana_Nt`)
  aren't on by default. Enable them before the first segment if you need them.
