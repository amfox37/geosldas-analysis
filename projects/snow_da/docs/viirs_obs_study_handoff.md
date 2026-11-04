# VIIRS obs study: handoff (2026-11-04)

For whoever continues this study (human or AI session started in geosldas-analysis). Read
`viirs_obs_study_plan.md` first: it has the questions, data, periods, methods, decision rule and log. This note
adds the current state, the next steps, and context that is not in the plan.

## State

| Item | Status |
|---|---|
| Plan, periods (P1–P7), decision rule for QA class 1 | confirmed by amfox 2026-11-04 |
| `lib/scf_qc.py` (Python mirror of the GEOSldas VIIRS + MODIS SCF reader QC) | done; 7 unit tests pass; bit-identical to the Fortran readers on 9 files (GridComp `6d0ed32`) |
| Obs download (`scripts/viirs_obs_study/get_obs.py`) | done: 228/231 files, 1.2 GB in `data/viirs_obs_study/cmg/` (git-ignored); `manifest.txt` there. Missing VNP10C1 2024-07-13..15 (not at NSIDC) |
| Obs-to-tile mapping | worked out from the GEOSldas code; documented in the plan ("Obs-to-tile mapping") |
| Part A (QA class 1), Part B (J1 vs J2) | **not started** |
| Part C (DA runs on Discover) | after A, B and the QA decision |
| Git | commits `f52c48e`, `a1d079b`, `b1c4448` (+ this note) are **local only**; amfox pushes (`git push origin main`) |

## Next steps

1. **Part A1, coverage** (`scripts/viirs_obs_study/coverage.py`, to write): per day and period, VIIRS cells kept
   with `Basic_QA` {0,3} vs {0,1,3}, by latitude band and SCF class; share of class 1 among kept and snowy cells;
   tiles that get obs only through class-1 cells; monthly maps of class-1 cells. Use
   `scf_qc.viirs_qc_without_qa()` (all reader tests except QA) and the QA grid.
   - **M36 tiles:** local files are enough — tile coords and grids from
     `/Users/amfox/Desktop/LS_DAv8_M36_as_test3/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_DAv8_M36_as_test3.ldas_tilecoord.bin`
     and `.ldas_tilegrids.bin`. Readers: `read_tilecoord`, `read_tilegrids` in
     `GEOSldas_GridComp/GEOSldas_App/util/shared/python/read_GEOSldas.py` (little-endian Fortran sequential).
   - **CF0360 tiles:** tile coords are local
     (`/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/M21C_test_CF0360.ldas_tilecoord.bin`);
     the obs-mapping grid is the hard-wired 1440 × 1080 lat/lon pert grid (see plan), so no tilegrids file is
     strictly needed. Ask amfox for `ldas_tilegrids.bin` from a `M21C_test_CF0360_jpss_*` run on Discover if you want
     to confirm.
2. **Part A2–A3, agreement with MODIS** (Aqua main, Terra second) by QA class, latitude band, month, SCF class; snow /
   no-snow agreement at the snow-update thresholds (`SCF_ANA_ALPHA`, `SCF_ANA_BETA`, in the GEOSldas
   `LDASsa_DEFAULT_inputs_ensupd.nml` / snow analysis code); then IMS for cells where VIIRS class 1 and MODIS disagree
   (reuse `projects/IMS/scripts/download_ims_ascii_to_nc.py`; new nearest regrid IMS → 0.05° CMG).
3. **Part B, J1 vs J2** on P3–P7 with the same data.
4. Write results into a report `docs/viirs_obs_study_report.md` and add a line to the plan's log for each step.

## Context not in the plan

- **Background:** the VIIRS reader is on GEOSldas_GridComp branch `feature/amfox/jpss_hdf5_reader` (now `6d0ed32`;
  local checkout `/Users/amfox/Desktop/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp`). Species
  59 VJ110C1, 60 VJ210C1, 61 VNP10C1 (CYGNSS L1 will move to 62 later). Reader-test report:
  `viirs_scf_reader_test_report.md`.
- **Decided (amfox):** VIIRS `errstd` stays -9999 (rule-based snow update ignores it); SCF capped at 1 (VIIRS
  Snow_Cover can exceed Clear_Index by 1; dropping would lose ~23% of full-snow cells). **Open:** QA class 1 — the
  question of this study; amfox decides with the rule in the plan.
- **Window fix:** SCF obs of tiles whose center is exactly on a 3-hourly window boundary (multiples of 45°) were
  assimilated twice. Fixed for VIIRS (`8c778e4`) and MODIS (PR #240, tested on Discover M36: 100 tiles on
  2 Jan 2020, all on ±45/±135°). Not relevant for obs-only analysis, but explains stripes at ±45°/±135° in older
  M36 SCF DA runs (e.g. Lauren Andrews's).
- **Snow update averages obs per tile:** with two SCF species in the same window, GEOSldas averages them per tile
  before the rule-based update (relevant for Part B / C).
- **If the reader QC changes:** update `lib/scf_qc.py` and rerun `scripts/viirs_obs_study/check_qc_vs_fortran.py`
  (instructions in that folder's README).

## Conventions

- Python: `/Users/amfox/mamba/envs/regrid/bin/python`. Data under `data/` (git-ignored); never commit data or
  large outputs.
- Script-first; figures from scripts; clear notebook outputs before committing (repo `CLAUDE.md`).
- Commit messages: `snow_da: <imperative summary>`. amfox pushes; ask before any push, `rm`, `mv`, `cp`.
- Document as you go: plan log, folder README, report. amfox asked for everything to be well documented.
