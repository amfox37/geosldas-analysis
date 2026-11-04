# VIIRS snow cover obs study: plan

Status: **plan confirmed by amfox 2026-11-04 (periods and decision rule); QC mirror done and checked; downloads started.** Started 2026-11-04 (amfox). Code: `scripts/viirs_obs_study/`, `lib/scf_qc.py`.

This study supports two decisions for the GEOSldas VIIRS snow cover fraction (SCF) reader
(`read_obs_JPSS_SCF` / `read_JPSS_SCF_hdf`, GEOSldas_GridComp branch `feature/amfox/jpss_hdf5_reader`, species 59
VJ110C1 = JPSS-1, 60 VJ210C1 = JPSS-2, 61 VNP10C1 = S-NPP). Background: `viirs_scf_reader_test_report.md`.

## Questions

**Q1. Keep or drop VIIRS `Basic_QA` class 1?** The reader keeps classes 0 ("good"), 1 ("poor") and 3 ("other", low
sun) and drops 2 ("bad"). Class 3 must be kept: it is most of the winter obs north of 45°N. Class 1 is a small share
but looked biased on one day (15 Jan 2020, VJ110C1 vs MODIS Aqua: correlation 0.41, bias +0.26, RMSD 0.43, against
correlation 0.90 and RMSD 0.10 for class 3; ~1% of kept cells, ~2% of snowy cells). Does class 1 add enough coverage to be worth
its error, across seasons and years? **Decision: amfox.**

**Q2. JPSS-1 vs JPSS-2, and is assimilating both better than one?** Both cross at ~13:30 local, ~50 min apart on
the same orbit, so they see nearly the same scene with different swath gaps and clouds. (a) How well do they agree,
and is either closer to MODIS? (b) How much coverage does the second satellite add? (c) Does assimilating both
improve the snow analysis over one? (a) and (b) are obs-only (this study, on the Mac); (c) needs GEOSldas runs on
Discover (Part C).

Note on (c): when two SCF species are assimilated in the same window, GEOSldas averages the obs per tile before the
rule-based snow update (`clsm_ensupd_upd_routines.F90`, update_type 11: "average in case there are multiple asnow
obs"). The update does not weight obs by their number, so a second satellite adds coverage and some noise
averaging, not weight.

## Data

| Product | Satellite | Overpass (local) | Available | Role |
|---|---|---|---|---|
| VJ110C1 v2 | JPSS-1 / NOAA-20 | 13:30 | 2018-01-05 → | studied |
| VJ210C1 v2 | JPSS-2 / NOAA-21 | 13:30 | 2023-02-10 → | studied |
| VNP10C1 v2 | S-NPP | 13:30 | 2012-01-19 → 2026 | studied (Q1) |
| MYD10C1 v6.1 | Aqua | 13:30 | | main reference (same overpass time) |
| MOD10C1 v6.1 | Terra | 10:30 | | second reference |
| IMS 4 km daily | (multi-sensor, analyst) | daily | | independent snow / no-snow reference |

All SCF products are daily on the 0.05° CMG grid. Download: `scripts/viirs_scf_reader_test/get_snow_cmg_obs.py`
(NSIDC via CMR, ~/.netrc). VIIRS files are ~7 MB/day each. IMS: from `projects/IMS/scripts/download_ims_ascii_to_nc.py`,
plus a new nearest-neighbour regrid to the CMG grid.

Data go under `data/` (git-ignored), not in the repo.

## QC used in the study

`lib/scf_qc.py` is a Python mirror of the two readers' QC (thresholds, SCF = Snow_Cover/Clear_Index in real*4, cap
at 1 for VIIRS). It is checked against the Fortran readers themselves (`check_qc_vs_fortran.py`: same cells, same
order, bit-identical SCF; last run GridComp `6d0ed32`, 9 files, all OK). For Q1, VIIRS cells are kept by all reader
tests except `Basic_QA`, and the QA class is kept as a label (`viirs_qc_without_qa`).

## Periods (confirmed 2026-11-04)

Seven consecutive days per period, chosen to cover accumulation, mid-winter, melt and the Southern Hemisphere winter,
and two years. 2024 has all three VIIRS satellites, so it serves Q1 and Q2 with the same downloads.

| Period | Dates | Why |
|---|---|---|
| P1 | 2020-01-13 .. 2020-01-19 | mid-winter; overlaps the CF0360 reader-test runs (Jan 2020) |
| P2 | 2020-04-11 .. 2020-04-17 | spring melt |
| P3 | 2024-01-13 .. 2024-01-19 | mid-winter, J1 + J2 + S-NPP |
| P4 | 2024-03-23 .. 2024-03-29 | early melt |
| P5 | 2024-05-18 .. 2024-05-24 | late melt (high latitudes, mountains) |
| P6 | 2024-07-13 .. 2024-07-19 | Southern Hemisphere winter (Andes, Australian Alps, New Zealand) |
| P7 | 2024-11-09 .. 2024-11-15 | early accumulation |

About 49 days × 5 products ≈ 1.5 GB, plus IMS.

## Part A — Q1, QA class 1 (Mac)

1. **Coverage.** Per day and period: cells kept with QA {0,3} vs {0,1,3}, by latitude band (0–30, 30–45, 45–60,
   60–90 N and S) and SCF class (0, 0–1, 1). Share of class 1 among kept and among snowy cells. **Tile coverage:**
   map cells to CF0360 and M36 tiles as GEOSldas does (nearest tile center, `ldas_tilecoord.bin`), and count tiles
   that have obs only through class-1 cells. Monthly maps of where class 1 occurs.
2. **Agreement with MODIS Aqua** (same day and cell, both passing their reader QC), per VIIRS QA class, latitude band,
   month and SCF class: N, bias, RMSD, correlation; snow / no-snow agreement at the thresholds of the snow update
   (`SCF_ANA_ALPHA`, `SCF_ANA_BETA`). Repeat with MODIS Terra.
3. **IMS check** for cells where VIIRS class 1 and MODIS disagree: which one IMS supports.
4. **Effect on the analysis.** Tile-mean SCF with and without class 1; how often dropping class 1 changes the snow
   update between add / remove / no change, using the forecast `asnow` of the January 2020 CF0360 ObsFcstAna files
   (Discover, `M21C_test_CF0360_jpss_*`; copy a few days).

**Decision rule** (confirmed 2026-11-04, before any results were looked at): drop class 1 if its RMSD vs MODIS Aqua is more
than 1.5× that of class 3 in most periods, unless dropping it removes obs from more than 2% of snowy tiles.

## Part B — Q2 (a, b), J1 vs J2 (Mac)

Periods P3–P7 (J2 data start 2023-02-10).

1. **Agreement:** J1 vs J2, same day and cell, by QA class, latitude band and SCF class: bias, RMSD, correlation;
   both vs MODIS Aqua.
2. **Added coverage:** cells and tiles (CF0360, M36) with obs from J1 only, J2 only, both. Per tile and week: days
   with at least one obs, J1 alone vs J1 or J2.
3. **Effect on the snow update:** tiles where the J1+J2 tile mean falls on a different side of the add / remove
   thresholds than J1 alone.

If B2 shows that J2 adds few tile-days over J1, Part C matters less.

## Part C — Q2 (c), DA impact (Discover; brief to be written after A and B)

CF0360 (M21C), update_type 11 (SCF only), 1 member as in the reader tests, a full snow season after J2 started (e.g.
2023-10-01 .. 2024-06-30). Runs: OL, DA-J1, DA-J2, DA-J1+J2, optional DA-MODIS-Aqua. Evaluate against obs that were
not assimilated: MODIS O–F, IMS snow / no-snow skill (`projects/IMS`, needs a CF0360 version), SNOTEL and GHCN snow
depth (`projects/SNOTEL`, `projects/GHCN_snwd`), plus SWE increments and number of updates. Report by region and
season, with day-to-day spread.

**Order:** (1) Discover test of reader commit `6d0ed32` passes (`GEOSldas/discover_brief_granule_id_logging.md`);
(2) Q1 decided — the QA choice changes the obs, so decide it before the season-long runs; (3) Parts A and B;
(4) Part C.

## Log

| Date | What |
|---|---|
| 2026-11-04 | Periods and decision rule confirmed (amfox); obs downloads started. |
| 2026-11-04 | Plan written. `lib/scf_qc.py` + unit tests; `check_qc_vs_fortran.py` against GridComp `6d0ed32`: 9/9 files identical. |
