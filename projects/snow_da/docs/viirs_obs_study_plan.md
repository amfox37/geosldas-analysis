# VIIRS snow cover obs study: plan

Status: **plan confirmed by amfox 2026-11-04 (periods and decision rule); QC mirror done and checked; obs downloaded; Parts A1 and A2 done (results: `viirs_obs_study_report.md`).** Started 2026-11-04 (amfox). Code: `scripts/viirs_obs_study/`, `lib/scf_qc.py`.

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

## Obs-to-tile mapping used for tile coverage

As in GEOSldas (`get_tile_num_for_obs` → `get_tile_num_from_latlon`, `LDAS_TileCoordRoutines.F90`; called from
`read_obs_*_SCF` with the grid built in `clsm_ensupd_enkf_update.F90`): the search grid is the perturbation grid
`pert_grid_f` (`get_pert_grid`, `LDAS_PertRoutines.F90`).

- **EASEv2 M36:** the pert grid is the M36 grid itself; an obs goes to the single tile of the M36 cell that contains
  it (none if the cell has no land tile).
- **CF0360:** the pert grid is a global lat/lon grid of 4·360 × 3·360 = 1440 × 1080 cells ("DE001440xPE001080",
  0.25° × 1/6°). Each tile is listed in the cell that contains its center (`pert_i_indg`, `pert_j_indg`); an obs goes
  to the tile of its cell with the smallest |Δlon| + |Δlat| to the tile center.
- In both cases the obs is dropped if it lies outside the chosen tile's bounding box (`min/max_lon/lat`), because
  the SCF species have FOV = 0.

Tile coords and grids come from the runs' `rc_out/*.ldas_tilecoord.bin` and `*.ldas_tilegrids.bin`. Implemented in
`lib/obs2tile.py`; checks and the CF0360 cell-loss finding in the report ("Obs-to-tile mapping").

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
than 1.5× that of class 3 in most periods, unless dropping it removes obs from more than 2% of snowy tiles. The 2% is
per period, i.e. snowy tile-days pooled over the period's 7 days (amfox, 2026-11-05).

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
| 2026-11-04 | **Obs-to-tile mapping confirmed against GEOSldas output.** CF0360 ObsFcstAna of `M21C_test_CF0360_jpss_5d4d636_OL` (2020-01-09..10, MYD10C1, MOD10C1, VJ110C1, VNP10C1) and M36 reproduced exactly (all but 1 obs per product). Needed: the runs' own tilecoord (copied from run A2; the local `M21C_test_CF0360` file has other tile boxes), land tiles only, and real*4 cell centers with fused multiply-add. CF0360 cell-loss finding confirmed: 54.5% of land CMG cells get a tile, 17.5% of land tiles never get an obs. A1 rerun with the runs' tilecoord: results unchanged to ≤0.02 points. |
| 2026-11-04 | **A1 coverage done** (`coverage.py`, `coverage_report.py`; report section A1). Class 1 = 0.24–1.1% of kept cells, 0.8–2.8% of snowy cells. Snowy tile-days with obs only via class 1: CF0360 0.9–2.0%, M36 0.5–1.0% per period (max single day 2.45%, CF0360). Coverage leg of the decision rule not met in any period (pooled over 7 days); 2% read per period (confirmed by amfox 2026-11-05). |
| 2026-11-05 | **A2 agreement with MODIS done** (`agreement.py`, `agreement_report.py`; report section A2). Class-1 RMSD vs MODIS Aqua 0.34–0.44, 1.5–29× class 3 in 19/19 period-satellites (snowy subset 15/19), 4.6–8× class 0; Terra the same. Bias +0.15..+0.25: class 1 sees snow MODIS does not (30% on the snowy side of 0.55 where MODIS is not, 1% the reverse). Only 11–25% of class-1 cells have a usable MODIS cell. **Both legs of the decision rule point to dropping class 1**; amfox decides. Next: A3 (IMS check of the disagreements). |
| 2026-11-05 | amfox: the 2% coverage threshold is per period (7 days pooled), not per day. With it, the coverage leg is not met in any period (max 1.96%, S-NPP P1, CF0360), so the rule's outcome is to drop class 1. |
| 2026-11-04 | `lib/obs2tile.py` (CMG cell → tile as GEOSldas) + `check_tile_map.py`. M36 checked end to end against `LS_DAv8_M36_as_test3` ObsFcstAna (MYD10C1, 2020-05-06): 17,783/17,827 tile obs reproduced. Finding: for CF0360 with FOV = 0 only 56% of land CMG cells get a tile (42% north of 60°N); not yet checked against CF0360 ObsFcstAna. |
| 2026-11-04 | Obs downloaded (`get_obs.py`, 1.2 GB): 228 of 231 files; manifest in `data/viirs_obs_study/cmg/manifest.txt`. Missing: VNP10C1 2024-07-13..15 (P6) — not in the NSIDC catalog (CMR); VNP10C1 2024-07-16 is about half the usual size (likely a partial day). |
| 2026-11-04 | Periods and decision rule confirmed (amfox); obs downloads started. |
| 2026-11-04 | Plan written. `lib/scf_qc.py` + unit tests; `check_qc_vs_fortran.py` against GridComp `6d0ed32`: 9/9 files identical. |
