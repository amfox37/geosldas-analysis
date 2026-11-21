# MODIS-to-VIIRS transition and Terra vs Terra+Aqua: analysis plan

Status: **draft, 2026-11-21 (amfox + Claude). Nothing below has been run.** Items marked **[CHECK]** are not verified.
Runs: `/discover/nobackup/projects/land_da/VIIRS_testing/main_run/` (README there has the run history).
Related: `viirs_obs_study_plan.md` (reader QC, QA class 1, J1 vs J2 obs agreement). This plan is its Part C,
on M36 and for three years instead of one CF0360 snow season.

## Questions

1. **Transition.** Does switching the assimilated snow cover from MODIS to VIIRS change (a) snow DA skill and
   (b) the water balance (SWE, melt, runoff, ET, storage)?
2. **Constellation.** Does assimilating Terra only instead of Terra+Aqua (MODIS), or J2 only instead of J1+J2
   (VIIRS), change the same things? This is the "could we drop to one satellite" question.
3. **Update structure.** MODIS Terra (10:30 local) and Aqua (13:30 local) fall in different 3-hourly analysis
   windows, so a tile can receive two increments per day. J1 and J2 (both 13:30) fall in the same window and are
   averaged into one. How much does this matter for the number, size and sign of increments?

## Runs (single member, snow only, M36, update_type 11)

All five arms start 2020-08-01 from the same restart (`restarts/`, `RESTART_ID LS_DAv8_M36`) and use the same
forcing and the same build (`install-viirs-f87aa3f`, GridComp `f87aa3f`). Everything except the assimilated species
is identical. All four SCF species have `getinnov = .true.` in every segment, so every arm has innovations for
every sensor; only `assim` differs.

| Arm | 2020-08-01 → 2023-08-01 (namelist) | 2023-08-01 → 2026-08-01 (namelist) |
|---|---|---|
| `viirs_OL` | none (`nml_OLa`) | none |
| `viirs_T_to_J2` | Terra (`nml_E1a`) | VIIRS J2 only (`nml_E1b`) |
| `viirs_TA_to_J1J2` | Terra + Aqua (`nml_E2a`) | VIIRS J1 + J2 (`nml_E2b`) |
| `viirs_T_control` | (is `viirs_T_to_J2`) | Terra (`nml_E1a_cont`) |
| `viirs_TA_control` | (is `viirs_TA_to_J1J2`) | Terra + Aqua (`nml_E2a_cont`) |

Species: 52 Aqua `MYD10C1`, 53 Terra `MOD10C1`, 59 J1 `VJ110C1`, 60 J2 `VJ210C1`; 61 S-NPP `VNP10C1` is not
assimilated in any arm (confirmed from each run's namelist, 2026-11-21). The `_cont` namelists differ from
`E1a`/`E2a` only in the MODIS `path` (merged `MODIS_obs/`, which adds the 2026 files that `gmao_ops` lacks).

**Continuous chains.** `T_to_J2` (pre) + `T_control` (post) is one uninterrupted Terra-only MODIS run, and
`TA_to_J1J2` + `TA_control` one uninterrupted Terra+Aqua MODIS run, because the controls restart from the pre-switch
arms' own 2023-08-01 restarts. They are the MODIS reference chains for the whole 6 years.

**Analysis years.** Use "snow years" 1 Aug – 31 Jul. They line up with the start (2020-08-01), the switch
(2023-08-01) and the end (2026-08-01): SY1–SY3 (2020-21 to 2022-23) before the switch, SY4–SY6 (2023-24 to
2025-26) after. Southern Hemisphere snow is minor; report NH results.

## Contrasts

| ID | Arms | Period | Same start state? | What it isolates |
|---|---|---|---|---|
| C1 | `T_control` vs `T_to_J2` | SY4–SY6 | **yes** (2023-08-01) | Terra vs J2 (sensor + overpass time) |
| C2 | `TA_control` vs `TA_to_J1J2` | SY4–SY6 | **yes** | Terra+Aqua vs J1+J2 |
| C3 | `T_to_J2→T_control` vs `TA_to_J1J2→TA_control` | SY1–SY6 | **yes** (2020-08-01) | Terra vs Terra+Aqua, MODIS, 6 years |
| C4 | `T_to_J2` vs `TA_to_J1J2` | SY4–SY6 | no (3 years of different history) | J2 vs J1+J2; treat as coverage comparison |
| C0 | each arm vs `viirs_OL` | SY1–SY6 | yes | DA impact (DA − OL) |

C1, C2 and C3 are the cleanest: observations are the only difference. Pre-switch years of C1/C2 are not a contrast
(the arms are identical there, since `T_control`/`TA_control` start at the switch); use them instead as the
before-switch baseline of each chain.

Confounds to state with every result: Terra is a morning overpass and J1, J2 and Aqua are afternoon overpasses
(the readers set 10:30 and 13:30 local), so T → J2 changes sensor **and** time of day; and with one member there is no
ensemble spread.

## Analyses

Keep to five. Headline outcomes (chosen now, before looking): **H1** fraction of tile-days with two increments;
**H2** net-to-absolute SWE increment ratio; **H3** DA − OL skill against IMS snow cover; **H4** DA − OL snow-year
runoff and end-of-melt SWE. Everything else supports these.

### A. Observation availability
Per arm and day, from `ana/ens_avg/*ldas_ObsFcstAna*` (3-hourly; has `species`, `assim_flag`, `tilenum`, `obs`):
obs per species, tiles with ≥1 obs, tiles with obs from both satellites, by latitude band and month. Compare Terra
vs Aqua overlap and J1 vs J2 overlap. This is the likely explanation for any T vs T+A difference, so do it first.

### B. Update structure and increments (Q3)
From `cat/ens0000/Y*/M*/<exp_id>.catch_progn_incr.<yyyymmdd>.nc4` (one file per day, 8 instantaneous 3-hourly steps;
single member, so there is no `ens_avg`) per tile-day:
- number of analysis steps with a non-zero snow increment (H1); fraction with two, by arm and latitude/longitude
  (windows are set by local time, so the UTC pattern varies with longitude);
- for double days: size of first vs second increment, sign agreement (reinforce vs undo);
- annual net and absolute SWE increments and net/abs (H2); fraction of increments at the cap
  `SCF_ANA_MAXINCRSWE = 5 kg m-2` per update (alpha 0.60, beta 0.55);
- from `getinnov`: innovation of the second sensor in the same arm vs in the OL (how much the first update already
  moved the state toward the second obs); the same for J1 vs J2.
Mechanism (read in the code, both checkouts): obs in the same update window are averaged per tile
(`clsm_ensupd_upd_routines.F90`, "average in case there are multiple asnow obs"); the increment depends on the
forecast SCF, so a second update acts on the already-updated state.

### C. Independent snow skill (Q1a, Q2)
Against IMS snow cover (`projects/IMS`), SNOTEL SWE (`projects/SNOTEL`) and GHCN snow depth (`projects/GHCN_snwd`):
skill of each arm and of the OL; report DA − OL, then arm contrasts C1–C4. Use the difference between post- and
pre-switch skill within each chain as a baseline (a "difference of differences"). NH seasonal-snow mask as in the
M21C snow analysis. Obs are not comparable across sensors, so do not score an arm against the sensor it assimilates
as a primary result; scoring each arm against the sensors it does not assimilate is a useful secondary check.
**[CHECK]** IMS, SNOTEL and GHCN coverage through 2026-07 in the existing project scripts.

### D. Water balance (Q1b, Q2)
Reuse the M21C snow-DA budget definitions (`projects/M21C_ls`, water-year differential budget, including the PEATCLSM
free-standing-water term; do not re-derive them). Per snow year and arm: DA − OL for SWE (peak and end of melt),
snowmelt, runoff, ET, `TWLAND`, RZMC; the snow-DA input as the summed `WESNN*_INCR`; and the partition of added snow
water. Then A − B maps and domain means for C1–C4.
**[CHECK]** whether these runs write `PEATCLSM_FSWCHANGE` and enough restarts for exact endpoints (restarts are
monthly under `rs/`; the M21C budget used 1 October endpoints, here 1 August is natural).

### E. Transition behaviour
The first 4–8 weeks after 2023-08-01 in C1/C2: step in SWE, SCF and increments (spin-up of the new sensor). Decide
whether to exclude a short window from the contrasts.

## How results will be read

- With one member and 3 snow years per contrast, use paired differences and 5° spatial-block bootstrap (as in M21C)
  over tiles, and report year-to-year spread explicitly. Where a difference is not distinguishable from zero, say
  "no resolvable difference" and do not interpret maps.
- C3 (6 years, same start) is the most reliable estimate of what Aqua adds; C4 is only a consistency check on it.
- A difference between MODIS and VIIRS arms is a statement about the whole configuration (sensor, overpass time,
  QC, number of updates per day), not about the sensors alone. Part of the aim of B is to separate the update
  structure from the rest.

## Code and data facts relied on (verified 2026-11-21)

- Build `install-viirs-f87aa3f`, all five arms; GridComp `HEAD` = `f87aa3f`, no uncommitted source changes.
- VIIRS reader keeps `Basic_QA` 0 and 3 only (class 1 dropped; matches Part A of the obs study).
- Window-boundary double assimilation fixes are in: `8c778e4` (VIIRS), `86afd34` + `8fbcc61` (MODIS, PR #240).
  The ±45°/±135° striping artefact should therefore not appear.
- Overpass hours in the readers: Terra 10.5, Aqua 13.5, J1 13.5, J2 13.5 (local).
- Analysis step is 3 h, inferred from the 3-hourly `catch_progn_incr` / `ObsFcstAna` output; `LANDASSIM_DT` itself was
  not found in the run files **[CHECK]**.
- Output volume: `ObsFcstAna` is ~160 MB and ~250 files per arm-month; `catch_progn_incr` is ~34 MB per daily file
  (~1 GB and 31 files per arm-month). `land_da` was at ~91% of its file quota on
  2026-11-19, so extract compact per-day or per-month summaries rather than copying or unpacking files, and check
  `showquota` (inodes) before large jobs.

## Planned layout

`config/arms.json` (arms, namelists, species, dates, contrasts, paths: scripts read it instead of hard-coding),
`scripts/viirs_transition/` (one script per analysis A–E, each with a `_discover.sh`/`.j` wrapper), outputs under
`output/` (not versioned), results in `viirs_transition_report.md`.

## Open items

Resolved 2026-11-21:

1. ~~Species and `assim` flags per arm.~~ Confirmed from each arm's current `run/LDASsa_SPECIAL_inputs_ensupd.nml`
   (as in the runs table); species 61 (S-NPP) has `assim = .false.` in every arm.
2. ~~Final run state.~~ All five arms have a restart at `20260801_0000` and the queue is empty, so every chain reached
   `END_DATE 20260801`.
3. ~~Does `ObsFcstAna` carry forecasts for getinnov-only species?~~ Yes. `ana/ens_avg/*ldas_ObsFcstAna*` lists all four
   species with `assim_flag` and finite `fcst` and `ana` for every obs (checked 2024-01-15 12z for
   `TA_to_J1J2`: 2741 / 510 non-assimilated Aqua / Terra obs and 3964 / 3765 assimilated J1 / J2 obs).
   **Species numbering inside these files is local, not the namelist index:** 1 = `MYD10C1` (Aqua), 2 = `MOD10C1`
   (Terra), 3 = `VJ110C1` (J1), 4 = `VJ210C1` (J2); `assim_flag` differs by arm and period. Scripts should map
   through `obsparam_descr` (or the table in `config/arms.json`), not the namelist indices 52/53/59/60. Because
   `assim_flag` is the right filter for "assimilated", `fcst` and `ana` of non-assimilated obs are usable for the
   cross-sensor checks in analysis B and C.

Still open:

1. `LANDASSIM_DT` / `LANDASSIM_T0`: no uncommented setting found in the exeinp or `run/` files. The 3 h step is
   inferred from the 3-hourly `ObsFcstAna` file names and the `catch_progn_incr` time axis (`time_increment` 03:00:00),
   which HISTORY.rc says must match `LANDASSIM_DT`. Confirm from the resolved resource file or a job log.
2. IMS / SNOTEL / GHCN data coverage to 2026-07 (analysis C).
3. Availability of `PEATCLSM_FSWCHANGE` in the daily output and of restarts at the budget endpoints (analysis D).
   Restarts exist at least at 2023-08-01 and 2026-08-01; check that monthly restarts are complete.

## Log

| Date | What |
|---|---|
| 2026-11-21 | Plan drafted from the run configuration and the GEOSldas source (reader overpass hours, per-tile obs averaging, SCF update parameters). |
| 2026-11-21 | `config/arms.json` written. Open items updated: species/assim flags confirmed, all five arms reached `END_DATE 20260801`, `ObsFcstAna` carries fcst/ana for non-assimilated species (local species ids 1-4). Increment path corrected to `cat/ens0000/` (daily files, 8 steps). |
