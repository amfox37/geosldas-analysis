# MODIS-to-VIIRS transition and Terra vs Terra+Aqua: analysis plan

Status: **2026-11-21 (amfox + Claude). Analyses A, B, C, D (+D2) and E have all been run (`scripts/viirs_transition/`, figures and tables in `docs/viirs_transition_report_figures/`); synthesis in `viirs_transition_report.md`. The decision rule below was CONFIRMED by amfox on 2026-11-21 (after the results were seen; see there) and the IMS SY5-6 behaviour is open.** Items marked **[CHECK]** are not verified.
Runs: `/discover/nobackup/projects/land_da/VIIRS_testing/main_run/` (README there has the run history).
Related: `viirs_obs_study_plan.md` (reader QC, QA class 1, J1 vs J2 obs agreement). This plan is its Part C,
on M36 and for three years instead of one CF0360 snow season.

## Questions

1. **Transition.** Does switching the assimilated snow cover from MODIS to VIIRS change (a) snow DA skill and
   (b) the water balance (SWE, melt, runoff, ET, storage)?
2. **Constellation.** Does assimilating Terra only instead of Terra+Aqua (MODIS), or J2 only instead of J1+J2
   (VIIRS), change the same things? This is the "could we drop to one satellite" question.
3. **Update structure.** MODIS Terra (10:30 local) and Aqua (13:30 local) fall in different 3-hourly analysis
   windows (always consecutive steps, by construction: the readers use fixed local overpass hours exactly 3 h apart and
   the windows are 3 h long; seen in `ObsFcstAna`), so a tile can receive two increments per day. J1 and J2 (both
   13:30) fall in the same window and are
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
| C4 | `T_to_J2` vs `TA_to_J1J2` | SY4–SY6 | no (3 years of different history) | J2 vs J1+J2, as a difference of differences (below) |
| C0 | each arm vs `viirs_OL` | SY1–SY6 | yes | DA impact (DA − OL) |

C1, C2 and C3 are the cleanest: observations are the only difference. Pre-switch years of C1/C2 are not a contrast
(the arms are identical there, since `T_control`/`TA_control` start at the switch); use them instead as the
before-switch baseline of each chain.

**VIIRS constellation (J2 vs J1+J2).** There is no same-start pair for this (decided 2026-11-21: no sixth arm, to keep
the design simple). Answer it as a difference of differences: the C4 difference in SY4–SY6 minus the same T-vs-T+A
difference in the pre-switch years SY1–SY3 (both arms there are MODIS, from the same start state). This applies to the
water balance (analysis D) as well as to skill (analysis C). It assumes the T-vs-T+A gap would have stayed the same
without the switch; state that assumption with the result, and compare with C3 (MODIS, six years) as a check.

Confounds to state with every result: Terra is a morning overpass and J1, J2 and Aqua are afternoon overpasses
(the readers set 10:30 and 13:30 local), so T → J2 changes sensor **and** time of day; the readers use fixed hours for
all six years, but I believe the real Terra and Aqua equator-crossing times drifted from about 2020–2022 (not verified,
and not visible in the CMG daily files), so the time-of-day difference may have grown over the record; and with one
member there is no ensemble spread.

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
Analysis A found that the VIIRS union covers about 6–7% fewer tiles than MODIS T+A, and that some sensors had days
with no obs at all (see "Code and data facts"). Interpret differences in increment counts and sizes against these
coverage differences, and exclude or flag the no-obs days in every contrast (e.g. `T_control` assimilated nothing on
2025-05-11 to 05-14).
Mechanism (read in the code, both checkouts): obs in the same update window are averaged per tile
(`clsm_ensupd_upd_routines.F90`, "average in case there are multiple asnow obs"); the increment depends on the
forecast SCF, so a second update acts on the already-updated state.

### C. Independent snow skill (Q1a, Q2)
Against IMS snow cover (`projects/IMS`), SNOTEL SWE (`projects/SNOTEL`) and GHCN snow depth (`projects/GHCN_snwd`):
skill of each arm and of the OL; report DA − OL, then arm contrasts C1–C4. Use the difference between post- and
pre-switch skill within each chain as a baseline (a "difference of differences"). NH seasonal-snow mask as in the
M21C snow analysis. Obs are not comparable across sensors, so do not score an arm against the sensor it assimilates
as a primary result; scoring each arm against the sensors it does not assimilate is a useful secondary check.

**IMS is used, but it is not independent.** IMS is a multi-sensor analysis that itself uses MODIS and VIIRS imagery, so
IMS skill can favour arms that assimilate the same imagery. Keep it (it is the only global check of snow cover
fraction, which is what is constrained), but always show it next to SNOTEL SWE and GHCN snow depth, say so in every
IMS figure and table, and do not rest a conclusion on IMS alone. SNOTEL (western US) and GHCN (sparse, mostly
mid-latitude) are the genuinely independent checks; they do not cover most of the snow domain.

Data status (checked 2026-11-21; nothing downloaded yet), all in the repo clone under `land_da`:
SNOTEL `projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet` runs 1999-10-01 to 2025-09-30 (needs extending to
2026-07; downloader has `--end`); GHCN only has validation outputs for 2000 to 2024, the raw build products are not
there (rebuild with `notebooks/ghcn_snwd_global_1998_present_build.ipynb`); IMS raw 24 km files exist for 2024 only
(`projects/IMS/output/gz/`; `scripts/download_ims_ascii_to_nc.py` fetches the rest from NSIDC G02156). The existing
OL/DA metrics are for `LS_OLv8_M36`/`LS_DAv8_M36` (2000 to 2024), not our arms, and the scripts are OL/DA-pair
specific, so they need generalising to N arms (driven by `arms.json`). Store extractions as consolidated
parquet/NetCDF under `land_da`, not as thousands of small files (file quota).

### D. Water balance (Q1b, Q2)
Reuse the M21C snow-DA budget definitions (`projects/M21C_ls`, water-year differential budget, including the PEATCLSM
free-standing-water term; do not re-derive them). Per snow year and arm: DA − OL for SWE (peak and end of melt),
snowmelt, runoff, ET, `TWLAND`, RZMC; the snow-DA input as the summed `WESNN*_INCR`; and the partition of added snow
water. Then A − B maps and domain means for C1–C4.
Checked 2026-11-21: `PEATCLSM_FSWCHANGE` (with `TWLAND`, `RUNSURFLAND`) is in the daily `tavg24_1d_lnd_Nt` files, and
all monthly restarts exist (73 per six-year arm, 37 per control), including every 1 August endpoint. Only the variable
list of one file was checked, not that every day has valid values; the budget closes to ~0.01 kg m-2, which supports it.
Run 2026-11-21: `wb_extract.py`/`wb_analyze.py` (budget, C0-C4) and `wb2_extract.py`/`wb2_analyze.py` (partition
fractions, RZMC/SFMC, significance maps, SY1 and gap-day sensitivities); results `D_*` and `D2_*`.

### E. Transition behaviour
The first 4–8 weeks after 2023-08-01 in C1/C2: step in SWE, SCF and increments (spin-up of the new sensor). Decide
whether to exclude a short window from the contrasts.
Run 2026-11-21 (`tr_*.py`, results `E_*`): no spin-up detected; no exclusion window needed (see the report).

## Results at a glance (2026-11-21; details and caveats in `viirs_transition_report.md`)

| Analysis | Files (`docs/viirs_transition_report_figures/`) | Main result |
|---|---|---|
| A obs availability | `obs_availability_*.md`, `figA_*` | VIIRS union ~6-7% fewer tiles than MODIS T+A; input gap days (below). |
| B increments | `B_summary.md`, `B_findings.md`, `B_1..4` | Double increments only in MODIS Terra+Aqua (51-60% of increment tile-days); T+A -> J1+J2 cuts annual absolute increments ~35%. |
| C skill | `C_skill_results.md`, `figC_1..4` | DA beats OL slightly; MODIS-vs-VIIRS and Terra-vs-T+A differences are 1-2 orders smaller than the RMSE; IMS DA benefit turns negative in SY5-6. |
| D water balance | `D_water_balance_tables.md`, `D_1..5` | Budget closes; DA adds ~65-134 kg m-2 per snow year, ~half to runoff, ~half to ET; Aqua adds ~46; MODIS T+A has ~36 more than J1+J2. |
| D2 additions | `D2_water_balance_additions.md`, `D2_1..4` | Partition fractions, RZMC response, significance maps; SY1 storage is a real transient. |
| E transition | `E_findings.md`, `E_summary.md`, `E_1..5` | No spin-up; first 8 weeks hold ~1% of the SWE integral. |

## How results will be read

- With one member and 3 snow years per contrast (n = 3 in time for C1, C2 and C4), the 5° spatial-block bootstrap
  captures spatial uncertainty only, not year-to-year variability, which is the larger uncertainty here. Report each
  snow year separately next to the pooled result, and use paired differences. Where a difference is not
  distinguishable from zero, say "no resolvable difference" and do not interpret maps.
- **Decision rule (CONFIRMED by amfox, 2026-11-21).** *Honest provenance:* the original version (|A - B| >= 10% of the
  DA - OL effect, same sign in all snow years) was proposed before any results were looked at. The skill and
  descriptive-term parts below were changed **after** seeing analysis C and D, because C showed the 10% reference is
  too weak for station skill (DA - OL is only about 1 mm RMSE at SNOTEL/GHCN and about 0 for IMS in SY4-6) and D showed
  spurious flags for residual and dStorage (reference near 0). The change was justified on those merits, but it is not
  pre-registered. Rule as confirmed:
  1. **Water-balance flux terms** (`I_snow`, ET, runoff, peak SWE, RZMC): a difference only if |A - B| is at least 10% of
     the DA - OL reference in every snow year and A - B has the same sign in every snow year (3; 6 for C3).
  2. **Skill metrics** (IMS accuracy, SNOTEL SWE/depth RMSE, GHCN depth RMSE): a difference only if A - B has the same
     sign in every snow year (3; 6 for C3) **and** the pooled paired block-bootstrap CI excludes 0. The 10%-of-DA-OL
     magnitude test is dropped. The C4 difference of differences has no per-year values, so it is CI-only and labelled so.
  3. **Residual, dStorage, dFSW, dP, snowmelt, end-of-melt SWE, SFMC, TWLAND:** reported descriptively, no flag.
  Anything not flagged is "no resolvable difference". Code: `scripts/viirs_transition/rule2.py`; output
  `rule2_skill.csv`, `rule2_water_balance.csv`, `rule2_changes.csv` (calls that changed versus the preliminary rule).
  `C_decision_rule_preliminary.csv` and `D_decision_rule.csv` keep the earlier all-terms 10% application.
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
- Analysis step is 3 h. `LANDASSIM_DT` itself was not found in the run files, but it is confirmed empirically: in the
  3-hourly `ObsFcstAna` the obs of each species at one step span 45° of longitude (3 h at 15°/h), and for a given
  longitude band Terra and Aqua sit in consecutive steps (`TA_control`, 2024-01-15 09z/12z/15z).
- **MODIS monitoring gap in 2026.** The OL and the VIIRS arms read monitor-only MODIS (species 52/53) from `gmao_ops`,
  which ends 2025-12-31, so their `ObsFcstAna` has no Aqua or Terra obs from 2026-01-01 (seen on 2026-01-15 and
  2026-07-15). The controls read `MODIS_obs/` and do have them. Assimilation is not affected, but any analysis that
  compares MODIS obs or innovations across arms (A, B and the cross-sensor checks) must take MODIS obs from the
  control arms for 2026.
- The MODIS obs seen by the arms are the same files before and after the switch: `MODIS_obs/` 2023 to 2025 are
  symlinks to the `gmao_ops` year directories (identical by construction); only the 214 files of 2026 are new
  (downloaded, `.061.hdf`). `gmao_ops` itself lacks 4 days in 2025.
- Species ids inside `ObsFcstAna` (1 Aqua, 2 Terra, 3 J1, 4 J2) were checked in 35 files: all five arms, before and
  after the switch, 2021 to 2026. Only `assim_flag` changes at 2023-08-01, as intended.
- **Obs coverage (analysis A, domain totals per day, SY4):** Terra about 21.9k tiles, Aqua about 18.2k; both on the
  same day at about 15.5k tiles (about 63% of the 24.5k tiles with any MODIS obs); Terra only 6.4k, Aqua only 2.6k.
  J1 about 21.2k, J2 about 20.9k; both at about 19.2k tiles (about 84% of the 22.8k with any VIIRS obs); J2 adds about
  1.6k tiles beyond J1 (about 7%). The VIIRS union (22.8k) is about 6–7% smaller than MODIS T+A (24.5k; SY5 23.6k vs
  25.2k), and J2 alone is about 5% below Terra alone. Coverage is almost all north of 30°N. J2 only starts 2023-02-10,
  so in SY3 it has about 10k tiles per day. Source: `docs/viirs_transition_report_figures/obs_availability_summary.md`
  and `output/obs_availability/` (git-ignored).
- **Input data gaps (days with no obs from a sensor that was being assimilated; no assimilation from it on those
  days, `zero_obs_days.csv`):** Aqua 2022-04-02 to 04-16 (15 days); Terra 2022-10-11 to 10-22 (12 days); Terra
  2025-05-11 to 05-14 (4 days; this leaves `T_control`, which assimilates Terra only, with no assimilation at all);
  J2 2024-02-28 (1 day). **Cause (checked 2026-11-21): the source files are missing from the `gmao_ops` archive**, not a reader or QC effect:
  `MYD10C1_V61/2022` has no files for 2022-04-01..04-16 (the Aqua gap in the archive starts 04-01, one day before the
  04-02 seen in the obs, probably a window offset), `MOD10C1_V61/2022` none for 2022-10-11..22 and `MOD10C1_V61/2025`
  none for 2025-05-11..14, and `VIIRS/VJ210C1/Y2024/M02/D28/` exists but is empty. Whether NSIDC has the MODIS days
  upstream was not checked; J1 on 2024-02-28 was not checked. The controls read the same 2025 gap.
  Also, 2023-08-01 is skipped in analysis A because its 00z file still carries the old segment's final analysis
  (it assimilates Terra, e.g. 1221 Terra and 4144 J2 obs in `T_to_J2`).
- Output volume: `ObsFcstAna` is ~160 MB and ~250 files per arm-month; `catch_progn_incr` is ~34 MB per daily file
  (~1 GB and 31 files per arm-month). `land_da` was at ~91% of its file quota on
  2026-11-19, so extract compact per-day or per-month summaries rather than copying or unpacking files, and check
  `showquota` (inodes) before large jobs.

## Planned layout

`config/arms.json` (arms, namelists, species, dates, contrasts, paths: scripts read it instead of hard-coding),
`scripts/viirs_transition/` (one script per analysis A–E, each with a `_discover.sh`/`.j` wrapper), outputs under
`output/` (not versioned), results in `viirs_transition_report.md`.

## Where things run

Everything that reads the model output (`ObsFcstAna`, `catch_progn_incr`, restarts) runs on Discover, with the GEOSpyD
py3.13 environment (`/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3`) and a SLURM wrapper per
script. Only compact summaries (per day or per month, per arm) are copied to the Mac for plotting notebooks; the
Mac-side notebooks must not assume the raw output is available. `arms.json` holds Discover paths and is read by the
Discover scripts.

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

1. `LANDASSIM_DT` / `LANDASSIM_T0` as an explicit setting (the 3 h step itself is confirmed from the data, see above;
   optional).
2. ~~Analysis C scoring.~~ **Closed 2026-11-21:** scoring generalised to the five arms and run (`skill_*.py`; GHCN
   `min_valid_days` rescaled 1500 -> 360 for the 2191-day window; SNOTEL inches converted to mm); results `C_*`.
3. ~~Confirm the decision rule.~~ **Closed 2026-11-21:** confirmed in the form in "How results will be read" (skill: same sign + bootstrap CI; water-balance flux terms: 10% + same sign; residual/storage descriptive).
4. Whether the real Terra and Aqua overpass times drifted (affects interpretation of the time-of-day confound only).
5. ~~Availability of `PEATCLSM_FSWCHANGE` and budget-endpoint restarts.~~ **Closed 2026-11-21** (see analysis D; caveat:
   file counts from a glob also match monthly means, count dailies by date).
6. **IMS DA benefit turns negative in SY5-SY6** (DA - OL accuracy +0.013 in SY1-3, +0.006 in SY4, -0.002..-0.004 in
   SY5-6), also in the MODIS control arms, so it is not caused by the VIIRS switch. Cause not investigated (IMS or OL
   change, MODIS obs quality in 2025-26, or a real effect). Check before resting any C1/C2 conclusion on IMS.
7. Upstream status of the input gap days (below): not checked whether NSIDC has the missing MODIS days (CMR approach in
   the project memory); an optional download would let the MODIS arms use them (needs amfox's go-ahead).
8. J1/J2 obs counts in week 8 of 2023 (26 Sep - 2 Oct) are ~25-30% below 2024-25 (analysis E); not investigated.
9. Reconcile the units of B and D increment totals (B: per tile over 52,732 tiles >= 30N; D: area-weighted over the
   46,602 seasonal-snow mask tiles): D's net `I_snow` (TA MODIS ~105-134) is larger than B's net per-tile mean (~79-82).
   Probably the different denominators, not verified.

## Log

| Date | What |
|---|---|
| 2026-11-21 | Plan drafted from the run configuration and the GEOSldas source (reader overpass hours, per-tile obs averaging, SCF update parameters). |
| 2026-11-21 | `config/arms.json` written. Open items updated: species/assim flags confirmed, all five arms reached `END_DATE 20260801`, `ObsFcstAna` carries fcst/ana for non-assimilated species (local species ids 1-4). Increment path corrected to `cat/ens0000/` (daily files, 8 steps). |
| 2026-11-21 | After review: no sixth arm (VIIRS constellation answered as a difference of differences); IMS kept but flagged as not independent; per-year reporting and a proposed decision rule added; 3 h step confirmed from `ObsFcstAna` longitude spans; MODIS monitoring gap in 2026 for the OL and VIIRS arms recorded; species ids verified in 35 files; data status and run locations added. |
| 2026-11-21 | **Analysis A run** (SLURM array 59362332, all five arms; 11 min for the 6-year arms): obs coverage per sensor and constellation overlap, input data-gap days found, 2026 MODIS gap confirmed in the OL and VIIRS arms. Plan updated with the coverage numbers and gap days. |
| 2026-11-21 | `.gitignore` now un-ignores `docs/viirs_transition_report_figures/*.png`; `obs_availability.j` walltime set to 15 min. |
| 2026-11-21 | **Validation data fetched and verified** (`scripts/validation_data/`, data in `land_da/VIIRS_testing/validation_data/`, 28 files, 145 MB). SNOTEL 2025-10-01 to 2026-07-31 (304/304 days, 913 stations, inches; the existing parquet covers to 2025-09-30; recent values are NRCS provisional); GHCN SNWD per-year parquet 2020-08-01 to 2026-07-31 (2191/2191 days, 29,284 stations, mm); IMS 2025 and 2026 (365/365 and 212/212 days) plus M36 nearest-neighbour versions, 2020 to 2024 already in `projects/IMS/output/`. Spot checks: IMS byte-identical to the existing 2024 raw files; SNOTEL overlap differences small (10 values, consistent with NRCS revisions); GHCN vs the old OL/DA product 0.004% of shared station-days differ, with some station-days present in only one product (likely NOAA revisions; not investigated). |
| 2026-11-21 | **Analyses B, C, D, D2 and E run** (SLURM arrays 59364204 B, 59364397/59364576 C, 59364182 D, 59364691/59364715 D2, 59364708 E). Results summarised in "Results at a glance" and `viirs_transition_report.md`. Open items 2 and 5 and the D [CHECK] closed; input-gap cause recorded; decision rule kept PROPOSED with the C finding and options; IMS SY5-6, upstream MODIS gaps, week-8 J1/J2 dip and B/D unit reconciliation added as open items. |
| 2026-11-21 | **Decision rule confirmed by amfox** (water-balance flux terms keep 10% + same sign; skill = same sign in every year + bootstrap CI excludes 0; residual/dStorage descriptive). Set after the results were seen; flags re-run with `rule2.py`. Changed calls: C2 SNOTEL depth (difference -> none), C4 post SNOTEL SWE (none -> difference); 14 water-balance terms now descriptive. |
