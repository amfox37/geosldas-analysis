# MODIS-to-VIIRS transition and constellation size: synthesis (draft, 2026-11-21)

Plan: `viirs_transition_plan.md` (design, contrasts, caveats). Detailed results: `viirs_transition_report_figures/`
(`B_findings.md`, `C_skill_results.md`, `D_water_balance_tables.md`, `D2_water_balance_additions.md`, `E_findings.md`,
`obs_availability_summary.md`). All numbers below are copied from those files.

**Status of this report:** single member, M36, snow only, 2020-08-01 to 2026-08-01. The decision rule of the plan is
**not confirmed**, so no contrast is declared "resolvable" here by rule; statements say what the tables show and which
differences have the same sign in every snow year. Read the caveats at the end first.

## Design in one paragraph

Five arms from one restart: OL, `T_to_J2` (Terra, then J2 from 2023-08-01), `TA_to_J1J2` (Terra+Aqua, then J1+J2) and two
MODIS controls that continue Terra / Terra+Aqua after 2023-08-01 from the pre-switch arms' own restarts. Snow years
(SY) run 1 Aug - 31 Jul, SY1-3 before and SY4-6 after the switch. C1 = Terra vs J2 and C2 = Terra+Aqua vs J1+J2
(same start state, SY4-6), C3 = Terra vs Terra+Aqua (MODIS, six years), C4 = J2 vs J1+J2 (difference of differences),
C0 = DA minus OL. A MODIS-vs-VIIRS difference is a difference of the **whole configuration** (sensor, overpass time,
QC, number of updates per day), not of the sensors alone.

## Q3. Update structure (analysis B; obs coverage from A)

- **Double increments are a MODIS Terra+Aqua feature.** 51-60% of tile-days with an increment get two (never three) in the
  MODIS Terra+Aqua configuration; 0% in any Terra-only or VIIRS configuration, because J1 and J2 share one window and are
  averaged into one obs per tile. 97% of doubles are in consecutive 3-hourly steps.
- **The pair mostly reinforces but partly cancels:** 69-75% have the same sign; pair net/abs is 61-64%, so about a third of
  the absolute adjustment cancels between Terra and Aqua. The second increment is about 10% smaller than the first.
- **Volume.** Mean annual absolute increment (kg m-2 yr-1, per tile >= 30N): Terra-only about 65-76, Terra+Aqua MODIS
  about 120-138, `TA_to_J1J2` after the switch 79-85 (a ~35% drop at the switch). The VIIRS configuration is closer to
  Terra-only in volume than to Terra+Aqua MODIS. Net/abs is 64-71% in all arms; 1.7-3.3% of increments hit the
  5 kg m-2 cap, so no arm is cap-dominated.
- **The first update does most of the work:** Aqua RMS O-F (SCF) goes 0.32 (OL) -> 0.244 (Terra assimilated) -> 0.237
  (Terra+Aqua); J1 0.320 -> 0.2725 (J2 only) -> 0.2689 (J1+J2).
- **Coverage (A):** the VIIRS union covers about 6-7% fewer tiles than MODIS Terra+Aqua; J2 adds about 7% over J1.

## Q1. Does the MODIS -> VIIRS switch change skill and water balance?

**Skill (analysis C; IMS is not independent of MODIS/VIIRS imagery, so SNOTEL SWE and GHCN depth carry the conclusion).**
DA beats OL in all three data sets over SY1-6, but the effect is small: IMS accuracy 82.2 -> 82.9%, SNOTEL SWE RMSE
-1.0..-1.3 mm (of ~208 mm), GHCN depth RMSE -1.2..-1.9 mm (of ~114 mm); the model has a large low bias at SNOTEL (SWE
about -134 mm). For the switch (C1, C2) IMS shows no resolvable difference; against SNOTEL SWE the MODIS controls are
slightly worse than VIIRS (RMSE +0.39 mm for Terra vs J2, same sign all three years; +0.26 mm for T+A vs J1+J2, SY6
reverses) and against GHCN worse by +0.93 / +0.23 mm (same sign all years). These differences are one to two orders of
magnitude smaller than the RMSE.

**Water balance (D, D2).** The differential budget closes (residual about 0.01 kg m-2 against 65-134 kg m-2 of snow
increment). DA adds substantial water: net snow increment `I_snow` (kg m-2 per snow year, domain mean over 46,602 seasonal-snow
tiles) is 65-81 for Terra-only and 105-134 for Terra+Aqua MODIS. It is partitioned about half to ET and half to runoff
(Terra+Aqua - OL: 44% ET, 52% runoff, 3% storage; Terra-only - OL: 51% / 46% / 4%), raises peak SWE only by about 1-2 kg
m-2 and end-of-melt SWE hardly at all, and raises root-zone moisture by about 0.012 (Terra-only) and 0.019 (Terra+Aqua)
m3 m-3 as a year mean.

- **C1 (Terra vs J2): no resolvable difference** in water: `I_snow` -3.2, -1.3, +8.7 across SY4-6 (mean about 2% of the
  DA effect), sign not consistent; partition fractions unstable; RZMC 0.0004 [-0.0001, 0.0010].
- **C2 (T+A MODIS vs J1+J2): clear difference, same sign every year.** MODIS has more water input: `I_snow` +36
  (27.7, 31.6, 47.8), runoff +22, ET +11, lower peak SWE (-0.85), RZMC +0.005. The extra goes 63% to runoff and 31% to ET.
  It grows from SY4 to SY6. About 70% of mask tiles are marked for `I_snow` and runoff on a same-sign-in-3-years criterion
  (no t-test with n = 3); end-of-melt SWE is almost never marked (< 2%).
- **Reading:** the water-balance difference between the MODIS and VIIRS configurations is driven by the Terra+Aqua
  two-updates-per-day structure (B), not by Terra vs a single VIIRS sensor (C1 ~ 0), consistent with the ~35% drop in
  absolute increments at the switch. Skill differences are small and mixed in sign by data set.

## Q2. Can we drop to one satellite?

- **MODIS, Terra vs Terra+Aqua (C3, six years, same start; the most reliable estimate).** Terra-only gets about 46 less
  `I_snow` per snow year (-43.5 to -53.2, same sign every year; about 49% of the mean DA effect in the rule's terms),
  28.5 less runoff, 16 less ET; the missing Aqua water goes 62% to runoff and 35% to ET. RZMC is lower by 0.0067 [-0.0076,
  -0.0058], resolvable in every season. About 74% of mask tiles are marked for `I_snow` (same sign in all 6 years and
  FDR-corrected t-test); end-of-melt SWE almost none. Against independent data Aqua helps slightly: SNOTEL SWE RMSE
  +0.23 mm and GHCN +0.60 mm for Terra-only (both same sign all years), but Terra-only is better against SNOTEL depth
  (-0.53 mm) and IMS is ~0 (+0.0003, CI includes 0). So dropping Aqua changes the water budget a lot and the skill very
  little.
- **VIIRS, J2 vs J1+J2 (C4).** Raw J2-only - J1+J2: `I_snow` about -14 (-14.3, -13.6, -14.1), runoff about -8, same sign
  in all three years. As a difference of differences (subtracting the Terra - T+A gap of SY1-3) the sign flips to
  about +29 to +30, because the pre-switch gap was -43 to -48; the VIIRS second sensor therefore matters much less than
  Aqua did in MODIS. Skill C4 is inconclusive (the GHCN baseline has the opposite sign pre and post, +0.62 vs -0.16 mm),
  so use C3 as the primary constellation estimate.
- **Partition:** the missing second-sensor water is mostly runoff (C3, C2) and for C4 storage differs (-0.14 of the
  difference), so interpret C4 partition with care.

## Transition behaviour (analysis E) and sensitivities

- **No spin-up at the 2023-08-01 switch.** The VIIRS-arm / MODIS-control ratio of weekly NH absolute increments is the same
  from week 0 in 2023 as in 2024 and 2025; J1/J2 counts and O-F are within the 2024-25 range (except an unexplained 25-30%
  dip in J1/J2 counts in week 8 of 2023). The first 8 weeks hold about 1.1-1.2% of the NH SWE integral. **No exclusion
  window is needed**; masking 1-28 Aug for C1/C2 is a conservative option (< 0.5% of the integral).
- **SY1 storage is a real transient, not a restart artefact:** DA - OL storage builds to +18 (T+A) and +11 (Terra) kg m-2
  in SY1 and plateaus. Excluding SY1 lowers the pooled T+A dStorage from 4.1 to 1.3 but barely changes `I_snow`,
  runoff or C3.
- **Input gaps** (days with no obs from an assimilated sensor) lower DA - OL `I_snow` by 1-2% when removed from every
  arm (C3 -45.98 -> -44.98; C1 1.38 -> 1.91, still not resolvable). Cause: files missing in the `gmao_ops` archive.

## Open points that affect interpretation

1. **IMS DA benefit turns negative in SY5-6** (-0.002..-0.004 accuracy), also in the control arms. Not a VIIRS effect,
   cause unknown; IMS-based C1/C2 statements should wait for it.
2. **Decision rule unconfirmed.** 10% of DA - OL works for water-balance terms but is too weak for station skill
   (DA - OL ~ 1 mm RMSE, ~0 IMS in SY4-6); options in the plan.
3. **Increment totals in B and D differ in units** (per tile >= 30N vs area-weighted mask mean): B net about 79-82 for
   T+A, D `I_snow` about 105-134; probably denominators, not reconciled.

## Caveats

- One member and one (2023) transition; three snow years per contrast (C1, C2, C4) and six for C3. The 5-degree block
  bootstrap captures spatial sampling only, not year-to-year variability; look at the per-year values.
- IMS is not independent; SNOTEL (western US/Alaska) and GHCN (sparse) cover little of the snow domain; SNOTEL recent
  values are provisional; model grid-cell means are compared with point data.
- Water-balance endpoints use single-member restart TWLAND at 1 August; the gap-day correction cannot adjust storage.
- Terra and Aqua overpass drift not checked; "MODIS vs VIIRS" includes the time of day (Terra 10:30 vs 13:30 local).
- Counting in B is per UTC-day file (double share slightly underestimated in the 180W-135W sector); MODIS O-F rows end
  2025-12-31 in all arms but the controls.

## Figure index (`viirs_transition_report_figures/`)

| Topic | Files |
|---|---|
| A obs availability | [figA_1](viirs_transition_report_figures/figA_1_tiles_per_day_by_species.png), [figA_2](viirs_transition_report_figures/figA_2_constellation_overlap.png), [figA_3](viirs_transition_report_figures/figA_3_tiles_by_latitude_band.png) |
| B increments | [B_1 double fraction](viirs_transition_report_figures/B_1_double_fraction.png), [B_2 by longitude](viirs_transition_report_figures/B_2_step_by_longitude.png), [B_3 net/abs, cap](viirs_transition_report_figures/B_3_netabs_cap.png), [B_4 second-sensor innovations](viirs_transition_report_figures/B_4_second_sensor_innov.png) |
| C skill | [figC_1 IMS](viirs_transition_report_figures/figC_1_ims_contrasts.png), [figC_2 SNOTEL SWE](viirs_transition_report_figures/figC_2_snotel_swe_contrasts.png), [figC_3 SNOTEL depth](viirs_transition_report_figures/figC_3_snotel_depth_contrasts.png), [figC_4 GHCN depth](viirs_transition_report_figures/figC_4_ghcn_depth_contrasts.png) |
| D water balance | [D_1 C0 MODIS](viirs_transition_report_figures/D_1_C0_budget_modis.png), [D_2 C0 VIIRS arms](viirs_transition_report_figures/D_2_C0_budget_viirs_arms.png), [D_3 contrasts by SY](viirs_transition_report_figures/D_3_contrasts_by_snow_year.png), [D_4 seasonal cycle](viirs_transition_report_figures/D_4_seasonal_cycle_SY4_6.png), [D_5 maps](viirs_transition_report_figures/D_5_maps_pooled.png) |
| D2 additions | [D2_1 partition](viirs_transition_report_figures/D2_1_partition.png), [D2_2 soil moisture](viirs_transition_report_figures/D2_2_soil_moisture.png), [D2_3 significance maps](viirs_transition_report_figures/D2_3_significance_maps_C3_C2.png), [D2_4 SY1 storage](viirs_transition_report_figures/D2_4_SY1_storage_spinup.png) |
| E transition | [E_1 snow state](viirs_transition_report_figures/E_1_snow_state_2023.png), [E_2 divergence](viirs_transition_report_figures/E_2_divergence_2023.png), [E_3 increments by week](viirs_transition_report_figures/E_3_increments_by_week.png), [E_4 new-sensor innovations](viirs_transition_report_figures/E_4_new_sensor_innov.png), [E_5 increment ratio](viirs_transition_report_figures/E_5_increment_ratio.png) |
