# E2: J1 obs dip (2023 week 8) and 1-28 Aug mask for C1/C2

## 1. The 2023-09-26..10-02 J1 obs dip is real missing input, not a reader/QC or counting effect
- Input: `/discover/nobackup/projects/gmao/polarm/JPSS/VJ110C1/2023/` (the default `obs_param_nml(59)%path`) has no file for
  2023-09-30 (DOY 273) and a half-size file for 2023-10-01 (3.30 MB vs ~6.5 MB). The J2 files for the same days are normal.
- ObsFcstAna obs per day (TA_to_J1J2): J1 14632 (09-29), **542 (09-30)**, **7661 (10-01)**, 15838 (10-02); J2 normal every day.
  The 542 on the missing day are window-edge obs from the neighbouring days' files (3 h windows across UTC midnight).
  Obs counts follow the file size; nothing points to a QC or reader change. n_steps = 8 every day.
- Week 26 Sep-2 Oct mean J1 obs/day over the mean of the two neighbouring weeks: 2023 **0.73**, 2024 0.98, 2025 0.98. J2: 0.95,
  0.98, 0.98. Excluding the two bad days the 2023 week mean is 15223/day, normal. So the dip is J1 only (E quoted "J1/J2").
- The same scan of 2023-08..2026-07 (`E2_dip_files.csv`, `E2_dip_obs_vs_files.png`): J1 missing 2024-02-04 and half-size 2024-02-05;
  J2 missing 2024-02-27 AND 2024-02-28 (the earlier gap list had only 02-28); J2 half-size 2025-09-18 (obs count barely lower,
  10306 vs ~10.8-13.9k). The J1 2024-02-04/05 gap and J2 2024-02-27 are not in the 34-day gap list used by B and D2.
- Effect on contrasts: all VIIRS arms lose the same J1 days; the J1+J2 arm (TA_to_J1J2) keeps J2 on 09-30; the J2-only arm
  (T_to_J2) is unaffected. Not a transition spin-up effect.

## 2. 1-28 Aug 2023 mask on C1/C2: no change (it is a no-op)
- The NH season mask already excludes all JJA days (and so all of August) for NH cells and stations; only SH cells/stations are
  affected. SNOTEL: 0 of 609 (SWE) / 633 (depth) stations are in the SH; GHCN: 1 of 8857. IMS: the per-arm valid-cell counts of
  SY4 are identical with and without removal of the 28 days (N = 12,560,708 for every arm), i.e. no eligible SH cell has
  contributing days in the window.
- Point estimates and verdicts (rule2 logic) are identical to 5 decimals in all 8 C1/C2 metric rows (`augmask_delta.csv`,
  `augmask_flags.csv`). The unmasked pass reproduces the C numbers (e.g. C1 GHCN RMSE +0.929, C1 SNOTEL SWE +0.390).
- Implication: the 1-28 Aug mask is only meaningful for the water-balance/SWE contrasts (D, E), not for the skill scores.
