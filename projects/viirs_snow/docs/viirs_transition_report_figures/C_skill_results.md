# Analysis C: independent snow skill of the five arms (draft, 2026-11-21)

> **IMS correction (2026-11-22).** The IMS numbers originally in this file used 2024+ IMS arrays that are flipped N-S on the M36 grid (see `IMS_findings.md`).
> All IMS statements below now use the corrected scores (`IMS_skill_by_arm_fixed.csv`, `IMS_skill_contrasts_fixed.csv`, flags in `rule2_skill_imsfix.csv`).
> `figC_1_ims_contrasts.png` and `C_pooled_contrasts_table.csv` were regenerated from the corrected scores (`skill_figures.py` now swaps in `IMS_skill_contrasts_fixed.csv`).
> `C_skill_by_arm.csv`, `C_skill_contrasts.csv`, `C_decision_rule_preliminary.csv` (superseded preliminary rule) and `rule2_skill.csv` still hold the original flipped-IMS values; SNOTEL and GHCN rows are unaffected (they do not use IMS).

Scripts (`scripts/viirs_transition/`): `skill_common.py`, `skill_extract.py` + `skill_extract.j` (per-arm extraction, SLURM array, 3-6 min per arm),
`skill_score.py` + `skill_score.j` (metrics, contrasts, paired 5-degree block bootstrap, 1000 draws, 3 min), `skill_figures.py`, `skill_rule.py`.
Compact intermediates: `land_da/VIIRS_testing/skill_data/skill_<arm>.npz` (not versioned). Tables here: `C_skill_by_arm.csv`,
`C_skill_contrasts.csv`, `C_pooled_contrasts_table.csv`, `C_decision_rule_preliminary.csv`; figures `figC_1..4`.

## Setup
- Window 2020-08-01..2026-07-31 (SY1-SY6), model = daily `tavg24_1d_lnd_Nt` (single member, `cat/ens0000`), same 112573 tiles for all arms.
- **IMS** (not independent of MODIS/VIIRS imagery): IMS 24 km on M36 (existing 2020-2024 files + the 2025/2026 files fetched 2026-11-21),
  model snow = `FRLANDSNO >= 0.5`, per-cell A/B/C/D counts, NH seasonal-snow cells (>= 10 IMS snow days, <= 90% of valid days; 45842 cells), "ALL" = NH
  without JJA / SH without DJF. Metrics: accuracy, CSI, hit rate, false-alarm ratio.
- **SNOTEL** (inches -> mm: SWE x25.4, depth x25.4; model SWE kg m-2 = mm, model `SNODPLAND` m x1000): nearest tile <= 40 km, |dz| <= 500 m, >= 360 valid obs days,
  >= 5 snow days per snow year; 609 (SWE) / 633 (depth) stations. **GHCN** (mm already; model depth = `SNODPLAND` x `FRLANDSNO`): same rules; 8857 stations.
  `min_valid_days` 1500 (25-y record, 9131 d) rescaled to 360 for the 2191-d window; snow days per year unchanged (5). Metrics: station RMSE, ubRMSE, |bias|, averaged over
  stations; CI by resampling 5-degree blocks of stations. Non-JJA days, same station-days for all arms of a window.
- delta = first minus second term (C0 = DA - OL, C1 = MODIS Terra - J2, C2 = MODIS T+A - J1+J2, C3 = Terra - T+A (MODIS series, 6 y), C4 = difference of differences
  (also C4a post, C4b pre shown separately)). For RMSE negative delta = first term better; for IMS accuracy positive = first term better.

## Findings
1. **DA vs OL (C0).** All DA arms beat OL against all three data sets over SY1-6 (pooled): IMS accuracy +0.014..+0.016 (OL 90.6% -> DA 92.1..92.3%; corrected IMS), SNOTEL SWE RMSE -1.0..-1.3 mm (of ~208 mm, ~0.5%),
   SNOTEL depth -4 mm (of ~582), GHCN depth -1.2..-1.9 mm (of ~114). Effects are small relative to the errors. The model has a large low bias at SNOTEL (SWE bias about -134 mm).
2. **IMS DA impact is stable in time (after the IMS orientation fix).** DA - OL accuracy is +0.015 in SY1-3 and +0.010..+0.019 in SY4-6, positive with CI excluding 0 in every snow year and every DA arm
   (T_to_J2: +0.015/+0.015/+0.015/+0.017/+0.019/+0.013). The earlier "negative in SY5-SY6, also for the MODIS control arms" was a scoring artefact: IMS 2024+ was flipped N-S (`IMS_findings.md`). Not a VIIRS effect and not a real drift.
3. **C1, C2 (same start state, SY4-6): IMS favours VIIRS slightly** (MODIS - VIIRS accuracy -0.0026 [-0.0030, -0.0023] for C1 and -0.0026 [-0.0032, -0.0020] for C2, same sign in SY4, 5, 6; about 17-18% of the DA - OL effect;
   a flag under the confirmed rule, but IMS is not independent of MODIS/VIIRS imagery, so this may partly reflect shared inputs).
   Against SNOTEL SWE, MODIS (control) is worse than VIIRS: RMSE +0.39 mm (Terra vs J2, all three years same sign) and +0.26 mm (T+A vs J1+J2, SY6 reverses); GHCN depth RMSE +0.93 / +0.23 mm
   (MODIS worse, same sign all years). SNOTEL depth: no resolvable difference for C1, +0.48 mm (CI spans 0) for C2. Sizes are 1-2 orders smaller than the RMSE.
4. **C3 (Terra vs T+A, MODIS, 6 y, same start):** Aqua helps slightly against SNOTEL SWE (Terra-only RMSE +0.23 mm, same sign all 6 y) and GHCN (+0.60 mm, same sign in all years),
   but Terra-only is better against SNOTEL depth (-0.53 mm); IMS accuracy -0.0005 [-0.0010, -0.0001] (Terra-only slightly worse, but the sign is not the same in all years: SY1-5 -0.0001..-0.0010, SY6 0.0000; no flag).
5. **C4 (J2 vs J1+J2) as a difference of differences:** IMS -0.0002 [-0.0007, 0.0003], GHCN -0.78 mm [-0.99, -0.60], SNOTEL SWE -0.10 and depth +0.27 mm (CIs span 0).
   The baseline (Terra vs T+A in SY1-3) is of the opposite sign to the post window for GHCN (+0.62 pre, -0.16 post), so the DoD is dominated by the baseline shift; treat as inconclusive and
   use C3 as the primary estimate of what a second daily sensor adds (as in the plan).
6. **Preliminary decision rule** (superseded by the confirmed rule, `rule2_skill.csv` / `rule2_skill_imsfix.csv`; `C_decision_rule_preliminary.csv` keeps the original, uncorrected-IMS application): reference effects (DA - OL) are only ~1 mm RMSE for station data and ~0.014 for IMS after the fix,
   so the preliminary 10% threshold was a very low bar for station data (the IMS ratio of 120% in that file used the flipped-IMS reference of ~0.0002 and is obsolete). Needs a decision on the reference (pooled SY1-6 DA-OL?) before it can be applied.

## Caveats
- One member, three snow years per contrast: CIs capture spatial sampling only; year-to-year spread is shown by the per-year dots. Block bootstrap treats stations/cells in a 5-degree block as dependent only.
- IMS is not independent (shown on every IMS table/figure); SNOTEL (western US/Alaska) and GHCN (sparse) cover little of the snow domain, SNOTEL provisional for recent months.
- SNOTEL/GHCN depth: model grid-cell mean vs point; SNOTEL tiles are biased low (elevation). Daily model means are 24 h averages centred on 12 UTC; station reports are daily values, not matched in time of day.
- Not done: scoring each arm against sensors it does not assimilate (secondary check in the plan, needs ObsFcstAna), per-season station scores, IMS-by-snow-fraction or lat-band breakdowns, Southern Hemisphere split.
- `C_skill_*.csv` rows with blank `season` are station data (season "ALL" = non-JJA/non-DJF-SH only); IMS rows also carry DJF/MAM/SON for NH cells.
