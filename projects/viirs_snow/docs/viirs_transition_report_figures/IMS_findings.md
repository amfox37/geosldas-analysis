# IMS SY5-6 "DA benefit turns negative": cause and corrected scores (2026-11-21)

**Cause: scoring artefact. The 24 km IMS arrays for 2024-01-01 onward are vertically flipped relative to 2000-2023 on the M36 grid.**

Evidence
- `IMS_product_monthly.csv` / `IMS_1_break_2024-01-01_maps.png`: on 2023-12-31 -> 2024-01-01 the share of M36 land cells labelled IMS code 1 ("sea") goes from 6.2% to 43% (5,194 -> 35,977 cells); NH snow fraction drops 0.62 -> 0.43 in one day. Central Asia, Tibet, India, Mexico and N South America turn into "sea". Same level in 2025 and 2026 (both fetched separately), so it is the 2024+ files, not a download problem.
- `ims_diag_orient.py`: with identical source-index maps (verified equal for 2023 and 2024), flipping the raw 2024/2025 arrays vertically (`raw[::-1]`) restores the sea fraction to 6.2-8.1% (as-is 43-46%; other transforms 45-59%). The raw files state "(1,1) starts at: lower left corner" (bottom-up rows); the 2000-2023 files were evidently built top-down. The raw 2024 arrays are internally coherent (`IMS_2_raw_2024_orientation_check.png`).
- Effect on scoring: the flipped field puts "sea" (valid, no-snow) over snowy Eurasia, so model snow there counts as false alarm, for OL, DA and controls alike. It starts 2024-01 (SY4, Jan-Jul only) and covers all of SY5-6, matching the pattern (also in the control arms; VIIRS is not involved).
- Other candidates (MODIS 2026 input, forcing change, model drift) are not needed to explain it; monthly DA-OL after the correction is positive every winter in every year (`IMS_3_monthly_accuracy.png`, `IMS_monthly_DA_minus_OL.csv`).

Corrected scores (`IMS_skill_by_arm_fixed.csv`, `IMS_skill_contrasts_fixed.csv`; scoring code unchanged, flipped 2024+ IMS, 1000 block bootstrap draws). Accuracy, all seasons:

| Contrast | Window | As scored in C | Corrected |
|---|---|---|---|
| DA-OL T_to_J2 | SY1-6 | +0.0070 | +0.0157 [0.0132, 0.0184] |
| DA-OL T_to_J2 | SY4 / SY5 / SY6 | +0.0059 / -0.0031 / -0.0018 | +0.0173 / +0.0187 / +0.0126 |
| DA-OL TA_to_J1J2 | SY5 / SY6 | -0.0044 / -0.0027 | +0.0194 / +0.0128 |
| DA-OL T_control, TA_control | SY4-6 | +0.0001 / -0.0005 | +0.0136 / +0.0143 |
| C1 T_control - T_to_J2 | SY4-6 | -0.0002 (CI incl. 0) | -0.0026 [-0.0030, -0.0023], same sign SY4, 5, 6 |
| C2 TA_control - TA_to_J1J2 | SY4-6 | 0.0000 (CI incl. 0) | -0.0026 [-0.0032, -0.0020], same sign SY4, 5, 6 |
| C3 MODIS_T - MODIS_TA | SY1-6 | +0.0003 (CI incl. 0) | -0.0005 [-0.0010, -0.0001]; SY6 0.0000, SY1/2 CI incl. 0 |
| C4 DoD | | +0.0008 [0.0003, 0.0012] | -0.0002 [-0.0007, 0.0003] (inconclusive) |

Absolute accuracy (OL): SY1-3 0.902, SY4-6 0.911 (was 0.878 / 0.767).

Remaining features
- A small DA-OL deficit of -0.001..-0.004 in Sep-Oct (snow onset) persists in some arms and years (e.g. 2023-09, 2024-09, 2025-10 for TA arms; also before the IMS break); worth a look but not a drift.
- IMS is not independent of VIIRS/MODIS imagery, so a "VIIRS arms score higher than MODIS controls" (C1, C2: ~17% of the DA-OL effect, same sign every year) may partly reflect the IMS inputs.
- IMS 2024+ is also probably wrong in any earlier analysis using `projects/IMS/output/ims_snowcover_24km_2024_on_m36_nearest.nc4` (e.g. the M21C 2000-2024 IMS comparison, year 2024): not checked.
- The regridded 2024/2025/2026 nc4 files themselves were not modified or deleted; the fix is applied in `ims_extract.py` from the raw files.
