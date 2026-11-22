# B vs D snow-increment totals: reconciliation

Script: `scripts/viirs_transition/bd_reconcile.py` (raw numbers: `BD_reconcile_table.csv`, gitignored).

**Conclusion: no bug. The two totals are the same quantity (sum of `WESNN1-3_INCR` over the 8 daily steps) averaged
over different tile sets with different weights.** B is an unweighted mean over all 52,732 tiles at or north of 30N;
D is an area-weighted mean over the 46,602-tile seasonal-snow mask. Area weighting is the dominant factor, the mask
the second, and everything else (gap days, the 0.01 kg m-2 threshold, the 365/n_days scaling) moves the number by
at most 2-3 kg m-2.

## Definitions

| | B (`incr_summary.py`, table B3) | D (`wb_analyze.py`, I_snow) |
|---|---|---|
| Quantity | sum over steps of the 3-layer increment | same (`wb_extract.py`) |
| Threshold | only steps with abs(incr) >= 0.01 kg m-2 | none |
| Days | 33 gap/boundary days dropped, then x 365/n_days | all days of the snow year |
| Tiles | lat >= 30N (52,732) | mask: lat > 20N, snow possible, JJA SCF < 0.2 (46,602; 45,051 of them are >= 30N) |
| Weighting | unweighted mean per tile | weighted by `catchment.def` area (column 7) |
| Period | per year (365/n_days) | per snow year (1 Aug - 31 Jul) |

## Decomposition (kg m-2 per snow year, one arm = one row)

Steps, in order: B headline -> D definition on the same tiles with no weighting (removes the 365/n_days scaling, the gap-day
treatment and the threshold) -> area weighting on the same tiles -> D's mask (the number in the D tables).

| arm | SY | B headline | + scaling / gap days / threshold | + area weighting (N30 tiles) | + mask instead of N30 | = D table |
|---|---|---|---|---|---|---|
| TA_to_J1J2 | SY1 | 81.2 | +0.2 | +18.4 | +9.0 | 108.7 |
| TA_to_J1J2 | SY2 | 81.6 | -2.0 | +16.1 | +9.0 | 104.6 |
| TA_to_J1J2 | SY3 | 79.3 | +1.0 | +27.6 | +10.8 | 118.7 |
| TA_control | SY4 | 81.2 | +0.0 | +27.0 | +9.3 | 117.5 |
| TA_control | SY5 | 89.8 | -0.1 | +18.8 | +10.8 | 119.3 |
| TA_control | SY6 | 93.0 | +0.0 | +28.5 | +12.3 | 133.8 |
| T_to_J2 | SY1 | 45.0 | +0.1 | +15.5 | +4.7 | 65.2 |
| T_to_J2 | SY2 | 47.0 | -0.7 | +14.7 | +4.7 | 65.6 |
| T_to_J2 | SY3 | 45.0 | -1.5 | +21.4 | +5.4 | 70.3 |
| T_control | SY4 | 47.0 | +0.0 | +20.6 | +4.7 | 72.3 |
| T_control | SY5 | 53.0 | -0.6 | +14.9 | +5.5 | 72.7 |
| T_control | SY6 | 54.1 | +0.0 | +20.5 | +5.9 | 80.5 |

The D column equals the "C0 MODIS_TA - OL" (TA_to_J1J2 SY1-3, TA_control SY4-6) and "MODIS_T - OL" rows of
`D_water_balance_tables.md`; OL has zero increments, so D's I_snow for C0 is just the arm's own total. The VIIRS-period
SY4-6 rows of TA_to_J1J2 and T_to_J2 (B 55-61 and 42-47; D 86-90 and 72-76) reconcile the same way
(see `BD_reconcile_table.csv`).

## Why area weighting matters so much

- Catchment tiles are very unequal in size (N of 30N: median 312, mean 631, 99th percentile 4,630 in the units of
  `catchment.def` column 7). The 3,154 tiles larger than 2,000 carry 33% of the area weight and have a mean increment of
  154 kg m-2 (SY1, TA_to_J1J2), against 80 for tiles as a whole. Large boreal/Arctic catchments have large increments, so
  the area-weighted mean is about 20 kg m-2 higher than the per-tile mean.
- The seasonal-snow mask removes snow-free and low-latitude tiles (about 4% of N30 tiles have no increment at all) and adds
  tiles between 20N and 30N, which raises the mean by a further 5-12 kg m-2.
- Both choices are legitimate; they answer different questions. Per-tile means describe how often and how hard the
  filter acts on a typical catchment; area-weighted means give kg m-2 over land and are what the water budget needs.

## What the small terms are

- **365/n_days scaling** (B only): extrapolates the days that remain to a full year. It over- or under-states by up to
  2-3 kg m-2 when the dropped gap days are not typical (SY2: Aqua gap in April, scaled B is 2.0 above the unscaled all-day value).
- **Threshold**: skipping increments below 0.01 kg m-2 changes the total by 0.01-0.2 kg m-2.
- **Gap days**: B's tile-level `tile_net` includes gap days; only the headline step sums exclude them.

## Side observation (no effect)

478 tiles (438 of them at or north of 30N) have a negative or zero area in `catchment.def` column 7; the area-weighted
mean changes by 0.007 kg m-2 with or without them, so it is harmless, but the weight vector is not strictly positive.
The same file is used by the M21C scripts.

## Consequences

- No edits to the B or D scripts or tables are needed.
- The report quotes B (about 79-82 for Terra+Aqua, per tile north of 30N) and D (105-134, area-weighted over the mask). Both
  are correct; the report should say which weighting each number uses (suggested wording in the hand-off).
- B's ratios (net/abs, share of double increments, VIIRS/MODIS ratio of totals) are unaffected, since they use the
  same weighting on both sides. The ratio of totals depends on weighting only moderately (TA_to_J1J2 / TA_control, SY4-6 mean: B 0.65, D 0.71).
