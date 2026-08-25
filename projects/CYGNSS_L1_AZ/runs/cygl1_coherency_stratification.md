# Does the CYGNSS L1 `coherency_state`/`coherency_ratio` flag predict L1-operator fit?

Follow-up to `cygl1_operator_diagnosis.md` (L1-vs-L3 operator diagnosis). Scripts:
`scripts/extract_cygl1_coherency_join.py`, `scripts/analyze_cygl1_coherency_stratification.py`
(primary, OL run) and `scripts/extract_cygl1_coherency_join_da_secondary.py`,
`scripts/analyze_cygl1_coherency_stratification_da_secondary.py` (secondary, DA run).

## Pre-registered prediction (stated before / independent of the result)

The CYGNSS L1 operator computes reflectivity via `mwRTM_get_lr_reflectivity` -- flat
Fresnel reflectivity, with an explicit in-code comment that no roughness or vegetation
correction is applied (confirmed during the earlier L1-vs-L3 diagnosis). This is a
coherent-scattering formulation. **Physics therefore predicts: coherent-flagged returns
(`coherency_state==1`) should show tighter/better-fit O-F than incoherent-flagged
returns (`coherency_state==0`).** If the opposite is found, that is equally
informative -- it would mean either the coherency flag isn't identifying what we think,
or the DDM3x5-crop-average scalar isn't actually isolating the coherent component the
flag describes.

**Bottom line up front:** the prediction is **partially supported and partially
confounded**. A naive pooled comparison looks consistent with the prediction, but a
within-tile control (comparing coherent vs. incoherent obs from the *same* tile) shows
that most of the pooled signal is a tile-selection artifact, not a per-observation
physical effect. What survives the within-tile control is a small (~0.27-0.31 dB),
statistically significant *reduction in O-F spread* for coherent obs (mild support for
the prediction), bundled with a much larger (~+2.0-2.3 dB), highly significant,
*unpredicted systematic positive bias shift* for coherent obs that dominates any naive
absolute-error comparison and makes that comparison net non-significant. See Step 5.

## Step 0: does the staging file already carry coherency fields?

**No.** `ncdump -h` on
`/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/Y2020/M0{1,2}/cygnss_l1_ddm3x5_crop_scalar_m36_<date>_all_cyg.nc4`
(the file GEOSldas actually reads) shows only: `obs_id, selection_order, sample_id,
ch_id, year, day, sc_num, status, tile_start/tile_count, active/mapped/unmatched pixel
counts, sp_nearest_tile_ig/jg/lat/lon/distance_km, observed_y_linear/db,
ddm_timestamp_*, sp_lat/lon/inc_angle, ddm_snr_db, incidence_angle_deg, frequency_hz,
unmatched_coefficient_fraction, vegetation_geometry_factor,
coefficient_weighted_opacity/attenuation, build/timing diagnostics, error_message,
product_json, obs_key, ddm_time_utc` on the `obs` dimension, and
`tile_index0/tile_ig/tile_jg/coefficient/coefficient_weight/active_pixel_count_by_tile`
on the `support` dimension. No `coherency_state`/`coherency_ratio` anywhere. Also
checked: `product_json` is an empty string for every `status==1` row (not a hidden
carrier of the raw granule payload) and the upstream per-spacecraft preprocessor files
(e.g. `.../preprocessor_groups/20200101_doy001_cyg01_selected100/cygnss_l1_ddm3x5_crop_scalar_m36_20200101_cyg01.nc4`)
have the identical schema -- coherency information is dropped somewhere between the raw
SDS read and this preprocessor's output.

This is **not** a one-step join. The coherency fields survive only in the earlier-stage
per-(day, spacecraft) QC-pass CSVs:
`/discover/nobackup/projects/land_da/CYGNSS_operator/artifacts/out_images/cygnss_qc_m36_window_counts_<YYYYMMDD>_cyg<NN>/cygnss_l1_qc_pass_<YYYYMMDD>_cyg<NN>.csv`
(re-verified present for all 60 Jan-Feb 2020 days x 8 spacecraft, with one legitimate
exception: `20200131` has only 7 QC-pass CSV dirs (missing `cyg04`) and the matching
staging file's own `spacecraft` global attribute independently confirms `cyg04` (and
`cyg02`) contributed zero rows to that day's merge -- a genuine data gap upstream, not a
join defect). These CSVs carry `sample_id, ch_id, ..., coherency_state, coherency_ratio,
...` and are keyed uniquely on `(sample_id, ch_id)` within each file (checked: 0
duplicate keys in every file touched by the join).

## Step 1: many-to-one check on the staging-file side

**Confirmed many-to-one, and GEOSldas's own reader already resolves it -- we
reproduced that exact resolution rather than inventing a new one.**

Reading `read_obs_cygnss_l1_scalar()` in
`GEOSldas_cygnss_operator/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90`
(lines ~3310-3660) shows the reader explicitly handles this: it assigns each raw obs to
an "owner tile" by matching `sp_nearest_tile_ig/jg` against the experiment's local
`tile_coord` (i_indg/j_indg), then **keeps only the candidate with the smallest
`sp_nearest_tile_distance_km` when more than one raw obs maps to the same owner tile in
the same assimilation window**, discarding the rest (there's a literal `N_duplicate`
counter in the log). The window itself is `(date_time - dtstep_assim/2, date_time +
dtstep_assim/2]` with `dtstep_assim = 10800 s` (3 h, matching the observed 3-hourly OFA
cadence and the `dtstep_assim(4)=10800` entry in `clsm_bias_routines.F90`).

Empirical counts (script: `extract_cygl1_coherency_join.py`, function
`step1_duplicate_check`):

- Coarse check over the **full raw candidate pool** (all tiles/spacecraft touched by any
  Jan-Feb 2020 staging file, not just the ones actually used in this AZ-domain OFA
  output; 38,902 status-ok raw rows binned into 3-h `(tile_ig, tile_jg, window)` bins):
  **5,059 of 33,843 bins (14.95%) have >1 competing raw candidate.** Distribution: 28,784
  bins with exactly 1 candidate, 5,059 with exactly 2, 0 with 3+.
- When the join is restricted to the exact `(owner tile, window)` pairs that actually
  appear in the OL-run OFA output (20,731 obs, see Step 2), only **5 of 20,731 (0.02%)**
  had more than one candidate competing for the same owner tile in the same window,
  correctly resolved by minimum `sp_nearest_tile_distance_km` per the reproduced
  algorithm.

The gap between 14.95% (global raw-candidate pool) and 0.02% (actual AZ-domain assimilated
obs) is expected: the staging files are built from all spacecraft passes globally-ish
merged per day, while the 909-tile AZ-limited experiment domain only pulls a small,
mostly non-overlapping-in-time subset of those candidates into its 3-h windows. Both
numbers are reported since they answer different questions (raw-pipeline duplication
rate vs. this-experiment's actual dedup rate); only the second one matters for the join
below, and it is small enough that mis-resolving the reader's tie-break would have
negligible effect if it happened at all here.

## Step 2: join (OL run, primary)

**Method:** exact-key joins only, no value-matching anywhere in the chain (per the
brief's warning about ~10% repeated OFA values at 4 decimal places):
1. OFA `(datetime, tilenum)` for species `CYGNSS_L1_DDM3X5_CROP_SCALAR`, resolved by name
   from each file's own `obsparam_descr`/`obsparam_species_id` (never a hardcoded
   index) -- run: `OLv8_M36_all_sensors_AZ_scaled`. `tilenum-1` indexes the tilecoord
   arrays (same convention as `extract_cygl1_l3_pairs.py`), giving `(i_indg, j_indg)`.
2. Reproduce the reader's exact owner-tile + 3-h-window + nearest-distance tie-break
   (Step 1) against the pre-loaded Jan-Feb 2020 staging-file candidate index, to
   identify the single raw `(year, day, sc_num, sample_id, ch_id)` GEOSldas actually
   used for that OFA row.
3. Look up `(sample_id, ch_id)` in the matching `(day, spacecraft)` QC-pass CSV for
   `coherency_state`/`coherency_ratio` -- an exact integer-key lookup, not fuzzy.

**Validation (strong, multiple independent cross-checks against numbers the user had
already computed independently before this script was written):**

| Quantity | User's independently-known value | This join's result |
|---|---|---|
| N obs, Jan-Feb 2020 | 20,731 | **20,731** (exact) |
| N distinct tiles | 498 | **498** (exact) |
| O-F sd | ~3.19 dB | **3.1863 dB** |

All 20,731 OFA rows matched a raw candidate (0 dropped for bad tilenum, 0 dropped for
no candidate in window) and all 20,731 found a coherency match in the QC-pass CSVs (0
missing). This three-way exact agreement is strong evidence the tile/window/owner-tile
reproduction of the Fortran reader is correct, not a coincidence.

Sanity-only (not used in the join): `obs_db (OFA) - raw_obs_db (staging)` has
mean 4.34 dB, sd 2.87 dB -- a real, non-constant offset, consistent with the
per-tile/per-time z-score obs scaling (`scale_obs_cygl1scal_zscore`) rather than a bug.

## Step 2b: DA-run secondary check

Run: `DAv8_M36_all_sensors_AZ_scaled_cygl1assim`. Species index **re-verified directly
from this experiment's own obsparam** (not assumed from the OL run): `CYGNSS_L1_DDM3X5_CROP_SCALAR`
= species **13** here as well (coincidentally the same index as the OL run -- confirmed
independently, per the brief's caution not to assume), with `obsparam_assim=1` (this run
actually assimilates the species, unlike the OL run) and `obsparam_errstd=4.4 dB` (the
value the user had already noted; the OL run's own errstd is 3.0 dB -- these are two
different experiment configs, not the same R). Same join method, same Jan-Feb 2020
window: 20,756 obs matched (vs. 20,731 in the OL run -- a very close but not identical
count, consistent with a slightly different tile-domain/getinnov footprint), 498
distinct tiles, all rows found a coherency match, `assim_flag==1` for all 20,756 rows
(confirming this is genuinely the assimilating run).

**This is treated as secondary evidence only** -- O-F in an assimilating run reflects
whatever the analysis increment did on prior cycles, not pure operator fit. Result
(Step 4/5 numbers below): **essentially identical to the OL run** (signed within-tile
diff +2.26 vs +2.25 dB; sd diff -0.287 vs -0.274 dB; both p<0.0001; |O-F| diff
non-significant in both, p=0.86 vs p=0.82). The DA run's own increments have not
materially altered the pattern, which strengthens confidence that the OL-run result
below isn't an artifact of that run's particular scale/errstd configuration.

## Step 3: does `obsvar` vary?

**No, it varies substantially -- both raw and normalized stratifications are reported.**
`obsvar` (OL run) ranges from 0.118 to 173.96 with mean 8.54, sd 11.75, and has 4,362
distinct values (rounded to 4 dp) across the 20,731 matched obs -- not remotely a
constant, despite the fixed `errstd` metadata scalar (3.0 dB OL / 4.4 dB DA), because
`obsvar` here reflects the z-score-rescaled obs error (`scale_obs_cygl1scal_zscore`
adjusts obs error to the local model climatology's dB variance, so it varies by
tile/season). Both raw O-F and normalized `z = O-F/sqrt(obsvar)` are reported below;
qualitatively they tell the same story (see per-bucket `mean z` values), so the
verdict does not depend on which one is used.

## Step 4: stratify by `coherency_state` and `coherency_ratio` (OL run, primary)

`coherency_state` legend: 0=not coherent, 1=coherent, 2=mixed, 3=indeterminate.
**Only states 0, 1, 2 occur in the Jan-Feb 2020 AZ-domain data; state 3 (indeterminate)
was never observed in this window.**

| state | N | mean O-F (dB) | median O-F (dB) | sd O-F (dB) | mean z(O-F) | mean fcst (dB) | mean co-loc. L3 SM6hr fcst (m3/m3, N with match) | n tiles | top tile share |
|---|---|---|---|---|---|---|---|---|---|
| 0 not_coherent | 17,020 | -0.510 | -0.580 | 3.213 | -0.260 | -11.227 | 0.1205 (N=6,197) | 496 | 0.6% (tile 5706) |
| 1 coherent | 2,913 | **+1.041** | +0.600 | 2.966 | +0.481 | -9.659 | 0.1147 (N=1,024) | 317 | 2.0% (tile 5383) |
| 2 mixed | 798 | +0.510 | +0.245 | 1.853 | +0.190 | -6.665 | 0.1335 (N=300) | 153 | 3.3% (tile 5923) |

No bucket has any tile >20% of its count (max is 3.3%, for state=2), so tile
dominance is not a first-order confound of the pooled numbers above (it *is*, however,
a first-order confound for the coherent/incoherent *comparison*, resolved in Step 5).

**coherency_ratio** (continuous, 0.068-8.843, pooled mean 1.16):
- Full-population Spearman rho(ratio, O-F) = **+0.260** (p~0, N=20,731): higher ratio
  associates with a *more positive* O-F -- same direction as the state-level bias shift.
- Full-population Spearman rho(ratio, |O-F|) = **-0.183** (p=4.4e-156): higher ratio
  associates with *smaller absolute error* -- naively consistent with the prediction, but
  see the within-tile caveat below; this pooled correlation is not tile-controlled.
- Quantile bins **within the coherent (state==1) subset only** (N=2,913, comfortably
  above the ~450/bin threshold even at 3 bins, so not capped at 2 as the brief's more
  conservative ~1,000-obs coherent-fraction estimate assumed -- actual coherent N here is
  2,913, about 3x that estimate):

  | ratio bin | N | mean O-F | median O-F | sd O-F |
  |---|---|---|---|---|
  | (2.00, 2.19] | 971 | 0.596 | 0.326 | 2.564 |
  | (2.19, 2.46] | 971 | 1.034 | 0.600 | 2.908 |
  | (2.46, 8.84] | 971 | 1.491 | 1.084 | 3.312 |

  **Notable and counter to the naive pooled correlation:** *within* the coherent subset,
  higher `coherency_ratio` (i.e. "more strongly coherent") is associated with a
  *larger*, more positive bias AND a *larger* sd -- the opposite of "more coherent =
  better fit." This is a genuine inconsistency worth flagging: the binary
  `coherency_state==1` flag and the continuous `coherency_ratio` do not point the same
  direction once you're already inside the coherent population.
- Quantile bins within the not-coherent (state==0) subset (N=17,020, 3 bins) for
  contrast: O-F rises from -1.140 (low ratio) to -0.437 to +0.047 (high ratio) while sd
  falls from 3.733 to 2.990 to 2.716 -- i.e. within *this* subset, higher ratio *does*
  track a tighter (lower-sd) fit, the reverse pattern from the coherent subset.

## Step 5: within-tile control (the important refinement)

Tiles with `coherency_state==1` obs: 317. Tiles with `coherency_state==0` obs: 496.
**Tiles with BOTH: 315** -- ample for this comparison.

For each of these 315 tiles, computed (coherent mean) - (incoherent mean), paired by
tile, for three different statistics:

**(a) Signed bias: mean(O-F | coherent) - mean(O-F | incoherent), per tile**
- mean = **+2.252 dB**, median = +1.692 dB, sd = 2.925 dB
- only 46/315 tiles (14.6%) have a *negative* diff (i.e. only 14.6% of tiles show
  coherent obs *less* positively biased than incoherent obs in that same tile)
- Wilcoxon signed-rank p = 1.2e-36 -- extremely robust
- Robust to per-tile sample-size cutoffs: mean stays +2.0 to +2.3 dB whether requiring
  >=1, >=3, >=5, or >=10 obs per category per tile.
- **This ~2 dB systematic positive shift for coherent obs, present in essentially every
  tile, is NOT predicted by the flat-Fresnel/coherent-scattering theory and dominates
  the naive pooled bias comparison in Step 4 (+1.041 vs -0.510 dB pooled, i.e. most of
  that 1.55 dB pooled gap is this same-tile bias effect, not a tile-selection
  artifact).**

**(b) Error magnitude ("tightness"): mean(|O-F| | coherent) - mean(|O-F| | incoherent), per tile**
- mean = +0.294 dB (median = -0.142 dB) -- mean and median disagree in sign, a sign the
  mean is skewed by a few large-bias tiles from (a)
- only 175/315 tiles (55.6%) show coherent *tighter* (smaller mean |O-F|) -- barely
  above chance
- Wilcoxon signed-rank p = 0.82 -- **not significant**
- Sample-size robustness check flips sign: mean diff is +0.29 dB at min-N>=1, falling to
  -0.10 dB at min-N>=10 -- i.e. this statistic is NOT robust and should not be treated as
  either supporting or refuting the prediction.
- **Verdict: once tile identity is controlled for, coherent obs are NOT reliably
  more accurate in absolute terms than incoherent obs from the same tile.** The naive
  pooled mean|O-F| gap in Step 4 (state 0: 2.217 dB vs state 1: 1.998 dB) looks like it
  supports the prediction, but that gap does not survive the within-tile control and is
  better explained by which tiles happen to produce coherent returns.

**(c) Spread: sd(O-F | coherent) - sd(O-F | incoherent), per tile** (258 tiles with >=2
obs in both categories; 57 tiles dropped, sd undefined with a single obs)
- mean = **-0.274 dB**, median = -0.271 dB
- 166/258 tiles (64.3%) show coherent obs with *lower* spread than incoherent obs in
  the same tile
- Wilcoxon signed-rank p = 2.7e-5 -- **significant**
- Robust to sample-size cutoffs: mean diff stays -0.26 to -0.31 dB at min-N>=1/3/5/10
- **This is the one statistic that survives the within-tile control cleanly in the
  predicted direction: coherent-flagged obs really do have modestly, robustly lower
  O-F spread (precision) than incoherent obs from the same tile, independent of which
  tiles happen to produce coherent returns.**

## Honest verdict

The pre-registered prediction -- coherent returns fit tighter because the operator is a
flat-Fresnel coherent-scattering model -- is **not simply confirmed or refuted; it
fractures into two different, opposite-sign effects once tile identity is controlled
for**:

1. **A naive pooled comparison (Step 4) looks like it supports the prediction**
   (smaller sd, smaller mean|O-F| for coherent vs. not-coherent), but Step 5 shows this
   is substantially a **tile-selection artifact**: coherent returns come disproportionately
   from a specific, smaller set of tiles (317 of 498, likely smoother/wetter/water-adjacent
   per the brief's own confound warning), and those tiles' *intrinsic* O-F behavior
   differs from the broader tile population regardless of coherency.
2. **What survives the within-tile control**: a small but real and statistically robust
   **precision improvement** for coherent obs (~0.27-0.31 dB lower sd, same tile,
   p=2.7e-5) -- mild, genuine support for the physical prediction on the *spread*
   dimension.
3. **What also survives, unpredicted, and dominates any absolute-error comparison**: a
   **~+2.0-2.3 dB systematic positive bias shift** for coherent obs relative to
   incoherent obs in the *same tile* (p=1.2e-36, present in 85% of qualifying tiles,
   reproduced almost exactly in the independent DA-run secondary check). This bias is
   large enough (nearly the full pooled O-F sd of ~3.2 dB) that it, not precision, is
   the dominant within-tile difference between coherent and incoherent CYGNSS L1
   observations of this operator.
4. The continuous `coherency_ratio` field does not cleanly reinforce even the modest
   precision story: *within* the coherent-flagged population, higher ratio (more
   strongly coherent) associates with *worse* fit (larger bias, larger sd), the opposite
   of what "more coherent" should mean if the categorical result were being driven by a
   single coherent-scattering mechanism whose strength scales with `coherency_ratio`.

Taken together, this is closer to the brief's third possibility than to a clean
confirmation: **the coherency flag is doing *something* real (it predicts a small,
robust precision difference and a much larger, robust bias difference, both surviving
tile-identity control), but it is not simply flagging "the operator's Fresnel model
applies well here."** The dominant within-tile effect (bias, not spread) and the
non-monotonic behavior of `coherency_ratio` within the coherent subset both suggest the
DDM3x5 crop-average scalar and/or the coherency SDS field are picking up a mix of
mechanisms (e.g. surface type / wetness / possibly a genuine coherent-scattering
contribution) rather than isolating a single well-behaved "coherent vs. incoherent
Fresnel-model-fit" axis. This does not by itself explain the broader L1-operator
underperformance found in `cygl1_operator_diagnosis.md`, but it does show that
`coherency_state`/`coherency_ratio`, used naively and pooled, would give a misleadingly
clean-looking "coherent fits better" answer that a tile-identity control shows is not
the real, dominant effect.

## Records dropped / caveats

- 0 OFA rows dropped for bad `tilenum`, 0 for no matching raw candidate in the
  reconstructed window, 0 for missing coherency lookup (OL run and DA run alike).
- 1 legitimate upstream data gap noted (20200131, spacecraft 02 and 04 contributed 0
  rows to that day's merge) -- does not affect the join, just reduces that day's N.
- 57 of 315 within-tile-control tiles excluded from the spread (sd) comparison only,
  because at least one category had a single observation (sd undefined); included in
  the bias and |O-F| comparisons.
- All numbers in this note are from Jan-Feb 2020 in the AZ-limited domain only, per the
  brief's scope; no claim is made about seasons/domains outside that window.
