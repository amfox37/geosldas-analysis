---
name: cygl1-operator-diagnosis
description: "L1-vs-L3 observation-error diagnosis: intrinsic noise (A) vs crude forward operator (B). Tasks 2-4 of the 2026-08-20 brief."
metadata:
  node_type: memory
  type: project
---

# CYGNSS L1 vs L3: is the O-F degradation intrinsic noise or a crude operator?

**Provenance note (open gap, reported as required by Section 9):** the brief states this
analysis "follows the CYGNSS L1 observation-error sweep analysed in
`projects/CYGNSS_L1_AZ/runs/cygl1_assim_R_sweep.md`." That file **does not exist** — searched
the full `CYGNSS_L1_AZ` project tree (`runs/` contains only `OLv8_M36_all_sensors_AZ.md`) and
broader filesystem locations reachable from this session; no match anywhere. The R-sweep
experiments themselves (`DAv8_M36_all_sensors_AZ_scaled_cygl1assim` + `_halfR` + `_quarterR`)
do exist on disk under `/discover/nobackup/projects/land_da/cygl1_operator_test/` and their
`errstd` values match the brief's stated 4.4/2.2/1.1 dB exactly (see Task 4 below), so the
*sweep itself* is real — only the named write-up file is missing. Proceeding on the strength of
the on-disk experiment configs and the independently-reproduced Section 8 numbers below, per
the brief's own instruction to not substitute a different run to paper over a gap, but simply
to flag it.

## Which open-loop run was used, and why

**Used:** `OLv8_M36_all_sensors_AZ_scaled` at
`/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_all_sensors_AZ_scaled`.

Resolution chain: the brief's candidate tag "OL_scaled_baseline" was traced to
`/gpfsm/dnb06/projects/p284/hsaf_cdr_test/omf_compare_sums/OL_scaled_baseline_allspecies/OL_scaled_baseline_allspecies_config.yaml`,
which records `expdir: /discover/nobackup/projects/land_da/cygl1_operator_test/`,
`expid: OLv8_M36_all_sensors_AZ_scaled`. All 13 species in this run's obsparam have `assim=.false.`
(confirmed directly from `output/SMAP_EASEv2_M36_GLOBAL/rc_out/Y2020/M01/OLv8_M36_all_sensors_AZ_scaled.ldas_obsparam.20200101_0000z.txt`)
— a genuine open loop, and `assim_flag==0` for all 368,992 extracted rows confirms this held
throughout 2020. This choice was then validated empirically: it reproduces the Section 8 target
numbers (below) to <0.1%, which a wrong run/species pairing would not do by chance.

A sibling unscaled run, `OLv8_M36_all_sensors_AZ` (same project dir, no `_scaled` suffix), also
exists with 2020 OFA output present but was **not** used for Tasks 1-3 — its obsparam has
`scale=.false.` for both CYGNSS species (see Task 4), which does not match the Section 8 target
numbers' implied correlations, and scaled O to F removes the climatology-matching that the
targets assume.

## Section 8 sanity-check gate: PASSED

Computed as the N_data-weighted mean of per-tile (population, ddof=0) O_stdv/F_stdv/OmF_stdv,
keeping tiles with N≥20 — this convention reproduces the targets; a pooled/grand-stdv convention
does not (e.g. gives L1 O_stdv=5.30, clearly wrong). Script:
`scripts/sanity_check_temporal_stats.py`.

| Species | Metric | Target | Computed | Diff |
|---|---|---|---|---|
| CYGNSS_L1_DDM3X5_CROP_SCALAR | O_stdv | 2.701 | 2.7007 | 0.01% |
| | F_stdv | 2.738 | 2.7375 | 0.02% |
| | OmF_stdv | 3.095 | 3.0948 | 0.01% |
| | N_data | 171483 | 171485 | 2 rows |
| | implied corr | 0.352 | 0.3523 | 0.09% |
| CYGNSS_SM_6hr | O_stdv | 0.026 | 0.0263 | 1.2% |
| | F_stdv | 0.029 | 0.0289 | 0.3% |
| | OmF_stdv | 0.023 | 0.0230 | 0.0% |
| | N_data | 197211 | 197333 | 122 rows |
| | implied corr | 0.658 | 0.6575 | 0.08% |

562/573 (L1) and 673/677 (L3) tiles had N≥20 and were kept. Differences are within reproduction
tolerance (independent re-extraction, rounding in the brief's 3-decimal targets, minor tile-count
differences). **Gate passed — proceeding to Tasks 1-4 as instructed.**

## Task 1: extraction

Script: `scripts/extract_cygl1_l3_pairs.py`. Species resolved by name from each OFA file's own
`obsparam_descr`/`obsparam_species_id` variables (never a hardcoded index) — resolved as
`CYGNSS_L1_DDM3X5_CROP_SCALAR -> 13`, `CYGNSS_SM_6hr -> 12` in this run's 13-species AZ obsparam
(differs from the generic full-domain indices 56/54 documented elsewhere — confirms the
index-instability warning). `tile_id`/`lon`/`lat` taken from the binary tilecoord file by row
index (`tilenum - 1`), not from OFA lon/lat. 2,927 OFA files found for 2020 (~1 short of the
theoretical 2,928 = 366×8); 368,992 total rows assembled, 0 dropped for bad tile index.

**Deliverable:** `output/cygnss_l1_l3_paired_ofa_2020.csv.gz` (10,570,411 bytes; 368,992 data
rows; columns `datetime,tile_id,lon,lat,species,species_name,obs,fcst,ana,assim_flag`), with
`#`-prefixed provenance header (experiment path/id, OFA filename pattern, tilecoord file,
resolved species mapping, code commit, lon/lat caveat) and a `.sha256` checksum file
(`cygnss_l1_l3_paired_ofa_2020.csv.gz.sha256`). Row counts: CYGNSS_L1_DDM3X5_CROP_SCALAR
171,612; CYGNSS_SM_6hr 197,380.

## Task 2: is the transfer function tight (A) or a broad cloud (B)?

Script: `scripts/task2_task3_analysis.py`. Forecast-only merge on matched `(tile_id, datetime)`
with non-null `fcst` for both species: **53,272 matched tile-times across 566 tiles**.

**Pooled result: a broad, essentially flat cloud, not a tight curve.**
- Pooled Pearson r = 0.0576 (R² = 0.0033); pooled Spearman ρ = 0.1123.
- Pooled linear fit F(L1) ~ F(L3): slope = 7.45 dB per m³/m³, intercept = -12.75, R² = 0.0033.
- See `output/figures/task2_fig1_pooled_hexbin_F_L1_vs_F_L3.png`: F(L1) spans roughly 0 to -70 dB
  across the entire F(L3) range (0.02-0.30 m³/m³) with no visible pooled trend.

**But per-tile, the relationship is much stronger** (`output/task2_per_tile_spearman.csv`, 545
tiles with N≥10 matched):
- Per-tile Spearman ρ: mean 0.4733, std 0.2028; percentiles 10%=0.2187, 25%=0.3461,
  median=0.4809, 75%=0.6136, 90%=0.7333, max=0.9235.
- Median per-tile Pearson R² = 0.1658.
- 47.0% of tiles have ρ > 0.5; only 2.0% have ρ < 0 (i.e. a negative sign is rare).
- `output/figures/task2_fig2_individual_tiles_scatter.png` (top-6 highest-N tiles) shows
  clearly visible, mostly monotonic increasing (less-negative-dB with higher SM) relationships
  at individual tiles — e.g. tile 6037: rho=0.82, dB range -11 to -3 over SM 0.05-0.23; tile
  5596: rho=0.79 (mostly flat ~-10 dB, with one large outlier at -54 dB); tiles 5707/5708/6159/
  6045: visible but noisier upward trends, rho 0.47-0.61.
- `output/figures/task2_fig3_spearman_distribution.png` — histogram, median marked.

**Interpretation:** this is the pattern the brief calls out as most informative for
distinguishing (A) from (B): F(L1) is **not** near-constant (which would be the most damning,
cleanly-fixable operator symptom) and it is **not** simply unrelated to F(L3) (which would argue
for pure intrinsic noise). Instead, each tile has its own real, fairly strong, monotonic
SM-to-dB relationship (median per-tile ρ≈0.48, up to 0.92), but the **mapping's offset/scale
differs substantially tile-to-tile** — consistent with per-tile static ancillary/coefficient
inputs (footprint geometry, clay/porosity, and especially the statically-baked vegetation
opacity factor — see Task 4) dominating the pooled scatter's spread. This is a signature that
leans toward **explanation (B)** (a real but crudely-calibrated per-tile operator, not simply an
intrinsically noisy observable) — a tight universal SM-dB curve is not physically expected given
the nonlinear Fresnel/Mironov relationship (see Task 4), but the *tile-to-tile* dispersion in
offset is a plausible instrument-agnostic operator-calibration artifact, not something obviously
attributable to observation noise alone.

## Task 3: does the raw observable itself carry the signal?

Matched observation pairs (both species have an `obs` in the same assimilation window):
**53,272 pairs across 566 tiles** — identical count to Task 2's forecast-merge, a structural
fact of this OFA logging (fcst and obs co-occur per observation event, so the obs-merge and
fcst-merge select the same rows). 545 tiles had ≥10 pairs.

`output/task3_per_tile_correlations.csv`, `output/figures/task3_fig1_correlation_comparison.png`:

| Comparison | mean | median | 10% | 90% |
|---|---|---|---|---|
| O(L1) vs O(L3), same tile-time | 0.4608 | 0.4681 | 0.2201 | 0.7175 |
| O(L1) vs F(L1), same samples | 0.3867 | 0.4154 | 0.0822 | 0.6416 |
| O(L3) vs F(L3), same samples | 0.7056 | 0.7302 | 0.5551 | 0.8416 |

**Interpretation:** the raw L1 observable's cross-correlation with the independent L3 SM
retrieval (median 0.468) is comparable to — in fact very slightly *higher* than — L1's own
correlation with its model forecast on the identical sample set (median 0.415). That is: **the
CYGNSS L1 observable does carry real soil-moisture information**, about as much of it as its own
forward model captures, and that information is at least as visible against an independent
retrieval as against the forecast it's meant to be compared to. Meanwhile L3's own O-F coherence
on the same tiles/times is far higher (median 0.730). This pattern says the *deficiency is not
primarily "the L1 signal is buried in noise relative to what a good operator could extract"* —
if it were, O(L1)-O(L3) would sit well *below* O(L1)-F(L1) (a noisy proxy should correlate worse
with an independent product than with the very model it's regressed against). Instead they're
about equal, meaning the operator is capturing roughly the same fraction of the true signal that
an independent instrument sees, but that fraction/ceiling itself (~0.4-0.5, vs L3's ~0.7) is
capped well below L3's — consistent with a real, moderately informative raw observable, whose
current forward-model realization (Task 2/4) has room to improve, rather than the observable
being fundamentally noise. Net: evidence leans toward **explanation (B)** — a crude but
fixable operator — over (A) — irreducible observational noise — though it does not rule out (A)
contributing as well (median cross-instrument correlation of 0.47 is real but far from "clean").

## Task 4: configuration confirmation

**Scaling-path configuration — the brief's hypothesis is confirmed, but on the *sibling*
unscaled run, not the run used for Tasks 1-3:**

- `OLv8_M36_all_sensors_AZ_scaled` (used above): both species have `scale=T`. CYGNSS_SM_6hr
  scalepath=`.../OLv8_M36_all_sensors_AZ/output/SMAP_EASEv2_M36_GLOBAL/stats/python_z_score_clim_quarter_degree`,
  scalename=`AZ_CYGNSS_SM_zscore_2020_doy1_2022_doy365_W_75d_Nmin_20_sp_ALL_all_pentads`,
  errstd=0.04 m³/m³. CYGNSS_L1 scalepath=`.../OLv8_M36_all_sensors_AZ/output/SMAP_EASEv2_M36_GLOBAL/stats/cygnss_l1_z_score_clim`,
  scalename=`AZ_CYGNSS_L1_zscore_all_pentads`, errstd=3.0 dB. Both non-empty — this run does
  **not** match the "CYGNSS_L1 has empty scaling paths" hypothesis.
- `OLv8_M36_all_sensors_AZ` (unscaled sibling, not used for Tasks 1-3): both species have
  `scale=F`. CYGNSS_SM_6hr scalepath=`/discover/nobackup/projects/land_da/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/stats`,
  scalename=`OLv8_M36_cd_CYGNSS_zscore_stats_all_pentads` — an exact match to the brief's
  hypothesized L3 scaling path. CYGNSS_L1 has `scalepath=''`, `scalename=''` — an exact match to
  the brief's hypothesized "L1 has empty scaling paths." So the brief's Task-4 hypothesis
  describes the unscaled run's obsparam precisely; it does not describe the scaled run
  (validated against Section 8) used for the rest of this analysis. Reporting both rather than
  picking one, per instructions not to paper over discrepancies.

**R-sweep errstd — confirmed, at obs_param index 13 in every case (not assumed):**
`DAv8_M36_all_sensors_AZ_scaled_cygl1assim` (full R): errstd=4.400000 dB. `..._halfR`:
errstd=2.200000 dB. `..._quarterR`: errstd=1.100000 dB. All three have `assim=T, scale=T` for
species 13 (CYGNSS_L1_DDM3X5_CROP_SCALAR), confirmed directly from each experiment's own
`ldas_obsparam.txt` for January 2020, not inferred from the base run's index. This exactly
matches the brief's stated 4.4/2.2/1.1 dB, variance ratio 1 : 1/4 : 1/16.

**Forward operator documentation** (source-traced in the personal checkout,
`/gpfsm/dnb34/amfox/GEOSldas_cygnss_operator/GEOSldas/src/Components/@GEOSldas_GridComp/`):

- Dispatch: `clsm_ensupd_read_obs.F90:10192` (`case ('CYGNSS_L1_DDM3X5_CROP_SCALAR')`) and
  `clsm_ensupd_upd_routines.F90:1730-1740` (dispatched on `varname=='cygl1scal'`).
- Forward operator: `cygnss_preprocessed_obs.F90::cygnss_preproc_get_obs_pred` (lines 491-600).
  Reads **only `sfmc` (surface soil moisture)** as the dynamic model state; ancillary static
  inputs are `mwp_clay`, `mwp_poros` (clay fraction, porosity — `clsm_ensupd_upd_routines.F90:1502-1505`),
  `sp_inc_angle` (incidence angle), and per-support-tile `coefficient(k)` + `freq` (GPS L1,
  1.57542 GHz). Functional form is a coefficient-weighted linear combination of reflectivity
  over DEM-footprint support tiles: `hx_linear = hx_linear + coefficient(k)*refl_lr` (line 588),
  then `obs_pred = 10*log10(hx_linear)` (line 593) if positive.
- Reflectivity itself: `mwRTM_get_lr_reflectivity` (`mwRTM_routines.F90:254-296`) — Mironov soil
  dielectric mixing model (line 282, driven by SM, clay only) into flat Fresnel r_vv/r_hh
  coefficients; SM clipped to `[1e-4, porosity]` (line 280). **An explicit code comment at
  `mwRTM_routines.F90:256-260` states this reflectivity calculation deliberately applies NO
  roughness or vegetation correction.** Vegetation is instead handled upstream and statically,
  baked into the per-tile `coefficient(k)` weights at coefficient-generation time
  (`/gpfsm/dnb06/projects/p284/CYGNSS_operator/cygnss_operator/tile_coefficient_product.py:384-393,401-403`,
  `exp(-weighted_opacity*geometry_factor)` under `vegetation_mode="tile-opacity-sp"`, confirmed
  via the `vegetation_mode` global attribute on the coefficient files) — i.e. vegetation
  attenuation is fixed at whatever climatological opacity was used when coefficients were built,
  not responsive to the actual simulated vegetation state at each timestep. `coefficient(k)` is a
  DEM-pixel-count footprint/glistening-zone area weighting, not a fitted regression.
- Nonlinear, non-uniform sensitivity (from `obs_error_variance_diagnosis.md`, reimplemented and
  independently evaluated with this domain's actual clay≈0.19/porosity≈0.42): SM→dB slope ≈48-53
  dB per m³/m³ near SM=0.05 (dry), falling to ≈9-10 dB per m³/m³ near SM=0.30-0.35 (wet);
  essentially incidence-angle independent (<1% difference across 10°/25.6°/40°). At this domain's
  actual mean SM (≈0.10 m³/m³), the locally relevant slope is ≈35-36 dB/(m³/m³), converting
  errstd 4.4 dB / 2.2 dB to SM-equivalent std ≈0.12 / ≈0.06 m³/m³ — both large relative to
  typical land-DA target accuracy (~0.04 m³/m³) and to this domain's physical SM range
  (~0.03-0.30 m³/m³).
- **DA-divergence red flag** (documented independently in
  `/discover/nobackup/projects/land_da/cygl1_operator_test/obs_error_variance_diagnosis.md`,
  not rerun here per Section 9): the first real CYGNSS-L1-assimilating test (errstd=4.4 dB, job
  57640274) showed O-F and O-A both drifting from ≈+0.5 dB (day 1) to ≈-1.7 dB (day 29), with
  O-A more negative than O-F on almost every day — i.e. each analysis cycle pushed the state
  *away* from matching future obs, not toward it (only 47.3% of individual obs got closer to the
  analysis than the forecast). This is flagged there as a possible sign/gain bug or a
  scaling-climatology invalidation once assimilation begins shifting the model state away from
  the open-loop climatology the z-score scaling was built from — not simply an R-mistuning
  issue. A half-errstd probe (job 57640766) was launched to check whether stronger gain makes the
  drift proportionally worse (pointing to a structural bug) — result not yet available as of this
  writing and out of scope for this diagnosis (no rerun performed here).

## Synthesis: (A) intrinsic noise vs (B) crude operator

The evidence assembled here is not fully dispositive, but on balance leans toward **(B) — a
real, moderately-informative observable currently realized through a crude, non-dynamically-
vegetation/roughness-corrected, statically-weighted forward operator** rather than toward (A)
pure intrinsic noise:

1. Task 2: F(L1) is neither near-constant (ruling out the most damning operator failure mode)
   nor pooled-uncorrelated with F(L3) — it has a real, fairly strong (median ρ≈0.48, up to 0.92)
   per-tile monotonic relationship whose *offset* varies tile-to-tile, exactly the signature of a
   physically real but poorly-calibrated per-tile coefficient/vegetation-baking scheme, not of
   pure sensor noise (which would show no per-tile structure at all).
2. Task 3: the raw observable's correlation with an independent product (O(L1) vs O(L3),
   median 0.468) is comparable to (and marginally better than) its correlation with its own,
   current-generation forward-modeled forecast (O(L1) vs F(L1), median 0.415) — i.e. the operator
   is not yet extracting all the signal the observable demonstrably carries against an
   independent instrument. If (A) were the dominant explanation, cross-instrument coherence
   should sit clearly *below* own-forecast coherence, not roughly matched.
3. Task 4: three operator bugs were already found and fixed during development (documented in
   `runs/OLv8_M36_all_sensors_AZ.md`); the reflectivity model explicitly omits roughness/vegetation
   dynamics by design (not by oversight, but the omission is real and unaddressed); vegetation
   correction is frozen at build-time climatology rather than cycling with the simulated state;
   and the one completed real-assimilation test shows a structural divergence pattern (O-A
   trending away from O-F, not toward it) that is hard to explain by R-mistuning alone.

This does **not** rule out (A) contributing materially — median per-tile correlations of
0.4-0.5 (both L1-L3 and L1-F(L1)) are real but modest, well below L3's ~0.7 O-F ceiling on the
same samples, so some portion of the gap is plausibly genuine roughness/incidence/representativeness
noise in the raw DDM scalar itself, not solely fixable by a better operator. The honest reading is
a mixture, with the balance of evidence in this specific analysis pointing more toward a fixable
operator/calibration problem than toward "the observation should be shelved."

## Files produced

- `output/cygnss_l1_l3_paired_ofa_2020.csv.gz` (+ `.sha256`) — Task 1 primary deliverable.
- `output/_cygl1_l3_pairs_2020.pkl` — scratch convenience copy, not part of the deliverable set.
- `output/task2_per_tile_spearman.csv`, `output/task3_per_tile_correlations.csv` — per-tile
  supporting tables.
- `output/figures/task2_fig1_pooled_hexbin_F_L1_vs_F_L3.png`,
  `task2_fig2_individual_tiles_scatter.png`, `task2_fig3_spearman_distribution.png`,
  `task3_fig1_correlation_comparison.png`.
- `scripts/extract_cygl1_l3_pairs.py`, `scripts/sanity_check_temporal_stats.py`,
  `scripts/task2_task3_analysis.py`.
