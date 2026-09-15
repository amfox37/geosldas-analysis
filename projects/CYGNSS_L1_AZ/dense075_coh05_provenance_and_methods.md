# CYGNSS L1 `dense075_coh05`: end-to-end provenance and methods

Working note tracing every figure in
`notebooks/final_dense075_coh05_omf_figures.ipynb` back through the DA
experiment, the observation thinning chain, and the CYGNSS L1 preprocessing
pipeline. Written 2026-09-08, covering the 24-month (Jan 2020 - Dec 2021)
`dense075_coh05` result.

Two repositories are involved:

- `../CYGNSS_operator` (i.e. `/Users/amfox/Desktop/CYGNSS_operator`) - the
  observation operator, QC, and coefficient preprocessing.
- this repository - thinning, DA experiment drivers, postprocessing, figures.

Companion run notes, in the order the work happened:

1. `runs/cygl1_assim_R_sweep.md` - the failure this whole thread responds to.
2. `runs/cygl1_operator_diagnosis.md` - L1-vs-L3 operator diagnosis.
3. `runs/cygl1_coherency_stratification.md` - is the coherency flag useful.
4. `runs/cygl1_thinning_and_localization_summary.md` - why thinning, not
   localization.
5. `runs/weekly_update_2026-08-27.md` - 22-month intermediate arm, hard gate.
6. `runs/cygl1_coh05_density_spectrum_and_24mo_result.md` - the result this
   note documents.

---

## 1. Preprocessing pipeline

The assimilated quantity `CYGNSS_L1_DDM3X5_CROP_SCALAR` is not a soil-moisture
retrieval. It is a raw DDM-derived scalar in dB, with a physics-based forward
operator precomputed into fixed per-tile coefficients. Five stages.

### Stage 1: raw input

CYGNSS Level-1 v3.2 Science Data Record daily granules, `cyg01`-`cyg08`,
`power-brcs` files, read from `$CYGNSS_L1_PATH`
(`<CYGNSS_DATA_ROOT>/cygnss/L1/v3.2` on Discover). Each sample carries four
simultaneous DDM channels; DDM dimensions are `delay=17`, `doppler=11`.

### Stage 2: hard QC screen

Nothing is ranked, selected, or thinned until this passes
(`README.CYGNSSObservationQC.md`, `docs/cygnss_m36_window_thinning_notes.md`):

```text
(quality_flags_2 & 2048) == 0        # poor-land-quality bit
sp_land_valid == 1
sp_land_confidence >= 2
ddm_snr > 2
sp_rx_gain > 1
sp_inc_angle < 65
pekel_sp_water_percentage_5km < 1
pekel_sp_water_flag == 0
sp_in_ddm == 1                       # specular point inside the DDM array
finite location/geometry fields
usable positive full 3x5 BRCS crop around the DDM peak
```

This stage writes the per-(day, spacecraft) QC-pass CSVs at

```text
<CYGNSS_operator>/artifacts/out_images/cygnss_qc_m36_window_counts_<YYYYMMDD>_cyg<NN>/
    cygnss_l1_qc_pass_<YYYYMMDD>_cyg<NN>.csv
```

**These CSVs are load-bearing for the coherency filter in Section 2.**
`coherency_state` and `coherency_ratio` are dropped between the raw SDS read
and the daily staging file, and survive nowhere else in the chain - verified by
`ncdump -h` on both the staging files and the upstream per-spacecraft
preprocessor group files (`runs/cygl1_coherency_stratification.md`, Step 0).

### Stage 3: the observable

```text
brcs(delay,doppler) -> reflectivity(delay,doppler) -> fixed 3x5 crop about the
DDM peak -> summed to one scalar -> observed_y_db
```

`reflectivity_peak` and the NBRCS variables are kept as diagnostics only; they
are not the assimilated scalar.

### Stage 4: forward operator and its tile-coefficient factorization

From `docs/cygnss_operator_forward_model_explainer.md` and
`README.CYGNSSTileCoefficientOperator.md`. The IGOT DDM integral is

```text
<sigma(i,j)> = INT_Sigma0  Gamma_f(m_v, T, rho) * g(rho) * <|chi(dtau_i, df_j)|^2> drho  +  nu
```

where `g(rho)` is the DEM-derived topographic scattering factor (Copernicus
DEM, weighted-least-squares plane fit via the compiled `dem_gradient`/FFTW
extension) and `chi` is the Woodward Ambiguity Function kernel for the selected
3x5 bins.

The central approximation: `Gamma_f` is the only soil-moisture-dependent term
and varies smoothly over a narrow linear range, while the geometric terms span
many orders of magnitude but are fixed for a given observation geometry. Assume
soil moisture constant within a GEOSldas tile and pull it out of the per-tile
integral:

```text
H(x) = sum_t C_t * R_t(x),     C_t = INT_{pixels in t} g(rho) * |chi|^2 drho
```

`scripts/preprocess_cygnss_coefficients.py` computes and freezes the `C_t`
(tile-space vegetation opacity folded in, default mode `tile-opacity-sp`),
writing one ragged NetCDF per observation group with `tile_start`,
`tile_count`, `tile_ig`, `tile_jg`, `coefficient` (product schema 0.5).

Runtime evaluation needs only:

```text
sfmc_clipped = min(max(SFMC_t, 1e-4), MWRTM_POROS_t)
epsilon_t    = Mironov(1.57542 GHz, sfmc_clipped, MWRTM_CLAY_t)
R_t          = LR Fresnel reflectivity(epsilon_t, incidence angle)
H(x)         = 10 * log10( sum_t C_t * R_t )
```

Two properties matter downstream:

- `R_t` applies **no roughness and no vegetation correction** - deliberate, per
  the code's own comment in `mwRTM_get_lr_reflectivity`. It is a
  coherent-scattering formulation, which is what motivated the coherency
  hypothesis in Section 2.
- The tile contract is **global EASE `ig`/`jg`**, never positional tile
  numbers. `sp_nearest_tile_index0`, `tile_index0`, `owner_tilenum` and friends
  are only valid for the exact domain and tile ordering used when the
  preprocessor ran (`docs/geosldas_cygnss_preprocessed_operator_handoff.md`).

### Stage 5: daily M36 staging and the GEOSldas read

Per-spacecraft products are merged into one daily file:

```text
/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/Y<YYYY>/M<MM>/
    cygnss_l1_ddm3x5_crop_scalar_m36_<YYYYMMDD>_all_cyg.nc4
```

That path and filename template are literally the `obsparam` `path`/`name`
entries for species 13 (see
`OLv8_M36_all_sensors_AZ_scaled_describe/OLv8_M36_all_sensors_AZ_scaled.ldas_obsparam.txt`).
The sample month staged locally in `example_obs/M06/` is a copy of it, and is
what the notebook's example-observation figures read.

`read_obs_cygnss_l1_scalar()` in `clsm_ensupd_read_obs.F90` then:

1. assigns each obs an owner tile by matching `sp_nearest_tile_ig/jg` against
   the experiment's local `tile_coord`;
2. bins into the centered 3-hour window `(t - dtstep/2, t + dtstep/2]` with
   `dtstep_assim = 10800 s`;
3. keeps only the candidate with the smallest `sp_nearest_tile_distance_km`
   per (owner tile, window), discarding the rest (there is an `N_duplicate`
   counter in the log).

Duplicate rate: 14.95% over the full raw candidate pool, but only 0.02% (5 of
20,731) among obs that actually reach this experiment's AZ-domain OFA output.

Finally `scale_obs_cygl1scal_zscore` z-scores the obs against the pentad
climatology `AZ_CYGNSS_L1_zscore_all_pentads`, built from the OL run's own
stats directory. This also rescales the observation error, which is why
`obsvar` on disk ranges 0.118-173.96 with ~4,362 distinct values despite a
scalar `errstd` in the obsparam metadata.

---

## 2. Experiment setup

### Domain and species

AZ box on EASEv2 M36, **909 tiles**, roughly 117.4-106.2 W and 29.3-39.6 N.
Thirteen species; only CYGNSS L1 is assimilated, and the other twelve are
monitored and act as the independent verification set.

| stats index (0-based) | OFA `species` ID | species | units |
|---:|---:|---|---|
| 0-3 | 1-4 | SMOS Tbh/Tbv asc/desc | K |
| 4-7 | 5-8 | SMAP L1C Tbh/Tbv asc/desc | K |
| 8-10 | 9-11 | ASCAT H-SAF MetOp-A/B/C | m3/m3 (fcst) |
| 11 | 12 | `CYGNSS_SM_6hr` (L3) | m3/m3 |
| 12 | 13 | `CYGNSS_L1_DDM3X5_CROP_SCALAR` | dB |

The notebook's `GROUPS` tuple uses the 0-based stats-array index; its OFA
reader filters on `species == 13`. Both are correct and the offset is a real
trap - the same species is index 56 in the generic full-domain obsparam
(`runs/cygl1_operator_diagnosis.md`). Always resolve species by name from each
file's own `obsparam_descr`/`obsparam_species_id`, never by a hardcoded index.

### Why this is a thinning experiment

**The R-sweep failure (2026-08-20, `runs/cygl1_assim_R_sweep.md`).**
Full-density CYGNSS L1 assimilation at `errstd` 4.4 / 2.2 / 1.1 dB against a
scaled OL, calendar 2020. Change in O-F standard deviation vs OL:

| group | full R | half R | quarter R |
|---|---:|---:|---:|
| SMOS | +3.9% | +6.8% | +10.5% |
| SMAP | +4.5% | +8.0% | +12.6% |
| ASCAT | +5.8% | +13.3% | +23.1% |
| CYGNSS L3 | +2.1% | +5.9% | +11.2% |

Nothing improved at any R, and the degradation scaled cleanly with gain. But
O-A relative to O-F was monotonic, correctly signed, and scaling properly:
**the EnKF was working; the observation was the problem.** The operator
diagnosis quantified it - obs-forecast correlation r ~ 0.35 for L1 against
0.64-0.93 for every other sensor, so its increments are correctly-sized noise.
Per-tile the picture is better (median Spearman rho ~ 0.48, up to 0.92); it is
the tile-to-tile dispersion in offset and scale that dominates the pooled
scatter, consistent with per-tile static ancillary inputs rather than an
intrinsically unusable observable.

**Localization tightening was tested and rejected**
(`runs/cygl1_thinning_and_localization_summary.md`). `check_compact()`
hard-aborts if `xcompact < 2 x` the largest relevant correlation length. The
binding constraint is not CYGNSS L1's own obs-error correlation length - it is
`xcorr_force_pert%pcp/sw/lw = 0.5 deg`, flooring `xcompact` at 1.00 deg against
a current 1.25 deg. Sweeping the entire legal range moves the local-interaction
count only from median 60 to 42. Thinning's target was median 1, a ~60x
reduction. Lowering the floor further would mean changing the forcing
perturbation ensemble for every species project-wide - a materially different,
higher-blast-radius experiment. **Thinning is the only practical lever.**

### The thinning ladder

`scripts/thin_cygl1_nested_density_6mo.py` builds strictly **nested** tiers, so
that any change is attributable to added observations rather than a different
sample. Obs are binned by true centered 3-hour window anchor
(`round(t / 10800) * 10800`), then greedily admitted in fixed order if at least
`min_sep_deg` from every obs already kept in that same window.

| tier | rule | character |
|---|---|---|
| `sparse` | min_sep = 5.0 deg | one isolated obs per window; outside every other obs's 1.25 deg GC ellipse |
| `intermediate` | force-keep sparse, then min_sep = 2.40 deg | limited overlap; local-interaction median ~1 |
| **`dense075`** | force-keep intermediate, then min_sep = **0.75 deg** | `scripts/build_cygl1_dense075_thinning.py` |
| `dense` | unthinned | full stream |

`xcompact = ycompact = 1.25 deg` is the actual localization radius used
throughout (corrected 2026-08-21 from an earlier, wrong 2.5 deg figure in
project memory - the `sparse` 5.0 deg floor was chosen against the old number
but remains valid since 5.0 > 2 x 1.25).

The `dense075` tier was calibrated 2026-09-06 because the gap between
`intermediate` and full density was too wide to interpret:

| min_sep (deg) | N (pre-coherency) | vs intermediate | median NN dist | local-interaction med/p90 |
|---:|---:|---:|---:|---:|
| 1.25 | 12,756 | 2.4x | 1.336 deg | 4/7 |
| 1.00 | 17,842 | 3.4x | 1.091 deg | 7/11 |
| **0.75 (chosen)** | **26,340** | **5.0x** | **0.824 deg (~1.8-2.3 M36 cells)** | **12/18** |
| 0.60 | 34,284 | 6.5x | 0.681 deg | 16/25 |
| 0.50 | 42,044 | 8.0x | 0.578 deg | 21/32 |

`write_thinned_files()` rewrites the daily staging files with
`tile_start`/`tile_count`/support arrays remapped, and writes only the dates
passed in. That is how the record was extended in three non-overlapping passes
(Jan-Jun 2020, Jul-Dec 2020, Jan-Dec 2021) without touching earlier output.
`build_cygl1_dense075_thinning.py` rebuilds from the ORIGINAL unthinned stream
(not from already-thinned files) so the support remapping stays correct.

### The coherency filter layered on top

`scripts/filter_cygl1_by_coherency.py` keeps only `coherency_ratio >= 0.5`,
joining back to the Stage-2 QC CSVs on `(sc_num, sample_id, ch_id)`. Obs whose
join fails are excluded. The join must go through the CSVs directly - not
through `build_cygl1_coherency_screening_experiment.py`'s cached
`file_idx`/`obs_idx` log, which is keyed to its own source root and date range
and does not line up with a different tier's row indices.

Motivation is physical and was pre-registered
(`runs/cygl1_coherency_stratification.md`): since `R_t` is flat Fresnel with no
roughness or vegetation term, coherent returns should fit better. The honest
finding was mixed:

- A naive pooled comparison supports the prediction, but a within-tile control
  (315 tiles carrying both classes) shows most of it is tile selection.
- Absolute error does **not** survive the control: mean |O-F| difference
  Wilcoxon p = 0.82, and the sign flips with per-tile sample cutoff.
- **Spread does survive cleanly**: sd(O-F) difference = -0.274 dB, 166/258
  tiles (64.3%), p = 2.7e-5, robust at min-N >= 1/3/5/10.
- An unpredicted systematic **+2.25 dB positive bias** for coherent obs shows
  up in essentially every tile (Wilcoxon p = 1.2e-36) and dominates any naive
  absolute-error comparison.

So the filter went in as a precision screen, not an accuracy one.

Three arms were run as a density spectrum with the filter fixed at 0.5:

| arm | thinning | N kept (Jan-Jun 2020) |
|---|---|---:|
| `intermediate_coh05` | min_sep 2.40 deg + coh >= 0.5 | 3,940 |
| **`dense075_coh05`** | min_sep 0.75 deg + coh >= 0.5 | **19,802** |
| `dense_coh05` | none + coh >= 0.5 | 84,567 |

### Run configuration

- Experiment `DAv8_M36_AZ_paired_cygl1_dense075_coh05`; paired open loop
  `OLv8_M36_AZ_paired_monitor`.
- `errstd = 2.75 dB`; `xcorr = ycorr = 0.625 deg`;
  `xcompact = ycompact = 1.25 deg`.
- The 0.625 deg is **not** a free choice. The first submission at 1.25 deg (the
  template default) crashed both jobs immediately with
  `LDAS ERROR (3000) from check_compact` - the `>= 2x` rule is enforced **per
  species**, and with `xcompact` held at the project-standard 1.25 deg shared
  with every other arm, 0.625 deg is exactly the ceiling.
- Ungated binary (a separate pre-`406206a` build, surgically relinked from the
  gated binary's own object files). This arm carries **no hard gate** - one
  intervention per experiment.
- Jan 2020 - Dec 2021, restarted from the global `LS_OLv8_M36_v2` source. An
  Apr-1 restart from `OLv8_M36_all_sensors_AZ`'s own AZ-cropped output hit a
  real `ldas_setup` bug (domain-cropped `.til.domain` mishandled as netCDF), so
  all paired experiments restart from the standard global source and take three
  extra months of spin-up.
- The final 24-month run was redone from BEG_DATE specifically to add the
  `catch_progn_incr` and `inst3_1d_lndfcstana_Nt` output collections, completed
  cleanly 2026-09-08. That is why `example_cat/M06/` exists.

### Comparison method: the OL cross-mask

**The baseline is not the plain OL.**
`scripts/postproc_drivers/run_cygl1_paired_OL_xmask_coh05_24mo.py` runs

```python
ol_main = load_exp('OLv8_M36_AZ_paired_monitor', exptag=exptag)
da_sup  = load_exp(da_expid)
da_sup['use_obs'] = True
postproc_ObsFcstAna([ol_main, da_sup], start_time, exp_end_time, ...)
```

so the OL's O-F is recomputed over **exactly the observation population the DA
arm actually assimilated**. For species 13 this is essential - the DA arm sees
a thinned and coherency-filtered stream the plain OL never saw. For species
1-12 the populations coincide anyway, so the cross-mask is equivalent to the
plain OL there.

This is what `_xmask_` in the filenames means, and it is why the notebook reads
`temporal_stats_OL_paired_monitor_xmask_dense075_coh05_20200101_20211231.nc4`
rather than the unmasked `temporal_stats_OL_paired_monitor_*` file sitting
beside it in the same directory.

### Postprocessing products

The `postproc_ObsFcstAna` toolkit computes monthly sums (skipping months whose
sums file already exists, so extensions append rather than recompute), then:

- `spatial_stats_*.pkl` - 24 months x 13 species, domain-aggregated.
- `temporal_stats_*.nc4` - 909 tiles x 13 species, full-period.

Both carry `N_data`, `OmF_mean`, `OmF_stdv`, `O_mean`, `F_mean`. Staged for the
notebook under `output/thinning_expts/`.

---

## 3. Results

`notebooks/final_dense075_coh05_omf_figures.ipynb` reads the four files above
plus `OLv8_M36_all_sensors_AZ_describe/OLv8_M36_all_sensors_AZ.ldas_tilecoord.bin`.
Conventions: `NMIN = 10`, species pooled into five systems by `N_data`-weighted
mean, percent difference always `100 * (DA - OL) / OL` so **negative means DA
fits better**. Map extent (-118.6, -105.2, 28.6, 40.4), latitude sub-mask at
37.5 N.

### Headline

| system | area-wtd map % | monthly-mean % | months improved | tiles improved |
|---|---:|---:|---:|---:|
| **CYGL1** | **-9.84** | **-12.36** | **24 / 24** | **75%** |
| CYGL3 | +0.36 | +0.57 | 4 / 24 | 46% |
| ASCAT | +1.06 | +1.44 | 1 / 24 | 41% |
| SMOS | +0.63 | +0.39 | 6 / 24 | 39% |
| SMAP | +0.67 | +0.71 | 7 / 24 | 37% |

Monthly CYGL1 OmF StDev: **3.020 -> 2.635 dB**, improving in every one of 24
months.

### Reconciling the three numbers in circulation

All three are correct; they are different weightings of the same result. Quote
whichever, but say which:

- **-12.88%** - `N_data`-weighted pooled stdv over the whole 24 months, from
  `scripts/compare_cygl1_coh05_omf.py --period 24mo`. This is the run note's
  headline.
- **-12.36%** - the notebook's mean of the 24 individual monthly percentages.
- **-9.84%** - the notebook's area-weighted mean over the 544 common tiles.
  Tile-space, so low-count tiles pull it toward zero.

### Honest caveat on the monitor species

The monitors are the co-primary metric and are at or near the noise floor, but
not exactly zero. ASCAT in particular is **+1.4% monthly with only 1 of 24
months improved** - small, but consistently the wrong sign rather than random
scatter. "Neutral" is defensible at this magnitude; "no impact" is not.

The `OmF_mean` percentages for CYGL1 (-208%) are numerically meaningless: the
OL mean is +0.090 dB and the DA mean is -0.114 dB, so the denominator is near
zero. The physically real statement is a **bias sign flip of about 0.20 dB**.

### Figure inventory

All under `output/thinning_expts/figures/`.

| output | content |
|---|---|
| `dense075_coh05_omf_stdv_maps_5x3.png` | Primary figure. 5 systems x (OL, DA, % diff). Shared per-row viridis for columns 1-2; shared symmetric RdBu_r capped at +/-30% for column 3. CYGL1's diff panel is broadly blue; the four monitor rows are near-white. |
| `dense075_coh05_omf_stdv_percent_map_summary.csv` | Per-system tile counts, tile mean/median %, area-weighted % (full and lat < 37.5), fraction of tiles improved. |
| `dense075_coh05_monthly_omf_stdv_5x2.png` | 24-month series, OL dashed grey vs DA colored, with a percent panel beside each. This is the durability evidence: no drift, no seasonal blow-up. |
| `dense075_coh05_monthly_omf_mean_5x2.png` | Same layout for OmF mean. |
| `dense075_coh05_nobs_omean_fmean_maps_5x3.png` | Observation support. CYGL1 gets 109 obs/tile over 24 months across 545 tiles, vs CYGL3 577, SMAP 272, ASCAT 289. `O_mean` -13.72 dB vs `F_mean` -13.80 dB confirms the z-score scaling is working. Nobs is from the DA run; obs/forecast means are from the paired OL. |
| `dense075_coh05_example_obs_june2020_maps_2x2.png` | Count, mean `observed_y_db`, mean DDM SNR, mean incidence angle, 0.25 deg binned from `example_obs/M06/`. |
| `dense075_coh05_example_obs_assim_windows_<date>_2x4.png` | Specular points by 3-hour window. Quietly informative: on 2020-06-01 the AZ box sees obs in only three of eight windows (03z n=10, 06z n=55, 09z n=52; five empty). CYGNSS's low-inclination sampling is bursty, not uniform. |
| `dense075_coh05_example_obs_coefficients_20200601_06z_*.png` | A single observation's normalized `C_t` field on the true EASEv2 M36 grid via `pcolormesh`, plus an all-obs overlay. |
| `dense075_coh05_example_slide_{obs,omf,ana,sfmc_ana_minus_fcst,rzmc_ana_minus_fcst,srfexc_incr,rzexc_incr}_20200601_06z.png` | Single-window walkthrough joining OFA rows to staging-file footprints by exact specular-point match (1e-5 deg tolerance, greedy, no reuse), alongside the `catch_progn_incr` state increments. |

### Where this sits against every other intervention tried

| approach | CYGL1 own skill vs OL | monitor Tb/SM |
|---|---:|---|
| thinning alone (`intermediate`, 22mo) | -9.2% | neutral |
| **thinning + coherency (`dense075_coh05`, 24mo)** | **-12.9%** | **neutral** |
| hard gate, full density (`dense_gated`, 22mo) | -0.4% | neutral (harm removed, no skill) |
| coherency alone, full density (`dense_coh05`, 6mo) | -3.3% | **Tb +16%, SM +23%** |

The last row is the load-bearing control: coherency screening at full density
still does large cross-species damage. **Density reduction, not per-obs quality
screening, is what prevents the harm.** The filter adds ~3.7 percentage points
of own-skill on top of thinning; it does not substitute for it.

Stability with record length: CYGL1 -8.71% (6mo) -> -11.80% (12mo) -> -12.88%
(24mo), monitors within +/-1.2% throughout. Not a short-window artifact.
`intermediate_coh05` was deliberately stopped at 12 months and essentially
matches (-13.18%, Tb -0.10%, SM -0.23%).

---

## 4. Caveats to carry onto any slide

1. **The unresolved thinning discrepancy.** 45-57% of same-window
   nearest-neighbour pairs in the `intermediate` output are closer than its own
   nominal 2.40 deg admission floor; the closest pair is 0.029 deg (~3 km) with
   timestamps 0.5 s apart. The greedy algorithm as written makes that
   impossible, and root cause is still unidentified. `dense075` is built by
   re-running the same `build_intermediate_candidate()` calls, so it inherits
   whatever this is. The paired result stands on its own evidence, but the
   "clean isolation" mechanism story does not.
2. One region (AZ), one 24-month window, one seed.
3. **Sign conventions differ between figure families.** The older thinning
   figures use `100 * (OL - DA) / OL` (positive = improvement); this notebook
   and the coherency figures use `100 * (DA - OL) / OL` (negative =
   improvement). Never compare colors or signs across them without checking the
   axis label.
4. Two divergent copies of the shared `postproc_ObsFcstAna` toolkit exist in
   this environment. These stats were built with the git-unregistered
   `hsaf_cdr_test` worktree copy, not the "official" `dnb34` checkout - that
   one is older, lacks NC4 ObsFcstAna support, and crashes. Unreconciled.
5. `dense_coh05` was never extended past 6 months, and is not planned to be.

---

## 5. Source index

### This repository

| path | role |
|---|---|
| `notebooks/final_dense075_coh05_omf_figures.ipynb` | the figures this note traces |
| `scripts/thin_cygl1_nested_density_6mo.py` | sparse/intermediate tiers, `write_thinned_files()` |
| `scripts/build_cygl1_dense075_thinning.py` | the 0.75 deg tier |
| `scripts/filter_cygl1_by_coherency.py` | `coherency_ratio >= 0.5`, any tier/date range |
| `scripts/compare_cygl1_coh05_omf.py` | the `--period {6mo,12mo,24mo}` comparison table |
| `scripts/postproc_drivers/run_cygl1_paired_density_coh05_24mo.py` | DA-arm stats |
| `scripts/postproc_drivers/run_cygl1_paired_OL_xmask_coh05_24mo.py` | OL cross-mask stats |
| `scripts/cygl1_localization_radius_sweep.py` | the localization dead-end check |
| `scripts/extract_cygl1_coherency_join.py` | OFA-to-QC-CSV exact-key join |
| `scripts/analyze_cygl1_coherency_stratification.py` | within-tile control |
| `OLv8_M36_all_sensors_AZ_scaled_describe/*.ldas_obsparam.txt` | species 13 config of record |
| `OLv8_M36_all_sensors_AZ_describe/*.ldas_tilecoord.bin` | tile lon/lat/area for all maps |
| `example_obs/M06/`, `example_ofa/M06/`, `example_cat/M06/` | June 2020 sample data for the example figures |

### `../CYGNSS_operator`

| path | role |
|---|---|
| `README.CYGNSSObservationQC.md` | the QC screen, variable inventory |
| `docs/cygnss_m36_window_thinning_notes.md` | window convention, best-obs ranking, QC CSV products |
| `README.CYGNSSCoefficientPreprocessor.md` | `preprocess_cygnss_coefficients.py` operator guide |
| `README.CYGNSSTileCoefficientOperator.md` | the `H(x) = sum C_t R_t` factorization |
| `docs/cygnss_operator_forward_model_explainer.md` | forward model walkthrough with figures |
| `docs/geosldas_cygnss_preprocessed_operator_handoff.md` | schema 0.5, the `ig`/`jg` tile contract |
| `docs/obs_error_variance_diagnosis.md` | R vs `fcstvar` decomposition, SM-equivalent error |
| `docs/discover_cygnss_preprocessor_quickstart.md` | Discover environment setup |

### Data paths (Discover, gitignored)

```text
staging (full stream) : /discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1/
thinned tiers         : .../cygl1_operator_test/CYGNSS_L1_thinned_{sparse,intermediate,dense075}_6mo/
coherency-filtered    : .../cygl1_operator_test/CYGNSS_L1_thinned_{intermediate,dense075,dense}_coh05/
QC-pass CSVs          : /gpfsm/dnb06/projects/p284/CYGNSS_operator/artifacts/out_images/
                            cygnss_qc_m36_window_counts_<date>_cyg<NN>/
stats output          : output/postproc_paired_density/stats_output/
```
