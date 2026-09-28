# CYGNSS L1 observation-error correlations: evidence from the fixed-operator open loop

*Arizona domain (MINLON -118, MAXLON -106, MINLAT 29, MAXLAT 40), species `CYGNSS_L1_DDM3X5_CROP_SCALAR`, 2020-01-01 to 2022-12-31. Written 2026-09-26.*

## Summary

- GEOSldas models the CYGNSS L1 observation-error correlation as an isotropic Gaussian in distance, `exp(-0.5 d²/xcorr²)`. The runs used `xcorr = 0.625°` (about 70 km) until now.
- The innovations from the fixed-operator open loop do not support that. Two L1 observations 3–5 km apart have an innovation correlation of only **0.57**. Once the forecast-error part is removed, the implied **observation-error correlation is about 0.5**, falling to about **0.2 at 0.13°** and **below 0.1 beyond 0.25°**. The Gaussian with `xcorr = 0.625°` implies 0.99, 0.96 and 0.89 at those separations.
- The correlated part is mostly **along-track**. At the same distance and footprint overlap, pairs from the same spacecraft and channel within 30 s are about **twice as correlated** as pairs from different spacecraft. Along a track the correlation also has a **long tail**: about 0.24–0.26 in innovation at 0.3–0.6° (5–12 s apart), where pairs from different spacecraft are at 0.09–0.12.
- For pairs from **different spacecraft**, the implied observation-error correlation is about **0.28 below 0.05°**, **0.05–0.08 at 0.05–0.15°**, and **zero beyond 0.15°**. Footprint overlap raises it only modestly. Differences in incidence angle lower it somewhat, but that effect is weak and noisy.
- **Consequence for the DA.** With `xcorr = 0.625°`, dense full-stream L1 obs make R numerically near-singular. The local condition number has a median of 8.8e6, and the analysis blew up: root-zone increments reached 582 mm, about 45% of analyses moved away from the obs, and the independent monitors degraded by 20–40%. Reducing `xcorr` to 0.15° fixed the blow-up. January 2020 then showed Tb monitors −4.75%, SM monitors −0.95% and L1 −0.74% vs the OL.
- **Recommendation.** No single isotropic Gaussian fits the data. A short-range, partly uncorrelated (nugget) structure fits better, with any longer-range correlation confined to along-track pairs. The practical options are in the last section.

## 1. Data and method

### 1.1 Innovations

| Item | Value |
|---|---|
| Run | `OLv8_M36_AZ_fixedop`: unscaled, monitor-only open loop; fixed L1 operator (GEOSldas_GridComp `cdfaaf4`/`11dfdb1`); build md5 `6a70b496`; full `CYGNSS_L1` stream |
| Forecast | ensemble mean (`ens_avg`), 24 members; 3-hourly analysis windows; all 8,768 `ObsFcstAna` files 2020–2022 |
| Obs | raw L1 dB, scaled offline with the new owner-tile z-score climatology (`scaling_params/cygnss_l1_z_score_clim/`, built from this OL). The offline scaling reproduces the model's own scaled values exactly (493/493 obs checked against `DA_L1_full_fixedop`). |
| Innovation | `d = scaled O − F`: mean −0.12 dB, variance 8.0 dB² (std 2.83) over all paired obs |

### 1.2 Pairs

- **Selection:** two L1 obs in the **same analysis window** within **0.6°** of each other, using obs (specular-point) lon/lat. This gives **1,053,216 pairs** over 1,095 days.
- **Never the same owner tile:** 100% of pairs are in different owner tiles, because the reader keeps at most one L1 obs per tile per window.
- **Matching:** each obs was matched to its preprocessed record in `CYGNSS_L1/Y*/M*/*_all_cyg.nc4` by location and raw dB value. The match rate was **100%**.

### 1.3 Pair attributes

- **Footprint overlap:** each obs' sparse footprint is `coefficient_weight` over (`tile_ig`, `tile_jg`), normalized to sum to 1. Overlap = Σ_t min(w₁,t, w₂,t). It is 0 for disjoint footprints and 1 for identical ones.
- **Acquisition type:**
  - *same-track*: same `sc_num` and `ch_id`, |Δt| < 30 s. 39% of pairs.
  - *same-sc other*: same spacecraft but a different channel or later. 4%.
  - *different-sc*: 56%.
- **Other attributes:** |Δt| from `ddm_timestamp_utc_sec`, and the incidence-angle difference from `sp_inc_angle`.

### 1.4 Statistics

- **Correlation:** Pearson correlation of the symmetrized pairs (d₁, d₂) ∪ (d₂, d₁).
- **Confidence intervals:** 90% **day-block bootstrap** CIs (100 resamples of whole days), because pairs within a day are not independent.

### 1.5 Separating obs-error from forecast-error correlation

This follows the Hollingsworth–Lönnberg approach. The innovation correlation at separation d is

`c(d) ≈ a·ρ_o(d) + b·ρ_f(d)`, where a + b = 1 are the obs- and forecast-error shares of innovation variance.

- **Assumption:** pairs from different spacecraft beyond 0.15° have negligible obs-error correlation.
- **Fit:** on those pairs, `b·exp(−d/L)` gives **b = 0.18** and **L = 0.86°**. That is a forecast-error share of 18% and an obs-error share of **a = 0.82**.
- **Cross-check:** the ensemble gives a = 0.875 (forecast spread about 1 dB² out of 8 dB²). The two are close. The ensemble is probably slightly under-dispersive.
- **Implied obs-error correlation:** ρ_o(d) = (c(d) − b·exp(−d/L)) / a.

These implied values depend on the fit's extrapolation to short range. Treat them as approximate, ±0.05 or so.

## 2. Results

### 2.1 Correlation vs distance, all pairs

| Distance (°) | N pairs | Innovation corr [90% CI] | Median overlap | % same-track | Implied obs-error corr |
|---|---|---|---|---|---|
| 0–0.05 | 30,296 | **0.566** [0.542, 0.589] | 0.91 | 95 | ~0.48 |
| 0.05–0.10 | 49,719 | **0.396** [0.374, 0.414] | 0.78 | 87 | ~0.28 |
| 0.10–0.15 | 32,599 | 0.317 [0.288, 0.347] | 0.59 | 62 | ~0.20 |
| 0.15–0.20 | 48,354 | 0.236 [0.217, 0.257] | 0.39 | 60 | ~0.11 |
| 0.20–0.30 | 112,844 | 0.202 [0.186, 0.215] | 0.21 | 43 | ~0.08 |
| 0.30–0.40 | 191,117 | 0.157 [0.144, 0.170] | 0.10 | 38 | ~0.05 |
| 0.40–0.50 | 278,862 | 0.162 [0.149, 0.176] | 0.04 | 33 | ~0.07 |
| 0.50–0.60 | 309,425 | 0.126 [0.118, 0.137] | 0.00 | 26 | ~0.04 |

The residual 0.04–0.07 beyond 0.3° comes entirely from same-track pairs (see 2.2).

### 2.2 Along-track vs other acquisitions

Innovation correlation, with 90% CI:

| Distance (°) | Same-track | Same spacecraft, other | Different spacecraft |
|---|---|---|---|
| 0–0.05 | **0.576** [0.556, 0.606] (N 28,767) | 0.378 [0.21, 0.58] (N 153) | **0.397** [0.323, 0.449] (N 1,376) |
| 0.05–0.10 | **0.435** [0.418, 0.457] | 0.103 [0.05, 0.29] | **0.209** [0.170, 0.253] |
| 0.10–0.15 | 0.385 [0.347, 0.428] | 0.162 | 0.218 [0.193, 0.252] |
| 0.15–0.20 | 0.290 [0.267, 0.322] | 0.183 | 0.144 [0.121, 0.174] |
| 0.20–0.30 | 0.293 [0.265, 0.315] | 0.158 | 0.133 [0.114, 0.150] |
| 0.30–0.40 | 0.240 [0.219, 0.255] | 0.110 | 0.112 [0.093, 0.127] |
| 0.40–0.50 | 0.264 [0.245, 0.281] | 0.099 | 0.117 [0.103, 0.131] |
| 0.50–0.60 | 0.240 [0.217, 0.264] | 0.097 | 0.091 [0.082, 0.102] |

Implied obs-error correlation, using the section 1.5 split:

| Distance (°) | Same-track | Different spacecraft |
|---|---|---|
| 0–0.05 | 0.49 | 0.28 |
| 0.05–0.10 | 0.33 | 0.06 |
| 0.10–0.15 | 0.28 | 0.08 |
| 0.15–0.20 | 0.17 | ~0 |
| 0.20–0.30 | 0.19 | ~0 |

Same-track pairs beyond 0.3° are 0.24–0.26 against 0.09–0.12 for different spacecraft, which implies an along-track obs-error correlation of about 0.12–0.17 still present 5–12 s and 30–60 km apart.

Observations:

- **Same-track pairs dominate the short-range correlation.** They make up 95% of pairs under 0.05° and 87% at 0.05–0.10°.
- **Different spacecraft:** obs-error correlation is modest only at the very shortest range (<0.05°, about 5 km) and is essentially zero beyond about 0.15°.
- **Same spacecraft, different channel or later pass:** behaves like different spacecraft (few pairs, wide CIs). The correlation is tied to the **track**, not the spacecraft.
- **Along-track, the correlation has two scales:**
  - a short one, decaying over about 2–3 s (section 2.3);
  - a long, weak one (~0.12–0.17 obs-error correlation) persisting to at least 12 s / 0.6°.

  The long one is consistent with slowly varying per-track errors such as calibration or antenna-gain drift, or a per-track geometry bias. The z-score scaling removes only the climatological per-tile bias, not per-track offsets.

### 2.3 Along-track decay with time separation

Same-track pairs only:

| \|Δt\| (s) | N | Innovation corr [90% CI] | Median overlap | Median distance (°) |
|---|---|---|---|---|
| ≤1.5 | 65,142 | 0.516 [0.498, 0.535] | 0.87 | 0.062 |
| 1.5–2.5 | 34,976 | 0.347 [0.326, 0.370] | 0.62 | 0.128 |
| 2.5–3.5 | 30,705 | 0.286 [0.256, 0.311] | 0.36 | 0.192 |
| 3.5–5.5 | 73,814 | 0.270 [0.243, 0.298] | 0.16 | 0.294 |
| 5.5–8.5 | 167,371 | 0.252 [0.239, 0.266] | 0.04 | 0.445 |
| 8.5–12.5 | 41,396 | 0.244 [0.214, 0.278] | 0.00 | 0.572 |

The short-range part decays within about 3 s, roughly 0.2° along track, as footprint overlap falls from 0.87 to 0.36. After that the correlation levels off at about 0.25 even when footprints no longer overlap at all (0.00 at 8.5–12.5 s). So the long-range along-track correlation is **not** footprint sharing. It is something the track carries with it.

### 2.4 Footprint overlap at fixed distance

Innovation correlation (N in parentheses) by overlap bin, for distances up to 0.3°.

**Same-track:**

| Distance (°) | 0–0.1 | 0.1–0.3 | 0.3–0.6 | 0.6–0.8 | 0.8–1.0 |
|---|---|---|---|---|---|
| 0–0.05 | – | 0.32 (36) | 0.43 (151) | 0.52 (1,206) | **0.60** (27,361) |
| 0.05–0.10 | 0.26 (103) | 0.19 (388) | 0.39 (5,500) | 0.46 (16,537) | 0.46 (20,554) |
| 0.10–0.15 | 0.55 (180) | 0.23 (1,447) | 0.32 (7,517) | 0.41 (6,060) | 0.48 (5,084) |
| 0.15–0.20 | 0.14 (2,005) | 0.27 (7,862) | 0.29 (11,262) | 0.36 (5,310) | 0.40 (2,759) |
| 0.20–0.30 | 0.31 (14,280) | 0.26 (16,424) | 0.27 (12,838) | 0.34 (3,951) | 0.36 (1,525) |

**Different spacecraft:**

| Distance (°) | 0–0.1 | 0.1–0.3 | 0.3–0.6 | 0.6–0.8 | 0.8–1.0 |
|---|---|---|---|---|---|
| 0–0.05 | – | – | – | 0.28 (304) | **0.41** (1,047) |
| 0.05–0.10 | – | −0.06 (80) | 0.25 (1,406) | 0.20 (2,733) | 0.23 (1,807) |
| 0.10–0.15 | −0.04 (131) | 0.18 (1,293) | 0.20 (5,771) | 0.27 (2,791) | 0.29 (1,276) |
| 0.15–0.20 | 0.02 (1,072) | 0.16 (5,862) | 0.16 (7,552) | 0.18 (2,312) | 0.17 (805) |
| 0.20–0.30 | 0.13 (16,083) | 0.13 (23,318) | 0.14 (15,101) | 0.12 (3,429) | 0.26 (843) |

Observations:

- **Within a track**, correlation rises steadily with overlap at fixed distance. At 0.05–0.15° it goes from about 0.2 to about 0.47.
- **Across spacecraft**, the rise is weaker. Nearly identical footprints (overlap > 0.8) from different spacecraft reach only 0.23–0.41 innovation correlation, an obs-error correlation of about 0.1–0.3.
- Footprint overlap is necessary for a strong correlation but not sufficient. The same footprint seen by a different receiver, at a different time and geometry, carries largely independent error. That fits the preprocessor's per-obs support: different incidence and azimuth, DEM-dependent weighting, and coherent vs incoherent scattering all make the effective sampling differ even where the tiles overlap.

### 2.5 Incidence-angle difference

Different-spacecraft pairs within 0.15°:

| Δθ_inc (°) | N | Innovation corr [90% CI] | Median overlap |
|---|---|---|---|
| 0–3 | 5,131 | 0.301 [0.258, 0.344] | 0.63 |
| 3–6 | 3,372 | 0.229 [0.173, 0.290] | 0.63 |
| 6–10 | 3,087 | 0.201 [0.154, 0.261] | 0.60 |
| 10–15 | 3,035 | 0.231 [0.188, 0.284] | 0.62 |
| 15–25 | 2,825 | 0.148 [0.078, 0.224] | 0.61 |
| >25 | 1,235 | 0.218 [0.147, 0.280] | 0.63 |

Restricted to overlap > 0.6, the correlations are 0.295 for 0–5°, 0.250 for 5–15° and 0.232 above 15°.

Similar geometry (Δθ < 3°) gives the highest correlation, and the overlap is the same in every bin. Beyond that the trend is weak and non-monotonic. **Correction to an earlier quick look:** the first 10-day pass suggested a sharp drop, 0.30 → 0.12 at Δθ > 15°. With 3 years and bootstrap CIs the drop is only about 0.30 → 0.15–0.23. Geometry matters somewhat, but less than track membership.

### 2.6 Seasonality

| Season | d < 0.1°, all | d < 0.1°, same-track | d < 0.1°, different-sc | d 0.3–0.6° |
|---|---|---|---|---|
| DJF | 0.465 | 0.522 | 0.191 | 0.154 |
| MAM | 0.451 | 0.478 | 0.225 | 0.139 |
| JJA | 0.481 | 0.508 | 0.247 | 0.148 |
| SON | 0.457 | 0.481 | 0.276 | 0.145 |

The seasonal variation is small. The different-spacecraft short-range correlation is somewhat higher in the second half of the year, which may relate to vegetation or roughness, but it isn't examined further here.

### 2.7 How the R models compare

These are model-implied innovation correlations, `a·ρ_R(d) + b·exp(−d/L)` with a = 0.82, b = 0.18, L = 0.86°, against the observed values:

| d (°) | Observed | Gaussian `xcorr` 0.625 | 0.25 | 0.15 | 0.10 | 0.15 + nugget 0.5 |
|---|---|---|---|---|---|---|
| 0.032 | 0.566 | 0.992 | 0.987 | 0.975 | 0.953 | 0.574 |
| 0.067 | 0.396 | 0.982 | 0.957 | 0.908 | 0.820 | 0.537 |
| 0.127 | 0.317 | 0.958 | 0.875 | 0.727 | 0.519 | 0.441 |
| 0.174 | 0.236 | 0.936 | 0.790 | 0.564 | 0.326 | 0.356 |
| 0.257 | 0.202 | 0.887 | 0.618 | 0.323 | 0.164 | 0.228 |
| 0.357 | 0.157 | 0.816 | 0.415 | 0.167 | 0.120 | 0.143 |
| 0.450 | 0.162 | 0.740 | 0.269 | 0.116 | 0.107 | 0.111 |
| 0.551 | 0.126 | 0.651 | 0.167 | 0.096 | 0.095 | 0.095 |

"Nugget 0.5" means ρ_R = 0.5·exp(−0.5 d²/xcorr²) for d > 0: half the obs-error variance is uncorrelated. This is the form tested in M21C_testing, `R = σ²[(1−a)C + aI]`.

- **Every pure Gaussian badly overstates the short-range correlation.** Even `xcorr = 0.10°` gives 0.95 at 3–4 km against 0.57 observed.
- The observations need a large **uncorrelated fraction** (about 50% of obs-error variance at d → 0) plus a short correlated part.
- **Nugget 0.5 with 0.15** matches at d < 0.05° and d > 0.25°, but still overstates the 0.05–0.2° range by 0.1–0.14. A shorter correlated scale fits better. An exponential with L ≈ 0.065° and amplitude about 0.8 fits the 0.03–0.07° implied values but under-predicts the 0.13–0.2° tail.
- The tail comes from the along-track long-range component, which **no isotropic function can represent**. It isn't a function of distance at all.

## 3. Why this matters for the DA

### 3.1 Conditioning

`assemble_obs_cov` (`clsm_ensupd_upd_routines.F90`) builds R from the isotropic Gaussian, in degrees, within each species. `enkf_increments` then solves `W = HPHᵀ∘GC + R` in **single precision** using Numerical Recipes `ludcmp`, which has no conditioning check. This is the same mechanism diagnosed for dense H SAF ASCAT in `M21C_testing` (notes section 12 and the CF0360 conditioning doc).

In the full L1 stream, obs in **adjacent owner tiles** are often 0.03–0.10° apart:
- a third of obs have a neighbour within 0.1°;
- the 10th percentile of nearest-neighbour distance is 0.033°.

With `xcorr = 0.625°`, R treats those pairs as 0.99 correlated, while the data give about 0.5 at most. That makes R near-singular. Local sets of obs within the ±1.25° localization box, January 2020:

| Run | cond(R), median / p90 / max | q = dᵀ(R+diag P_f)⁻¹d / n (≈1 if consistent) |
|---|---|---|
| Full stream, `xcorr` 0.625 | 8.8e6 / 2.8e8 / 3.5e10 | 2.75 (diagonal-only 0.70) |
| Full stream, `xcorr` 0.15 | 4.2e2 / 2.5e3 / 4.2e4 | 0.79 |
| Thinned (dense075_coh05; ≥0.75° spacing), `xcorr` 0.625 | 2.4e1 / 8.2e2 / 6.9e4 | 0.61 |

Offline, adding a nugget of 0.1 or 0.2 to `xcorr = 0.625` lowers the median condition number to 3e2 or 1.6e2 (q 1.44 or 1.13). So a nugget fixes the conditioning too, but as section 2.7 shows, the correlation length would still be far too long.

### 3.2 DA outcomes, January 2020

All runs use R = 2.75 dB, `xcompact` 1.25°, the fixed operator and the new L1 climatology. Scores are O-F std as (DA − OL)/OL against `OLv8_M36_AZ_fixedop`, cross-masked on each DA run's own scaled obs (`score_cygl1_arm.py`).

| Run | Tb monitors (8) | SM monitors (ASCAT ×3 + CYGNSS L3) | CYGNSS L1 | RZEXC incr p99 / max | Wrong-way analyses* |
|---|---|---|---|---|---|
| Full stream, `xcorr` 0.625 (`DA_L1_full_fixedop`) | **+20.6%** | **+39.6%** | +2.7% | 22 / 582 mm | 46.7% |
| Full stream, `xcorr` 0.15 (`DA_L1_full_xc015_fixedop`) | **−4.75%** | **−0.95%** | −0.74% | 0.8 / 8 mm | 35.8% |
| Thinned, `xcorr` 0.625 (`DA_L1_dense075_coh05_fixedop`) | −0.40% | +0.67% | −0.02% | 0.3 / 3 mm | 15.7% |

\* Share of assimilated L1 obs with |O−F| > 0.5 dB where sign(A−F) ≠ sign(O−F).

**Full stream at 0.625:**
- Root-zone increments grew over the month: p99 went from 2 mm on Jan 1 to 61 mm on Jan 19.
- The worst single event (2020-01-10 06z, tile at −108.11°, 30.31°) was an RZEXC of −279 mm. Nine of the 10 nearest obs had O−F between −1.7 and −7.5 dB, yet the analysis moved the other way in obs space.

**Full stream at 0.15:**
- This is the best independent-monitor result for any CYGNSS assimilation in this project to date. It is **one month**; the 6-month run is in progress.
- About 36% of analyses still move away from their own obs, against 16% in the thinned run. This is consistent with 0.15 still over-correlating close same-track pairs (section 2.7).

### 3.3 Earlier results this affects

- **"Dense/full-stream L1 hurts more than thinned."** Pre-fix, this was probably driven largely by R conditioning at the prevailing `xcorr`, not by obs density itself.
- **The Desroziers/R reassessment** planned after the 6-month runs should treat the error correlation (structure and length) as a free parameter alongside `errstd`. The two interact: lowering the assumed correlation raises the effective weight of dense clusters.

## 4. Options for representing L1 obs-error correlation

In order of effort:

1. **Short isotropic `xcorr`, namelist only.** 0.10–0.15° is what's being tested now. It fixes the conditioning and matches the data beyond about 0.2°, but overstates short-range correlation (0.9 modelled vs about 0.3–0.5 implied at 0.03–0.07°). The FOV ≥ `xcorr` check is skipped for `cygl1scal` (`clsm_ensupd_upd_routines.F90` around line 530). `xcompact ≥ 2·xcorr` must still hold.
2. **Short `xcorr` plus a nugget, small code change.** This is the M21C `obs_cov_nugget` pattern, `fac = (1−a)·fac`. A nugget of about 0.5 with a correlated scale of about 0.07–0.1° best matches sections 2.1 and 2.7. Caveat from M21C: the obs perturbations are still drawn with the full correlated field, so the perturbation covariance would not match the new R unless that is changed too.
3. **Along-track superobbing or thinning, preprocessing.** Average or keep one obs per about 2–3 s (about 0.15–0.2°) of each track, then use a near-diagonal R (small `xcorr`). This removes most of the correlated structure at its source and reduces the cost. The weak long-range along-track component would remain, but uncorrelated obs in R is a reasonable approximation for it.
4. **Track-aware R, larger code change.** Correlate obs only within the same `sc_num`/`ch_id` and a time window (short and long component), and set it to zero across tracks. This matches the physics best. It needs the track identifiers carried through the reader into `assemble_obs_cov`.

**Suggested next steps:**
- Finish the 6-month runs: full stream at `xcorr` 0.15, thinned, SMAP-only and L3-only.
- Estimate R per stream from them with Desroziers diagnostics.
- Then test option 2 (nugget) against option 3 (along-track superob) on the full stream, one intervention at a time.

## 5. Caveats

- **Domain and time:** one domain (Arizona, arid to semi-arid), and 3 years of one CYGNSS L1 product version (`cygnss_l1_version v3.2`, preprocessor schema 0.5). Other regimes, such as vegetated or wet ones, may differ.
- **Obs and forecast definitions:** innovations use the ensemble-mean forecast of an open loop. The obs-error correlations are inferred through a Hollingsworth–Lönnberg split that assumes different-spacecraft pairs beyond 0.15° carry no obs-error correlation, and that the forecast-error correlation is exponential in distance. The implied ρ_o values are approximate, ±0.05.
- **Pair selection:**
  - Same-window only (3-hour windows). Along-track pairs therefore span only up to about 12 s, so the long-range along-track tail is observed only to 0.6°.
  - Pairs are never in the same owner tile, because the reader keeps one obs per tile per window. The shortest-range statistics come from adjacent-tile pairs.
- **Footprint definition:** overlap uses the preprocessor's normalized `coefficient_weight` support. Other overlap measures, such as coefficient-weighted or cosine similarity, were not compared.
- **Sample sizes:** some cells in the overlap and acquisition tables have few pairs (e.g. same-spacecraft-other at < 0.05°: N = 153). Rely on the bootstrap CIs where given.

## 6. Reproducibility

| Item | Location |
|---|---|
| Pair builder | `projects/CYGNSS_L1_AZ/scripts/obs_error_correlation/build_l1_pairs.py`: reads the OL ObsFcstAna, the L1 climatology and the preprocessed files; runtime about 15 min |
| Pair dataset | `projects/CYGNSS_L1_AZ/output/obs_error_correlation/l1_pairs_3yr.parquet` (1,053,216 rows; gitignored) |
| Analysis | `projects/CYGNSS_L1_AZ/scripts/obs_error_correlation/analyze_l1_pairs.py`; output `output/obs_error_correlation/analyze_l1_pairs.out` |
| OL run | `/discover/nobackup/projects/land_da/cygl1_operator_test/OLv8_M36_AZ_fixedop` |
| DA runs | `DA_L1_full_fixedop`, `DA_L1_full_xc015_fixedop`, `DA_L1_dense075_coh05_fixedop` (same project dir); templates in `templates/DA_L1_full*`, `templates/DA_L1_dense075_coh05` |
| DA scoring | `projects/CYGNSS_L1_AZ/scripts/postproc_drivers/score_cygl1_arm.py --da-expid <EXP> --arm-tag <tag> --ol-expid OLv8_M36_AZ_fixedop` |
| DA update diagnostics (section 3.1 and the increment/wrong-way columns of 3.2) | `projects/CYGNSS_L1_AZ/scripts/obs_error_correlation/da_l1_update_diagnostics.py --run DA_L1_full_fixedop:0.625 --run DA_L1_full_xc015_fixedop:0.15 --run DA_L1_dense075_coh05_fixedop:0.625 --month 202001` (reproduces the tables exactly) |
| Related M21C analysis | `/gpfsm/dnb06/projects/p284/M21C_testing/NOTES_2026-09-21_CF0360_HSAF_negative_soil_moisture.md` §12, and `notes_2026-09-21_files/scripts/r_conditioning.py` |

The offline `xcorr`/nugget sensitivity numbers in section 3.1 and the day-by-day increment growth in 3.2 came from one-off scratch scripts and are recorded here only.
