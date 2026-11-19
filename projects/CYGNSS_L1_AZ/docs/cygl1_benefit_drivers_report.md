# Where CYGNSS L1 DA helps and hurts (Arizona domain): drivers, and what they say about the fixed roughness term

*2026-11-19. Arizona limited domain (lon −118..−106, lat 29..40), EASEv2 M36 (909 tiles, 552 in the L1 mask), fixed-operator
runs, 2020-01-01 to 2022-12-31. Primary arm `DA_L1_full_xc015_coh040216_fixedop` (coherency filter, errstd 2.75) against
`OLv8_M36_AZ_fixedop`, with `DA_L3_fixedop` as control. **Everything below is exploratory:** the covariate list was fixed before
looking, but no outcome-specific hypothesis was. Code: [Reproducing](#reproducing). Companion reports:
[IV/TC](cygl1_iv_tc_az_report.md), [L1/L3/SMAP O−F](cygl1_fixedop_2020_2022_l1_l3_smap_report.md),
[operator skill](cygl1_operator_skill_report.md).*

## Summary

- **L1 has two separate stories, not one.** The outcomes split by reference system. The SMAP-reference measures (SMAP and SMOS O−F gain,
  IV against SMAP L3) agree with each other (r = 0.73). The ASCAT-reference measures (ASCAT O−F, ASCAT IV, TC ΔR²_mod) agree with each
  other (r = 0.71). Between the two groups the correlation is only 0.21 [0.08, 0.35]. A single combined score hides this.
- **TC ΔR²_mod is mostly an ASCAT quantity.** The TC triplet is SMAP L3 + model + ASCAT H121. The ASCAT leg carries about two thirds of
  the ΔR²_mod variance and the SMAP leg about one third (337 tiles where clamping does not interfere).
- **Where it works.** Against SMAP, L1 helps more where soil moisture varies a lot and the operator's SM term has skill, and less on
  barren tiles. Tiles that are better against both references are in the south-central and southeast domain and have more obs
  (median 1,314 against 708).
- **Where it hurts.** The ASCAT penalty is concentrated on the northeast high plateau (the Four Corners hotspot holds 71% of the ASCAT
  O−F penalty), but the family structure is domain-wide. Excluding the hotspot, TC against ASCAT IV still correlates at 0.72. The same
  tiles show no ASCAT penalty for L3 DA, so the penalty belongs to the L1 arm.
- **The fixed-roughness hypothesis is not supported by local data, and not tested.** Roughness-specific operator diagnostics
  (incidence-angle slope, coherency tail, outlier fraction) add nothing out of spatial block. Terrain proxies add a little, but their sign
  flips between the SMAP and ASCAT groups, which looks like confounding. Operator SM-term skill, obs count, baseline skill and a coarse
  elevation/aridity/vegetation gradient matter more. Only an IGOT perturbation experiment can isolate the roughness term
  ([Next](#next)).

## Data and method

| Item | Setting |
|---|---|
| Tiles | 909 AZ tiles; analysis on the L1 mask (`n_l1_obs ≥ 100`, 552 tiles). TC is valid in 365 of them, ASCAT IV in 357, SMAP L3 IV in 548 |
| Outcomes (positive = L1 better than OL) | TC ΔR²_mod (= −ΔfMSE) vs SMAP L3 / ASCAT H121; IV anomaly-correlation change vs each reference; percent O−F stddev reduction (sign flipped) for SMAP, SMOS, ASCAT, L3 obs; L1 minus L3 versions |
| Covariates | obs supply and quality; operator skill (per tile and season); signal/noise (H(x) spread, gain, innovation/spread); terrain, CTI, soil, canopy height; IGBP fractions; LAI/NDVI/opacity climatology; precipitation and aridity; OL soil-moisture climate; baseline OL skill |
| Uncertainty | 2° × 2° block bootstrap (28 blocks, 2000 draws), Benjamini–Hochberg FDR per outcome; spatial-block CV (5-fold, 10 repeats) for models |
| Sources | Local: O−F stats, operator tables, obs-quality bundle. Discover: per-tile TC/IV, operator skill, static covariates, vegetation, climate, OL soil-moisture climate. MODIS MCD12C1 v061 (2021) for IGBP fractions |

Caveats on the inputs:
- Vegetation is one repeating annual climatology from the boundary files, the same for all three years; there is no real VOD series.
- PET uses land skin temperature, not 2 m air temperature, so the aridity index is approximate.
- Outcomes exist per tile (full period or per year) or per domain-season. **Nothing is tile-season.** Seasonal swings cannot be tested
  against benefit within a tile.
- The `DA_SMAP_fixedop` arm is not included: its ASCAT IV columns do not exist yet (job drafted, not run).
- The IGBP dominant class agrees with the existing land-cover columns in only about 72% of tiles; the existing columns are per-obs samples.

## Results

### 1. The outcomes form two groups

- Within the SMAP group r = 0.73 [0.65, 0.80]; within the ASCAT group r = 0.71 [0.59, 0.79]; between groups 0.21 [0.08, 0.35].
- Two factors explain 83% of the variance (311 complete tiles). TC ΔR²_mod loads mainly on the ASCAT factor (0.77), with some SMAP loading (0.39).
- ASCAT O−F and ASCAT IV correlate at 0.86, so the ASCAT group is effectively two sources.
- Outside the hotspot the ASCAT-group correlations persist (TC vs ASCAT IV 0.72, TC vs ASCAT O−F 0.56). Inside it, ASCAT O−F is worse in 98%
  of tiles (median −7.4%) while SMAP O−F improves (+2.4%).

### 2. A combined score

The first combined score (mean of z(TC) and z(SMAP/SMOS O−F)) tracks the SMAP group at 0.85 and the ASCAT group at 0.54, so it is
dominated by SMAP-type measures. **Recommended primary:** the equal-weight average of the two groups' ranks. The SMAP group is IV SMAP L3,
SMAP O−F and SMOS O−F. The ASCAT group is ASCAT IV and ASCAT O−F. TC is kept as a sensitivity check, because it is the least spatially reliable
component and dropping it adds 30 tiles.

Always report the four categories, not only the scalar. Of 552 tiles (157 have no ASCAT):

| Category | Tiles | Where |
|---|---|---|
| Better on both | 120 | south-central and southeast; more obs |
| SMAP better, ASCAT worse | 175 | north-central and northeast plateau; grassland, +625 m, higher SM-term skill, more humid |
| ASCAT better only | 39 | western edge and northeast rim |
| Worse on both | 61 | western edge and northeast rim; many barren tiles |

On the 339 tiles where TC and the SMAP/SMOS O−F gain both exist, 149 are better on both measures, 134 on O−F only, 19 on TC only and 37 on neither.
The mean of the two measures is not a stable ordering of tiles: year-to-year consistency of the O−F measures is only 0.34–0.57.

### 3. What explains the pattern (univariate and multivariate)

Spatial-block CV R² (gradient boosting) against a lat/lon-only model and a permuted null:

| Outcome | GB | lat/lon only | null |
|---|---|---|---|
| SMAP/SMOS O−F gain | 0.50 | 0.34 | −0.07 |
| ASCAT O−F gain | 0.49 | 0.18 | −0.06 |
| IV vs SMAP L3 | 0.26 | 0.18 | −0.06 |
| IV vs ASCAT | 0.31 | −0.06 | −0.05 |
| TC ΔR²_mod | 0.14 | −0.03 | −0.10 |
| TC ΔR²_mod, L1 minus L3 | 0.01 | −0.03 | n/a |

- **Structure beyond smoothness:** gradient boosting beats lat/lon-only for every outcome except TC L1 minus L3, which has no demonstrable predictive structure.
- **Circularity:** some of the explained variance is by construction. The SMAP Tb correlation `op_r_SMAPh` and the OL SMAP O−F stddev (the denominator of the outcome) are the top covariates for the SMAP O−F gain. OL baseline skill comes from the same TC/IV triplets as the outcomes. Do not read those as drivers.
- **SMAP group, real signals (modest):** SM-term skill in JJA (ρ +0.41 [0.20, 0.56]), OL soil-moisture variability (≈ +0.38), barren fraction (−0.39). Top versus bottom quintile: about +2.5 to +3 points of O−F reduction.
- **ASCAT group:** the ASCAT O−F gain falls with elevation and slope (ρ ≈ −0.56 [−0.69, −0.36], top versus bottom quintile about −8.5 points) and rises with clean-flag fraction and CTI spread. IV and TC are baseline-limited: the OL skill is the top predictor.
- **The two groups have different covariate profiles** (profile correlation −0.14). Elevation and aridity covariates have opposite signs.
- **L1 minus L3:** elevation is the main driver of the SMAP-group difference (ρ +0.66 [0.53, 0.74]); L3 degrades at elevation more than L1 does.
- **Collinearity:** elevation, slope, DEM, precipitation, aridity, SM mean and spring greenness are one gradient, so the "top 12" covariates are only a few independent axes. Group-level results are more trustworthy than single-covariate ranks.

### 4. Time

- Domain level: the L1 O−F gain against SMAP falls from +2.7% (DJF) to +0.7% (MAM) and +0.5% (JJA). ASCAT is negative in every season (worst in SON, −4.1%).
  L1 minus L3 against SMAP is worst in MAM (−1.7).
- Operator SM-term skill rises slightly in JJA (+0.04) and dips in SON (−0.03); it does not dip in spring. SM-term variance, LAI, H(x) std and the outlier fraction rise in MAM/JJA.
- Tiles with a larger spring jump in vegetation opacity or LAI have worse benefit (opacity against ASCAT ρ −0.57 [−0.66, −0.39]). This is spatial and confounded, and it is consistent with spring greening as a degrader.
- Year-to-year variation is mostly noise (tile ICC 0.34–0.52) and no covariate varies by year, so year-anomaly tests were not possible.

### 5. The fixed-roughness hypothesis

How roughness enters, per the docs: it is fixed inside the IGOT pixel factor and only reshapes the geometry weights, which are frozen in preprocessing.
The reflectivity is flat Fresnel with no roughness or vegetation term, and H is purely coherent. The σ_S and σ_L values are not yet recovered from MiXIL.

What the data can say:
- The operator-residual diagnostics tied to a wrong roughness (incidence-angle slope, coherency tail, outlier fraction, dynamic ratio, SNR)
  add essentially nothing out of block (block-CV ΔR² ≈ 0 or negative) over SM-term skill and log obs count. None is in any top-five importance list.
- The static terrain block (terrain, canopy, ASCAT z0, incidence slope, geometry fraction) adds block-CV ΔR² of +0.08 to +0.11 for the SMAP and ASCAT
  O−F gains, and under +0.03 for TC, ASCAT IV and SMAP IV. Its sign flips between families (slope +0.24 for the SMAP O−F, −0.25 for the ASCAT O−F).
  The multivariate model gives a similar picture: terrain/roughness added to signal/noise gains +0.03 for the ASCAT O−F, none for the SMAP O−F and TC,
  and the ASCAT IV gain falls from +0.066 to +0.009 once baseline skill is included.
- Earlier OL residual diagnostics (O − scaled H(x)): the incidence-angle slope is +0.30 dB per 10° and positive in 90% of tile-seasons; it is
  0.24 in flat tiles and 0.42 in the roughest tercile. Returns with coherency below 0.5 sit about 1 dB below coherent ones in every season, and the
  diffuse fraction is 0.28 in rough terrain against 0.15 in flat. Both fit a wrong roughness damping, but weakly, and terrain, latitude, coherency and
  opacity covary. There is no seasonal bias in the tile-demeaned residual (within ±0.05 dB), which does not support a season-dependent roughness.
  The seasonal failure appears as lost SM-term skill, not a mean shift. Scaling may absorb a seasonal bias.

**Reading.** The hypothesis is not supported and not tested. The local data cannot separate a roughness effect from terrain and aridity, and the
fixed term can only be probed through proxies.

## Limits

- One domain, 28 blocks: the intervals are rough and the FDR is optimistic. Spatial confounding cannot be separated from mechanism.
- Outcome sample sizes differ (339–548 tiles), so the SMAP and ASCAT groups use different tiles, and the 157 tiles without ASCAT are not random
  (low elevation, barren, low SM-term skill).
- Spearman screens are univariate and monotone. Model importance is shared among collinear covariates.
- Several top covariates restate their outcome; the SMAP O−F R² is partly predictable by construction.
- The Four Corners hotspot was re-derived (130 tiles); the report's 73-tile mask is not on disk.
- TC is the least spatially reliable outcome (split-half reliability near 0.44, 0.18 once neighbours are excluded).

## Next

1. **IGOT roughness perturbation** (the decisive test). On about five tiles per terrain/vegetation stratum (200–500 OL obs each), recompute H(x) with σ_S and,
   separately, σ_L at 0.5×, 1× and 2×, using the same OL soil-moisture states. Apply the existing per-tile scaling and compare O−H skill and the O−F gain in MAM and DJF.
   If H only changes by a common factor that scaling absorbs, the hypothesis is rejected cheaply. If roughness is the weak link, SM-term skill and benefit should shift more on
   high-slope tiles and in spring.
2. **DA_SMAP_fixedop reference arm:** run the drafted ASCAT IV job on Discover, to separate L1-specific effects from generic DA effects on the plateau.
3. **Seasonal TC/IV:** needs a new Discover run with season-restricted dates; it is the only way to test seasonal swings per tile.
4. **A second domain** with a different climate, to break the spatial confounding.

## Reproducing

All scripts are in `projects/CYGNSS_L1_AZ/scripts/benefit_drivers/`; outputs are in `output/benefit_drivers/` (git-ignored) with one README per analysis folder.
Run with `/Users/amfox/mamba/envs/regrid/bin/python -I`.

1. `build_benefit_table.py`, `build_covariate_table.py`, `roughness_diagnostics.py`: local O−F, covariate and residual tables.
2. Discover extraction (per-tile TC/IV, operator skill, static, vegetation, climate, OL soil moisture): built on Discover; copy to `output/benefit_drivers/discover/`.
3. `build_igbp_landcover.py`: needs the MCD12C1 v061 files in `data/landcover_igbp/` (LP DAAC download; the script reads HDF4 via `hdp dumpsds` from the cygnss env).
4. `build_analysis_table.py`: joins all of the above into `analysis_table.csv` (909 × 350).
5. Analyses: `univariate_screen.py`, `multivariate_attribution.py`, `seasonal_roughness_analysis.py`, and `outcomes_{common,structure,combined,disagree}.py`.
