# Current Scientific State

> **Authority note:** This document records the accepted August 2026 scientific state of the manuscript. It takes precedence over older manuscript drafts, exploratory analysis notes, and superseded interpretations wherever they conflict, unless explicitly superseded by a later documented decision.
>
> The current GitHub analysis products and manuscript figures are the authoritative numerical and graphical sources. Older May 2026 manuscript drafts remain useful source material for methods, background, references, and wording, but they are not authoritative for the current scientific interpretation.

## 1. Experiment and scope

The manuscript examines a global GEOS Land Data Assimilation System (GEOS LDAS) experiment spanning **June 2000–May 2024**.

The principal paired experiments are:

- **OL:** open-loop simulation
- **DA:** multi-sensor data-assimilation simulation

The paired OL and DA experiments use the same meteorological forcing and land-model configuration; their principal difference is the application of land data assimilation in DA.

The experiment is performed on the **36-km EASEv2 M36 grid** using the GEOS Catchment Land Surface Model (CLSM) and a **24-member ensemble Kalman filter**.

The assimilated observing system evolves through the record and includes:

- MODIS snow-cover fraction (SCF)
- ASCAT soil-moisture retrievals
- SMOS brightness temperatures
- SMAP brightness temperatures
- CYGNSS soil-moisture retrievals

The manuscript should not claim that every temporal change in the analysis can be uniquely attributed to intrinsic changes in observation information content. Although the land-model and DA configuration are held fixed within the experiment, the system was not necessarily optimally retuned for every observing-system era.

## 2. Governing manuscript question

The manuscript is no longer primarily a generic multi-sensor validation study.

The central question is:

> **How does an evolving land observing system change the strength, spatial structure, physical consequences, and temporal consistency of a long-term multi-sensor land reanalysis?**

The manuscript should be organized around three connected scientific themes:

1. **Evolving observational influence**
   - How does the density, composition, and effectiveness of observational constraint change through time?

2. **Physical propagation of assimilation impacts**
   - To what extent do corrections to directly observed quantities propagate through modeled snow, soil-water, runoff, evapotranspiration, storage, and energy pathways?

3. **Temporal consistency**
   - Does the changing observing system alter long-term trends or introduce structural transitions into the analyzed record?

The governing scientific narrative is:

> **observing system evolves → direct analysis influence strengthens and changes character → corrections propagate through modeled land processes → the correction regime becomes demonstrably nonstationary → but broad long-term state behavior is generally preserved**

The conceptual endpoint is the distinction between:

- **evolving analysis influence**, and
- **preservation of long-term state behavior**.

These ideas are related but must not be conflated.

## 3. Relationship to Fox et al. (2025) and Fox et al. (2026)

For purposes of drafting this manuscript, the submitted/revised CYGNSS assimilation paper should be treated as a published paper and cited as:

> **Fox et al. (2026)**

The published ASCAT/SMAP paper remains:

> **Fox et al. (2025)**

These papers should be referenced extensively where appropriate, particularly in the Methods.

### Fox et al. (2025)

Use primarily for the established ASCAT + SMAP GEOS LDAS implementation, including:

- GEOS LDAS / CLSM soil-moisture assimilation framework
- 24-member EnKF implementation
- ASCAT soil-moisture assimilation
- SMAP brightness-temperature assimilation
- localization
- observation-space transformations / rescaling
- observation errors
- propagation of surface observational information into root-zone soil moisture
- interpretation of ASCAT versus SMAP and joint assimilation

### Fox et al. (2026)

Use primarily for the current multi-sensor soil-moisture configuration, including:

- ASCAT, SMOS, SMAP, and CYGNSS handling
- current multi-sensor EnKF configuration
- CYGNSS assimilation
- QC and observation handling
- bias correction / observation-space transformation terminology
- monitored-only observations and OmF diagnostics
- cross-sensor OmF interpretation
- practical multi-sensor assimilation caveats

The present manuscript should therefore avoid re-documenting every aspect of the established soil-moisture DA system. Its Methods should summarize only what is required to understand the 2000–2024 experiment and then refer to Fox et al. (2025, 2026) for implementation details.

Foundational methods should still cite their original sources where appropriate rather than replacing them entirely with self-citations.

## 4. Observing-system periods

**Figure 1 is the authoritative observing-system chronology.**

It defines **P1–P9** and includes, among other observing-system changes:

- the early MODIS Terra-only period
- the transition to Terra + Aqua
- ASCAT
- SMOS
- SMAP
- CYGNSS
- subsequent mission/configuration changes represented in the registry

P1–P9 must be used consistently wherever detailed observing-system periods are needed.

Broader periods may still be used where scientifically necessary, for example for validation sample size or the snow-process analysis, but they must be explicitly described as aggregations of the Fig. 1 chronology rather than as an independent competing period taxonomy.

Do not resurrect the old manuscript period definitions as the primary chronology.

The authoritative registry places:

- **P6 start: 1 April 2015**
- **P9 start: 1 December 2021**

The older November 2021 P9 date is superseded.

## 5. Evolution of observational influence

The broad observing-system evolution is:

- early record: predominantly MODIS SCF constraint
- 2007 onward: microwave soil-moisture observations become available
- SMOS introduces L-band brightness-temperature information
- **April 2015 / SMAP is the dominant system-level transition**
- CYGNSS later adds additional tropical/subtropical sampling

The central interpretation is not simply that more sensors are added, but that the **strength and character of observational constraint evolves substantially across the record**.

OmF diagnostics show that observational consistency generally improves as the observing system becomes denser and more capable.

The SMAP era is associated with a marked strengthening in soil-moisture observational influence.

Generic sensor-by-sensor comparisons should be kept concise because those questions overlap with Fox et al. (2025, 2026). The current manuscript should emphasize **temporal evolution** rather than re-running a sensor-comparison paper.

## 6. Independent evaluation

Independent soil-moisture evaluation supports the interpretation that later multi-sensor assimilation improves analyzed soil moisture, particularly near the surface.

Surface-soil-moisture responses are generally stronger than root-zone responses because the assimilated microwave observations constrain the near-surface state more directly.

Root-zone impacts are:

- physically meaningful,
- generally weaker than surface impacts,
- spatially heterogeneous,
- propagated indirectly through ensemble-estimated covariances and land-model dynamics.

This surface/root-zone distinction remains useful, but it should not be generalized into the older claim that unobserved hydrological quantities respond only weakly to DA. The snow-water-budget analysis demonstrates that assimilation-induced mass corrections can propagate strongly into unobserved hydrological fluxes and reservoirs.

ERA5-Land should be described as a **large-scale consistency comparison**, not as independent truth or fully independent validation.

## 7. Snow validation

MODIS SCF assimilation strongly constrains snow occurrence and snow-cover timing.

Independent snow evaluation shows clear improvements in snow occurrence / SCF behavior, whereas changes in SWE and snow depth are generally more modest or heterogeneous.

The correct interpretation is:

> SCF is the directly constrained snow quantity, whereas snow mass and depth are model-mediated quantities.

However, the limited direct improvement in SWE or snow depth does **not** imply that snow DA has little hydrological consequence.

The new process analysis demonstrates that assimilation-induced snow-water corrections propagate substantially into subsequent melt-season hydrology.

## 8. Snow-DA hydrologic propagation

The clean attribution period is the **six complete MODIS-only seasons, 2001–2006**, before microwave soil-moisture assimilation begins.

This period is particularly useful because snow DA is the only direct satellite constraint on modeled land-water storage.

The physical pathway is:

> **snow DA increment → snow storage / melt → runoff + infiltration / soil moisture → ET + storage**

Spatial, seasonal, and controlled analyses show that stronger snow-water corrections are associated with subsequent changes in:

- snowmelt
- root-zone soil moisture
- runoff
- evapotranspiration
- total water storage

These are **model-internal process diagnostics**, not independent validation.

### Controlled attribution

Within-tile controls remove or account for:

- persistent tile climatology
- common year effects
- open-loop snow amount

Spatial-block bootstrap uncertainty is used.

Earlier signed controlled coefficients include approximately:

- snowmelt: **0.751**
- infiltration: **0.394**
- RZMC: **0.000174 m³ m⁻³ per kg m⁻²**
- ET: **0.201**
- runoff: **0.301**
- total water/storage: **0.333**

These controlled regressions are supporting evidence rather than the primary water-partition result.

### Non-overlapping-window sensitivity

A stricter Oct–Mar predictor with subsequent non-overlapping response windows retains positive responses for:

- RZMC: **6.01 × 10⁻⁵**
- ET: **0.103**
- runoff: **0.0553**
- total land water: **0.147**

The corresponding 5° spatial-block confidence intervals remain above zero.

The important exception is delayed infiltration:

- coefficient approximately **0.0067**
- confidence interval approximately **−0.0049 to 0.0186**

Therefore, do **not** claim that delayed infiltration remains robust under the strict non-overlapping formulation.

Snowmelt remains a useful pathway diagnostic.

## 9. Common-window snow-DA water budget

This is the **main hydrological process result** and should be emphasized ahead of the regression analyses.

Across six complete water years, the mean signed snow-DA input is approximately:

> **58.24 kg m⁻² yr⁻¹**

The aggregate residual is approximately:

> **−1.21 kg m⁻² yr⁻¹**, or about **−2.1%**

For **snow-addition tile-years**, the directly accounted fate of assimilation-added snow water is approximately:

- **43.1% surface runoff**
- **21.2% baseflow**
- **35.9% evapotranspiration**
- **4.2% storage change**
- **−2.7% peatland free-standing water**
- **−1.7% residual**

Therefore:

> **total runoff = 64.3%**

The principal interpretation is:

> **Most snow water added by MODIS SCF assimilation subsequently leaves the land column as runoff; a substantial fraction leaves through evapotranspiration, while relatively little remains in storage by the end of the water year.**

The positive-input sample contains approximately:

> **247,545 tile-years**

The 5° spatial-block interval for total runoff is approximately:

> **61.1–67.4%**

Annual direct runoff fractions span approximately:

> **55.6–70.1%**

This indicates that the pooled runoff-dominated result is not produced by a single anomalous year.

### Independent controlled water-year estimator

A separate controlled regression gives approximately:

- runoff: **0.749**
- ET: **0.182**
- storage: **0.085**

With peak-snow control:

- runoff: **0.751**, CI **0.715–0.787**
- ET: **0.181**, CI **0.153–0.212**
- storage: **0.084**, CI **0.074–0.096**
- residual: **−0.016**, CI **−0.022 to −0.011**

The correct statement is that:

> direct water accounting and an independent controlled estimator both support a physically closed, runoff-dominated response.

Do **not** claim that the two estimators are numerically identical.

## 10. Water-budget audit and residual

A forcing audit confirms that OL and DA precipitation are effectively identical.

The float32 comparison gives approximately:

- maximum annual tile discrepancy: **0.102 kg m⁻²**
- domain-mean discrepancy: **0.000127 kg m⁻²**

This supports the interpretation that the water-budget differences arise from DA rather than unmatched precipitation forcing.

`WCHANGELAND` does not contain the analysis mass injection in the form needed for this differential budget. The accepted budget therefore uses the change in `TWLAND` between instantaneous 00Z 1 October restart endpoints, together with an explicit PEATCLSM free-standing water term for the store that `TWLAND` excludes.

The residual is small, systematically negative, and nearly constant across water years.

The likely explanation is the approximation involved in using monthly-mean endpoint storage rather than instantaneous same-date storage.

The manuscript should acknowledge the systematic residual rather than describing closure as exact.

## 11. Snowmelt interpretation

Snowmelt is a **pathway / transit term**, not a terminal sink in the added-water budget.

Controlled snowmelt coefficients are approximately **0.68–0.75**, depending on formulation.

A large coefficient is physically expected because assimilation-added snow water must subsequently leave the snowpack.

The scientific question is therefore not whether added snow melts, but **where that water goes after leaving the snowpack**.

The terminal partition is:

- runoff
- evapotranspiration
- end-of-period storage
- residual

Do not include snowmelt as a terminal water-budget fraction.

## 12. Snow-removal case

Snow-removal tile-years show a different partition, approximately:

- runoff: **38.3%**
- ET: **46.1%**
- storage: **10.5%**
- residual: **+5.2%**

This sample is smaller and lower magnitude than the snow-addition case.

The removal result should remain secondary or supplemental unless further work resolves whether the asymmetry is physically meaningful or primarily reflects sampling/noise differences.

Do not headline the removal case.

## 13. Root-zone soil-moisture response after snow DA

The principal defensible within-year RZMC quantities are:

- peak DA−OL RZMC: approximately **0.0189 m³ m⁻³**
- May is the most common month of peak response
- MJJ mean DA−OL RZMC: approximately **0.0108 m³ m⁻³**

These quantities can be used in the main text.

Persistence/residence-time diagnostics are more problematic because:

- approximately **67.2%** of cases are right-censored
- RZMC differences are already positive in October
- states inherit differences from earlier assimilation

Therefore the approximate **4.5-month persistence/residence** result should not be a main-text headline. If retained, it belongs in the Supplement with clear caveats.

## 14. Snow DA and later microwave soil-moisture DA

This investigation produced two distinct results.

### 14.1 Magnitude/activity hypothesis — not supported

The hypothesis:

> stronger snow DA causes greater subsequent soil-moisture DA activity

is not supported robustly.

Raw bins, within-tile controls, and tile-wise correlations do not provide a sufficiently interpretable positive magnitude relationship.

Do not make this a substantive manuscript claim.

### 14.2 Signed interaction — secondary supported result

Where snow DA adds water, later microwave soil-moisture analysis corrections tend on average to be more negative.

The preferred interpretation is:

> **Later microwave soil-moisture assimilation does not simply become more active where snow DA was stronger, but its signed corrections tend on average to oppose prior snow-added water.**

This is a population-level tendency, not a strong tile-year predictor.

Tile-wise correlations are weak.

Do not describe the interaction as full cancellation: the snow-to-hydrology response remains clearly detectable.

This result is secondary and can be summarized briefly in the main text, with detailed diagnostics in the Supplement.

## 15. Hydro-energy response

Under warm, snow-free conditions, locations with wetter DA−OL RZMC generally exhibit:

- greater evapotranspiration
- greater latent heat flux
- lower sensible heat flux
- greater evaporative fraction

The latent-heat components respond coherently.

This provides **internal physical-consistency evidence** that assimilation-induced soil-water changes propagate into surface-energy partitioning.

It is not independent validation of surface flux skill.

A short main-text synthesis is sufficient; detailed figures may remain supplemental.

## 16. Long-term trends

The trend analysis uses paired OL, DA, and DA−OL monthly fields for June 2000–May 2024.

The primary interpretation is:

> **DA generally preserves the broad long-term behavior of OL rather than introducing widespread secular drift, while modifying regional trends in root-zone soil moisture and total land-water storage and, over smaller areas, evapotranspiration and runoff.**

### Precipitation

OL and DA precipitation are identical by construction apart from negligible numerical precision differences.

Trend behavior is therefore essentially identical:

- OL and DA trend maps strongly agree
- OL/DA slope correlation approximately **0.9998**
- paired DA−OL has **zero FDR-significant tiles**

This is a **null-behavior / forcing-control sanity check**, not independent evidence of trend-estimator skill.

### Snow cover

SCF declines over the record in both simulations at nearly identical rates:

- OL: **−0.000554 yr⁻¹**
- DA: **−0.000549 yr⁻¹**
- DA−OL: **−0.000006 yr⁻¹**

The DA−OL trend is essentially null.

Interpretation:

> SCF assimilation largely preserves rather than creates the long-term SCF decline.

### Snow mass and snow depth

There is no important DA-induced long-term trend.

FDR-significant DA−OL trends are essentially absent.

### Root-zone soil moisture

RZMC exhibits a clear regional response to DA.

The area with significant regional trends increases from approximately:

- **4.36% in OL**
- to **9.18% in DA**

Significant paired DA−OL trends occur at approximately:

> **7.01% of valid-land tiles**

and are predominantly positive.

There is **no significant global land-mean RZMC trend**.

The correct interpretation is therefore regional modification of trends, not a global DA-induced wetting trend.

### Evapotranspiration, total runoff, and total land-water storage

A focused extension applied the same paired trend framework to evapotranspiration, total runoff, and total land-water storage. None of the nine area-weighted OL, DA, or DA−OL domain series has a significant whole-record trend after autocorrelation correction and false-discovery-rate control.

At the tile scale, significant paired DA−OL trends occur at **4,121 ET tiles (3.66%; 3.75% of land area)**, **5,434 total-runoff tiles (4.83%; 5.00% of land area)**, and **8,590 total-storage tiles (7.63%; 8.15% of land area)**. The changes are predominantly positive: 3,598 of the significant ET tiles, 4,443 runoff tiles, and 8,030 storage tiles have positive DA−OL trends.

OL and DA tile slopes remain strongly correlated, with correlations of **0.948 for ET, 0.968 for runoff, and 0.932 for total storage**. Thus, as for RZMC, DA largely preserves the background spatial trend structure while superimposing smaller regional trend modifications. The effect is broadest for total land-water storage.

The accepted long-term interpretation is therefore:

> **The evolving DA system modifies regional water-cycle trends without producing a resolved global secular trend in RZMC, ET, runoff, or total land-water storage.**

## 17. Observing-system transitions and April 2015

Whole-record trend behavior and structural observing-system transitions are separate questions.

The clearest transition occurs at:

> **April 2015 — P5 to P6 — introduction of SMAP brightness-temperature assimilation**

The interrupted-series estimate of the RZMC DA−OL level shift is approximately:

- valid land: **+0.00102 m³ m⁻³**
- warm snow-free: approximately **+0.00127 m³ m⁻³**

The approximate bootstrap interval for the valid-land response is:

> **0.00049–0.00156 m³ m⁻³**

The transition is accompanied by marked increases in soil-moisture analysis activity.

For example, absolute prognostic soil-water correction activity increases by approximately:

> **+10.61 kg m⁻²**

The corresponding focused analysis of ET, runoff, and total land-water storage provides additional physical context for the April 2015 transition. Total land-water storage DA−OL increases by **2.144 kg m⁻²** at P6, with a **95% bootstrap interval of 0.932–3.266 kg m⁻² and q = 0.009**. This is a formally significant positive storage shift.

ET DA−OL also increases at P6 by **0.502 kg m⁻² month⁻¹**, with an individual bootstrap interval of **0.094–0.883 kg m⁻² month⁻¹**, and PELT independently places an ET changepoint exactly in April 2015. However, the prescribed-date coefficient has **q = 0.072** and therefore does **not** survive the separate nine-series boundary-family FDR correction. It should be described as convergent or suggestive rather than formally FDR-significant.

Total runoff has no comparable P6 response: its estimated level change is **+0.230 kg m⁻² month⁻¹**, with an interval of **−0.170 to +0.623 kg m⁻² month⁻¹ and q = 0.478**.

The original soil-water transition analysis shows:

- **10 primary diagnostic series** have accepted independent changepoints exactly in April 2015
- **9 of those 10** also have significant known-date P6 level changes
- no accepted changepoint occurs in paired OL or DA state-control series

The focused flux/storage analysis independently identifies April 2015 changepoints in both ET and total storage. This strengthens the interpretation of April 2015 as a system-level land-water transition while retaining the distinction between formal and suggestive individual responses.

A second physically relevant transition occurs near P2 during the MODIS-only period. Total runoff DA−OL increases by **0.790 kg m⁻² month⁻¹** at P2, with a **95% interval of 0.379–1.176 kg m⁻² month⁻¹ and q = 0.0045**. The P2 runoff slope is also positive at **0.108 kg m⁻² month⁻¹ yr⁻¹** with an interval of **0.021–0.194 and q = 0.0349**. PELT places an accepted runoff changepoint in June 2002, one month before the registered July 2002 P2 boundary.

Because snow-cover assimilation is the only satellite land DA operating during this early period, the runoff transition is physically consistent with the independent MODIS-only snow-water-budget result showing a runoff-dominated fate of assimilation-added snow water. The temporal agreement nevertheless remains an observing-system association and should not be presented as proof that the introduction of Aqua uniquely caused the runoff change.

The accepted independent-break inventory from the original state/correction analysis contains:

- **37 total accepted breaks**
- **20 paired DA−OL**
- **17 correction diagnostics**
- **20** matching known observing-system boundaries within ±3 months
- **2 additional** matching within ±6 months
- **15 unmatched**

The manuscript should say that the April 2015 change:

- **coincides with** the introduction of SMAP
- is **consistent with** a SMAP-driven change in analysis behavior

Do not claim that statistical coincidence uniquely proves SMAP is the sole cause of every contemporaneous change.

The overall transition interpretation is:

> **The evolving observing system produces clear changes in correction activity and discrete changes in the analyzed land-water state, most clearly in April 2015, without producing a corresponding global secular drift. The April 2015 transition is expressed most strongly through soil-water correction activity, root-zone soil moisture, and total land-water storage, whereas the clearest runoff transition occurs earlier during the MODIS-only period.**

## 18. Current main-text figure sequence

The intended figure sequence is:

1. **Fig. 1** — observing-system timeline / P1–P9
2. **Fig. 2** — mean assimilated observations per day
3. **Fig. 3** — monthly assimilated observation counts
4. **Fig. 4** — full-period OmF standard-deviation improvement
5. **Fig. 5** — monthly evolution of OmF improvement
6. **Fig. 6** — period/sensor OmF improvement maps
7. **Fig. 7** — ISMN soil-moisture skill by validation period
8. **Fig. 8** — IMS snow-cover skill maps
9. **Fig. 9** — SNOTEL SWE skill
10. **Fig. 10** — GHCN snow-depth skill
11. **Fig. 11** — ERA5-Land soil-moisture comparison bars
12. **Fig. 12** — ERA5-Land soil-moisture improvement maps
13. **Fig. 13** — ERA5-Land snow comparison
14. **Fig. 14** — snow-DA water budget
15. **Fig. 15** — monthly snow-DA hydrologic pathway
16. **Fig. 16** — long-term RZMC and SCF trends
17. **Fig. 17** — observing-system transitions

### Fig. 16 plotting decision

For Fig. 16:

- RZMC DA−OL should retain a separately optimized color scale.
- SCF DA−OL should use the same scale as the corresponding OL/DA panels so that the near-null assimilation-induced SCF trend is not visually exaggerated.

### Fig. 17 plotting decision

For Fig. 17:

- panel (a) should explicitly identify the plotted standardized quantities as z scores / standardized seasonally adjusted values
- P6 / SMAP should be immediately identifiable
- panels (b–d) remain in native units with bootstrap confidence intervals
- the independent breakpoint-agreement matrix remains supplemental rather than being added to Fig. 17

## 19. Main text versus Supplement

### Main text

The main manuscript should contain:

- P1–P9 observing-system chronology
- concise observation-density and OmF evolution
- independent soil-moisture validation
- snow occurrence / SWE / depth validation
- concise ERA5-Land comparison
- common-window snow-water budget
- monthly snow-to-hydrology pathway
- concise controlled-attribution corroboration
- within-year RZMC response
- concise hydro-energy physical-consistency result
- long-term RZMC / SCF trends
- April 2015 observing-system transition
- optional short signed snow-to-SM counteraction synthesis

### Supplement

The Supplement should contain most of the diagnostic machinery, including:

- detailed controlled regression hierarchy
- strict Oct–Mar non-overlap analysis
- infiltration exception
- percent-change companion analyses
- alternative block sizes / trimming
- snow-removal partition
- individual annual budget details
- budget technical audits
- RZMC persistence/residence-time diagnostics
- signed snow-to-SM detailed analysis
- precipitation trend control
- SFMC trends
- snow-mass and snow-depth trends
- ET / total-runoff / total-land-water-storage trend maps
- ET / total-runoff / total-land-water-storage transition-series diagnostics
- breakpoint-boundary agreement
- changepoint sensitivity analyses

## 20. Revised interpretation of directness of constraint

The earlier manuscript interpretation that DA impacts become weak or negligible as variables become less directly observed is too strong.

The accepted interpretation is:

> Observational constraints remain strongest and most immediate for quantities closely related to the assimilated observations, but assimilation-induced state and mass corrections can propagate through model dynamics into unobserved hydrological reservoirs, fluxes, and surface-energy partitioning.

The snow-water budget is the clearest evidence for this revised interpretation.

## 21. Claims that should no longer be made

Do not state or imply that:

- changing observing-system effects can be summarized simply by the old pre-ASCAT / pre-SMAP / SMAP taxonomy
- stronger snow DA necessarily causes greater later soil-moisture DA activity
- the approximately 4.5-month RZMC residence time is a secure main-text result
- delayed infiltration remains robust under the strict non-overlapping-window analysis
- snowmelt is a terminal sink in the snow-water budget
- the snow-removal partition is as well constrained as the snow-addition result
- the two water-budget estimators are numerically identical
- precipitation provides an independent trend validation
- ERA5-Land is independent truth
- DA creates a significant global mean RZMC trend
- the April 2015 coincidence proves unique causation by SMAP
- snow DA has little hydrological effect merely because SWE/depth validation responses are modest
- all unobserved hydrological quantities respond only weakly to assimilation
- fixing the DA configuration across the record means all temporal differences uniquely represent changes in intrinsic sensor information content

## 22. Remaining factual/editorial checks

These are cleanup items rather than unresolved scientific questions:

- verify final corrected-precipitation description for the actual experiment
- verify exact experiment/spin-up dates and naming
- verify MODIS cloud / surface-temperature / QC description
- verify snow-albedo climatology description
- verify SMOS alias-free-zone / retrieval-screening terminology
- verify SMOS and SMAP observation-count definition, including H/V handling
- verify the surface soil-moisture frozen-temperature threshold
- define precisely what is counted as a MODIS observation in Figs. 2–3
- ensure OL is not described as using localization or observation-error settings that apply only to DA
- provide a concise rationale for ERA5-Land rather than ERA5
- remove redundant old observing-system table if Fig. 1 fully replaces it
- keep terminology and sign conventions consistent throughout figures and text

## 23. Manuscript-writing rules

The rebuild should be **destructive rather than additive**.

The May manuscript is source material, not the structure that must be preserved.

New analyses should displace redundant generic system and sensor detail rather than simply increasing manuscript length.

The preferred architecture is:

### 1. Introduction

### 2. Methods and Data
- 2.1 GEOS LDAS and Catchment model
- 2.2 Assimilated observing system
- 2.3 Land data assimilation
- 2.4 Experiments and observing-system periods
- 2.5 Evaluation datasets and metrics
- 2.6 Snow-DA hydrologic pathway and water-budget analysis
- 2.7 Long-term trends and observing-system transitions

### 3. Results
- 3.1 Evolution of the assimilated observing system
- 3.2 Evolution of observational influence
- 3.3 Independent soil-moisture evaluation
- 3.4 Snow-cover, SWE, and snow-depth evaluation
- 3.5 Comparison with ERA5-Land
- 3.6 Hydrologic consequences of snow-cover assimilation
- 3.7 Long-term trends and observing-system transitions

### 4. Discussion
- 4.1 Evolving observational constraint in long-term land analysis
- 4.2 From direct observational constraint to model-mediated hydrologic response
- 4.3 Observing-system transitions and long-term trend fidelity
- 4.4 Limitations and implications for future multi-sensor reanalysis

### 5. Conclusions

The final manuscript should ideally be **shorter in prose than the May draft despite containing more science**.

## 24. Do-not-revert list

Unless explicitly superseded by a later documented decision:

1. **June 2000–May 2024** is the experiment period.
2. **Fig. 1 / P1–P9** is the authoritative observing-system chronology.
3. **April 2015 / P6 / SMAP** is the dominant identified observing-system transition.
4. The main snow-process result is the **common-window differential water budget**.
5. For snow-addition tile-years, approximately **64.3% of added snow water becomes total runoff**.
6. Snowmelt is a transit term, not a terminal water-budget sink.
7. The controlled estimators corroborate but do not numerically duplicate the direct budget.
8. The strict non-overlap delayed-infiltration result is not significant.
9. RZMC peak and MJJ response may be used; the residence-time result remains supplemental.
10. The snow-DA → later SM-DA magnitude/activity hypothesis is not supported.
11. The signed later SM correction has a weak but coherent opposing tendency.
12. Hydro-energy changes are physical-consistency diagnostics, not independent flux validation.
13. DA generally preserves broad long-term behavior while modifying regional RZMC and total-land-water-storage trends and, over smaller areas, ET and runoff trends.
14. No significant global land-mean whole-record trend is resolved in OL, DA, or DA−OL ET, total runoff, or total land-water storage.
15. At P6, total land-water storage DA−OL has a significant +2.144 kg m⁻² level shift; ET shows convergent but not FDR-significant evidence for a positive shift, and runoff has no resolved P6 level change.
16. The clearest total-runoff transition occurs near P2 in the MODIS-only period and is consistent with, but does not independently prove, the snow-DA runoff pathway.
17. There is no significant global land-mean RZMC trend.
18. SCF long-term decline is essentially unchanged by DA.
19. Precipitation is a null forcing-control check, not independent validation.
20. Fig. 16 shows RZMC and SCF trends.
21. Fig. 17 shows observing-system transitions.
22. Fox et al. (2025) and Fox et al. (2026) should carry much of the detailed soil-moisture DA methods burden.
23. The manuscript's core argument is **evolving observational influence + physical propagation + temporal consistency**, not simply whether DA improves validation metrics.
