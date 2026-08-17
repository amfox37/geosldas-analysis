# Supporting Information for

## Evolving Observational Influence and Temporal Consistency in a Multi-Sensor Land Reanalysis, 2000–2024

*Author list to be inserted from the final manuscript front matter.*

This Supporting Information provides additional methodological detail, statistical calibration, sensitivity analyses, and supporting diagnostics for the analyses described in Sections 2.6 and 2.7 of the main manuscript. The definitions and estimands are identical to those used in the main analysis; the material below documents their implementation and robustness rather than introducing alternative primary analyses.

---

# Text S1. Additional methods for snow-DA hydrologic attribution

## S1.1 Differential snow-DA water budget

The differential budget is evaluated over the six complete MODIS-only water years WY2001–WY2006, each running from October through September. This interval precedes microwave soil-moisture assimilation in June 2007, so DA−OL differences are attributable to snow-cover assimilation alone. The domain is the 48,067-tile Northern Hemisphere seasonal-snow mask, defined by a 20°N southern limit, snow-cover-possible thresholds of 0.05 snow-cover fraction and 5 kg m⁻² snow mass, and a maximum mean JJA snow-cover fraction of 0.20 to exclude permanent ice. All quantities are evaluated in native M36 tile space, and fluxes are converted to monthly water-equivalent totals before aggregation over the water year.

The snow-DA input, $I_{\mathrm{snow}}$, is the water-year sum of the native prognostic `snow_net` increments, with positive values denoting addition of snow water. The stored monthly increment product is already the sum of native submonthly increments within each month, so no temporal differencing is applied. The closing terms are the DA−OL changes in surface runoff, baseflow, evapotranspiration, and total land-water storage:

$$
I_{\mathrm{snow}} = \Delta R_{\mathrm{surf}} + \Delta R_{\mathrm{base}} + \Delta ET + \Delta S + \epsilon
$$

where $\epsilon$ is the residual. Precipitation does not enter the differential budget because the two experiments use nominally identical forcing. This is verified rather than assumed: the production audit checks the maximum absolute annual tile-level and area-weighted domain-mean DA−OL precipitation difference against a 0.2 kg m⁻² annual tile guard, and the six-year integrated area-weighted domain difference is −5.5 × 10⁻⁵ kg m⁻².

Two classes of quantity are excluded from the closing equation to avoid double counting. Snowmelt and infiltration redistribute water within the land column rather than removing it from the differential budget, and are retained as pathway diagnostics rather than terminal sinks. RZMC and SFMC are diagnostic states already contained in `TWLAND` and are retained as timing and impact diagnostics; the compressed monthly inputs contain no separate native prognostic soil-water storage state, so no independent soil-water storage term is available.

$\Delta S$ is calculated from the September-to-September change in the DA−OL monthly-mean `TWLAND` state, not from the integrated `WCHANGELAND` tendency. The tendency diagnostic excludes the discontinuous mass introduced by the analysis update and therefore cannot equal the change in the DA−OL state anomaly; it is retained only as a process-balance diagnostic. Monthly-mean endpoints introduce a small temporal approximation that is considered when interpreting the residual.

### Partition estimator

For a selected tile-year sample $G$ and response term $k$, the partition fraction is

$$
f_k = \frac{\sum_{(i,t) \in G} A_i\, Y_{k,i,t}}{\sum_{(i,t) \in G} A_i\, I_{\mathrm{snow},i,t}}
$$

where $A_i$ is the M36 tile area, $Y_{k,i,t}$ is the DA−OL response for term $k$ at tile $i$ in water year $t$, and $I_{\mathrm{snow},i,t}$ is the signed snow-water input.

Partition fractions are ratios of aggregate area-weighted response to aggregate area-weighted input rather than averages of tile-year ratios, avoiding unstable ratios where local snow input approaches zero.

The main manuscript emphasizes $G = \{I_{\mathrm{snow}} > 0\}$, the snow-addition sample, comprising 247,545 tile-years across 45,490 tiles. The all-tile sample (288,402 tile-years, 48,067 tiles) and the snow-removal sample (37,395 tile-years, 13,012 tiles) are retained as diagnostics and sensitivities.

## S1.2 Spatial-block uncertainty for the water-budget fractions

Partition uncertainty was estimated with 1,000 spatial-block bootstrap replicates. M36 tiles were assigned to approximately 5° × 5° latitude–longitude blocks, and complete blocks, including all associated tile-years, were resampled with replacement. For each realization, the area-weighted response and snow-input sums were recomputed and their ratio retained. The 2.5th and 97.5th percentiles define the reported 95% interval. A 10° × 10° block analysis was used as a sensitivity.

Blocks are formed by binning tile longitude and latitude into fixed-width cells and assigning each unique bin pair a block code; each realization draws the original number of blocks as a multinomial allocation of block counts under a fixed seed. Resampling is spatial, not temporal.

## S1.3 Controlled water-year regression

The controlled regression is fitted as

$$
\tilde{Y}_{i,t} = \alpha + \beta\, \tilde{I}_{i,t} + \delta\, \tilde{S}^{\mathrm{OL}}_{i,t} + \gamma_t + \varepsilon_{i,t}
$$

where a tilde denotes deviation from that tile's six-water-year mean, and

- $\tilde{Y}_{i,t}$ is the within-tile anomaly of the signed DA−OL terminal response (total runoff, ET, storage, or residual);
- $\tilde{I}_{i,t}$ is the within-tile anomaly of signed water-year snow-water input;
- $\tilde{S}^{\mathrm{OL}}_{i,t}$ is the within-tile anomaly of the day-weighted March–May mean OL snow mass (`SNOMASLAND`), weighted by days in month;
- $\gamma_t$ are water-year fixed effects;
- $\beta$ is the controlled response per unit anomalous snow-water input.

Within-tile anomalies remove persistent spatial differences among locations, year fixed effects remove domain-wide interannual anomalies, and the OL MAM snow-mass anomaly controls for interannual variations in the underlying spring snowpack. The production design contains an intercept, the snow-input anomaly, the OL MAM snow-mass anomaly, and five year indicators. The sample is all 288,402 complete signed tile-years, not only the positive-input subsample.

Coefficient uncertainty uses the same spatial-block resampling described in S1.2, with 5° blocks primary and 10° blocks as a sensitivity. Because the design matrix and sample are common to all responses, the bootstrap is implemented through per-block sufficient statistics, so every replicate refits the full design.

The response coefficients for total runoff, ET, storage, and the residual sum to exactly one. This is a property of the construction rather than an empirical finding: the four terms close the same differential budget by definition and are fitted on an identical design matrix and sample, so their coefficients on $\tilde{I}$ necessarily sum to the coefficient from regressing $\tilde{I}$ on itself.

## S1.4 Why the controlled coefficients are not identical to the direct fractions

The controlled coefficients and direct partition fractions estimate related but different quantities. The direct calculation is an area-weighted mass partition conditional on positive snow input, whereas the regression uses signed within-tile interannual variation with year and OL-snow controls. Numerical agreement is therefore not expected; agreement in the response ordering provides independent corroboration.

The accepted production estimates are:

| Term | Direct positive-input accounting | Controlled regression $\beta$ [95% CI] |
|---|---:|---:|
| Total runoff | 64.3% | 0.749 [0.711, 0.783] |
| Evapotranspiration | 35.9% | 0.182 [0.155, 0.213] |
| Storage | 3.9% | 0.085 [0.075, 0.097] |
| Residual | −4.1% | −0.016 [−0.022, −0.011] |

Direct fractions use the snow-addition sample (247,545 tile-years); regression coefficients use the production run with the OL MAM snow-mass control and 5° spatial blocks over all 288,402 signed tile-years. Both columns sum to unity once the residual is included.

Both approaches identify runoff as the dominant response.

## S1.5 Additional robustness tests

- **10° spatial blocks.** Coarser resampling widens all intervals but does not change any sign or ordering; the positive-input runoff fraction interval widens from 61.1–67.2% to 58.3–69.7%.
- **1st–99th percentile predictor-anomaly trimming.** Applied to the seasonal regression to confirm that the fitted slopes are not driven by extreme input anomalies. All four classified downstream responses retain intervals excluding zero.
- **Non-overlapping October–March analysis.** The predictor is restricted to October–March snow input and responses begin in April or May, giving exactly zero shared months. All four classified downstream responses survive; standardized effects range from 24% (runoff) to 67% (ET) of the corresponding signed-MAM values but remain positive.
- **Water-year accounting-boundary sensitivity.** Repeating the partition on a September–August boundary changes the runoff fraction from 64.3% to 64.5% and ET from 35.9% to 36.0%, both far below the pre-set 5 percentage-point reporting threshold. Storage changes by +0.3 and the residual by −0.6 percentage points. Annual residuals are negative in 6 of 6 years under both boundaries.
- **Snow-removal sample.** Tile-years with $I_{\mathrm{snow}} < 0$ partition differently, with a larger ET share (46.1%) and a positive residual.
- **Soil-moisture persistence diagnostic.** Time from the within-year RZMC peak to a 0.001 m³ m⁻³ threshold is reported for the snow-addition sample.

Three interpretive limitations apply:

1. The strict non-overlap test cannot resolve delayed infiltration; restricting the predictor to October–March discards input reaching the soil column later in the same water year, so the reduced coefficients bound rather than measure the seasonal response.
2. Snow-removal partitioning is qualitatively different from snow-addition partitioning and is not used for the headline result.
3. RZMC persistence is strongly right-censored — 67.2% of snow-addition tile-years never fall below the threshold after their within-year peak before September — and the mean DA−OL RZMC anomaly is already positive in October, so these counts include inherited state differences from prior assimilation. Residence-time estimates are not used as headline results.

---

# Text S2. Additional statistical methods for trends and observing-system transitions

## S2.1 Distinct statistical questions

The temporal analysis addresses three distinct questions: whole-record trends, changes at prescribed observing-system boundaries, and independently detected abrupt changepoints. These are estimated separately because a discrete level or slope change need not produce a significant whole-record trend, and a long-term trend need not correspond to an abrupt transition.

## S2.2 Paired DA−OL trend estimand

For paired model variables, the primary assimilation-impact quantity is the monthly paired series

$$
D_i(t) = \mathrm{DA}_i(t) - \mathrm{OL}_i(t)
$$

The primary assimilation-impact trend is estimated directly from the paired DA−OL series. OL and DA are restricted to common finite monthly support before differencing; their separate trends are retained as contextual control fields. OL, DA, and DA−OL maps retain separate complete-field FDR families.

## S2.3 Trend-preserving seasonal adjustment

Seasonality is removed with a trend-preserving calendar-month adjustment implemented in four steps:

1. estimate a preliminary exact Theil–Sen slope on the raw series;
2. temporarily remove that slope;
3. estimate each calendar month's climatology from the detrended values;
4. subtract **only** that estimated seasonal climatology from the original, untrended series.

This construction removes the calendar-month climatology without absorbing part of the underlying linear trend into the seasonal adjustment. The trend removed at step 2 is restored at step 4 and is not removed before the final test.

Calendar months failing the configured minimum climatology sample count (three) are omitted from the adjusted series. Production outputs expose a `calendar_month_used` matrix and per-tile `n_calendar_months_used` so that seasonally selective dropping is visible.

## S2.4 Theil–Sen slope and modified Mann–Kendall test

Trend magnitude was estimated with the exact Theil–Sen median slope, the median of all pairwise slopes between observations, expressed in source units per year. Significance was assessed separately using the Mann–Kendall test with the stated Hamed–Rao autocorrelation adjustment.

The Hamed–Rao implementation is a conservative adaptation of the modified Mann–Kendall variance correction. Rank autocorrelation is calculated from Sen-detrended residuals at actual monthly lags 1–24, requiring at least 24 pairs per lag and using a 0.05 autocorrelation significance level. Only significantly *positive* lag correlations inflate the variance, and the resulting variance factor is floored at one. With intermittent eligibility, lag correlations pair values that remain separated by the requested calendar-month lag while the Hamed–Rao weight uses the valid observation count; this is a documented gap-aware approximation to the regularly sampled formula.

By construction the correction can only make a result harder to declare significant. The uncorrected independent-sample Mann–Kendall p-value is retained in the output to show the size of the correction.

## S2.5 Spatial multiple testing and trend confidence intervals

Benjamini–Hochberg false-discovery-rate control at $q = 0.05$ is applied across all finite tiles in each output field. The family structure is explicit:

- each complete mapped field is its own FDR family;
- OL, DA, and DA−OL fields are separate families, even where they share an underlying sample;
- FDR is computed from the complete field, not independently by spatial chunk. Spatial subset runs are written to a distinct diagnostic path and explicitly labeled as non-global FDR; only a complete-mask run supplies production FDR values;
- maps use the `significant_fdr` flag.

Separately, the reported Theil–Sen confidence interval is a **first-order autocorrelation adjustment**, not the exact inversion of the modified Mann–Kendall test. SciPy's nominal 95% Sen limits are retained as `slope_ci_*_nominal`, and the primary interval expands each nominal half-width by the square root of the same Hamed–Rao variance factor used for significance.

The Sen confidence interval and modified-Mann–Kendall/FDR decision are separate inferential products and need not agree exactly; mapped significance is defined by the FDR-controlled test. Production outputs expose an explicit `fdr_ci_disagreement` flag wherever the FDR test is significant but the adjusted Sen interval contains zero. Disagreement is absent for the paired state fields and the snow sensitivities, and reaches 22.1% of significant tiles for absolute soil-water activity and 8.9–13.6% for the SFMC correction diagnostics. The first-order interval is not plotted as inferential whiskers for the activity fields.

## S2.6 Known-date interrupted time-series model

The response is an area-weighted monthly domain series over the complete 288-month record. Tile weights are M36 land-tile areas from the GEOS LDAS tile-coordinate file; missing OL or DA values are cross-masked before any paired mean or difference is formed. The segmented design contains:


| Component | Count |
|---|---:|
| Intercept | 1 |
| P1 baseline slope | 1 |
| Calendar-month fixed effects | 11 |
| Level changes at P2–P9 | 8 |
| Slope hinges at P2–P6, P8, P9 | 7 |
| **Total fitted parameters** | **28** |

The design is full rank over 288 months. P7 spans only 15 months: it receives a level-change term but no independently fitted slope hinge, and the P6 slope is carried until the P8 boundary. This constraint is inherited from the period registry and fixed before estimation.

Three coefficient types appear in the output:

- **Level change** — immediate vertical displacement at a boundary, after accounting for seasonality and the piecewise trend; source units.
- **Slope change** — change in the rate of evolution across a boundary; source units per year.
- **Period slope** — the resulting slope within a defined period; source units per year.

Period slopes are not assimilation impacts: several significant control period slopes appear in both OL and DA and are consistent with shared forcing or climate variability.

## S2.7 Prais–Winsten AR(1)

Serial dependence is fitted with iterative Prais–Winsten AR(1): the first-order residual persistence is estimated, the regression is transformed accordingly, the first observation is retained rather than discarded, and the procedure is refitted until the AR(1) estimate stabilizes. The production configuration allows at most 25 iterations with a convergence tolerance of 10⁻⁷ and caps the absolute AR(1) at 0.98.

Prais–Winsten supplies the fitted serial-dependence model used by the production bootstrap and is not itself the final uncertainty estimator. Both OLS with Newey–West standard errors and Prais–Winsten with Newey–West standard errors were anti-conservative in the fixed-seed 288-month AR(1) no-transition test, whereas the bootstrap passed the BH false-discovery guard.

## S2.8 Innovation-resampling bootstrap

Primary transition p-values and confidence intervals come from a fixed-seed fitted-AR(1) innovation-resampling bootstrap:

1. fit the complete Prais–Winsten segmented model;
2. obtain the whitened innovations from the fitted residuals;
3. center the innovations;
4. resample the centered innovations with replacement;
5. propagate the resampled innovations through an AR(1) process, discarding a 120-month burn-in;
6. add the simulated noise to the fitted deterministic model;
7. refit the entire segmented model to the reconstructed series;
8. repeat for 1,999 replicates.

The complete model is refitted on every replicate, so the bootstrap distribution reflects estimation uncertainty in all 28 parameters jointly.

The simulation AR(1) parameter is the upper 95% large-sample confidence bound on the fitted AR(1), capped at 0.98; both the fitted and simulation values are reported. A 28-parameter segmented fit can absorb part of the true residual persistence into the deterministic component and bias the estimated AR(1) downward, so simulating at the upper bound makes the intervals conditional on a conservative persistence estimate.

Inference uses centered two-sided empirical bootstrap p-values and basic bootstrap confidence intervals. Newey–West/Bartlett and independent-sample standard errors remain in the output as diagnostics and do not determine transition significance.

## S2.9 FDR families for the known-boundary analysis

FDR control is applied separately within each predeclared coefficient family, where a family is one named boundary or period across all registered domain series. For the original state and correction-diagnostic analysis this gives:

| Family type | Count |
|---|---:|
| Level-change boundaries (P2–P9) | 8 |
| Slope-change boundaries | 7 |
| Independently estimated period slopes | 8 |
| **Total predeclared families** | **23** |

P7 has no period-slope family, consistent with the P7 constraint in S2.6.

Families were predeclared by boundary and coefficient type so that each FDR correction corresponds to one scientific transition question.

The focused ET, runoff, and storage extension retains its own separate nine-series boundary families and is not combined with the original state and correction-diagnostic families. The extension's families are smaller, so its FDR threshold is less stringent than pooling would produce; the P6 ET level change fails FDR at $q = 0.072$ under this less stringent structure.

## S2.10 Independent changepoint detection with PELT

Independent detection uses `ruptures` 1.1.10 and consumes the same complete area-weighted monthly series written by the known-transition stage, preserving its run matrix, masks, units, paired support, and provenance.

**Seasonal handling.** Calendar-month effects are estimated jointly with one global linear trend, and only the month effects are removed. Detection is therefore performed on a trend-preserving seasonally adjusted response, consistent with S2.3.

**Two methods.**

- *Primary:* piecewise AR(1), intercept, and linear-trend Gaussian profile-likelihood cost.
- *Sensitivity:* Gaussian profile-likelihood linear cost after Prais–Winsten prewhitening.

**Prewhitening persistence.** The prewhitening AR(1) is the median of estimates from overlapping 60-month locally detrended windows stepped every 30 months, so that a long level or slope change is not absorbed as apparent persistence. The ordinary global-trend residual AR(1) is retained as a diagnostic.

**Segment floor.** Each segment must contain at least 24 months.

**Penalties.** The penalty basis is the BIC parameter-count penalty. Multipliers 0.5, 0.75, 1.0, 1.25, 1.5, and 2.0 are all retained in the output; the primary penalty is 1.25 × BIC.

**Acceptance rule.** A changepoint is accepted only when all three conditions hold:

1. it is present in the primary method at 1.25 × BIC;
2. it recurs within three months in at least half of the primary-method penalty grid;
3. a break from the independent prewhitened formulation occurs within three months.

An accepted date is therefore both penalty-stable and cross-method.

**Known-boundary comparison.** Accepted dates are compared one-to-one with P2–P9 using ±3 months as the primary tolerance and ±6 months as a sensitivity. Unmatched accepted breaks are retained as independently detected structural changes and are not forced onto the nearest observing-system boundary. P7 can be displayed and matched but is exempt from agreement scoring, because its 15-month duration is structurally incompatible with the 24-month segment floor.

Per the calibration in S2.11, this configuration recovers abrupt changes well and localizes gradual slope hinges poorly. PELT results are used as corroboration of abrupt transitions; known-date slope inference comes from the interrupted-series model.

## S2.11 Synthetic calibration

Calibration results are tabulated in Table S2. These tests characterize the behavior of the implemented pipeline under the predeclared synthetic scenarios and are not claims of universal operating characteristics.

Scenario configurations are: trend calibration on 288-month synthetic fields with a fixed seasonal cycle, seed 20260812, 100 series per scenario; interrupted-series calibration with 24 AR(1) = 0.8 no-transition series and 24 series for each planted P6 effect, 499 bootstrap replicates; and PELT calibration with 24 series per scenario, seed 20260814.

Under the null scenarios, BH-controlled trend detection and accepted PELT breaks produce no false positives, and the interrupted-series analysis yields one null BH discovery across the 15 tested transition families. Abrupt P6 level shifts are recovered by PELT in 85.4–91.7% of cases within six months. Interval coverage is 91.7% for both planted P6 level and slope effects, and the autocorrelation-adjusted Sen interval covers at 87–97% against 43–56% for the nominal interval. FDR power is approximately 25% for planted level shifts and approximately 8.3% for planted slope shifts, so a non-significant slope coefficient is weak evidence about the absence of a slope change. Gradual P6 slope hinges are recovered by PELT in only 4.2% of cases, with a median localization error of approximately 24.5 months; the absence of a PELT break is therefore not evidence that a gradual slope transition is absent.


---

# Table S1. Statistical estimands

| Scientific question | Primary estimand | Estimator | Dependence handling | Multiple-testing treatment | What the result supports |
|---|---|---|---|---|---|
| Where does assimilation-added snow water go? | Area-weighted fraction of $I_{\mathrm{snow}}$ appearing in each terminal budget term | Ratio of area-weighted aggregate response to area-weighted aggregate input over the snow-addition sample | 5° spatial-block bootstrap, 1,000 replicates (10° sensitivity) | None; a small predeclared set of partition terms | Average mass fate of added water over WY2001–WY2006 |
| Does anomalous snow input predict an anomalous response at the same location? | Controlled marginal response per unit input, $\beta$ | OLS on within-tile anomalies with year fixed effects and OL MAM snow-mass control | 5° spatial-block bootstrap on the full design (10° sensitivity) | None; four closing responses on a common design | Marginal within-location partition; corroboration, not causal proof |
| Does a field drift systematically over 24 years? | Theil–Sen slope of the paired DA−OL monthly series, per year | Exact Theil–Sen median pairwise slope after trend-preserving seasonal adjustment | Hamed–Rao modified Mann–Kendall variance inflation, lags 1–24, factor floored at 1 | BH FDR at 0.05 within each complete field; OL, DA, DA−OL separate families | Presence and sign of a monotonic long-term tendency at each tile |
| Does the series step at a known observing-system date? | Level-change coefficient at that boundary, in source units | Segmented Prais–Winsten AR(1) regression, 28 parameters | Fitted-AR(1) innovation bootstrap, 1,999 replicates, conservative simulation AR(1) | BH FDR within each of 8 level-change boundary families | An abrupt offset in the domain series coincident with the boundary |
| Does the rate of change alter at a known date? | Slope-change coefficient at that boundary, in source units per year | Same segmented model; hinges only where segment length permits | Same innovation bootstrap | BH FDR within each of 7 slope-change boundary families | A change in rate of evolution; low power by construction |
| Does the data alone place a break near that date? | Accepted changepoint date | PELT with piecewise AR(1)+linear cost, confirmed against a prewhitened linear formulation | Piecewise AR(1) cost; locally estimated prewhitening AR(1) | None; acceptance is by penalty-stability and cross-method consensus | Independent corroboration of an abrupt structural change |

These six quantities carry different units and interpretations and are not compared directly to one another.

---

# Table S2. Synthetic calibration

| Method | Null / effect scenario | Serial dependence | Replicates / series | False-positive behavior | Recovery / coverage | Interpretive limitation |
|---|---|---|---:|---|---|---|
| Trend (modified MK + BH) | No trend | White noise | 100 series | 3% pointwise, 0% BH | — | Pointwise stippling would overstate significance |
| Trend (modified MK + BH) | No trend | AR(1) = 0.8 | 100 series | 12% pointwise, 0% BH | — | FDR, not pointwise tests, controls the map |
| Trend (modified MK + BH) | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 100% BH detection | Benchmark effect size; not a general power claim |
| Theil–Sen nominal CI | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 43–56% coverage | Nominal interval is badly undercovered |
| Theil–Sen adjusted CI | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 87–97% coverage | First-order adjustment only; not a test inversion |
| Interrupted series | No transition | AR(1) = 0.8 | 24 series, 499 bootstrap | 1 null BH discovery across 15 families | — | Guard permits ≤2; not a zero-false-positive claim |
| Interrupted series | Planted P6 level change | AR(1) = 0.8 | 24 series, 499 bootstrap | — | 91.7% coverage; ~25% FDR power | Low power; non-detection is not absence |
| Interrupted series | Planted P6 slope change | AR(1) = 0.8 | 24 series, 499 bootstrap | — | 91.7% coverage; ~8.3% FDR power | Very low power for slope hinges |
| PELT | No changepoint | White noise | 24 series | 0 accepted breaks | — | — |
| PELT | No changepoint | AR(1) = 0.8 | 24 series | 0 accepted breaks | — | Guard permits ≤0.25 mean false breaks |
| PELT | Isolated abrupt P6 level shift | AR(1) = 0.8 | 24 series | — | 91.7% within 6 months | — |
| PELT | Two opposing level shifts | AR(1) = 0.8 | 24 series | — | 85.4% within 6 months | Multiple breaks reduce recovery |
| PELT | Level shift near segment edge | AR(1) = 0.8 | 24 series | — | 91.7% within 6 months | — |
| PELT | Gradual P6 slope hinge | AR(1) = 0.8 | 24 series | — | 4.2% within 6 months; median localization error ≈24.5 months | Absence of a PELT break is not evidence against a gradual change |

Cells marked "—" are not applicable to that scenario. No cell is inferred; every populated value appears in a production validation report.

---

# Supplementary figure inventory

Final S-numbering is pending until the main Results are complete. Filenames are current tracked review products under `projects/M21C_ls/docs/paper_figures/` and `projects/M21C_ls/docs/trends_breakpoints_report_figures/`.

| Current filename | Purpose |
|---|---|
| `fig08_supp_ims_terra_aqua_scope_bars.png` | IMS snow-cover categorical skill for the Terra-only versus Terra+Aqua observing-system scopes, as domain-aggregate bars. |
| `fig08_supp_ims_terra_aqua_scope_maps.png` | Spatial distribution of the same Terra-only versus Terra+Aqua IMS skill contrast. |
| `fig_supp_snow_da_octmar_attribution.png` | Non-overlapping October–March snow-DA attribution, showing the sequential control ladder for all six responses in native units. |
| `fig_supp_snow_da_boundary_sensitivity.png` | Water-year accounting-boundary sensitivity of the snow-water partition, comparing October–September with September–August. |
| `figS05_precipitation_trends.png` | Precipitation trend control demonstrating that the common forcing produces no paired DA−OL trend. |
| `figS06_sfmc_trends.png` | Surface soil-moisture trend control for OL, DA, and paired DA−OL. |
| `figS07_snow_mass_depth_trends.png` | Snow-mass and snow-depth trend controls over the seasonal-snow domain. |
| `figS08_breakpoint_boundary_agreement.png` | Agreement between independently detected changepoints and the P2–P9 observing-system boundaries. |
| `p6_level_changes.png` | Detailed P6 level-change estimates across the primary transition series. |
| `p6_soil_water_convergence.png` | Convergence of known-date and independently detected evidence for the April 2015 soil-water transition. |
| `boundary_agreement_matrix.png` | Full accepted-changepoint versus known-boundary agreement matrix. |
| `phase2_flux_storage_trend_maps.png` | Whole-record trend maps for ET, total runoff, and total land-water storage (OL, DA, DA−OL). |
| `phase2_flux_storage_transition_series.png` | DA−OL ET, runoff, and storage domain series with accepted independent changepoints marked. |

Ordering and numbering will be assigned once the main-text figure set is frozen.
