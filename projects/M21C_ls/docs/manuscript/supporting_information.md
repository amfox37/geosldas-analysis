# Supporting Information for

## Evolving Observational Influence and Temporal Consistency in a Multi-Sensor Land Reanalysis, 2000–2024

*Author list to be inserted from the final manuscript front matter.*

This Supporting Information provides additional methodological detail, statistical calibration, sensitivity analyses, and supporting diagnostics for the analyses described in Sections 2.6 and 2.7 of the main manuscript. The definitions and estimands are identical to those used in the main analysis; the material below documents their implementation and robustness rather than introducing alternative primary analyses.

---

# Text S1. Additional methods for snow-DA hydrologic attribution

## S1.1 Differential snow-DA water budget

The differential budget is evaluated over the six complete MODIS-only water years WY2001–WY2006, each running from October through September. This interval precedes the introduction of microwave soil-moisture assimilation in June 2007, so differences between the DA and OL experiments over this period are attributable to snow-cover assimilation rather than to a mixture of snow and soil-moisture analysis increments. The domain is the 48,067-tile Northern Hemisphere seasonal-snow mask used throughout the snow-process analysis, defined by a 20°N southern limit, snow-cover-possible thresholds of 0.05 snow-cover fraction and 5 kg m⁻² snow mass, and a maximum mean JJA snow-cover fraction of 0.20 to exclude permanent ice. All quantities are evaluated in native M36 tile space, and all fluxes are converted to monthly water-equivalent totals before aggregation over the water year.

The snow-DA input, $I_{\mathrm{snow}}$, is the water-year sum of the native prognostic `snow_net` increments, with positive values denoting addition of snow water. The stored monthly increment product is already the sum of native submonthly increments within each month, so no temporal differencing is applied. The closing terms are the DA−OL changes in surface runoff, baseflow, evapotranspiration, and total land-water storage:

$$
I_{\mathrm{snow}} = \Delta R_{\mathrm{surf}} + \Delta R_{\mathrm{base}} + \Delta ET + \Delta S + \epsilon
$$

where $\epsilon$ is the residual. Precipitation does not enter the differential budget because the two experiments use nominally identical forcing. This is verified numerically rather than assumed: the production audit reports the maximum absolute annual tile-level and area-weighted domain-mean DA−OL precipitation difference against a 0.2 kg m⁻² annual tile guard. The observed six-year integrated area-weighted domain difference is −5.5 × 10⁻⁵ kg m⁻², i.e. numerically negligible.

Two classes of quantity are deliberately excluded from the closing equation:

- **Snowmelt and infiltration are transit terms, not sinks.** They describe water moving between reservoirs within the land column rather than leaving the differential budget, and including them alongside runoff and ET would double-count the same water. Both are retained as pathway diagnostics.
- **RZMC and SFMC are diagnostic states already contained in `TWLAND`.** Adding them to $\Delta S$ would double-count soil water. They are retained as timing and impact diagnostics. The compressed monthly inputs contain no separate native prognostic soil-water storage state, so no independent soil-water storage term is available.

$\Delta S$ is calculated from the September-to-September change in the DA−OL monthly-mean `TWLAND` state, not from the integrated `WCHANGELAND` tendency. The tendency diagnostic does not include the discontinuous mass introduced by the analysis update and therefore cannot equal the change in the DA−OL state anomaly; it is retained only as a process-balance diagnostic. Using monthly-mean endpoints introduces a small temporal approximation that is considered when interpreting the residual.

### Partition estimator

For a selected tile-year sample $G$ and response term $k$, the partition fraction is

$$
f_k = \frac{\sum_{(i,t) \in G} A_i\, Y_{k,i,t}}{\sum_{(i,t) \in G} A_i\, I_{\mathrm{snow},i,t}}
$$

where $A_i$ is the M36 tile area, $Y_{k,i,t}$ is the DA−OL response for term $k$ at tile $i$ in water year $t$, and $I_{\mathrm{snow},i,t}$ is the signed snow-water input.

This is a ratio of the area-weighted aggregate response to the area-weighted aggregate input. It is **not** the mean of individual tile-year response-to-input ratios. The distinction is material: individual ratios have $I_{\mathrm{snow},i,t}$ in the denominator and become unstable — and in the limit unbounded — wherever local snow input is small, so their average would be dominated by tile-years carrying almost none of the water. The aggregate ratio weights each tile-year by the water it actually contributes.

The main manuscript emphasizes $G = \{I_{\mathrm{snow}} > 0\}$, the snow-addition sample, comprising 247,545 tile-years across 45,490 tiles. The all-tile sample (288,402 tile-years, 48,067 tiles) and the snow-removal sample (37,395 tile-years, 13,012 tiles) are retained as diagnostics and sensitivities.

## S1.2 Spatial-block uncertainty for the water-budget fractions

Uncertainty in the partition fractions is estimated with a spatial-block bootstrap. Nearby M36 tiles share meteorological and hydrological variability and cannot reasonably be treated as independent observations, so resampling individual tile-years would substantially understate uncertainty.

The implementation is as follows:

- Primary blocks are approximately 5° × 5°, formed by binning tile longitude and latitude into fixed-width cells and assigning each unique bin pair a block code.
- Every tile-year belonging to the same spatial block remains grouped; blocks, not individual tile-years, are the resampling unit.
- Each bootstrap realization draws the original number of spatial blocks with replacement, implemented as a multinomial allocation of block counts.
- The area-weighted numerator and denominator are both recomputed from the resampled block sums, so each realization yields a complete re-estimate of the ratio $f_k$ rather than a perturbation of the point estimate.
- The production analysis uses 1,000 replicates with a fixed seed.
- Reported 95% intervals are the 2.5th and 97.5th percentiles of the replicate distribution.
- A 10° block analysis is retained as a sensitivity.

This is a **spatial** block bootstrap. It resamples geographic regions, not time windows, and should not be described as a temporal block bootstrap.

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

The production design matrix contains an intercept, the $I_{\mathrm{snow}}$ anomaly, the OL MAM snow-mass anomaly, and five year dummy variables relative to the reference water year. The sample is all 288,402 complete signed tile-years, not only the positive-input subsample.

The three controls serve distinct purposes:

**Within-tile anomalies** remove permanent spatial contrasts. Without them, a pooled regression could report that larger snow-DA input accompanies larger runoff simply because mountainous tiles are both snowier and generate more runoff than lowland tiles, with no within-location relationship at all.

**Year fixed effects** remove anomalies shared broadly across the domain within a given water year — for example a generally snowy year in which both snow corrections and runoff responses are elevated nearly everywhere simultaneously.

**OL MAM snow mass** asks whether anomalously large snow-DA input still predicts the response after accounting for whether the unassimilated model itself had an anomalously large spring snowpack at that tile in that year.

Coefficient uncertainty uses the same spatial-block resampling described in S1.2, with 5° blocks primary and 10° blocks as a sensitivity. Because the design matrix and sample are common to all responses, the bootstrap is implemented through per-block sufficient statistics, so every replicate refits the full design.

The response coefficients for total runoff, ET, storage, and the residual **sum to exactly one** in the fitted production calculation. This is a property of the construction, not an empirical finding: those four terms close the same differential budget by definition, and they are fitted on an identical design matrix and sample, so their fitted coefficients on $\tilde{I}$ necessarily sum to the coefficient obtained by regressing $\tilde{I}$ on itself.

## S1.4 Why the controlled coefficients are not identical to the direct fractions

The controlled regression coefficients are not expected to reproduce the direct partition fractions numerically. The direct partition is an area-weighted mass accounting conditional on tile-years for which snow DA adds water. The regression instead uses signed within-tile interannual variation, removes each tile's mean, includes all complete signed tile-years, and controls for common year effects and OL MAM snow mass. Agreement in the qualitative partitioning—particularly the dominance of runoff—therefore provides independent corroboration rather than a duplicate calculation.

Stated compactly: the direct fraction describes the *average fate* of the water that snow assimilation added, while the regression coefficient describes the *marginal partition* of one additional unit of anomalous input at a given location.

The accepted production estimates are:

| Term | Direct positive-input accounting | Controlled regression $\beta$ [95% CI] |
|---|---:|---:|
| Total runoff | 64.3% | 0.749 [0.711, 0.783] |
| Evapotranspiration | 35.9% | 0.182 [0.155, 0.213] |
| Storage | 3.9% | 0.085 [0.075, 0.097] |
| Residual | −4.1% | −0.016 [−0.022, −0.011] |

Direct fractions are from the snow-addition sample (247,545 tile-years); regression coefficients are from the production run with the OL MAM snow-mass control and 5° spatial blocks over all 288,402 signed tile-years. Both columns sum to unity once the residual is included.

Both approaches independently identify runoff as the dominant fate of assimilation-added snow water. They should not be presented as two estimates of the same parameter.

## S1.5 Additional robustness tests

- **10° spatial blocks.** Coarser resampling widens all intervals but does not change any sign or ordering; the positive-input runoff fraction interval widens from 61.1–67.2% to 58.3–69.7%.
- **1st–99th percentile predictor-anomaly trimming.** Applied to the seasonal regression to confirm that the fitted slopes are not driven by extreme input anomalies. All four classified downstream responses retain intervals excluding zero.
- **Non-overlapping October–March analysis.** The predictor is restricted to October–March snow input and responses begin in April or May, giving exactly zero shared months. All four classified downstream responses survive; standardized effects range from 24% (runoff) to 67% (ET) of the corresponding signed-MAM values but remain positive.
- **Water-year accounting-boundary sensitivity.** Repeating the partition on a September–August boundary changes the runoff fraction from 64.3% to 64.5% and ET from 35.9% to 36.0%, both far below the pre-set 5 percentage-point reporting threshold. Storage changes by +0.3 and the residual by −0.6 percentage points. Annual residuals are negative in 6 of 6 years under both boundaries.
- **Snow-removal sample.** Tile-years with $I_{\mathrm{snow}} < 0$ partition differently, with a larger ET share (46.1%) and a positive residual.
- **Soil-moisture persistence diagnostic.** Time from the within-year RZMC peak to a 0.001 m³ m⁻³ threshold is reported for the snow-addition sample.

Three interpretive limitations are stated explicitly:

1. The strict non-overlap test cannot resolve delayed infiltration. Restricting the predictor to October–March necessarily discards input that would have reached the soil column later in the same water year, so the reduced coefficients bound rather than measure the seasonal response.
2. Snow-removal partitioning is qualitatively different from snow-addition partitioning and is not used for the headline result.
3. RZMC persistence is strongly right-censored — 67.2% of snow-addition tile-years never fall below the threshold after their within-year peak before September — and the mean DA−OL RZMC anomaly is already positive in October, so these counts include inherited state differences from prior assimilation. Residence-time estimates are therefore not used as headline results.

---

# Text S2. Additional statistical methods for trends and observing-system transitions

## S2.1 Distinct statistical questions

The temporal-consistency analysis separates three statistical questions. Whole-record trend estimation asks whether a variable exhibits a systematic long-term tendency over June 2000–May 2024. Known-date interrupted-series analysis asks whether the level or slope of an area-weighted monthly series changes at a pre-specified observing-system boundary. Independent changepoint detection asks whether an abrupt structural change can be identified without supplying the observing-system dates. These are complementary rather than interchangeable tests.

The distinction is not merely procedural. A record can contain a large, well-resolved level change at a known date and no significant whole-record trend, because a sequence of discrete offsets is not a monotonic tendency. Conversely, a resolved trend need not imply any identifiable transition. Each result below is reported against the question its estimator actually answers.

## S2.2 Paired DA−OL trend estimand

For paired model variables, the primary assimilation-impact quantity is the monthly paired series

$$
D_i(t) = \mathrm{DA}_i(t) - \mathrm{OL}_i(t)
$$

and the primary trend is the trend in $D_i(t)$. The impact is **not** defined as the difference between separately estimated OL and DA trends.

The paired formulation preserves month-by-month common forcing and identical sample support, so common meteorological variability cancels before any slope is estimated. It also asks the scientifically intended question directly: whether the effect of assimilation itself changes systematically through time. Differencing two independently estimated Theil–Sen slopes would instead compare two robust summaries of two different trajectories, would not cancel common forcing at the monthly scale, and has no straightforward paired significance test.

OL and DA trends are retained as context and control fields. Each complete field keeps its own spatial FDR family; they are fitted on the same underlying paired samples and masks as the corresponding DA−OL field.

## S2.3 Trend-preserving seasonal adjustment

Seasonality is removed with a trend-preserving calendar-month adjustment implemented in four steps:

1. estimate a preliminary exact Theil–Sen slope on the raw series;
2. temporarily remove that slope;
3. estimate each calendar month's climatology from the detrended values;
4. subtract **only** that estimated seasonal climatology from the original, untrended series.

The trend is removed at step 2 solely to obtain an unbiased seasonal climatology, and is restored at step 4. It is not removed before the final test.

The reason for the construction is that subtracting ordinary calendar-month means from a trending series distorts the trend. Because each calendar month is sampled at a different mean position along the trend, raw monthly means absorb part of the long-term tendency, and subtracting them converts a smooth linear trend into an artificial annual step function. The trend-preserving form removes the seasonal cycle while restoring the original long-term tendency.

Calendar months failing the configured minimum climatology sample count (three) are omitted from the adjusted series. Production outputs expose a `calendar_month_used` matrix and per-tile `n_calendar_months_used` so that seasonally selective dropping is visible.

## S2.4 Theil–Sen slope and modified Mann–Kendall test

These are two different tools with two different jobs, and are reported separately.

**Theil–Sen** is the effect-size estimator. It is the exact median of all pairwise slopes between observations, expressed in source units per year. Relative to ordinary least squares it is robust to outliers and makes no Gaussian-residual assumption, which is appropriate for monthly land fields and for zero-heavy increment diagnostics. It supplies magnitude, not significance.

**Mann–Kendall** supplies the significance test. It is a rank-based test of monotonic tendency and does not require Gaussian residuals. It tests whether a monotonic tendency is distinguishable from noise, but does not itself estimate the size of that tendency.

**Hamed–Rao adaptation.** Significance uses a conservative adaptation of the Hamed–Rao modified Mann–Kendall variance correction. Rank autocorrelation is calculated from Sen-detrended residuals at actual monthly lags 1–24, requiring at least 24 pairs per lag and using a 0.05 autocorrelation significance level. Only significantly *positive* lag correlations inflate the variance, and the resulting variance factor is floored at one. With intermittent eligibility, lag correlations pair values that remain separated by the requested calendar-month lag while the Hamed–Rao weight uses the valid observation count; this is a documented gap-aware approximation to the regularly sampled formula.

Positive serial persistence means that 288 monthly observations carry substantially less independent information than 288 independent samples would. The correction accounts for this by inflating the test variance, and by construction can only make a result harder to declare significant, never easier. The uncorrected independent-sample Mann–Kendall p-value is retained in the output solely to show the size of the correction.

## S2.5 Spatial multiple testing and trend confidence intervals

Benjamini–Hochberg false-discovery-rate control at $q = 0.05$ is applied across all finite tiles in each output field. The family structure is explicit:

- each complete mapped field is its own FDR family;
- OL, DA, and DA−OL fields are separate families, even where they share an underlying sample;
- FDR is computed from the complete field, not independently by spatial chunk. Spatial subset runs are written to a distinct diagnostic path and explicitly labeled as non-global FDR; only a complete-mask run supplies production FDR values;
- maps use the `significant_fdr` flag.

Separately, the reported Theil–Sen confidence interval is a **first-order autocorrelation adjustment**, not the exact inversion of the modified Mann–Kendall test. SciPy's nominal 95% Sen limits are retained as `slope_ci_*_nominal`, and the primary interval expands each nominal half-width by the square root of the same Hamed–Rao variance factor used for significance.

Because the interval and the test are related but not mathematically identical procedures, exclusion of zero by the interval and FDR significance can disagree. This occurs most often for tied or zero-heavy increment diagnostics. Production outputs expose an explicit `fdr_ci_disagreement` flag wherever the FDR test is significant but the adjusted Sen interval contains zero. In the Phase 1 production batch, disagreement is absent for the paired state fields and the snow sensitivities, but reaches 22.1% of significant tiles for absolute soil-water activity and 8.9–13.6% for the SFMC correction diagnostics.

This disagreement is reported rather than removed. Neither inferential product is adjusted to mimic the other, and the first-order interval is not plotted as inferential whiskers for the activity fields. Mapped inferential significance is based on the FDR-controlled test. A test-inverted or block-bootstrap interval remains the appropriate later sensitivity if interval-based mapped inference is required.

## S2.6 Known-date interrupted time-series model

The response is an area-weighted monthly domain series over the complete 288-month record. Tile weights are M36 land-tile areas from the GEOS LDAS tile-coordinate file; missing OL or DA values are cross-masked before any paired mean or difference is formed.

The segmented design contains:

| Component | Count |
|---|---:|
| Intercept | 1 |
| P1 baseline slope | 1 |
| Calendar-month fixed effects | 11 |
| Level changes at P2–P9 | 8 |
| Slope hinges at P2–P6, P8, P9 | 7 |
| **Total fitted parameters** | **28** |

The design is full rank over 288 months. P7 spans only 15 months: it receives a level-change term but no independently fitted slope hinge, and the P6 slope is carried until the P8 boundary. This constraint is inherited from the period registry and fixed before estimation; it is a design decision, not an inference from the data.

Three coefficient types appear in the output and must not be conflated:

- **Level change** — the estimated immediate vertical displacement at a boundary, after accounting for seasonality and the piecewise trend. Reported in source units.
- **Slope change** — the estimated change in the rate of evolution across a boundary. Reported in source units per year.
- **Period slope** — the resulting slope within a defined period. Reported in source units per year.

A period slope can be significant because of shared forcing or climate variability with no assimilation impact at all; several significant control period slopes appear in both OL and DA and should not be read as DA effects.

## S2.7 Prais–Winsten AR(1)

Monthly residuals from a segmented fit are serially correlated. Ordinary regression assumes residual independence and consequently underestimates coefficient uncertainty when that assumption fails.

Serial dependence is fitted with iterative Prais–Winsten AR(1): the first-order residual persistence is estimated, the regression is transformed accordingly, the first observation is retained rather than discarded, and the procedure is refitted until the AR(1) estimate stabilizes. The production configuration allows at most 25 iterations with a convergence tolerance of 10⁻⁷ and caps the absolute AR(1) at 0.98.

Prais–Winsten is **not** the final uncertainty estimator in this analysis. It supplies the fitted serial-dependence model that the production bootstrap then uses. This distinction matters: both OLS with Newey–West standard errors and Prais–Winsten with Newey–West standard errors were anti-conservative in the fixed-seed 288-month AR(1) no-transition test, whereas the bootstrap passed the BH false-discovery guard.

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

The complete model is refitted on every replicate, so the bootstrap distribution reflects estimation uncertainty in all 28 parameters jointly rather than in one coefficient conditioned on the others.

The simulation AR(1) parameter is deliberately conservative: it is the **upper 95% large-sample confidence bound** on the fitted AR(1), capped at 0.98. Both the fitted value and the simulation value are reported. A 28-parameter segmented model has enough flexibility to absorb part of the true residual persistence into the deterministic fit, biasing the estimated AR(1) downward; simulating at the upper bound makes the resulting intervals conditional on a conservative persistence estimate and reduces the risk of anti-conservative transition inference.

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

The rationale is that each family answers one specific scientific question — "does anything change at this boundary, across the registered response series?" — and correction should be applied within that question. Pooling all 1,591 fitted coefficient rows into a single omnibus FDR family would mix unrelated questions, combining calendar-month nuisance terms and baseline slopes with the transition coefficients of interest, and would not correspond to any question the analysis asks.

The focused ET, runoff, and storage extension retains its **own** separate nine-series boundary families and is not combined with the original state and correction-diagnostic families. This separation is stated explicitly because it affects interpretation: the extension's families are smaller, so its FDR threshold is less stringent than pooling would produce. Where a coefficient in the extension fails FDR — as the P6 ET level change does at $q = 0.072$ — it fails under the *more* permissive of the two available family structures.

## S2.10 Independent changepoint detection with PELT

Independent detection uses `ruptures` 1.1.10 and does not reimplement PELT. It consumes the same complete area-weighted monthly series written by the known-transition stage, preserving its run matrix, masks, units, paired support, and provenance.

**Seasonal handling.** Calendar-month effects are estimated jointly with one global linear trend, and only the month effects are removed. Detection is therefore performed on a trend-preserving seasonally adjusted response, consistent with S2.3.

**Two methods.**

- *Primary:* piecewise AR(1), intercept, and linear-trend Gaussian profile-likelihood cost.
- *Sensitivity:* Gaussian profile-likelihood linear cost after Prais–Winsten prewhitening.

**Prewhitening persistence.** The prewhitening AR(1) is the median of estimates from overlapping 60-month locally detrended windows stepped every 30 months. Estimating persistence locally prevents a long level or slope change from being absorbed as apparent persistence. The ordinary global-trend residual AR(1) is retained as a diagnostic.

**Segment floor.** Each segment must contain at least 24 months.

**Penalties.** The penalty basis is the BIC parameter-count penalty. Multipliers 0.5, 0.75, 1.0, 1.25, 1.5, and 2.0 are all retained in the output; the primary penalty is 1.25 × BIC.

**Acceptance rule.** A changepoint is accepted only when all three conditions hold:

1. it is present in the primary method at 1.25 × BIC;
2. it recurs within three months in at least half of the primary-method penalty grid;
3. a break from the independent prewhitened formulation occurs within three months.

An accepted date is therefore both penalty-stable and cross-method, not the output of one favorable tuning.

**Known-boundary comparison.** Accepted dates are compared one-to-one with P2–P9 using ±3 months as the primary tolerance and ±6 months as a sensitivity. Unmatched accepted breaks are retained as independently detected structural changes and are not forced onto the nearest observing-system boundary. P7 can be displayed and matched but is exempt from agreement scoring, because its 15-month duration is structurally incompatible with the 24-month segment floor.

As the calibration in S2.11 shows, this configuration corroborates strong abrupt changes well and localizes gradual slope hinges poorly. PELT results are accordingly used as corroboration of abrupt transitions; known-date slope inference remains the job of the interrupted-series model.

## S2.11 Synthetic calibration

The values below are fixed-seed regression benchmarks from the production validation runs. They characterize the configurations actually used and are not claims of universal operating characteristics.

**Trend calibration** (288-month synthetic fields, fixed seasonal cycle, seed 20260812, 100 series per scenario):

- white-noise no-trend: 3% pointwise modified-MK positives, 0% BH detections;
- AR(1) = 0.8 no-trend: 12% pointwise positives, 0% BH detections;
- AR(1) = 0.8 with imposed ±0.02 units yr⁻¹ trend: 100% BH detection;
- nominal Sen confidence-interval coverage 43–56% in the AR(1) scenarios;
- first-order autocorrelation-adjusted interval coverage 87–97%.

The contrast between 12% pointwise and 0% BH false positives under strong autocorrelation is the practical reason mapped significance uses FDR rather than pointwise stippling. The contrast between 43–56% and 87–97% coverage is the reason the adjusted interval replaces the nominal one.

**Interrupted-series calibration** (24 AR(1) = 0.8 no-transition series and 24 series for each planted P6 effect, 499 bootstrap replicates):

- one null BH discovery across the 15 tested transition families;
- 91.7% interval coverage for both planted P6 level and slope effects;
- correct mean-effect direction and relative mean bias below 50%;
- FDR power approximately 25% for the planted level shift;
- FDR power approximately 8.3% for the planted slope shift.

The inference is therefore deliberately conservative and is particularly low-powered for slope changes. A non-significant slope coefficient is weak evidence about the absence of a slope change.

**PELT calibration** (24 series per scenario, seed 20260814):

- zero accepted breaks in the seasonal white-noise null;
- zero accepted breaks in the AR(1) = 0.8 null;
- 91.7% recovery within six months for an isolated abrupt P6 level shift;
- 85.4% recovery across two opposing level shifts;
- 91.7% recovery for a level shift near the minimum-segment edge;
- 4.2% recovery within six months for a gradual P6 slope hinge, with a median localization error of approximately 24.5 months.

Independent changepoint evidence is therefore meaningful for abrupt breaks. The absence of a PELT break is **not** evidence that a gradual slope transition is absent.

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

A partition fraction, a regression slope, a Theil–Sen trend, an interrupted-series level coefficient, a slope hinge, and a PELT date are six different quantities with six different units and interpretations. They are never compared directly to one another.

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

Final S-numbering is **pending** until the main Results are complete. Filenames below are current tracked review products under `projects/M21C_ls/docs/paper_figures/` and `projects/M21C_ls/docs/trends_breakpoints_report_figures/`.

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

Ordering and final numbering will be assigned once the main Results section is complete and the main-text figure set is frozen.
