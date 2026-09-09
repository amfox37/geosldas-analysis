# Supporting Information for

## Evolving Observational Influence and Temporal Consistency in a Multi-Sensor Land Reanalysis, 2000–2024

*Author list to be inserted from the final manuscript front matter.*

This Supporting Information provides additional methodological detail, statistical calibration, sensitivity analyses, and supporting diagnostics for the analyses described in Sections 2.5–2.7 of the main manuscript, covering the evaluation datasets and metrics, the snow-DA hydrologic attribution, and the long-term trends and observing-system period differences. The definitions and estimands are identical to those used in the main analysis; the material below documents their implementation and robustness rather than introducing alternative primary analyses.

---

# Text S1. Additional evaluation datasets and metrics

This section documents the station matching, quality control, regridding, uncertainty estimation, and unit conversions summarized in Section 2.5. Each evaluation compares OL and DA on a common sample, so the reported differences reflect assimilation rather than differing support.

## S1.1 ISMN soil-moisture validation

In-situ soil moisture was taken from a local archive of the International Soil Moisture Network (ISMN), read with the Python `ismn` package version 1.5.2. Six networks carry sufficient coverage for the analysis windows: SNOTEL, SCAN, USCRN, SMOSMANIA, OzNet, and ARM. SNOTEL, SCAN, and ARM are evaluated over all three validation periods; USCRN, SMOSMANIA, and OzNet over the latter two. The validation periods are the SCF-only period P1–P2 (1 June 2000 to 31 May 2007), the pre-SMAP microwave period P3–P5 (1 June 2007 to 31 March 2015), and the SMAP-era microwave period P6–P9 (1 April 2015 to 1 June 2024).

Quality control retained only soil-moisture records carrying the ISMN quality flag `G`. Because ISMN geophysical quality flags for frozen conditions, including in-situ or GLDAS soil and air temperature below 0 °C, are excluded by that requirement, no separate model-based snow-cover or frozen-soil mask is applied in this station workflow. Duplicate timestamps were averaged before the records were averaged to daily means, and the daily observation timestamps were shifted by 12 h to align with the GEOS LDAS daily files.

Surface soil moisture was taken from the available sensor depth nearest 0.05 m. Root-zone soil moisture used network-specific strict multilayer composites, which require the full expected set of layer observations at each timestamp; timestamps with any required layer missing yield a missing root-zone value. ARM contributes surface results only, as no root-zone composite is defined for that network.

Each station was matched to the nearest M36 tile by squared latitude/longitude degree distance with a maximum of 0.1 degree². The workflow also stores the corresponding great-circle distance, but that quantity is not the selection criterion. A station was retained only where at least 1,000 valid paired daily samples were available in each required analysis window, which enforces common station membership across the periods compared for a given network, layer, and metric so that period comparisons are not driven by station turnover.

Skill was summarized with Pearson correlation (R), anomaly correlation (anomR), and unbiased RMSE (ubRMSE). Anomalies were formed by removing a day-of-year climatology estimated with a 31-day moving window, applied circularly across the year end, and requiring at least 30 valid samples for each daily climatological estimate; days failing that minimum yield missing anomalies. Plotted values are paired network-mean differences, DA−OL for R and anomR and OL−DA for ubRMSE, so that positive values indicate DA improvement for all three metrics. The plotted 95% intervals are Student-*t* intervals formed from the across-station standard error of the paired station differences.

The common-site counts underlying the manuscript figure are 34 surface and 14 root-zone sites for SNOTEL, 70 and 39 for SCAN, 88 and 55 for USCRN, 20 and 20 for SMOSMANIA, 16 and 10 for OzNet, and 4 and 0 for ARM. One point departs from these counts: the SNOTEL root-zone anomaly correlation in the SCF-only period rests on 13 rather than 14 stations, because the day-of-year climatology is not estimable at one station within that window.

## S1.2 IMS snow-cover validation

Daily snow-cover validation used the NOAA/National Ice Center Interactive Multisensor Snow and Ice Mapping System (IMS) Daily Northern Hemisphere Snow and Ice Analysis, distributed by NSIDC as data set G02156 Version 1. The 24-km daily product was used; its ASCII files carry the IMS file-version tag `v1.3`.

IMS categories were mapped to the `SMAP_EASEv2_M36_GLOBAL` grid by categorical-safe nearest-neighbor assignment, with no interpolation, smoothing, or averaging of categories. One deterministic representative GEOS LDAS tile was selected per M36 cell using the largest `frac_cell`, with ties broken by the smallest tile ID, and the nearest IMS source pixel to that representative tile was accepted only within a maximum source-target distance of 60 km.

Categorical values were converted to binary snow/no-snow observations before scoring: IMS code 4 denotes snow, finite non-fill codes 0–3 denote no snow, and fill values are excluded. A model cell was classified as snow-covered where snow-cover fraction was at least 0.5. Scores were accumulated on a common daily valid mask requiring finite IMS binary snow cover and finite OL and DA snow-cover fraction at a valid land cell, so OL and DA use the same daily sample. Cells with fewer than 10 IMS observed-snow days over the analysis record were excluded from the map and table domain.

For each experiment the analysis accumulates hits, false alarms, misses, and correct rejections, and computes accuracy, hit rate, miss rate, false alarm ratio, and correct rejection rate. In the manuscript figure the difference for each metric is oriented so that positive indicates improvement: DA−OL for accuracy, hit rate, and correct rejection rate, and OL−DA for miss rate and false alarm ratio.

Confidence intervals for the summary table were generated in the production postprocessing mode, which rebuilds the tables from stored per-cell contingency counts and therefore resamples cells rather than days. Representative cells were resampled with replacement over 2,000 replicates, contingency counts were reaggregated, metrics recomputed, and the 2.5th and 97.5th percentiles reported. This is recorded in the metadata of the archived per-cell counts product used for the manuscript table.

The Terra-only versus Terra+Aqua comparison supporting the Section 3.4 hit-rate discussion uses the same IMS product and scoring, evaluated over the P1 scope (1 June 2000 to 30 June 2002) and the P2 scope (1 July 2002 to 31 May 2007) and restricted to the metric-specific cells valid in both periods, with the snow-domain mask inherited from the upstream IMS product. Comparing the two scopes on that common support characterizes the change in categorical snow skill associated with expanding from Terra-only to Terra+Aqua SCF assimilation; it does not isolate a causal Aqua effect, because the two periods also differ in their meteorology and in the state of the model.

## S1.3 SNOTEL SWE validation

Snow water equivalent was evaluated against daily `WTEQ` observations from the USDA NRCS SNOTEL network, obtained through the Air and Water Database (AWDB). Observations were converted from inches of water equivalent to kg m⁻² by multiplying by 25.4, and compared with model `SNOMASLAND` by exact daily date.

Each station was collocated to the nearest M36 tile by great-circle distance with a maximum station-tile separation of 40 km. Rows lacking station identifier, date, or coordinates were removed, duplicate station-day records were averaged, model fill values greater than 1 × 10¹⁴ and negative model snow values were set to missing, and statistics were computed from finite paired samples only. No frozen-ground or snow-threshold filter was applied, because the target variable is snow itself; instead the all-season statistic excludes JJA to prevent summer no-snow conditions from dominating the aggregate. The plotted seasons are ALL (excluding JJA), SON, DJF, and MAM.

Station-level RMSE, ubRMSE, and absolute bias were computed for each experiment and season. The top-row bars are station-mean values with 95% confidence intervals from a 1,000-replicate station bootstrap using percentile limits and a fixed seed. Those bars use all stations with a finite station-level metric for the relevant season and experiment.

A tighter representativeness filter applies only to the mapped station-difference panels, which are restricted to stations where the absolute difference between tile and station elevation is less than 500 m and both OL and DA station metrics are finite; 611 stations meet that requirement. This elevation filter is not applied to the top-row station-mean bars. In the manuscript figure the mapped and summarized error differences are shown as OL−DA, so positive values indicate lower DA error; earlier SNOTEL output files used the opposite DA−OL convention.

## S1.4 GHCN-Daily snow-depth validation

Snow depth was evaluated against the `SNWD` element of the NOAA/NCEI Global Historical Climatology Network Daily (GHCN-Daily) archive. A local station-level snapshot was built from the GHCN-Daily `by_year` files on 19 March 2026, retaining only rows with element `SNWD`, a non-missing data value, and a blank quality flag. The validation period is 1 January 2000 to 31 December 2024.

Stations were required to have finite coordinates, at least 1,500 valid SNWD observation days, and at least 5 observed-snow days per year on average. Each station was collocated to the nearest M36 tile by great-circle distance and retained only where that tile lay within 40 km and the absolute tile-minus-station elevation difference was no greater than 500 m, with common station membership required for OL and DA.

Model snow depth was formed as the grid-cell-average depth `SNODPLAND` × SCF and converted to the millimetre units of the observations. Model fill values and negative depths were set to missing, snow-cover fraction was clipped to [0, 1], and metrics were computed from finite paired samples with observed and modelled depth in the physical range 0 to 5000 mm.

The manuscript uses the baseline-core analysis, which scores only days actually reported by GHCN and does not fill missing observation days with zero. A common OL/DA filter applied before the metrics required at least 1,500 paired days and at least 5 observed-snow days per year on average, retaining 9,451 Northern Hemisphere stations. Station-level RMSE, ubRMSE, and absolute bias were averaged across stations, with 95% confidence intervals from a 1,000-replicate station bootstrap using percentile limits and a fixed seed. The all-season statistic excludes JJA and the plotted seasons are ALL, DJF, MAM, and SON. The manuscript figure shows OL−DA, so positive values indicate lower DA error.

## S1.5 ERA5-Land consistency comparison

The reanalysis comparison uses monthly averaged ERA5-Land fields from the Copernicus Climate Data Store. The workflow can also generate an ERA5 comparison, but the manuscript figures use ERA5-Land throughout. This comparison is a large-scale consistency check against an independent reanalysis, not independent validation: ERA5-Land is a model product and is not treated as observational truth.

ERA5-Land fields were regridded to the M36 grid with xESMF `conservative_normed` remapping and periodic longitude handling, and the model and reference records were intersected by monthly period, giving 288 aligned monthly values from June 2000 through May 2024. GEOS LDAS tile values were aggregated to M36 cells as in the strict production workflow, and cells with no finite tile value were masked.

Surface soil moisture compares GEOS LDAS `SFMC` with ERA5-Land `swvl1`. Root-zone soil moisture compares `RZMC` with a 0–100 cm estimate formed as 0.07·`swvl1` + 0.21·`swvl2` + 0.72·`swvl3`. Soil-moisture metrics are computed under a joint warm, snow-free mask that requires, in both the model and ERA5-Land, a surface soil temperature above 275.15 K and a snow-cover fraction below 0.01. Soil-moisture metrics are reported over the same three broader validation periods used for ISMN, P1–P2, P3–P5, and P6–P9, and comprise R, anomaly correlation after removing the month-of-year climatology within each period, and ubRMSE.

Snow variables compare ERA5-Land `snowc`, `sd`, and `sde` with model `FRLANDSNO`, `SNOMASLAND` divided by 1000, and `SNODPLAND` × `FRLANDSNO` as appropriate. The snow comparison does not use the seasonal-snow mask defined for the snow-budget analysis; it uses the production ever-snow domain, which includes an M36 cell where either GEOS LDAS or ERA5-Land has a snow-cover fraction above 1 × 10⁻⁶ in at least one month of the record. ERA5-Land no-snow missing values in SWE and snow depth are treated as zero in the final plotting workflow. Snow metrics are RMSE, ubRMSE, and bias, and a grid-cell metric is retained only where at least 24 valid paired monthly samples are available.

The bar summaries are spatial means over the finite M36 cells contributing to each metric. Land-fraction weighting from the tile-coordinate `frac_cell` field is applied specifically to the map-panel annotations, and is not applied to every spatial mean reported from this comparison.

---

# Text S2. Additional methods for snow-DA hydrologic attribution

## S2.1 Differential snow-DA water budget

The differential budget is evaluated over the six complete MODIS-only water years WY2001–WY2006, each running from October through September. This interval precedes microwave soil-moisture assimilation in June 2007, so DA−OL differences are attributable to snow-cover assimilation alone. The domain is the 48,067-tile Northern Hemisphere seasonal-snow mask, defined by a 20°N southern limit, snow-cover-possible thresholds of 0.05 snow-cover fraction and 5 kg m⁻² snow mass, and a maximum mean JJA snow-cover fraction of 0.20 to exclude permanent ice. All quantities are evaluated in native M36 tile space, and fluxes are converted to monthly water-equivalent totals before aggregation over the water year.

The snow-DA input, $I_{\mathrm{snow}}$, is the water-year sum of the native prognostic `snow_net` increments, with positive values denoting addition of snow water. The stored monthly increment product is already the sum of native submonthly increments within each month, so no temporal differencing is applied. The closing terms are the DA−OL changes in surface runoff, baseflow, evapotranspiration, and total land-water storage:

$$
I_{\mathrm{snow}} = \Delta R_{\mathrm{surf}} + \Delta R_{\mathrm{base}} + \Delta ET + \Delta S + \Delta W_{\mathrm{peat}} + \epsilon
$$

where $\epsilon$ is the residual. Precipitation does not enter the differential budget because the two experiments use nominally identical forcing. This is verified rather than assumed: the production audit checks the maximum absolute annual tile-level and area-weighted domain-mean DA−OL precipitation difference against a 0.2 kg m⁻² annual tile guard, and the six-year integrated area-weighted domain difference is −5.5 × 10⁻⁵ kg m⁻².

Two classes of quantity are excluded from the closing equation to avoid double counting. Snowmelt and infiltration redistribute water within the land column rather than removing it from the differential budget, and are retained as pathway diagnostics rather than terminal sinks. RZMC and SFMC are diagnostic states already contained in `TWLAND` and are retained as timing and impact diagnostics; the compressed monthly inputs contain no separate native prognostic soil-water storage state, so no independent soil-water storage term is available.

$\Delta S$ is calculated from the change in the DA−OL `TWLAND` state between instantaneous 00Z 1 October endpoints, not from the integrated `WCHANGELAND` tendency. The tendency diagnostic excludes the discontinuous mass introduced by the analysis update and therefore cannot equal the change in the DA−OL state anomaly; it is retained only as a process-balance diagnostic. The endpoints are reconstructed from `catch_internal_rst` restart files as the 24-member ensemble mean of the prognostic storage components. The earlier September monthly-mean proxy is retained as a comparison diagnostic: it agrees to about 2% in six-year aggregate but differs by up to about 20% in individual water years and inverts sign where the true change is small.

$\Delta W_{\mathrm{peat}}$ is the DA−OL change in PEATCLSM free-standing surface water. `TWLAND` is assembled from the soil, canopy-interception, and snow stores and by construction excludes this store, so on peat tiles (defined by a porosity of at least 0.90) water moving into or out of free-standing water leaves the `TWLAND`-based budget entirely and would otherwise be absorbed into the residual. The term is identically zero on non-peat tiles, matching the model's own initialization. Both closing terms are verified against an independently constructed absolute budget for each experiment: the production run asserts that the differential budget equals DA minus OL of the per-run balances, term by term, to within 5 × 10⁻³ kg m⁻².

### Partition estimator

For a selected tile-year sample $G$ and response term $k$, the partition fraction is

$$
f_k = \frac{\sum_{(i,t) \in G} A_i\, Y_{k,i,t}}{\sum_{(i,t) \in G} A_i\, I_{\mathrm{snow},i,t}}
$$

where $A_i$ is the M36 tile area, $Y_{k,i,t}$ is the DA−OL response for term $k$ at tile $i$ in water year $t$, and $I_{\mathrm{snow},i,t}$ is the signed snow-water input.

Partition fractions are ratios of aggregate area-weighted response to aggregate area-weighted input rather than averages of tile-year ratios, avoiding unstable ratios where local snow input approaches zero.

The main manuscript emphasizes $G = \{I_{\mathrm{snow}} > 0\}$, the snow-addition sample, comprising 247,545 tile-years across 45,490 tiles. The all-tile sample (288,402 tile-years, 48,067 tiles) and the snow-removal sample (37,395 tile-years, 13,012 tiles) are retained as diagnostics and sensitivities.

## S2.2 Spatial-block uncertainty for the water-budget fractions

Partition uncertainty was estimated with 1,000 spatial-block bootstrap replicates. M36 tiles were assigned to approximately 5° × 5° latitude–longitude blocks, and complete blocks, including all associated tile-years, were resampled with replacement. For each realization, the area-weighted response and snow-input sums were recomputed and their ratio retained. The 2.5th and 97.5th percentiles define the reported 95% interval. A 10° × 10° block analysis was used as a sensitivity.

Blocks are formed by binning tile longitude and latitude into fixed-width cells and assigning each unique bin pair a block code; each realization draws the original number of blocks as a multinomial allocation of block counts under a fixed seed. Resampling is spatial, not temporal.

## S2.3 Controlled water-year regression

The controlled regression is fitted as

$$
\tilde{Y}_{i,t} = \alpha + \beta\, \tilde{I}_{i,t} + \delta\, \tilde{S}^{\mathrm{OL}}_{i,t} + \gamma_t + \varepsilon_{i,t}
$$

where a tilde denotes deviation from that tile's six-water-year mean, and

- $\tilde{Y}_{i,t}$ is the within-tile anomaly of the signed DA−OL terminal response (total runoff, ET, storage, peatland free-standing water, or residual);
- $\tilde{I}_{i,t}$ is the within-tile anomaly of signed water-year snow-water input;
- $\tilde{S}^{\mathrm{OL}}_{i,t}$ is the within-tile anomaly of the day-weighted March–May mean OL snow mass (`SNOMASLAND`), weighted by days in month;
- $\gamma_t$ are water-year fixed effects;
- $\beta$ is the controlled response per unit anomalous snow-water input.

Within-tile anomalies remove persistent spatial differences among locations, year fixed effects remove domain-wide interannual anomalies, and the OL MAM snow-mass anomaly controls for interannual variations in the underlying spring snowpack. The production design contains an intercept, the snow-input anomaly, the OL MAM snow-mass anomaly, and five year indicators. The sample is all 288,402 complete signed tile-years, not only the positive-input subsample.

Coefficient uncertainty uses the same spatial-block resampling described in S2.2, with 5° blocks primary and 10° blocks as a sensitivity. Because the design matrix and sample are common to all responses, the bootstrap is implemented through per-block sufficient statistics, so every replicate refits the full design.

The response coefficients for total runoff, ET, storage, peatland free-standing water, and the residual sum to exactly one. This is a property of the construction rather than an empirical finding: these five terms close the same differential budget by definition and are fitted on an identical design matrix and sample, so their coefficients on $\tilde{I}$ necessarily sum to the coefficient from regressing $\tilde{I}$ on itself.

## S2.4 Why the controlled coefficients are not identical to the direct fractions

The controlled coefficients and direct partition fractions estimate related but different quantities. The direct calculation is an area-weighted mass partition conditional on positive snow input, whereas the regression uses signed within-tile interannual variation with year and OL-snow controls. Numerical agreement is therefore not expected; agreement in the response ordering provides independent corroboration.

The accepted production estimates are:

| Term | Direct positive-input accounting | Controlled regression $\beta$ [95% CI] |
|---|---:|---:|
| Total runoff | 64.3% | 0.749 [0.711, 0.783] |
| Evapotranspiration | 35.9% | 0.182 [0.155, 0.213] |
| Storage | 4.2% | 0.085 [0.074, 0.097] |
| Peatland free-standing water | −2.7% | −0.017 [−0.022, −0.012] |
| Residual | −1.7% | 0.0007 [0.0002, 0.0012] |

Direct fractions use the snow-addition sample (247,545 tile-years); regression coefficients use the production run with the OL MAM snow-mass control and 5° spatial blocks over all 288,402 signed tile-years. Both columns sum to unity once the residual is included.

Both approaches identify runoff as the dominant response.

## S2.5 Additional robustness tests

- **10° spatial blocks.** Coarser resampling widens all intervals but does not change any sign or ordering; the positive-input runoff fraction interval widens from 61.1–67.4% to 58.3–69.7%.
- **1st–99th percentile predictor-anomaly trimming.** Applied to the seasonal regression to confirm that the fitted slopes are not driven by extreme input anomalies. The classified downstream responses retain intervals excluding zero.
- **Non-overlapping October–March analysis.** The predictor is restricted to October–March snow input and responses begin in April or May, giving exactly zero shared months. The classified downstream responses survive; standardized effects range from 24% (runoff) to 67% (ET) of the corresponding signed-MAM values but remain positive.
- **Snow-removal sample.** Tile-years with $I_{\mathrm{snow}} < 0$ partition differently: total runoff 38.3%, ET 46.1%, storage 9.9%, peatland free-standing water 0.7%, and residual 5.1%. These fractions are signed response divided by signed input, so the positive percentages describe same-direction water removal. This sample is reported as a sensitivity and is not used for the headline positive-input result.
- **Soil-moisture persistence diagnostic.** Time from the within-year RZMC peak to a 0.001 m³ m⁻³ threshold is reported for the snow-addition sample.

Three interpretive limitations apply:

1. The strict non-overlap test cannot resolve delayed infiltration; restricting the predictor to October–March discards input reaching the soil column later in the same water year, so the reduced coefficients bound rather than measure the seasonal response.
2. Snow-removal partitioning is qualitatively different from snow-addition partitioning and is not used for the headline result.
3. RZMC persistence is strongly right-censored — 67.2% of snow-addition tile-years never fall below the threshold after their within-year peak before September — and the mean DA−OL RZMC anomaly is already positive in October, so these counts include inherited state differences from prior assimilation. Residence-time estimates are not used as headline results.

---

# Text S3. Additional statistical methods for long-term trends and observing-system period differences

## S3.1 Distinct statistical questions

The temporal analysis addresses two distinct questions: whether a variable exhibits a systematic whole-record trend, and when the regional differences identified by that trend analysis emerge within the P1–P9 observing-system chronology. These are estimated separately because a change confined to part of the record need not produce a significant whole-record trend, and a whole-record trend need not be concentrated in any single period.

## S3.2 Paired DA−OL trend estimand

For paired model variables, the primary assimilation-impact quantity is the monthly paired series

$$
D_i(t) = \mathrm{DA}_i(t) - \mathrm{OL}_i(t)
$$

The primary assimilation-impact trend is estimated directly from the paired DA−OL series. OL and DA are restricted to common finite monthly support before differencing; their separate trends are retained as contextual control fields. OL, DA, and DA−OL maps retain separate complete-field FDR families.

## S3.3 Trend-preserving seasonal adjustment

Seasonality is removed with a trend-preserving calendar-month adjustment implemented in four steps:

1. estimate a preliminary exact Theil–Sen slope on the raw series;
2. temporarily remove that slope;
3. estimate each calendar month's climatology from the detrended values;
4. subtract **only** that estimated seasonal climatology from the original, untrended series.

This construction removes the calendar-month climatology without absorbing part of the underlying linear trend into the seasonal adjustment. The trend removed at step 2 is restored at step 4 and is not removed before the final test.

Calendar months failing the configured minimum climatology sample count (three) are omitted from the adjusted series. Production outputs expose a `calendar_month_used` matrix and per-tile `n_calendar_months_used` so that seasonally selective dropping is visible.

## S3.4 Theil–Sen slope and modified Mann–Kendall test

Trend magnitude was estimated with the exact Theil–Sen median slope, the median of all pairwise slopes between observations, expressed in source units per year. Significance was assessed separately using the Mann–Kendall test with the stated Hamed–Rao autocorrelation adjustment.

The Hamed–Rao implementation is a conservative adaptation of the modified Mann–Kendall variance correction. Rank autocorrelation is calculated from Sen-detrended residuals at actual monthly lags 1–24, requiring at least 24 pairs per lag and using a 0.05 autocorrelation significance level. Only significantly *positive* lag correlations inflate the variance, and the resulting variance factor is floored at one. With intermittent eligibility, lag correlations pair values that remain separated by the requested calendar-month lag while the Hamed–Rao weight uses the valid observation count; this is a documented gap-aware approximation to the regularly sampled formula.

By construction the correction can only make a result harder to declare significant. The uncorrected independent-sample Mann–Kendall p-value is retained in the output to show the size of the correction.

## S3.5 Spatial multiple testing and trend confidence intervals

Benjamini–Hochberg false-discovery-rate control at $q = 0.05$ is applied across all finite tiles in each output field. The family structure is explicit:

- each complete mapped field is its own FDR family;
- OL, DA, and DA−OL fields are separate families, even where they share an underlying sample;
- FDR is computed from the complete field, not independently by spatial chunk. Spatial subset runs are written to a distinct diagnostic path and explicitly labeled as non-global FDR; only a complete-mask run supplies production FDR values;
- maps use the `significant_fdr` flag.

Separately, the reported Theil–Sen confidence interval is a **first-order autocorrelation adjustment**, not the exact inversion of the modified Mann–Kendall test. SciPy's nominal 95% Sen limits are retained as `slope_ci_*_nominal`, and the primary interval expands each nominal half-width by the square root of the same Hamed–Rao variance factor used for significance.

The Sen confidence interval and modified-Mann–Kendall/FDR decision are separate inferential products and need not agree exactly; mapped significance is defined by the FDR-controlled test. Production outputs expose an explicit `fdr_ci_disagreement` flag wherever the FDR test is significant but the adjusted Sen interval contains zero. Disagreement is absent for the paired state fields and the snow sensitivities reported here.

## S3.6 Regional timing decomposition of RZMC trends

The whole-record trend analysis establishes where assimilation alters root-zone
soil moisture. This analysis asks when those regional differences emerge within
the P1–P9 chronology.

Six fixed domains are used: global valid land and five lat/lon boxes.

| Region | Latitude | Longitude | Tiles |
|---|---|---|---:|
| Global valid land | — | — | 112,573 |
| Australia | −45 to −10 | +112 to +154 | 6,417 |
| Southern Africa | −35 to −15 | +10 to +41 | 3,253 |
| North Africa / Middle East | +12 to +37 | −18 to +60 | 14,063 |
| Western North America | +30 to +60 | −125 to −100 | 4,555 |
| Northern Eurasia | +50 to +70 | +60 to +140 | 7,461 |

Each region is the set of `valid_land` tiles inside its box, intersected with
the tiles finite in both OL and DA in every month. All 112,573 valid-land tiles
are finite in all 288 months, so the contributing tile count and represented
area are constant throughout the record in every region and no apparent change
can arise from varying spatial support. Series are area-weighted by M36 tile
area. The five regional boxes together cover 31.8% of valid land and do not
tile the globe; global valid land is the complete domain rather than their sum.

The box bounds were drawn from broad coherent structures in the Fig. 16 RZMC
DA−OL trend field, not from significant-trend tiles. Because they are
nevertheless motivated by Fig. 16, the analysis is an exploratory timing
decomposition of that result rather than an independent test that these regions
are special. Any persistent shift contributes to a positive whole-record trend
regardless of when it occurs, so selection on the trend field can enrich for
differences at any boundary; this applies equally to the P2 and P6 comparisons.

## S3.7 Period means and adjacent-period differences

The paired DA−OL series at every tile is seasonally adjusted using the same
four-step Theil–Sen procedure described in S3.3 for the tile-scale trend fields.
The adjusted tile series $D_i^*(t)$ are then area-weighted within region $R$:

$$
D_R^*(t) = \frac{\sum_{i \in R} A_i D_i^*(t)}{\sum_{i \in R} A_i},
$$

where $A_i$ is M36 tile area. Seasonal adjustment therefore precedes spatial
aggregation, allowing each tile to retain its own calendar-month climatology.
Figures 16 and 17 use the same adjusted paired tile series; only the subsequent
summary differs. Figure 16 estimates a whole-record slope at each tile, whereas
Fig. 17 calculates regional period means.

For region $i$ and period $P_j$, the period mean is the arithmetic mean of the
adjusted monthly series over that period, and the reported quantity at each
boundary is the adjacent-period difference

$$
\Delta_{j} = \overline{\Delta \mathrm{RZMC}}_{P_{j+1}} - \overline{\Delta \mathrm{RZMC}}_{P_{j}} .
$$

This is a difference in average assimilation influence between two periods. It
is **not** an instantaneous discontinuity at the boundary: a gradual change
within either period produces the same quantity as an abrupt one at the
boundary, and the analysis does not distinguish them.

Serial dependence is handled with an AR(1) effective sample size. A single
first-order autocorrelation $\rho$ and residual standard deviation $\sigma$ are
estimated from residuals about the period means of that region. In production
$\rho$ is clipped to the interval $[0, 0.98]$, so negative estimates are treated as
zero and no single estimate can collapse the effective sample size entirely. Each
period mean of $n$ months is assigned the standard error $\sigma/\sqrt{n_{\mathrm{eff}}}$ with

$$
n_{\mathrm{eff}} = \max\left(n\,\frac{1-\rho}{1+\rho},\; 1\right).
$$

The standard error of an adjacent-period difference is the quadrature sum of
the two period-mean standard errors, and 95% intervals are formed as
$\Delta_j \pm 1.96$ standard errors. The effective sample size is floored at one. Estimated $\rho$ ranges from 0.63
to 0.80 across the six regions, so effective sample sizes are roughly a tenth of
the month counts. No P7-related adjacent-period comparison is resolved,
consistent with the short 15-month duration of P7.

Benjamini–Hochberg control at 0.05 is applied separately at each of the eight
P2–P9 adjacent-period boundaries across the six regional tests, matching the boundary-family
convention used elsewhere in this analysis. Nine of the 48 comparisons are
resolved.

## S3.8 Serial-correlation and uncertainty sensitivity

The effective-$n$ interval was compared against a moving-block bootstrap
(24-month blocks) and a fitted-AR(1) parametric bootstrap, both with 2,000
replicates.

| Method | Resolved of 48 | Median SE relative to effective-$n$ |
|---|---:|---:|
| Effective-$n$ AR(1) | 9 | 1.00 |
| Fitted-AR(1) bootstrap | 9 | 0.90 |
| Moving-block bootstrap | 12 | 0.70 |

Effective-$n$ and the parametric bootstrap give the identical significance
classification. The block bootstrap produces smaller standard errors and three
additional marginal detections (global valid land P5−P4; southern Africa and
western North America P9−P8). The headline P2−P1 and P6−P5 findings are
resolved under all three methods. The effective-$n$ interval is retained
because it is the most conservative of the three and agrees exactly with the
fitted-AR(1) bootstrap.

## S3.9 Surface soil-moisture depth sensitivity

The analysis was repeated for surface soil moisture using identical regions,
seasonal adjustment, period definitions, uncertainty calculation and
false-discovery-rate structure; only the variable differs, and both are computed
by the same code path. SFMC is a secondary depth-sensitivity check and is not
pooled into the RZMC families, which remain six regional tests per boundary.

At P2−P1 the same three regions are resolved in both layers, with SFMC
differences 1.14–1.23 times the corresponding RZMC differences. At P6−P5, four
regions are resolved for RZMC and none for SFMC:

| Region | RZMC | SFMC | SFMC/RZMC |
|---|---:|---:|---:|
| Global valid land | +1.52 | +1.01 | 0.67 |
| Australia | +3.46 | +1.96 | 0.57 |
| Southern Africa | +3.57 | +2.46 | 0.69 |
| North Africa / Middle East | +1.62 | +0.14 | 0.09 |

Units ×10⁻³ m³ m⁻³. SFMC residual standard deviation is 1.22–1.59 times that of
RZMC and SFMC standard errors on the P6−P5 difference are 1.05–1.31 times
larger, while the effects themselves are 31–91% smaller. The loss of
significance is therefore driven predominantly by smaller effect sizes rather
than by inflated uncertainty. SFMC resolves 8 of 48 comparisons overall.

---

# Table S1. Statistical estimands

| Scientific question | Primary estimand | Estimator | Dependence handling | Multiple-testing treatment | What the result supports |
|---|---|---|---|---|---|
| Where does assimilation-added snow water go? | Area-weighted fraction of $I_{\mathrm{snow}}$ appearing in each terminal budget term | Ratio of area-weighted aggregate response to area-weighted aggregate input over the snow-addition sample | 5° spatial-block bootstrap, 1,000 replicates (10° sensitivity) | None; a small predeclared set of partition terms | Average mass fate of added water over WY2001–WY2006 |
| Does anomalous snow input predict an anomalous response at the same location? | Controlled marginal response per unit input, $\beta$ | OLS on within-tile anomalies with year fixed effects and OL MAM snow-mass control | 5° spatial-block bootstrap on the full design (10° sensitivity) | None; five closing responses on a common design | Marginal within-location partition; corroboration, not causal proof |
| Does a field drift systematically over 24 years? | Theil–Sen slope of the paired DA−OL monthly series, per year | Exact Theil–Sen median pairwise slope after trend-preserving seasonal adjustment | Hamed–Rao modified Mann–Kendall variance inflation, lags 1–24, factor floored at 1 | BH FDR at 0.05 within each complete field; OL, DA, DA−OL separate families | Presence and sign of a monotonic long-term tendency at each tile |
| When do the regional RZMC differences identified in Fig. 16 emerge? | Difference between adjacent P1–P9 regional mean DA−OL RZMC | The same tile-adjusted paired series used by the trend analysis, area-weighted before calculating period means | AR(1) effective sample size; fitted-AR(1) bootstrap sensitivity | BH FDR separately at each boundary across six regions | A change in average assimilation influence between observing-system periods; not an instantaneous discontinuity |

These quantities carry different units and interpretations and are not compared directly to one another.

---

# Table S2. Synthetic calibration of the trend analysis

| Method | Null / effect scenario | Serial dependence | Replicates / series | False-positive behavior | Recovery / coverage | Interpretive limitation |
|---|---|---|---:|---|---|---|
| Trend (modified MK + BH) | No trend | White noise | 100 series | 3% pointwise, 0% BH | — | Pointwise stippling would overstate significance |
| Trend (modified MK + BH) | No trend | AR(1) = 0.8 | 100 series | 12% pointwise, 0% BH | — | FDR, not pointwise tests, controls the map |
| Trend (modified MK + BH) | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 100% BH detection | Benchmark effect size; not a general power claim |
| Theil–Sen nominal CI | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 43–56% coverage | Nominal interval is badly undercovered |
| Theil–Sen adjusted CI | ±0.02 units yr⁻¹ | AR(1) = 0.8 | 100 series | — | 87–97% coverage | First-order adjustment only; not a test inversion |

Cells marked "—" are not applicable to that scenario. No cell is inferred; every populated value appears in a production validation report. The regional uncertainty comparison in S3.8 is a sensitivity analysis on real series rather than synthetic calibration and is reported there.

---

# Supplementary figure inventory

Final S-numbering is pending until the main Results are complete. Filenames are current tracked review products under `projects/M21C_ls/docs/paper_figures/`, `projects/M21C_ls/docs/trends_breakpoints_report_figures/`, and `projects/M21C_ls/docs/regional_rzmc_diagnostic_figures/`.

| Current filename | Purpose |
|---|---|
| `fig08_supp_ims_terra_aqua_scope_bars.png` | IMS snow-cover categorical skill for the Terra-only versus Terra+Aqua observing-system scopes, as domain-aggregate bars. |
| `fig08_supp_ims_terra_aqua_scope_maps.png` | Spatial distribution of the same Terra-only versus Terra+Aqua IMS skill contrast. |
| `fig_supp_snow_da_octmar_attribution.png` | Non-overlapping October–March snow-DA attribution, showing the sequential control ladder for all six responses in native units. |
| `figS05_precipitation_trends.png` | Precipitation trend control demonstrating that the common forcing produces no paired DA−OL trend. |
| `figS06_sfmc_trends.png` | Surface soil-moisture trend control for OL, DA, and paired DA−OL. |
| `figS07_snow_mass_depth_trends.png` | Snow-mass and snow-depth trend controls over the seasonal-snow domain. |
| `phase2_flux_storage_trend_maps.png` | Whole-record trend maps for ET, total runoff, and total land-water storage (OL, DA, DA−OL). |
| `regional_masks_on_trend_field.png` | The six fixed regional domains drawn on the Fig. 16 RZMC DA−OL trend field. |
| `regional_sfmc_period_means.png` | Regional surface-soil-moisture DA−OL by observing-system period, the depth-sensitivity counterpart to Fig. 17. |

Ordering and numbering will be assigned once the main-text figure set is frozen.
