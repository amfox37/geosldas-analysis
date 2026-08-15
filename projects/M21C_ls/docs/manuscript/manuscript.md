# Evolving Observational Influence and Temporal Consistency in a Multi-Sensor Land Reanalysis, 2000–2024

## Abstract

## 1. Introduction

## 2. Methods and Data

### 2.1 GEOS LDAS and Catchment model

### 2.2 Assimilated observing system

### 2.3 Land data assimilation

### 2.4 Experiments and observing-system periods

### 2.5 Evaluation datasets and metrics

### 2.6 Snow-DA hydrologic pathway and differential water-budget analysis

The hydrologic consequences of snow-cover assimilation were examined during the six complete MODIS-only water years, WY2001–WY2006, where each water year extends from October through September. This period predates the introduction of microwave soil-moisture assimilation in June 2007 and therefore provides the cleanest interval in which to attribute differences between DA and OL to snow-cover assimilation. The analysis was restricted to the Northern Hemisphere seasonal-snow domain used throughout the snow-process analysis. All quantities were evaluated in native M36 tile space, and fluxes were converted to monthly water-equivalent totals before aggregation over the water year.

The primary analysis was a differential water budget constructed from the signed snow-water analysis increments and the corresponding DA−OL hydrologic responses over a common October–September window. The snow-DA input, $I_{\mathrm{snow}}$, was calculated as the water-year sum of the native prognostic `snow_net` increments, with positive values denoting addition of snow water and negative values denoting removal. The closing terms were the DA−OL changes in total runoff, evapotranspiration, and total land-water storage,

\[
I_{\mathrm{snow}} =
\Delta R_{\mathrm{surf}} +
\Delta R_{\mathrm{base}} +
\Delta ET +
\Delta S +
\epsilon,
\]

where total runoff is the sum of surface runoff and baseflow and $\epsilon$ is the residual. Because OL and DA use the same precipitation forcing, precipitation does not enter the differential budget; a separate numerical check was used to verify that the annual DA−OL precipitation difference remained negligible. Snowmelt and infiltration were retained as pathway diagnostics but were not included as terminal budget terms because they redistribute water within the land column rather than remove it from the budget.

The storage term required particular care because the standard land-water tendency diagnostic does not include the discontinuous mass introduced by the analysis update. We therefore calculated $\Delta S$ directly from the change in DA−OL total land-water storage (`TWLAND`) between the September monthly mean preceding each water year and the September monthly mean at its end. The integrated DA−OL `WCHANGELAND` tendency was retained only as a diagnostic and was not used to close the budget. Surface and root-zone soil moisture were likewise treated as diagnostic state responses rather than independent storage terms because their water is already contained within total land-water storage. The use of monthly-mean storage endpoints introduces a small temporal approximation into the budget, which is considered when interpreting the residual.

Water-budget partitions were evaluated both for the full sample and separately for tile-years in which snow assimilation added ($I_{\mathrm{snow}}>0$) or removed ($I_{\mathrm{snow}}<0$) water. The main analysis emphasizes the snow-addition sample. Area-weighted partition fractions were calculated for surface runoff, baseflow, evapotranspiration, storage change, and the residual. Uncertainty in these fractions was estimated with 1,000 spatial-block bootstrap realizations using approximately $5^\circ \times 5^\circ$ blocks, thereby accounting for spatial dependence among neighboring tiles. Annual domain budgets were also retained to assess the consistency of the partitioning across the six individual water years. Figure 14 uses the same accepted positive-input sample and spatial-block uncertainty calculation throughout.

The direct accounting was complemented by an independent controlled regression analysis. Signed water-year snow input and each terminal budget response were first expressed as within-tile anomalies to remove persistent spatial differences among locations; common year effects were included, and the OL water-year peak snow mass was used as a covariate to reduce confounding by interannual variations in the underlying snow climatology. The resulting coefficients provide an independent estimate of the runoff, evapotranspiration, storage, and residual response per unit snow-water input. These regressions are used as corroboration of the direct mass accounting rather than as the primary estimate of water partitioning.

Finally, the timing of the hydrologic response was examined from monthly DA−OL fields across the same October–September water-year cycle. Monthly snow-water input was related to subsequent snowmelt, root-zone soil moisture, runoff, and evapotranspiration for snow-addition tile-years. Snowmelt is interpreted as a transit term marking the release of assimilation-added water from the snowpack rather than as a terminal sink. A stricter sensitivity analysis used only October–March snow input and non-overlapping subsequent response windows; detailed regression variants, alternative bootstrap configurations, the snow-removal partition, storage diagnostics, and soil-moisture persistence analyses are reported in the Supplement.

### 2.7 Long-term trends and observing-system transitions

Long-term changes in the land states were evaluated separately from discrete changes associated with the evolving observing system. The analysis used the 288 monthly values from June 2000 through May 2024. For variables available in both experiments, OL and DA were first restricted to identical finite monthly samples, and the primary measure of assimilation impact was the trend of the paired DA−OL series rather than the difference between trends estimated independently for OL and DA. Trends were evaluated for precipitation, surface and root-zone soil moisture, snow-cover fraction, snow mass, and snow depth. Soil-moisture and precipitation analyses used the valid-land domain, whereas snow-state analyses used the Northern Hemisphere seasonal-snow domain. All spatial aggregates were weighted by M36 tile area.

At each tile, the calendar-month climatology was removed using a trend-preserving seasonal adjustment, and long-term change was estimated using the exact Theil–Sen median slope (Helsel et al., 2020). Trend significance was assessed with a Mann–Kendall test incorporating a conservative adaptation of the autocorrelation variance correction of Hamed and Rao (1998), followed by Benjamini–Hochberg false-discovery-rate control at 0.05 across each complete spatial field (Benjamini and Hochberg, 1995; Wilks, 2016). OL, DA, and paired DA−OL trends were evaluated on the same underlying samples but retained separate spatial FDR families. Area-weighted domain-mean trends were evaluated separately from the tile-scale fields. Significance shown in the spatial maps is based on the FDR-controlled test rather than on pointwise confidence intervals.

Potential observing-system discontinuities were then assessed using the P1–P9 chronology defined in Fig. 1. This analysis used area-weighted monthly domain-mean series and was deliberately separated from the whole-record trend analysis because a discrete level change or change in slope need not imply a secular trend over the complete record. Known observing-system boundaries were evaluated using a segmented regression containing an intercept, a baseline linear trend, calendar-month effects, level changes at each P2–P9 boundary, and slope-change terms where the duration of the adjacent periods permitted independent estimation. Serial dependence was represented using iterative Prais–Winsten AR(1) fitting (Prais and Winsten, 1954), and the transition-coefficient uncertainty reported here was estimated using an innovation-resampling bootstrap based on the fitted AR(1) model, with a conservative estimate of residual persistence used when generating the bootstrap realizations (Wilks, 2019). False-discovery-rate control was applied separately within each observing-system-boundary family. Because P7 spans only 15 months, it was assigned a level-change term but no independent slope change; the P6 slope was retained until the P8 boundary.

As an independent check on the prescribed observing-system dates, abrupt changepoints were identified without supplying the P1–P9 boundaries to the detector. Two formulations of the Pruned Exact Linear Time algorithm (PELT; Killick et al., 2012), implemented in the ruptures package (Truong et al., 2020), were applied to trend-preserving, seasonally adjusted monthly series across a predeclared range of penalties. A changepoint was accepted only if it was retained by the primary piecewise AR(1)-plus-linear formulation, remained stable across the penalty range, and was reproduced within three months by a separately prewhitened linear formulation. Accepted changepoints were subsequently compared with the known P2–P9 dates using ±3 months as the primary matching tolerance and ±6 months as a sensitivity. Unmatched changepoints were retained rather than assigned to the nearest observing-system event. This conservative procedure is most sensitive to abrupt structural changes and has limited ability to localize gradual slope transitions; consequently, the known-date segmented analysis remains the primary test of slope changes.

Known-date and independently detected changes were interpreted together. Agreement between an independently detected discontinuity and an observing-system boundary strengthens the evidence for a change in analysis regime, but temporal correspondence alone is not considered evidence that the newly introduced observing system uniquely caused every contemporaneous change.

## 3. Results

### 3.1 Evolution of the assimilated observing system

### 3.2 Evolution of observational influence

### 3.3 Independent soil-moisture evaluation

### 3.4 Snow-cover, SWE, and snow-depth evaluation

### 3.5 Comparison with ERA5-Land

### 3.6 Hydrologic consequences of snow-cover assimilation

During the six complete MODIS-only water years, snow-cover assimilation introduced substantial snow-water corrections that subsequently propagated through the modeled hydrologic cycle. Across WY2001–WY2006, the area-weighted mean signed snow-DA input was 58.2 kg m⁻² yr⁻¹ (Fig. 14a). The corresponding runoff response varied from 55.6% to 70.1% of the annual net input, indicating that runoff was consistently a major fate of the assimilation-induced snow water rather than the result being dominated by a single year. The six-year mean budget had a small negative residual of −2.67 kg m⁻² yr⁻¹, equivalent to −4.6% of the net snow-water input.

The partition is clearer when the analysis is restricted to the 247,545 tile-years in which snow assimilation added water (Fig. 14b). Of the added snow water, 43.1% subsequently left as surface runoff and 21.2% as baseflow, giving a total runoff fraction of 64.3%. The 5° spatial-block bootstrap interval for total runoff was 61.1–67.2%. Evapotranspiration accounted for a further 35.9%, whereas only 3.9% remained as additional land-water storage at the end of the water year. The residual was −4.1%. Thus, most snow water added by MODIS SCF assimilation subsequently left the land column as runoff, with a substantial additional loss through evapotranspiration and relatively little remaining in storage by the end of the accounting period.

![Figure 14](../paper_figures/fig14_snow_da_water_budget.png)

**Figure 14.** Differential water-budget response to MODIS snow-cover assimilation during the six complete MODIS-only water years, WY2001–WY2006. (a) Area-weighted annual accounting over the Northern Hemisphere seasonal-snow domain. Diamonds show the signed snow-water input from the analysis increments; bars show the corresponding DA−OL changes in total runoff, evapotranspiration (ET), September-to-September total land-water storage, and the residual. (b) Fate of positive snow-DA input for the 247,545 snow-addition tile-years. Error bars denote 95% intervals from the 5° spatial-block bootstrap. Surface runoff and baseflow together account for 64.3% of the added snow water, with a 61.1–67.2% bootstrap interval.

An independent controlled water-year analysis produced the same runoff-dominated interpretation. After removal of persistent within-tile differences and accounting for common year effects and OL snow amount, the estimated response per unit snow-water input was 0.751 for runoff (95% interval 0.715–0.787), 0.181 for ET (0.153–0.212), and 0.084 for end-of-year storage (0.074–0.096), with a small residual of −0.016 (−0.022 to −0.011). The controlled estimates are not numerically identical to the direct accounting, but independently support the same partitioning of assimilation-added snow water.

The monthly evolution provides the corresponding process pathway (Fig. 15). Positive snow-water increments accumulated during the snow season and were followed by increased DA−OL snowmelt and a sustained positive soil-moisture response into spring. Across snow-addition tile-years, the mean within-year peak RZMC difference was 0.0189 m³ m⁻³, with May the most common month of the peak; the mean May–July RZMC difference was 0.0108 m³ m⁻³. Increased runoff and ET accompanied the spring release of the added snow water. A stricter analysis using only October–March snow input and subsequent non-overlapping response windows retained positive responses in RZMC, ET, runoff, and total land-water storage. Delayed infiltration was the exception and was not distinguishable from zero under this stricter formulation. These results support the sequence of snow-water addition, subsequent melt, soil-water recharge, and eventual partitioning into runoff, ET, and storage, while treating snowmelt itself as a transit term rather than a terminal component of the water budget.

![Figure 15](../paper_figures/fig15_snow_da_monthly_pathway.png)

**Figure 15.** Monthly hydrologic response to snow-water addition during WY2001–WY2006 for the same 247,545 snow-addition tile-years used in Fig. 14b. Area-weighted October–September means are shown for (a) signed snow-DA input, (b) DA−OL snowmelt, (c) DA−OL surface (SFMC) and root-zone (RZMC) soil moisture, and (d) DA−OL total runoff and evapotranspiration. The mean tile-year peak RZMC response is 0.0189 m³ m⁻³, and May is the most common month of peak RZMC response. The analysis period ends in September 2006, before microwave soil-moisture assimilation begins.

The modeled energy response was consistent with this hydrologic propagation. Under warm, snow-free conditions, wetter DA−OL root-zone states were associated with increased ET and latent heat flux, reduced sensible heat flux, and increased evaporative fraction. These relationships provide model-internal evidence that assimilation-induced soil-water changes propagate into the surface-energy balance, but do not constitute independent validation of the simulated fluxes.

After microwave soil-moisture assimilation began, there was no robust evidence that locations experiencing stronger snow DA subsequently required larger soil-moisture analysis corrections. The signed corrections nevertheless showed a systematic tendency: where preceding snow assimilation had added water, later soil-moisture corrections were on average more negative. This population-level response suggests partial opposition to snow-induced wetting rather than a simple increase in subsequent assimilation activity; the tile-year relationship was weak and did not eliminate the snow-to-hydrology response described above.

### 3.7 Long-term trends and observing-system transitions

Despite the substantial evolution in observational constraint, DA generally preserved the broad long-term behavior of the OL simulation. The common precipitation forcing provides a useful null-behavior check: 3,719 OL and 3,726 DA tiles had significant precipitation trends, with 3,603 significant in both experiments and with the same sign, and the complete OL–DA slope correlation was 0.9998. No paired DA−OL precipitation trend survived field-significance control. Similarly, paired DA−OL trends were essentially absent for snow mass and snow depth. Surface soil moisture showed somewhat greater local sensitivity to DA, but significant paired trends were limited to 1.25% of the valid-land domain.

The principal long-term state differences are summarized in Fig. 16. For RZMC, the area with significant trends increased from 4.36% in OL to 9.18% in DA. Significant trends in the paired DA−OL series covered 7.01% of valid land; of these, 7,267 tiles had positive trends and 625 had negative trends. The positive changes are particularly evident across Australia, North Africa and the Middle East, southern Africa, and parts of the Americas and northern Eurasia. Nevertheless, the area-weighted OL, DA, and DA−OL RZMC trends were all nonsignificant. DA therefore modifies the regional evolution of root-zone soil moisture without producing a resolved global land-mean wetting trend.

The snow-cover result is markedly different. Area-weighted SCF declined significantly over the seasonal-snow domain in both experiments, at −0.000554 yr⁻¹ in OL and −0.000549 yr⁻¹ in DA, corresponding to a decline of about 0.013 over the 24-yr record. The paired DA−OL trend was −0.000006 yr⁻¹ and was not significant. Thus, although snow-cover assimilation strongly affects the instantaneous analyzed snow state, it leaves the broad long-term SCF decline essentially unchanged.

![Figure 16](../paper_figures/fig16_longterm_rzmc_scf_trends.png)

**Figure 16.** Long-term June 2000–May 2024 trends in (a–c) root-zone soil moisture (RZMC) and (d–f) snow-cover fraction (SCF) for the open-loop (OL), data-assimilation (DA), and paired DA−OL series. Trends are exact Theil–Sen slopes after trend-preserving removal of the calendar-month climatology. Black stippling denotes trends significant after Benjamini–Hochberg false-discovery-rate control at 0.05. RZMC uses the valid-land domain and SCF the Northern Hemisphere seasonal-snow domain. The DA−OL panels show the trend of the paired difference series rather than the difference between independently estimated OL and DA trends. The RZMC DA−OL panel uses a separately optimized color scale, whereas the SCF DA−OL panel retains the OL/DA scale to avoid visually amplifying the near-null assimilation-induced SCF trend.

The absence of widespread DA-induced long-term trends does not imply that the influence of the observing system is stationary. The monthly soil-water diagnostics show a pronounced change at the P5–P6 boundary in April 2015, coincident with the introduction of SMAP brightness-temperature assimilation (Fig. 17a). The interrupted-series estimate gives an immediate positive RZMC DA−OL level change of 0.00102 m³ m⁻³ over valid land (95% bootstrap interval 0.00049–0.00156 m³ m⁻³) and 0.00127 m³ m⁻³ under warm, snow-free conditions (Fig. 17b). The transition is accompanied by increased surface and root-zone analysis-correction RMS (Fig. 17c) and a 10.61 kg m⁻² increase in absolute prognostic soil-water correction activity (Fig. 17d). Signed net soil-water correction also shifts positively by 1.53 kg m⁻². All six P6 changes shown in Fig. 17b–d survive the boundary-specific false-discovery-rate test.

![Figure 17](../paper_figures/fig17_observing_system_transitions.png)

**Figure 17.** Changes in soil-water data-assimilation behavior across the P1–P9 observing-system periods. (a) Standardized area-weighted monthly RZMC DA−OL and soil-water analysis-correction diagnostics during June 2000–May 2024. Background shading denotes the P1–P9 periods defined in Fig. 1; the P5–P6 boundary in April 2015 marks the introduction of SMAP brightness-temperature assimilation. The seasonally adjusted series are standardized using their full-record mean and sample standard deviation for display only; statistical inference is performed in native units. (b) Estimated P6 level changes in RZMC DA−OL for valid-land and warm, snow-free conditions. (c) Corresponding changes in SFMC and RZMC analysis-correction RMS. (d) Changes in absolute and signed net prognostic soil-water correction activity. Symbols show the estimated level changes and horizontal bars the 95% fitted-AR(1) bootstrap intervals. Statistical significance uses boundary-family false-discovery-rate control at 0.05.

Independent changepoint detection strongly corroborates the April 2015 transition. Ten primary DA-impact or correction diagnostics have accepted changepoints exactly in April 2015, and nine of these also have significant level changes in the prescribed-date analysis. In contrast, no paired OL or DA state-control series has an accepted changepoint at that date. Across the complete analysis, 37 independent changepoints were accepted, of which 20 matched a known observing-system boundary within ±3 months and two additional breaks matched within ±6 months; 15 remained unmatched and were not assigned to an observing-system change. The convergence of the known-date and independently detected analyses identifies April 2015 as the clearest change in the soil-water correction regime. Its coincidence with the start of SMAP assimilation is consistent with a SMAP-driven change in analysis behavior, although the timing alone does not establish SMAP as the unique cause of every contemporaneous change.

Taken together, the trend and transition analyses distinguish two aspects of temporal consistency. The evolving observing system demonstrably changes the magnitude and character of the analysis corrections, most clearly after April 2015, while the long-term analyzed state remains much closer to the underlying OL trajectory: snow trends are largely preserved, and the principal long-term DA effect is a regional modification of root-zone soil-moisture trends rather than a widespread secular drift.

## 4. Discussion

### 4.1 Evolving observational constraint in long-term land analysis

### 4.2 From direct observational constraint to model-mediated hydrologic response

### 4.3 Observing-system transitions and long-term trend fidelity

### 4.4 Limitations and implications for future multi-sensor reanalysis

## 5. Conclusions

## References

Benjamini, Y., and Y. Hochberg, 1995: Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society, Series B (Methodological)*, **57**, 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x.

Hamed, K. H., and A. R. Rao, 1998: A modified Mann–Kendall trend test for autocorrelated data. *Journal of Hydrology*, **204**, 182–196. https://doi.org/10.1016/S0022-1694(97)00125-X.

Helsel, D. R., R. M. Hirsch, K. R. Ryberg, S. A. Archfield, and E. J. Gilroy, 2020: *Statistical Methods in Water Resources*. U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 pp. https://doi.org/10.3133/tm4A3.

Killick, R., P. Fearnhead, and I. A. Eckley, 2012: Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, **107**, 1590–1598. https://doi.org/10.1080/01621459.2012.737745.

Prais, S. J., and C. B. Winsten, 1954: Trend estimators and serial correlation. Cowles Commission Discussion Paper, Statistics No. 383, University of Chicago, Chicago, 26 pp. [Available online at https://cowles.yale.edu/sites/default/files/files/pub/cdp/s-0383.pdf.]

Truong, C., L. Oudre, and N. Vayatis, 2020: Selective review of offline change point detection methods. *Signal Processing*, **167**, 107299. https://doi.org/10.1016/j.sigpro.2019.107299.

Wilks, D. S., 2016: “The stippling shows statistically significant grid points”: How research results are routinely overstated and overinterpreted, and what to do about it. *Bulletin of the American Meteorological Society*, **97**, 2263–2273. https://doi.org/10.1175/BAMS-D-15-00267.1.

Wilks, D. S., 2019: *Statistical Methods in the Atmospheric Sciences*. 4th ed. Elsevier, 840 pp. ISBN 978-0-12-815823-4.
