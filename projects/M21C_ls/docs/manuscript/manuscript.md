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

The direct accounting was complemented by an independent controlled regression analysis. Signed water-year snow input and each terminal budget response were first expressed as within-tile anomalies to remove persistent spatial differences among locations; common year effects were included, and OL spring snow amount was used as a covariate to reduce confounding by interannual variations in the underlying snow climatology. The resulting coefficients provide an independent estimate of the runoff, evapotranspiration, storage, and residual response per unit snow-water input. These regressions are used as corroboration of the direct mass accounting rather than as the primary estimate of water partitioning.

Finally, the timing of the hydrologic response was examined from monthly DA−OL fields across the same October–September water-year cycle. Monthly snow-water input was related to subsequent snowmelt, root-zone soil moisture, runoff, and evapotranspiration for snow-addition tile-years. Snowmelt is interpreted as a transit term marking the release of assimilation-added water from the snowpack rather than as a terminal sink. A stricter sensitivity analysis used only October–March snow input and non-overlapping subsequent response windows; detailed regression variants, alternative bootstrap configurations, the snow-removal partition, storage diagnostics, and soil-moisture persistence analyses are reported in the Supplement.

### 2.7 Long-term trends and observing-system transitions

Long-term changes in the land states were evaluated separately from discrete changes associated with the evolving observing system. The analysis used the 288 monthly values from June 2000 through May 2024. For variables available in both experiments, OL and DA were first restricted to identical finite monthly samples, and the primary measure of assimilation impact was the trend of the paired DA−OL series rather than the difference between trends estimated independently for OL and DA. Trends were evaluated for precipitation, surface and root-zone soil moisture, snow-cover fraction, snow mass, and snow depth. Soil-moisture and precipitation analyses used the valid-land domain, whereas snow-state analyses used the Northern Hemisphere seasonal-snow domain. All spatial aggregates were weighted by M36 tile area.

At each tile, the calendar-month climatology was removed using a trend-preserving seasonal adjustment, and long-term change was estimated using the exact Theil–Sen median slope (Helsel et al., 2020). Trend significance was assessed with a Mann–Kendall test incorporating the autocorrelation variance correction of Hamed and Rao (1998), followed by Benjamini–Hochberg false-discovery-rate control at 0.05 across each complete spatial field (Wilks, 2016, 2019). OL, DA, and paired DA−OL trends were evaluated on the same underlying samples but retained separate spatial FDR families. Area-weighted domain-mean trends were evaluated separately from the tile-scale fields. Significance shown in the spatial maps is based on the FDR-controlled test rather than on pointwise confidence intervals.

Potential observing-system discontinuities were then assessed using the P1–P9 chronology defined in Fig. 1. This analysis used area-weighted monthly domain-mean series and was deliberately separated from the whole-record trend analysis because a discrete level change or change in slope need not imply a secular trend over the complete record. Known observing-system boundaries were evaluated using a segmented regression containing an intercept, a baseline linear trend, calendar-month effects, level changes at each P2–P9 boundary, and slope-change terms where the duration of the adjacent periods permitted independent estimation. Serial dependence was represented using iterative Prais–Winsten AR(1) fitting, and uncertainty in the transition coefficients was estimated by generating realizations from the fitted AR(1) error process and refitting the complete model (Wilks, 2019). False-discovery-rate control was applied separately within each observing-system-boundary family. Because P7 spans only 15 months, it was assigned a level-change term but no independent slope change; the P6 slope was retained until the P8 boundary.

As an independent check on the prescribed observing-system dates, abrupt changepoints were identified without supplying the P1–P9 boundaries to the detector. Two formulations of the Pruned Exact Linear Time algorithm (PELT; Killick et al., 2012) were applied to trend-preserving, seasonally adjusted monthly series across a predeclared range of penalties. A changepoint was accepted only if it was retained by the primary piecewise AR(1)-plus-linear formulation, remained stable across the penalty range, and was reproduced within three months by a separately prewhitened linear formulation. Accepted changepoints were subsequently compared with the known P2–P9 dates using ±3 months as the primary matching tolerance and ±6 months as a sensitivity. Unmatched changepoints were retained rather than assigned to the nearest observing-system event. This conservative procedure is most sensitive to abrupt structural changes and has limited ability to localize gradual slope transitions; consequently, the known-date segmented analysis remains the primary test of slope changes.

Known-date and independently detected changes were interpreted together. Agreement between an independently detected discontinuity and an observing-system boundary strengthens the evidence for a change in analysis regime, but temporal correspondence alone is not considered evidence that the newly introduced observing system uniquely caused every contemporaneous change.

## 3. Results

### 3.1 Evolution of the assimilated observing system

### 3.2 Evolution of observational influence

### 3.3 Independent soil-moisture evaluation

### 3.4 Snow-cover, SWE, and snow-depth evaluation

### 3.5 Comparison with ERA5-Land

### 3.6 Hydrologic consequences of snow-cover assimilation

### 3.7 Long-term trends and observing-system transitions

## 4. Discussion

### 4.1 Evolving observational constraint in long-term land analysis

### 4.2 From direct observational constraint to model-mediated hydrologic response

### 4.3 Observing-system transitions and long-term trend fidelity

### 4.4 Limitations and implications for future multi-sensor reanalysis

## 5. Conclusions

## References

Hamed, K. H., and A. R. Rao, 1998: A modified Mann–Kendall trend test for autocorrelated data. *Journal of Hydrology*, **204**, 182–196. https://doi.org/10.1016/S0022-1694(97)00125-X.

Helsel, D. R., R. M. Hirsch, K. R. Ryberg, S. A. Archfield, and E. J. Gilroy, 2020: *Statistical Methods in Water Resources*. U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 pp. https://doi.org/10.3133/tm4A3.

Killick, R., P. Fearnhead, and I. A. Eckley, 2012: Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, **107**, 1590–1598. https://doi.org/10.1080/01621459.2012.737745.

Wilks, D. S., 2016: “The stippling shows statistically significant grid points”: How research results are routinely overstated and overinterpreted, and what to do about it. *Bulletin of the American Meteorological Society*, **97**, 2263–2273. https://doi.org/10.1175/BAMS-D-15-00267.1.

Wilks, D. S., 2019: *Statistical Methods in the Atmospheric Sciences*. 4th ed. Elsevier, 840 pp. ISBN 978-0-12-815823-4.
