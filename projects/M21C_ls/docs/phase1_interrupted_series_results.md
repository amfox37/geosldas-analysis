# M21C Phase 1 known-transition results

## Production contract

The P1-P9 interrupted-series production run was completed on 2026-08-12 from
source commit `b6d6dc7`. It expands the 21 Phase 1 matrix rows to 43
area-weighted monthly domain series: OL, DA, and paired `DA - OL` for 11 paired
runs, plus ten DA-only diagnostics. All series cover the complete 288-month
June 2000-May 2024 record.

The output audit passes all 43 models, 1,591 coefficient rows, and 23
boundary-specific FDR families. Thirty estimates are FDR significant: 18
primary estimands and 12 paired controls. The primary set contains 12 level
changes, one slope change, and five period slopes. The control detections are
all period slopes; no OL or DA control has an FDR-significant known-boundary
level or slope change.

## Clearest boundary result: SMAP entry

P6 begins in April 2015, when SMAP brightness-temperature assimilation enters
the observing system. It contains nine FDR-significant primary level changes,
all in soil-water DA impacts or correction diagnostics:

| Primary estimand | P6 level change (95% bootstrap interval) |
|---|---:|
| RZMC DA - OL, valid land | +0.001021 (+0.000488, +0.001558) m3 m-3 |
| RZMC DA - OL, warm snow-free | +0.001267 (+0.000436, +0.002037) m3 m-3 |
| Soil-water absolute activity | +10.612 (+7.056, +14.715) kg m-2 |
| Soil-water signed net correction | +1.530 (+0.589, +2.467) kg m-2 |
| SFMC increment absolute mean | +0.000476 (+0.000311, +0.000638) m3 m-3 |
| SFMC increment RMS | +0.001107 (+0.000677, +0.001505) m3 m-3 |
| RZMC increment signed mean | +0.000007 (+0.000003, +0.000010) m3 m-3 |
| RZMC increment absolute mean | +0.000041 (+0.000027, +0.000055) m3 m-3 |
| RZMC increment RMS | +0.000120 (+0.000077, +0.000162) m3 m-3 |

The valid-land SFMC `DA - OL` level change is positive (+0.001121 m3 m-3,
pointwise p=0.018) but does not survive its 43-series boundary family
(q=0.077). The agreement between valid-land and warm snow-free RZMC, the
increase in soil-water correction activity, and the absence of a significant
P6 change in paired OL/DA controls make P6 the strongest candidate for a real
observing-system impact in this domain-mean analysis. The date remains a known
system transition rather than proof that SMAP alone caused every coincident
change.

## Other transitions

- **P2, MODIS Aqua SCF begins (July 2002):** signed snow-water correction rises
  by 2.594 kg m-2 and absolute snowpack-correction activity rises by 4.413
  kg m-2. Both are FDR significant.
- **P3, ASCAT-A begins (June 2007):** SFMC increment RMS has a +0.000818
  m3 m-3 level change. The valid-land SFMC `DA - OL` slope changes by
  -0.001409 m3 m-3 yr-1; its estimated period slope changes from +0.000437 in
  P2 to -0.000972 m3 m-3 yr-1 in P3.
- **P4, P5, P7, P8, and P9:** no primary level or slope change survives the
  boundary-specific FDR correction. This is failure to resolve a change, not
  evidence that the transition had no effect. P7 remains level-only because
  its 15 months cannot support an independent slope claim.

## Period slopes and controls

Period slopes answer a different question from transition coefficients. The 12
significant control slopes include changes shared by OL and DA: for example,
P4 precipitation declines by about 2.24 kg m-2 month-1 yr-1 in both runs, and
P9 layer-1 soil temperature increases by about 0.349 K yr-1 in both. These are
consistent with common forcing or climate variability and should not be read
as DA impacts. Maps and time-series figures must distinguish period slopes from
level/slope changes at observing-system boundaries.

## Calibration and limits

The fixed-seed validation uses 24 AR(1)=0.8 no-transition series and 24 series
for each planted P6 level and slope effect, with 499 bootstrap replicates. It
passes the predeclared guards: one null BH discovery across 15 transition
families, 91.7% interval coverage for both planted effects, correct mean-effect
directions, and relative mean bias below 50%. FDR power is only 25% for the
planted level effect and 8.3% for the planted slope effect. The conservative
persistence treatment is therefore useful for attribution but deliberately
low-powered, especially for slope hinges.

These are area-weighted domain means, not tile-level transition maps. Dynamic
mask sensitivities retain monthly tile count and represented area, but changing
support can still affect interpretation. The next stage is independent,
autocorrelation-aware changepoint detection and comparison against the known
registry dates; figures should follow only after that validation passes.

## Derived outputs

The following ignored products are under `output/trends_breakpoints/`:

- `phase1_interrupted_series_coefficients.csv`
- `phase1_interrupted_series_models.csv`
- `phase1_interrupted_series_monthly.nc`
- `phase1_interrupted_series_audit.csv`
- `interrupted_series_validation.csv`
