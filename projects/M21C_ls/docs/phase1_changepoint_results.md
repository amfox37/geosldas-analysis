# M21C Phase 1 independent changepoint results

## Production contract

The independent PELT analysis was completed on 2026-08-12 from source commit
`2c88fc8`. Its 43 input domain series are the unchanged area-weighted monthly
series produced by the known-transition stage at commit `b6d6dc7`. The audit
passes all 43 models, 152 primary-penalty candidate rows, 1,446 complete
method/penalty-grid rows, and 344 P2-P9 comparison rows.

A break is accepted only when it is present in the primary AR(1)+linear method
at 1.25 times BIC, is stable within three months across at least half of the
penalty grid, and agrees within three months with the independently prewhitened
linear method. This yields 37 accepted breaks. All 37 occur in primary
`DA - OL` or DA-only diagnostic series: 20 in paired deltas and 17 in diagnostic
values. No OL or DA control series has an accepted break.

## Agreement with known boundaries

Twenty accepted breaks match a scored P2-P9 boundary within the primary
+/-3-month tolerance. Two more match only under the +/-6-month sensitivity.

| Boundary | Observing-system transition | +/-3 months | Additional +/-6 months |
|---|---|---:|---:|
| P2, 2002-07 | MODIS Aqua SCF begins | 4 | 0 |
| P3, 2007-06 | ASCAT-A assimilation begins | 4 | 0 |
| P4, 2010-05 | SMOS assimilation begins | 2 | 0 |
| P5, 2013-04 | ASCAT-B assimilation begins | 0 | 1 |
| P6, 2015-04 | SMAP assimilation begins | 10 | 0 |
| P7, 2018-08 | CYGNSS assimilation begins | exempt, 0 | 0 |
| P8, 2019-11 | ASCAT-C assimilation begins | 0 | 0 |
| P9, 2021-12 | ASCAT-A assimilation ends | 0 | 1 |

Fifteen accepted breaks remain more than six months from a scored boundary.
They are retained as independently detected structural changes and should not
be attributed to the nearest observing-system event without further evidence.

## Strongest convergence: P6/SMAP

Ten primary series independently break exactly in April 2015:

- RZMC `DA - OL` over valid land and the warm snow-free sensitivity domain;
- signed, absolute-mean, and RMS RZMC correction diagnostics;
- signed, absolute-mean, and RMS SFMC correction diagnostics;
- signed net and absolute soil-water correction activity.

Nine of these ten were also FDR-significant P6 level changes in the known-date
interrupted-series analysis. The additional independently detected series is
SFMC signed increment mean, whose known-date level coefficient was not
significant. Exact date agreement, penalty and method stability, replication
across valid-land and warm snow-free RZMC, and the complete absence of accepted
OL/DA control breaks make P6 the clearest observing-system transition in the
analysis. It remains a controlled attribution result, not proof that every
April 2015 change is caused only by SMAP.

## Other boundary matches

- **P2:** snow signed and absolute correction activity, valid-land SFMC
  `DA - OL`, and locally snowy SNODPLAND `DA - OL` break within one month.
- **P3:** precipitation `DA - OL` and RZMC increment RMS break exactly at the
  boundary; SFMC increment RMS and locally snowy SNOMASLAND `DA - OL` follow by
  one month. The known-date SFMC `DA - OL` slope change is not independently
  localized, consistent with poor slope-hinge power in validation.
- **P4:** seasonal-snow SNOMASLAND `DA - OL` and RZMC increment RMS fall within
  one month, although neither corresponding known-date boundary coefficient
  survived FDR. These are candidates for visual inspection rather than a firm
  transition claim.
- **P5 and P9:** only one series each falls within the wider six-month
  sensitivity; neither has a primary-tolerance match.
- **P7 and P8:** no accepted break matches either boundary. P7 remains exempt
  from scoring because its 15-month period is structurally incompatible with
  the 24-month minimum segment.

## Calibration and limits

The fixed-seed calibration uses 24 series per scenario. White-noise and
AR(1)=0.8 nulls each produce zero accepted breaks. Recovery within six months
is 91.7% for an isolated P6 level shift, 85.4% across two opposing level shifts,
and 91.7% for a level shift near the minimum-segment edge. These pass the
predeclared false-positive and recovery guards.

A gradual P6 slope hinge is recovered within six months in only 4.2% of cases
and has a median localization error of 24.5 months. Independent PELT results
must therefore be interpreted as evidence about strong abrupt structural
changes. Failure to detect a gradual transition is not evidence of no slope
change; known-date slope inference comes from the interrupted-series model.

The results are area-weighted domain means. They do not establish spatial
uniformity or replace tile-level transition maps. Dynamic-mask series preserve
their monthly sample support in the upstream product, but changing support
remains an interpretation sensitivity.

## Derived outputs

The following ignored products are under `output/trends_breakpoints/`:

- `phase1_changepoint_detections.csv`
- `phase1_changepoint_penalty_grid.csv`
- `phase1_changepoint_boundary_comparison.csv`
- `phase1_changepoint_models.csv`
- `phase1_changepoint_monthly.nc`
- `phase1_changepoint_audit.csv`
- `changepoint_validation.csv`
