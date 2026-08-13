# M21C Phase 1 trend production results

## Status

The first complete production trend batch was generated on 2026-08-12 from
the audited OL v2 / DA v3 monthly products for June 2000-May 2024. All 21
declared fields pass `scripts/audit_phase1_trend_outputs.py`: 17 primary
estimands and four predeclared mask sensitivities. The NetCDFs and their CSV/log
manifests are derived products under `output/trends_breakpoints/`; their shared
source commit is `ac1cbc1438fb1eaedc5aacfe0ec2445e7fcea414`.

These are whole-record trends, not yet transition-attributed trends. A sequence
of observing-system level changes can appear as a secular slope, so the broad
activity signals below make the planned P1-P9 interrupted-series analysis more
important, not less.

## Primary fields

`Significant` means modified Mann-Kendall followed by global Benjamini-Hochberg
FDR at 0.05 within that complete field. Area percentages use M36 tile area.
Directions refer to the fitted Sen slope; for paired states they describe the
trend in `DA - OL`.

| Field | Supported tiles | Significant tiles | Supported area significant | Positive | Negative | Zero Sen slope |
|---|---:|---:|---:|---:|---:|---:|
| SFMC DA - OL | 112,573 | 1,412 (1.25%) | 1.21% | 993 | 419 | 0 |
| RZMC DA - OL | 112,573 | 7,892 (7.01%) | 7.47% | 7,267 | 625 | 0 |
| SNOMASLAND DA - OL, seasonal snow | 48,067 | 0 | 0% | 0 | 0 | 0 |
| SNODPLAND DA - OL, seasonal snow | 48,067 | 0 | 0% | 0 | 0 | 0 |
| FRLANDSNO DA - OL, seasonal snow | 48,067 | 0 | 0% | 0 | 0 | 0 |
| TSOIL1 DA - OL | 112,573 | 8,143 (7.23%) | 7.71% | 6,610 | 1,533 | 0 |
| PRECTOTCORRLAND DA - OL | 112,573 | 0 | 0% | 0 | 0 | 0 |
| Signed snow correction, seasonal snow | 48,067 | 0 | 0% | 0 | 0 | 0 |
| Absolute snow activity, seasonal snow | 48,067 | 0 | 0% | 0 | 0 | 0 |
| Signed soil-water correction | 112,573 | 1,957 (1.74%) | 1.80% | 486 | 1,471 | 0 |
| Absolute soil-water activity | 112,573 | 106,009 (94.17%) | 97.96% | 99,806 | 0 | 6,203 |
| SFMC correction mean | 112,573 | 12,785 (11.36%) | 12.07% | 1,189 | 11,596 | 0 |
| SFMC correction absolute mean | 112,573 | 101,791 (90.42%) | 92.28% | 97,685 | 4,106 | 0 |
| SFMC correction RMS | 112,573 | 102,004 (90.61%) | 92.66% | 98,682 | 3,322 | 0 |
| RZMC correction mean | 112,573 | 3,063 (2.72%) | 2.79% | 936 | 2,127 | 0 |
| RZMC correction absolute mean | 112,573 | 112,572 (100.00%) | 100.00% | 112,572 | 0 | 0 |
| RZMC correction RMS | 112,573 | 112,547 (99.98%) | 99.99% | 112,547 | 0 | 0 |

## Sensitivities

| Field | Supported tiles | Significant tiles | Supported area significant | Positive | Negative |
|---|---:|---:|---:|---:|---:|
| SFMC DA - OL, warm snow-free months | 104,975 | 879 (0.84%) | 0.87% | 531 | 348 |
| RZMC DA - OL, warm snow-free months | 104,975 | 6,795 (6.47%) | 6.96% | 6,269 | 526 |
| SNOMASLAND DA - OL, locally snowy months | 38,395 | 63 (0.16%) | 0.18% | 7 | 56 |
| SNODPLAND DA - OL, locally snowy months | 38,395 | 156 (0.41%) | 0.37% | 131 | 25 |

The warm/snow-free result shows that most of the RZMC signal remains when
frozen and snow-covered months are removed. The locally-snowy restriction
reveals small snow-mass and snow-depth fields that do not survive FDR in the
all-month seasonal-snow analysis. Their differing dominant directions and very
small spatial fractions require maps and transition attribution before physical
interpretation.

## Main reading

1. Correction magnitude/activity changes far more broadly than signed state
   corrections. The root-zone absolute-mean and RMS activity trends cover
   essentially the entire land area, while signed root-zone corrections are
   significant over about 3% of tiles.
2. The paired state impact changes most clearly for RZMC and TSOIL1. Positive
   slopes dominate both significant fields, meaning DA tends to become wetter
   in the root zone and warmer in layer 1 relative to OL over those tiles.
3. No tile survives global FDR for the primary snow-state, snow-activity, or
   precipitation-control trends. This weakens a simple explanation based on a
   common precipitation trend, but does not test transition-specific effects.
4. The near-global activity trends should not yet be described as gradual
   climate-era trends. The next analysis must estimate level changes at the
   known observing-system transitions and use OL as the control where it exists.

## Support and interval diagnostics

All global/static successful tiles use all 12 calendar months and have complete
valid support. Dynamic masks retain 3-12 calendar months. Seventy-two warm-mask
tiles and 136 locally-snowy tiles pass the initial support gate but fail after
calendar-month adjustment; they remain explicit status-4 outputs.

The first-order autocorrelation-adjusted Sen interval is not the inversion of
the modified-MK test. They agree for the paired SFMC/RZMC and snow sensitivities,
but disagreement is common for tied, zero-heavy activity diagnostics. For
example, 23,454 absolute soil-water-activity tiles are FDR-significant while the
approximate interval contains zero, and 6,203 significant tiles have an exact
zero median Sen slope. The NetCDFs expose `fdr_ci_disagreement`; maps must use
`significant_fdr`. A block-bootstrap or test-inverted interval is a later
sensitivity if interval-based mapped inference is required.

## Verification

- Input audit: 234/234 checks pass.
- Phase 1 output audit: 21/21 fields pass.
- Focused tests: 27 pass.
- All fields embed the same trend configuration and source commit.
- Production execution used eight joblib process workers and completed the
  final 21-field batch in 40.9 minutes.

## Next analysis

Build the known-transition interrupted time-series stage from the same monthly
loader and P1-P9 registry. Estimate level changes at every transition, estimate
slope changes only where segment length permits, treat P7 as level-only, and
carry the Phase 1 support/mask definitions into every comparison.
