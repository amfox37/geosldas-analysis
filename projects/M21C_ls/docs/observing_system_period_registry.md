# M21C observing-system period registry

The machine-readable source of truth is
`projects/M21C_ls/config/observing_system_registry.json`. It contains the
P1-P9 fine periods, V1-V3 validation periods, and sensor availability used by
the unified paper figures and trend/breakpoint work.

The registry is loaded and validated by
`projects/M21C_ls/scripts/m21c_periods.py`. Do not reproduce the dates in a
notebook or plotting script.

## Fine periods

| ID | Inclusive dates | Months | Observing-system change at start | Slope use |
| --- | --- | ---: | --- | --- |
| P1 | 2000-06 to 2002-06 | 25 | Start of analysis record; MODIS Terra SCF | Cautious |
| P2 | 2002-07 to 2007-05 | 59 | MODIS Aqua SCF begins | Yes |
| P3 | 2007-06 to 2010-04 | 35 | ASCAT-A SSM begins | Yes |
| P4 | 2010-05 to 2013-03 | 35 | SMOS Tb begins | Yes |
| P5 | 2013-04 to 2015-03 | 24 | ASCAT-B SSM begins | At 24-month floor |
| P6 | 2015-04 to 2018-07 | 40 | SMAP Tb begins | Yes |
| P7 | 2018-08 to 2019-10 | 15 | CYGNSS begins | No |
| P8 | 2019-11 to 2021-11 | 25 | ASCAT-C SSM begins | Cautious |
| P9 | 2021-12 to 2024-05 | 30 | ASCAT-A SSM ends | Yes |

The periods are contiguous, have no overlaps, and cover all 288 months from
June 2000 through May 2024.

## Segment constraints

The initial trend and changepoint design uses a 24-month minimum segment.
P7 is only 15 months long, so it must not receive a period-specific slope
claim. Interrupted time-series treatment of the P6-P7 and P7-P8 transitions is
limited to level changes unless P7 is combined with an adjacent period for a
specific sensitivity analysis.

A changepoint method with a 24-month minimum cannot independently recover both
boundaries of P7. P7 is therefore exempt from detection-agreement scoring; its
absence from an independently detected breakpoint list is not counted as a
miss. P5 lies exactly on the floor, while P1 and P8 are only one month above
it, so their slope estimates require explicit stability checks.

## Broader validation periods

- V1, SCF-only period: P1-P2, 84 months.
- V2, microwave pre-SMAP period: P3-P5, 94 months.
- V3, SMAP-era microwave period: P6-P9, 110 months.

These broader periods remain suitable for validation summaries and for trend
sensitivity analyses requiring longer segments.
