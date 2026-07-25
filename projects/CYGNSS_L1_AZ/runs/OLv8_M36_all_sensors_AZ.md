# OLv8_M36_all_sensors_AZ

GEOSldas land data assimilation experiment used as a long-running observation
reader validation and monitoring case.

## Summary

`OLv8_M36_all_sensors_AZ` is an open-loop-style land surface run over a
CONUS-limited box. The original 12-month window begins at `BEG_DATE=20200101`;
the run has now been extended to `END_DATE=20220101`.

The run monitors how 13 satellite soil-moisture or brightness-temperature
observation species compare against the model's simulated land state. It is not
assimilating any of these observations. Each species is configured with:

- `assim=.false.`
- `scale=.false.`
- `getinnov=.true.`

This means observations are read and observation-minus-forecast (`O-F`)
innovation statistics are logged, but nothing feeds back into the model state.
This is a validation/monitoring run, not an active DA run.

## Purpose

The run validates a newly developed set of observation-operator readers over a
long multi-year period rather than only a short test window:

- CYGNSS L1 DDM: a preprocessed scalar soil-moisture proxy derived from
  GNSS-reflectometry data.
- H SAF ASCAT: MetOp-A/B/C H121 CDR backscatter-based soil moisture.

These newer readers are monitored alongside the existing SMOS and SMAP
brightness-temperature readers so their behavior can be compared against
established observation streams.

## Monitored Species

The run monitors 13 observation species:

- SMOS: 4 brightness-temperature species.
- SMAP: 4 brightness-temperature species.
- H SAF ASCAT: 3 species, one each for MetOp-A, MetOp-B, and MetOp-C.
- CYGNSS: L1 DDM scalar and L3 soil moisture.

## Engineering Notes

Two new observation operators were developed and debugged on separate git
feature branches, then merged into `feature/amfox/cygnss-ascat-hsaf-v8` so
CYGNSS L1 and H SAF ASCAT could run together.

Several real bugs were found and fixed through smaller-scale tests before the
long run was trusted:

- CYGNSS reader tile-index resolution bugs.
- CYGNSS coefficient schema bump from `0.4` to `0.5`.
- False-positive grid detection where September coefficient file paths matched
  a 9 km grid pattern and aborted the run.

The run self-resubmits in monthly SLURM segments because the full integration is
too long for a single job. It has been extended twice:

- After the September grid-detection bug stalled the run.
- On 2026-07-23, to push the end date from November 2021 to January 2022 after
  CYGNSS coefficient generation caught up through December 2021.

MetOp-A has a hard documented end-of-record date, 2021-11-15, baked into the H
SAF ASCAT reader. The code has been verified to handle this gracefully with a
clean early return rather than crashing when the run passes that date.

## Current Status

As of 2026-07-23, SLURM job `57297680` has been submitted. It is
self-resubmitting monthly from a November 2021 restart point toward the new
January 2022 end date.
