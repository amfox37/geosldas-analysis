# Raw Increment Summaries: Snow DA → Soil-Moisture DA Work

Companion notes for:

- `projects/M21C_ls/scripts/build_raw_catch_progn_increment_summaries.py` (+ `_discover.sh`, `submit_raw_catch_progn_increment_summaries.j`)
- `projects/M21C_ls/scripts/build_raw_inst3_fcstana_diagnostic_summaries.py` (+ `_discover.sh`, `submit_raw_inst3_fcstana_diagnostic_summaries.j`)
- output directory (once run): `projects/M21C_ls/output/raw_increment_summaries/`

## Purpose

`docs/Coauthor_feedback.md` records Lauren's review comment: "Add a short
discussion of snow–soil moisture assimilation interactions: does soil moisture
DA behave differently where SCF assimilation is active?" This workstream
builds the data needed to answer that with something more defensible than the
DA-minus-OL state-difference proxies used in the first-pass monthly synthesis
notebook (`monthly_synthesis_snow_sm_diagnostics.md`).

The two notebooks/doc sets are complementary, not duplicates:

- **Monthly synthesis notebook** — uses already-built coarse monthly
  compressed files (`OLv8_land_variables_*`, `LS_monthly_increments_*`).
  Quick to run, available now, but `RZMC_INC` there is a difference of
  monthly means (see caveat below).
- **This raw-increment workstream** — reads native/submonthly GEOSldas output
  directly (`catch_progn_incr`, `inst3_1d_lndfcstana_Nt`) and produces
  cumulative water-budget and burden metrics that are closer to what a
  reviewer would expect for a quantitative snow→SM-DA interaction claim.

## Two kinds of "increment," two different aggregations

The core distinction driving both builder scripts:

| Increment type | Native units | Physically meaningful aggregation | Why |
| --- | --- | --- | --- |
| Prognostic mass increments (`catch_progn_incr`: `WESNN*_INCR`, `SRFEXC_INCR`, `RZEXC_INCR`, `CATDEF_INCR`) | `kg m-2` | **Sum** over time — cumulative water added/removed | These are mass increments; they add linearly into a real water budget. |
| Diagnostic state corrections (`inst3_1d_lndfcstana_Nt`: `*_ANA - *_FCST`) | `m3 m-3`, `K` | **Mean / mean-of-abs / RMS** — never summed | These are volumetric/temperature state corrections, not mass; summing them has no physical meaning. They're activity/burden metrics. |

## Which file answers which question

| Question | File family | Key output variables | Aggregation |
| --- | --- | --- | --- |
| How much net/gross snowpack mass did DA add or remove this month? | raw `catch_progn_incr` | `snow_net`, `snow_abs_netpack`, `snow_abs_layers` | cumulative sum (`kg m-2`) |
| How much net/gross soil water did DA add or remove this month? | raw `catch_progn_incr` | `soil_water_net_approx`, `soil_water_abs_activity` | cumulative sum, **CATDEF sign-flipped** for the net version |
| How much SM-state correction ("burden") did DA apply this month? | raw `inst3_1d_lndfcstana_Nt` | `SFMC_INC_MEAN/_ABS_MEAN/_RMS`, `RZMC_INC_*`, etc. | mean / mean-of-abs / RMS of 3-hourly ANA−FCST diffs |
| Quick-look SM increment, already built, coarser | monthly compressed `LS_monthly_increments_2000_2024.nc` | `SFMC_INC`, `RZMC_INC` | difference of **monthly-mean** ANA and FCST — see caveat below |
| Propagated hydrology response (needs a matched OL run) | monthly OL/DA land + flux-aux files | `SFMC`, `RZMC`, `FRLANDSNO`, `SNOMASLAND`, `EVLAND`, `LHLAND`, etc. | DA−OL seasonal means |

The last row only exists where a matched OL run covers the period. Everything
above it comes from the DA run alone, so it's available even for periods
(e.g. 2000–2007) where no full OL companion run exists.

## CATDEF sign caveat (the one to not get wrong)

`CATDEF` is **catchment deficit**: it runs opposite to water content. A
positive `CATDEF_INCR` means the analysis made the catchment *drier*, while a
positive `RZEXC_INCR`/`SRFEXC_INCR` means *wetter*. You cannot add the three
raw `*_net` fields together — `catdef_net` must have its sign flipped first:

```
soil_water_net_approx ≈ srfexc_net + rzexc_net − catdef_net
```

This is implemented in `build_raw_catch_progn_increment_summaries.py`, and
`catdef_net`/`soil_water_net_approx` both carry a `note` attribute in the
output NetCDF warning about this. For *activity* metrics
(`soil_water_abs_activity = sum(abs(srf) + abs(rz) + abs(catdef))`) the sign
issue doesn't apply, since magnitudes are summed regardless of sign. Treat
`soil_water_net_approx` as approximate even with the sign fixed — the three
Catchment moisture prognostics aren't independent additive reservoirs.

## RZMC_INC: two versions, do not conflate

- `LS_monthly_increments_2000_2024.nc` (`common/python/io/read_monthly_increments.py`):
  built from `inst3_1d_lndfcstana_Nt.monthly.YYYYMM.nc4`, i.e.
  `RZMC_INC = RZMC_ANA - RZMC_FCST` using **monthly-mean** ANA and FCST
  fields. This is a difference-of-means.
- `inst3_fcstana_raw_monthly_diagnostic_*.nc` (this workstream): computes
  `ANA - FCST` at each raw 3-hourly timestep, then averages
  (mean / mean-of-abs / RMS) over the month. This is the
  mean-of-3-hourly-differences, and is the version the "SM DA work" row above
  is meant to use.

These are not interchangeable. The monthly-synthesis notebook's Analysis 3
currently uses the difference-of-means version because it's the only one that
existed when that notebook was built; once the raw diagnostic file is
available, prefer it for any burden statement that needs to go in the paper.

## Status

- Both builder scripts, their `_discover.sh` wrappers, and SLURM submit jobs
  (`submit_raw_catch_progn_increment_summaries.j`,
  `submit_raw_inst3_fcstana_diagnostic_summaries.j`, account `g0610`, 8h
  walltime) are written and committed.
- Neither has been run yet — `projects/M21C_ls/output/raw_increment_summaries/`
  does not exist locally or (as far as we know) on Discover.
- Separately, the monthly OL/DA flux-aux files
  (`flux_core`, `latent_components`, `water_budget`, `energy_context`) have
  already been built and copied down into
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`,
  feeding Analysis 2 of the monthly synthesis notebook. That's unrelated to
  this raw-increment workstream but uses the same Discover-build-then-copy
  pattern.

## Planned first gating test

Once both raw outputs exist locally:

1. Predictor: cumulative MAM `snow_net` (or `snow_abs_netpack`) from
   `catch_progn_raw_monthly_cumulative_*.nc`, restricted to the seasonal-snow
   mask already defined in the monthly synthesis notebook.
2. Response: seasonal mean `RZMC_INC_ABS_MEAN` (AMJ/MJJ/JJA) from
   `inst3_fcstana_raw_monthly_diagnostic_*.nc`, same tiles.
3. One snow region, one or two years first, as a sanity/gating check before
   committing to the full 2000–2024 binned-scatter treatment used in
   Analysis 1/3 of the monthly synthesis notebook.

This directly targets Lauren's review comment and is intended to either
support or rule out a short main-text/supplemental paragraph on snow-SM DA
interaction, not to become a new headline figure.
