# M21C trends and breakpoints: analysis contract

This document defines the input and provenance gate for the replacement of the
legacy 2025 trend exploration. No new trend result should be interpreted until
the audit described below passes.

## Existing monthly pipeline is authoritative

The trend layer consumes the products already validated by
`notebooks/monthly_synthesis_snow_sm_diagnostics.ipynb`. It does not rebuild
monthly fields from raw GEOS LDAS output and does not introduce a parallel
definition of DA increments.

The existing monthly-synthesis records are:

- `output/monthly_synthesis_diagnostics/monthly_synthesis_time_coverage.csv`
- `output/monthly_synthesis_diagnostics/monthly_synthesis_input_inventory.csv`
- `output/monthly_synthesis_diagnostics/monthly_synthesis_variable_registry.csv`
- `output/monthly_synthesis_diagnostics/analysisD_domain_monthly_timeseries.csv`

The first three files define and audit twelve aligned 288-month by 112,573-tile
datasets. The Analysis D table supplies existing area-weighted regional series
where its domain and variable coverage match a trend question. Tile-level work
reads the same validated monthly NetCDF products listed by the inventory.

`config/trend_breakpoint_inputs.json` is the machine-readable file contract for
those products. It records filenames, experiment/family ownership, source-root
tokens, expected dimensions and dates, and the tile-coordinate file. Both the
audit and series loader read this contract so they cannot silently diverge.

## Increment contract

- `catch_raw_cumulative` contains raw prognostic water/mass increments summed
  over each month in `kg m-2`.
- `inst3_raw_diagnostic` contains signed mean, absolute mean, and RMS
  ANA-minus-FCST state corrections. These are diagnostic state quantities and
  are never cumulatively summed.
- `CATDEF` is a deficit. The existing `soil_water_net_approx` product applies
  the documented sign reversal before combining it with `SRFEXC` and `RZEXC`.
- The legacy `LS_monthly_increments_2000_2024.nc` file is not an input to the
  new trend workflow.

See `docs/raw_increment_summaries_methods.md` for the full derivation.

## Variable selection

`config/trend_breakpoint_variable_selection.json` selects variables from the
existing monthly datasets without redefining their units or aggregation. Phase
1 covers soil moisture, snow, precipitation context, and raw-derived increment
activity. Phase 2 adds water-budget and energy-partitioning responses after the
statistical implementation passes synthetic and real-data validation.

`config/phase1_trend_runs.json` is the production matrix. It contains exactly
one primary run for each of the 17 Phase 1 variables plus four predeclared
sensitivities. Soil-moisture and general context fields use `valid_land` as the
primary mask. Snow states and snow-water corrections use `seasonal_snow`.
`SFMC` and `RZMC` additionally use the month-specific
`warm_snowfree_monthly` mask; `SNOMASLAND` and `SNODPLAND` additionally use
`locally_snowy_monthly`. The first production batch estimates the registered
primary impact series (`delta` for paired OL/DA fields and `value` for DA-only
diagnostics). Separate OL and DA trend fields are controls to generate when a
primary impact result needs attribution; they are not silently included in the
Phase 1 FDR families.

For paired OL/DA quantities, the primary assimilation-impact estimand is the
trend of the paired `DA - OL` monthly series. The workflow will not define the
impact by subtracting two Theil-Sen slopes.

## Shared monthly-series loader

`scripts/trend_breakpoint_series.py` supplies the only new analysis-facing
loader. It validates each opened dataset against the input contract and returns:

- strictly paired OL, DA, and `DA - OL` tile fields for model variables;
- unchanged monthly values for DA-only cumulative prognostic increments and
  ANA-minus-FCST diagnostic statistics;
- monthly water totals only when a source field is explicitly a monthly mean
  rate in `kg m-2 s-1`;
- area-weighted domain means together with the contributing tile count and
  represented area for every month.

The tile weights are the M36 land-tile areas from the GEOS LDAS tile-coordinate
file, in `km2`. Missing OL or DA values are cross-masked before any paired mean
or difference is calculated. There is no unweighted fallback.
The production contract also requires the 112,573 tile areas to total between
100 and 160 million `km2`; the observed total is about 130.4 million `km2`.
This broad physical check prevents an `m2`/`km2` label error without constraining
the weighting numerically.

The loader also centralizes the masks already used by the monthly-synthesis
notebook: valid land, NH snow possible, NH seasonal snow, warm static,
month-specific warm snow-free, and month-specific locally snowy. Static mask
thresholds are unchanged. The dynamic warm mask requires both OL and DA snow
fraction below 0.05 and layer-1 soil temperature above 277.15 K. The dynamic
snow mask restricts the static NH seasonal-snow domain to months where either
run has snow fraction above 0.05 or snow mass above 5 kg m-2.

## Period and changepoint constraints

The only observing-system date source is
`config/observing_system_registry.json`. The shared loader validates date
coverage, contiguity, month counts, V1-V3 membership, and segment-reliability
flags.

P7 is 15 months long. It is unsuitable for a period-specific slope and exempt
from changepoint detection-agreement scoring under the planned 24-month minimum
segment. P5 is exactly at the floor; P1 and P8 are only one month above it.

Before applying PELT, monthly seasonality and autocorrelation must be handled
explicitly. The planned comparison is an autoregressive changepoint cost versus
PELT on appropriately prewhitened residuals, with block-bootstrap stability and
synthetic AR-series false-positive tests. Climatology-period and seasonal-dummy
formulations are required sensitivity checks.

## Known-transition interrupted series

`config/interrupted_time_series.json` and
`scripts/interrupted_time_series.py` define the first transition-attribution
stage. The response is an area-weighted monthly domain series. Every paired
matrix row is fitted three ways on identical support: OL, DA, and `DA - OL`.
DA-only increment diagnostics contribute their registered `value` series. The
21 Phase 1 matrix rows therefore expand to 43 models: 11 paired rows times three
series plus ten DA-only rows. `DA - OL` or `value` remains the primary estimand;
OL and DA are attribution controls.

The segmented design contains an intercept, a P1 baseline slope, 11
calendar-month fixed effects, level changes at every P2-P9 boundary, and slope
hinges for P2-P6, P8, and P9. It is full rank with 28 parameters over the
288-month record. P7 receives a level change but no slope hinge: its 15-month
slope is constrained to the P6 slope until P8 establishes a new slope. This is
a model constraint carried from the period registry, not an inference made from
the data.

Serial dependence is fitted with iterative Prais-Winsten AR(1). Primary
transition p-values and confidence intervals use a fixed-seed fitted-AR(1)
innovation-resampling bootstrap: innovations from the fitted whitened residuals are centered,
resampled, propagated through an AR(1), and the complete model is refitted.
Because a 28-parameter segmented fit can bias residual persistence downward,
the simulation AR(1) is the upper 95% large-sample confidence bound on the
fitted value (capped at 0.98); both values are reported. This makes uncertainty
conditional on a conservative persistence estimate rather than pretending the
estimated AR coefficient is known exactly. The default production analysis
uses 1,999 replicates. Centered
two-sided empirical p-values and basic bootstrap intervals are reported.
Newey-West/Bartlett and independent-sample standard errors remain in the output
as diagnostics; they do not determine transition significance. This choice is
deliberate: both OLS-HAC and Prais-Winsten-HAC alone were anti-conservative in
the fixed-seed 288-month AR(1) no-transition test, whereas the bootstrap passed
the BH false-discovery guard.

`scripts/validate_interrupted_time_series.py` provides the larger fixed-seed
calibration outside the unit suite. Its default 72-series experiment uses 24
AR(1)=0.8 no-transition series, 24 with a P6 level change, and 24 with a P6
slope change. The gate permits at most two null BH discoveries across the 15
transition families and requires at least 80% target-interval coverage, the
correct mean-effect direction, and relative mean bias no larger than 50%.
These are workflow regression limits, not universal power claims.

BH FDR at 0.05 is applied separately for each named boundary across all 43
domain series. Level changes have eight families, slope changes have seven, and
independently estimated period slopes have eight; P7 has no period-slope FDR
family. This preserves the scientific question represented by each transition
and gives 23 predeclared families rather than combining unrelated coefficients
in one omnibus correction.

`scripts/build_phase1_interrupted_series.py` writes three atomic products under
`output/trends_breakpoints/`: a coefficient table, a model-diagnostic table,
and a monthly NetCDF containing observed, fitted, residual, contributing-tile,
and represented-area series. `scripts/audit_phase1_interrupted_series.py`
independently reconstructs the 43-series contract, verifies all transition
terms and FDR families, enforces the P7 constraint, and confirms
`delta = DA - OL` for every paired domain series.
The coefficient table retains each response's source `units` and separately
labels `estimate_units`: level changes use the source units, while baseline,
slope-change, and period-slope estimates use source units per year.

## Independent changepoint detection

`config/changepoint_detection.json` and
`scripts/changepoint_detection.py` define the independent detection stage. The
implementation uses `ruptures` 1.1.10, pinned in
`requirements-changepoints.txt`; it does not reimplement PELT. It consumes the
same 43 complete area-weighted monthly series written by the known-transition
stage, preserving its run matrix, masks, units, paired support, and provenance.

Calendar-month effects are estimated jointly with one global linear trend and
only the month effects are removed. Detection is therefore performed on a
trend-preserving seasonally adjusted response. Two PELT methods are run across
the complete predeclared penalty grid:

- the primary method uses a piecewise AR(1), intercept, and linear-trend
  Gaussian profile-likelihood cost;
- the sensitivity method uses a Gaussian profile-likelihood linear cost after
  Prais-Winsten prewhitening.

The prewhitening AR(1) is the median of estimates from overlapping 60-month
locally detrended windows stepped by 30 months. This prevents a long level or
slope change from being mistaken for persistence. The ordinary global-trend
residual AR(1) is retained as a diagnostic. Each segment must contain at least
24 months.

The primary penalty is 1.25 times the BIC parameter-count penalty. Multipliers
0.5, 0.75, 1.0, 1.25, 1.5, and 2.0 are all retained in the output. A break is
accepted only when it appears in the primary method at 1.25 times BIC, recurs
within three months in at least half the primary-method penalty grid, and has a
prewhitened-method break within three months. Thus an accepted date is both
penalty-stable and cross-method, not merely the output of one favorable tuning.

Accepted dates are compared one-to-one with P2-P9. The primary agreement
tolerance is +/-3 months and +/-6 months is retained as a sensitivity. P7 can
be matched and displayed but is exempt from detection-agreement scoring because
its entry and exit are only 15 months apart under a 24-month segment floor.
Unmatched accepted breaks are retained; they are not forced onto the nearest
observing-system boundary.

`scripts/validate_changepoint_detection.py` tests seasonal white-noise and
AR(1)=0.8 nulls, an isolated P6 level shift, a gradual P6 slope hinge, two
opposing level shifts, and a level shift near the minimum-segment edge. The
fixed guard allows at most 0.25 accepted false breaks per null series and
requires at least 70% recovery within six months for the abrupt isolated,
two-break, and edge scenarios. Gradual slope recovery is reported but does not
control the gate: its poor date localization is an explicit method limitation,
and known-date slope inference remains the job of the interrupted-series model.

`scripts/build_phase1_changepoints.py` writes accepted-candidate diagnostics,
the complete penalty grid, P2-P9 comparisons, model diagnostics, and the
seasonally adjusted monthly fields. `scripts/audit_phase1_changepoints.py`
reconstructs the 43-series contract, verifies the 344 boundary rows and P7
exemption, checks every accepted break against the consensus/stability/segment
rules, and confirms the observed series are unchanged from the interrupted
series product.

## Core trend statistics

`config/trend_statistics.json` and `scripts/trend_statistics.py` implement the
first inferential stage. The primary field remains the strictly paired monthly
`DA - OL` series supplied by `trend_breakpoint_series.py`; OL and DA trends can
also be calculated as controls without changing their common sample.

The default support gate requires at least 60 valid months, at least 80% of the
months eligible under the selected mask, and at least 15 years between the first
and last valid values. The eligible-month denominator is explicit. Thus a
locally snowy mask is evaluated against its snowy months rather than all 288
calendar months. Every output retains `n_eligible`, `n_valid`, `n_trend`,
`valid_fraction`, `span_years`, and a status code.

Seasonality is removed with a trend-preserving calendar-month adjustment. A
preliminary exact Sen slope is removed before calculating each calendar month's
climatology; only that detrended seasonal climatology is then subtracted from
the original series. This avoids turning a perfectly linear monthly trend into
an artificial annual step function.

The reported slope is the exact Theil-Sen median pairwise slope in source units
per year. SciPy's nominal 95% Sen limits are retained as
`slope_ci_*_nominal`. The primary `slope_ci_*` interval expands each nominal
half-width by the square root of the same Hamed-Rao variance factor used for
significance. This is a first-order autocorrelation adjustment, not a replacement
for a block-bootstrap interval. `ci_excludes_zero` is a pointwise diagnostic;
mapped significance must use `significant_fdr`. The first-order interval is not
the inversion of the modified-MK test, particularly for zero-heavy or tied
increment diagnostics, so the two need not agree exactly. The output exposes
`fdr_ci_disagreement` wherever the FDR test is significant but the adjusted Sen
interval contains zero. This disagreement must be reported rather than forcing
one inferential product to mimic the other. A test-inverted or block-bootstrap
interval remains the appropriate later sensitivity if interval-based mapped
inference is required.

Significance uses a conservative adaptation of the Hamed-Rao modified
Mann-Kendall variance correction (Hamed and Rao, 1998,
doi:10.1016/S0022-1694(97)00125-X). Rank autocorrelation is calculated from
Sen-detrended residuals at actual monthly lags 1-24. Only significantly positive
lag correlations inflate the variance, and the factor cannot fall below one.
The ordinary independent-MK p-value is retained only to show the size of the
correction. Benjamini-Hochberg FDR at 0.05 is then applied across all finite
tiles in the output. Production FDR must therefore be calculated in a complete
field run, not independently in spatial chunks.

With intermittent eligibility or missing values, lag correlations pair values
that remain separated by the requested calendar-month lag, while the Hamed-Rao
weight uses the valid observation count. This is a documented gap-aware
approximation to the regularly sampled formula. `lag1_rank_autocorrelation`
reports the rank correlation used by the correction;
`lag1_residual_pearson_autocorrelation` is retained as a familiar residual
diagnostic and is not used in the variance factor.

`scripts/build_trend_statistics.py` writes the tile-level NetCDF contract. In
addition to support, slope, nominal and adjusted slope limits, corrected and
ordinary p-values, and BH q-values, it records both lag-1 autocorrelation
diagnostics, variance inflation, the number of positive autocorrelation lags,
and FDR significance.
Because BH correction depends on the complete family of tests, spatial subset
runs require an explicit diagnostic output path and are labeled as non-global
FDR. Only a complete-mask run supplies production FDR values.

Calendar months failing the configured climatology sample count are omitted
from the adjusted trend series. The output exposes a `calendar_month_used`
matrix and `n_calendar_months_used` per tile, and the builder prints their tile
frequency distribution. This makes seasonally selective dropping visible in
both NetCDF diagnostics and run logs.

`scripts/run_phase1_trends.py` executes the production matrix one complete
field at a time. Each NetCDF is written to an incomplete temporary path and
atomically renamed only after a successful close. Before each run, the runner
audits any existing output and skips it only if it passes, so an interrupted
batch can be resumed without recomputing completed fields. Diagnostic tile
subsets have distinct filenames and are explicitly labeled as non-global FDR.
The production runner defaults to eight joblib process workers. A 5,000-tile
real-data benchmark gave identical data variables and coordinates for serial,
four-process, and eight-process calculations; eight processes reduced elapsed
time from 49.5 seconds serial to 14.5 seconds. The backend and worker count do
not alter the FDR family and the output records its execution strategy.
`scripts/audit_phase1_trend_outputs.py` checks dimensions, required variables,
matrix provenance, status and finite-value support, p/q ranges, FDR flags,
explicit FDR/CI disagreement, and complete-field FDR scope. Derived batch manifests,
logs, and output audits live under `output/trends_breakpoints/`.

`scripts/validate_trend_statistics.py` measures false positives and power on
288-month synthetic fields with a fixed seasonal cycle and white or AR(1)
noise. With seed 20260812 and 100 series per scenario, the default configuration
gave 3% pointwise modified-MK positives and 0% BH detections for white-noise
no-trend series; 12% pointwise positives but 0% BH detections for AR(1)=0.8
no-trend series; and 100% BH detection for AR(1)=0.8 trends of +/-0.02 units per
year. Nominal Sen coverage was 43-56% in the AR(1) scenarios; first-order
adjusted coverage was 87-97%. These are regression benchmarks, not universal
operating characteristics.

## Audit gate

`scripts/audit_trend_breakpoint_inputs.py` will fail unless:

- all required monthly datasets cover June 2000-May 2024 with 288 months and
  112,573 tiles;
- state files identify OL v2 and DA v3 source experiments;
- raw increment and diagnostic products identify DA v3;
- every selected variable exists in the existing monthly-synthesis registry
  and input inventory;
- source units and dimensions agree with that inventory;
- the period registry validates, including the P7 constraint.

Selected source variables ending exactly in `_INC` or `_INCR` are rejected as
legacy monthly increment products. Raw diagnostics such as `_INC_MEAN`,
`_INC_ABS_MEAN`, and `_INC_RMS` remain valid because they are not cumulative
monthly increment fields.

The audit writes machine-readable reports under
`output/trends_breakpoints/`. Those reports are derived products and are not a
second source of scientific definitions.

## Production status

The complete Phase 1 matrix was generated on 2026-08-12. All 21 fields pass the
production audit, use a single source commit and embedded configuration, and
are summarized in `docs/phase1_trend_results.md`. This closes the first item in
the original roadmap.

The known-transition stage was also completed on 2026-08-12. All 43 expanded
area-weighted series and 23 transition FDR families pass the production audit;
the calibration and results are summarized in
`docs/phase1_interrupted_series_results.md`.

## Remaining roadmap

The next stages remain, in order:

1. Add autocorrelation-aware independent changepoint detection, using OL as the
   control and exempting P7 from agreement scoring.
2. Validate changepoint behavior on synthetic autocorrelated no-break
   and known-break series, including climatology and penalty sensitivities.
3. Produce global maps, area-weighted regional series, transition summaries,
   and a provenance-complete methods/results report only after those tests pass.
