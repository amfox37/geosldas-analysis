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
per year. SciPy's nominal 95% Sen slope limits are retained as descriptive
diagnostics. They are not autocorrelation-adjusted and must not be used as the
primary significance test; synthetic AR(1) validation confirms undercoverage.

Significance uses a conservative adaptation of the Hamed-Rao modified
Mann-Kendall variance correction (Hamed and Rao, 1998,
doi:10.1016/S0022-1694(97)00125-X). Rank autocorrelation is calculated from
Sen-detrended residuals at actual monthly lags 1-24. Only significantly positive
lag correlations inflate the variance, and the factor cannot fall below one.
The ordinary independent-MK p-value is retained only to show the size of the
correction. Benjamini-Hochberg FDR at 0.05 is then applied across all finite
tiles in the output. Production FDR must therefore be calculated in a complete
field run, not independently in spatial chunks.

`scripts/build_trend_statistics.py` writes the tile-level NetCDF contract. In
addition to support, slope, nominal slope limits, corrected and ordinary
p-values, and BH q-values, it records lag-1 residual autocorrelation, variance
inflation, the number of positive autocorrelation lags, and FDR significance.
Because BH correction depends on the complete family of tests, spatial subset
runs require an explicit diagnostic output path and are labeled as non-global
FDR. Only a complete-mask run supplies production FDR values.

`scripts/validate_trend_statistics.py` measures false positives and power on
288-month synthetic fields with a fixed seasonal cycle and white or AR(1)
noise. With seed 20260812 and 100 series per scenario, the default configuration
gave 3% pointwise modified-MK positives and 0% BH detections for white-noise
no-trend series; 12% pointwise positives but 0% BH detections for AR(1)=0.8
no-trend series; and 100% BH detection for AR(1)=0.8 trends of +/-0.02 units per
year. These are regression benchmarks, not universal operating characteristics.

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

## Remaining roadmap

The next stages remain, in order:

1. Run and review complete Phase 1 trend fields for paired soil moisture, snow,
   precipitation context, and the valid DA-only increment diagnostics.
2. Estimate known observing-system level and slope changes, treating P7 as
   level-only and reporting the short-segment cautions carried by the registry.
3. Add autocorrelation-aware independent changepoint detection, using OL as the
   control and exempting P7 from agreement scoring.
4. Validate trend and changepoint behavior on synthetic autocorrelated no-break
   and known-break series, including climatology and penalty sensitivities.
5. Produce global maps, area-weighted regional series, transition summaries,
   and a provenance-complete methods/results report only after those tests pass.
