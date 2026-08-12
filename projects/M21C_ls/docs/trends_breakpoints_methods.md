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

1. Implement paired Theil-Sen trends on `DA - OL`, with serial-correlation-aware
   significance and false-discovery-rate control across tiles.
2. Estimate known observing-system level and slope changes, treating P7 as
   level-only and reporting the short-segment cautions carried by the registry.
3. Add autocorrelation-aware independent changepoint detection, using OL as the
   control and exempting P7 from agreement scoring.
4. Validate trend and changepoint behavior on synthetic autocorrelated no-break
   and known-break series, including climatology and penalty sensitivities.
5. Produce global maps, area-weighted regional series, transition summaries,
   and a provenance-complete methods/results report only after those tests pass.
