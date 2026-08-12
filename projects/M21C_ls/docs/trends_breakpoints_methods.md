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

The audit writes machine-readable reports under
`output/trends_breakpoints/`. Those reports are derived products and are not a
second source of scientific definitions.
