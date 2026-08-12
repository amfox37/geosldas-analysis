# Legacy M21C trend exploration (2025)

These files preserve the trend exploration run from 30 October through
3 November 2025. They are retained for provenance only and are not the current
M21C trend workflow.

## Contents

- `trend_analysis.py`: initial linear-trend implementation.
- `theilsen_mk_states_and_increments.py`: Theil-Sen slopes and
  Mann-Kendall significance.
- `plot_trends.ipynb`: OL, DA, and DA-minus-OL trend maps and global summary
  statistics.

The files were originally tracked under `common/python/stats/` and
`projects/utils/notebooks/`. They are project-specific: both scripts contain
hard-coded M21C input paths.

## Original inputs

The calculations used the older local products under
`/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/`:

- archived `LS_OLv8_M36` monthly states;
- `LS_DAv8_M36_v2` monthly states;
- the legacy `LS_monthly_increments_2000_2024.nc` product.

The current M21C analysis uses consolidated products built from
`LS_OLv8_M36_v2` and `LS_DAv8_M36_v3`. The old and current state files are not
numerically identical.

## Known limitations

1. The scripts cumulatively sum `SFMC_INC` and `RZMC_INC` from the legacy
   increment file. Those fields are monthly diagnostic ANA-minus-FCST state
   differences in `m3 m-3`, not cumulative water-mass increments, and must not
   be interpreted as accumulated water.
2. Current cumulative mass diagnostics come from raw `catch_progn_incr`
   summaries. Current diagnostic state-correction activity comes from raw
   `inst3_1d_lndfcstana_Nt` mean, absolute-mean, and RMS summaries.
3. Full-record trends cross several observing-system changes. The old analysis
   did not model those known transitions or estimate changepoints.
4. The old Mann-Kendall significance maps did not apply the serial-correlation
   and field-significance controls planned for the replacement workflow.

Do not cite the old maps as current M21C results.

## Surviving outputs

The original untracked outputs remain in
`/Users/amfox/Desktop/GEOSldas_diagnostics/Jupyter/`, including:

- `LS_theilsen_mk_states_increments.nc`
- `LS_trends_masked.nc`
- `trend_significance_masks.nc`
- `trend_summary.csv`
- twelve `*_trend_map.png` files

See `projects/M21C_ls/docs/raw_increment_summaries_methods.md` and
`projects/M21C_ls/docs/monthly_synthesis_report_out.md` for the corrected
monthly-product definitions used by later work.
