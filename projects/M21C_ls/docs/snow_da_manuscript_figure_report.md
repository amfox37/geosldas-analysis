# Snow-DA Manuscript Figure Report

The four figures below are plotting-only products built from the accepted machine-readable analysis tables. No statistical analysis was recomputed.

## Outputs

| Figure | PNG | PDF | PNG dimensions | Figure size |
|---|---|---|---:|---:|
| Figure 14: snow-DA water-budget accounting | `projects/M21C_ls/output/paper_figures/fig14_snow_da_water_budget_4panel.png` | `projects/M21C_ls/output/paper_figures/fig14_snow_da_water_budget_4panel.pdf` | 4200 x 3120 px at 300 DPI | 14.0 x 10.4 in |
| Figure 15: monthly snow-DA pathway | `projects/M21C_ls/output/paper_figures/fig15_snow_da_monthly_pathway.png` | `projects/M21C_ls/output/paper_figures/fig15_snow_da_monthly_pathway.pdf` | 3484 x 2614 px at 300 DPI | 11.5 x 8.6 in |
| Supplemental Figure S3: non-overlapping attribution | `projects/M21C_ls/output/paper_figures/fig_supp_snow_da_octmar_attribution.png` | `projects/M21C_ls/output/paper_figures/fig_supp_snow_da_octmar_attribution.pdf` | 4084 x 2325 px at 300 DPI | 13.5 x 7.4 in |
| Supplemental Figure S4: accounting-boundary sensitivity | `projects/M21C_ls/output/paper_figures/fig_supp_snow_da_boundary_sensitivity.png` | `projects/M21C_ls/output/paper_figures/fig_supp_snow_da_boundary_sensitivity.pdf` | 3484 x 1594 px at 300 DPI | 11.5 x 5.2 in |

Tracked review PNGs:
- `projects/M21C_ls/docs/paper_figures/fig14_snow_da_water_budget_4panel.png`
- `projects/M21C_ls/docs/paper_figures/fig15_snow_da_monthly_pathway.png`
- `projects/M21C_ls/docs/paper_figures/fig_supp_snow_da_octmar_attribution.png`
- `projects/M21C_ls/docs/paper_figures/fig_supp_snow_da_boundary_sensitivity.png`

## Sources And Changes

### Figure 14: snow-DA water-budget accounting

Sources:
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/annual_domain_budgets.csv`
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/six_year_integrated_partitions.csv`
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness/water_year_boundary_sensitivity.csv`

Difference from the old figure: Combines the old annual and positive-partition diagnostics into one main-paper hierarchy: six all-tile annual budgets at left and the positive-input fate at right, with runoff components explicitly joined and the residual treated as closure.

### Figure 15: monthly snow-DA pathway

Sources:
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/monthly_climatology_snow_addition.csv`
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/soil_moisture_summary_snow_addition.csv`
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/water_year_snow_da_budget/soil_moisture_peak_timing.csv`

Difference from the old figure: Retains the old pathway concept but removes infiltration and persistence/residence-time claims, emphasizes RZMC, and aligns snow input, snowmelt, soil moisture, runoff, and ET on one Oct-Sep axis.

### Supplemental Figure S3: non-overlapping attribution

Sources:
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness/analysisA_octmar_signed_controls.csv`

Difference from the old figure: Expands the preliminary four-panel control sequence to all six responses, retains native units in separate panels, and makes the infiltration interval crossing zero explicit.

### Supplemental Figure S4: accounting-boundary sensitivity

Sources:
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness/water_year_boundary_sensitivity.csv`
- `projects/M21C_ls/output/monthly_synthesis_diagnostics/targeted_snow_hydrology_robustness/water_year_september_timing_diagnostic.csv`

Difference from the old figure: Replaces the duplicated bar chart with paired points and matched intervals, retaining total runoff and its surface/baseflow split while moving September timing facts to the caption notes.

## Validation

- Figure 14a reproduces the six accepted WY2001-WY2006 all-tile budget values and each year closes as `input = runoff + ET + storage + residual`.
- Figure 14b uses 247,545 positive-input tile-years. Surface runoff plus baseflow equals total runoff, and runoff + ET + storage + residual equals 100% before rounding.
- Figure 14b total-runoff interval is the positive-input 5-degree spatial-block interval: 61.1-67.2%.
- Figure 15 uses the same `I_snow > 0` sample (247,545 tile-years) as the accepted pathway analysis. Mean peak RZMC is 0.0189 m3 m-3, May is the most common peak month, and mean MJJ RZMC is 0.0108 m3 m-3.
- The Oct-Mar supplemental predictor contains exactly October-March and has zero overlap with the AMJ and MJJ response windows in all six years.
- Both boundary rows come from the same targeted output table and use the same area-weighted positive-input partition and 5-degree block-bootstrap machinery.
- All figure samples end by September 2006, before microwave soil-moisture DA begins in June 2007.

## Caption Notes

September signed input is 9.01 kg m-2 (15.5% of Oct-Sep annual input), September DA-OL snowmelt is 9.16 kg m-2, and the snow-mass difference is only 0.24 kg m-2. These belong in the boundary-sensitivity caption rather than in the graphic.

PDFs use embedded TrueType text (`pdf.fonttype = 42`) and vector lines, markers, bars, and error bars. PNGs are exported at 300 DPI.
