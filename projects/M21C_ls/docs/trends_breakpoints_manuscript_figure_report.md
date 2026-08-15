# Trends And Observing-System Manuscript Figure Report

These figures are plotting-only products generated from the accepted Phase 1 trend, interrupted-series, and changepoint outputs. No scientific analysis was refitted.

## Outputs And Captions

### Main Figure 16: long-term state trends

- PNG: `projects/M21C_ls/output/paper_figures/fig16_longterm_rzmc_scf_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/fig16_longterm_rzmc_scf_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/fig16_longterm_rzmc_scf_trends.png`
- Sources: `RZMC_{ol,da,delta}_valid_land_trend_statistics.nc`, `FRLANDSNO_{ol,da,delta}_seasonal_snow_trend_statistics.nc`
- Provenance: Exact production slope; mapped significance is significant_fdr; Robinson projection; 60 S cutoff.

**Draft caption.** Long-term June 2000-May 2024 trends in (a-c) root-zone soil moisture (RZMC) and (d-f) snow-cover fraction (SCF) for the open-loop (OL), data-assimilation (DA), and paired DA-OL series. Trends are exact Theil-Sen slopes after trend-preserving removal of the calendar-month climatology. Black stippling denotes trends significant after Benjamini-Hochberg false-discovery-rate control at 0.05. RZMC uses the valid-land domain and SCF the Northern Hemisphere seasonal-snow domain. The DA-OL panels show the trend of the paired DA-OL series rather than the difference between independently estimated OL and DA trends.

### Main Figure 17: observing-system transitions

- PNG: `projects/M21C_ls/output/paper_figures/fig17_observing_system_transitions.png`
- PDF: `projects/M21C_ls/output/paper_figures/fig17_observing_system_transitions.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/fig17_observing_system_transitions.png`
- Sources: `phase1_changepoint_monthly.nc`, `phase1_interrupted_series_coefficients.csv`, `observing_system_registry.json`
- Provenance: Monthly production seasonal_adjusted fields; display-only full-record z scores; native-unit P6 inference.

**Draft caption.** Changes in soil-water data-assimilation behavior across the P1-P9 observing-system periods. (a) Standardized area-weighted monthly RZMC DA-OL and soil-water analysis-correction diagnostics during June 2000-May 2024. Background shading denotes the P1-P9 periods defined in Fig. 1; the P5-P6 boundary in April 2015 marks the introduction of SMAP brightness-temperature assimilation. Each seasonally adjusted series is standardized by its full-record mean and sample standard deviation for visual comparison only. (b-d) Estimated P5-P6 level changes from the interrupted time-series analysis for RZMC DA-OL, soil-moisture analysis-correction RMS, and prognostic soil-water correction activity. Symbols show the estimate and horizontal bars the 95% fitted-AR(1) bootstrap interval. Statistical significance uses boundary-family false-discovery-rate control at 0.05.

### Supporting Figure: precipitation trends

- PNG: `projects/M21C_ls/output/paper_figures/figSXX_precipitation_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figSXX_precipitation_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figSXX_precipitation_trends.png`
- Sources: `PRECTOTCORRLAND_{ol,da,delta}_valid_land_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on valid land.

**Draft caption.** Long-term precipitation trends for OL, DA, and the paired DA-OL series on valid land. Black stippling denotes production FDR significance. The common OL/DA pattern and null DA-OL result provide a forcing-control check.

### Supporting Figure: SFMC trends

- PNG: `projects/M21C_ls/output/paper_figures/figSXX_sfmc_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figSXX_sfmc_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figSXX_sfmc_trends.png`
- Sources: `SFMC_{ol,da,delta}_valid_land_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on valid land.

**Draft caption.** Long-term surface soil-moisture trends for OL, DA, and the paired DA-OL series on valid land. Black stippling denotes production FDR significance.

### Supporting Figure: snow mass and depth trends

- PNG: `projects/M21C_ls/output/paper_figures/figSXX_snow_mass_depth_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figSXX_snow_mass_depth_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figSXX_snow_mass_depth_trends.png`
- Sources: `SNOMASLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc`, `SNODPLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on the Northern Hemisphere seasonal-snow mask.

**Draft caption.** Long-term (a-c) snow-mass and (d-f) snow-depth trends for OL, DA, and the paired DA-OL series on the production Northern Hemisphere seasonal-snow mask. Black stippling denotes production FDR significance.

### Supporting Figure: breakpoint-boundary agreement

- PNG: `projects/M21C_ls/output/paper_figures/figSXX_breakpoint_boundary_agreement.png`
- PDF: `projects/M21C_ls/output/paper_figures/figSXX_breakpoint_boundary_agreement.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figSXX_breakpoint_boundary_agreement.png`
- Sources: `phase1_changepoint_boundary_comparison.csv`, `phase1_changepoint_detections.csv`
- Provenance: Accepted consensus breaks only; +/-3-month primary and +/-6-month sensitivity definitions retained.

**Draft caption.** Accepted independent changepoints relative to known P2-P9 dates for the primary Phase 1 estimands. Values are detected-minus-known months; blue is early, red late, white exact, and grey indicates no accepted match within six months. P7 is hatched because its 15-month duration is detection-exempt under the predeclared minimum-segment rule.

## Validation

- RZMC FDR-significant counts reproduce 4,909 OL, 10,329 DA, and 7,892 paired DA-OL tiles; significant DA-OL signs reproduce 7,267 positive and 625 negative.
- SCF area-weighted slopes reproduce -0.000554 yr-1 (OL), -0.000549 yr-1 (DA), and -0.000006 yr-1 (DA-OL).
- Precipitation reproduces 3,719 OL and 3,726 DA significant tiles, 3,603 same-sign overlaps, slope correlation 0.9998, and zero DA-OL FDR tiles.
- SFMC reproduces 6,992 OL, 8,966 DA, and 1,412 DA-OL FDR-significant tiles.
- Snow mass reproduces 12 OL, 5 DA, and zero DA-OL FDR tiles; snow depth reproduces 7 OL, 6 DA, and zero DA-OL FDR tiles.
- P6 begins 2015-04-01. All six plotted native-unit estimates, fitted-AR(1) bootstrap intervals, and boundary-family FDR flags reproduce the production coefficient table.
- Ten primary series have accepted breaks exactly in April 2015; nine also have significant known-date P6 level changes; no accepted break occurs in paired OL or DA state controls.
- The accepted-break inventory remains 37 total: 20 paired DA-OL and 17 correction diagnostics; 20 match within +/-3 months, two additional within +/-6 months, and 15 remain unmatched.
- All maps use exact production Theil-Sen slopes and `significant_fdr`; pointwise confidence-interval exclusion is not used for mapped inference.

## Plotting Choices

- Maps follow the existing report convention: Robinson projection, 60 S cutoff, grey land, thin coastlines, segmented RdBu_r scales centered on a white zero bin, and black stippling.
- OL and DA share a symmetric color scale within each row. DA-OL uses a separately labeled symmetric scale where needed; snow mass/depth retain the OL/DA scale so negligible differences are not visually exaggerated.
- Figure 17 panel (a) shows unsmoothed monthly production series. Full-record z scoring is display-only; interrupted-series inference remains in native units.
- The requested optional all-boundary decorative summary was not produced; the accepted breakpoint-agreement matrix already provides the auditable all-boundary view.

## Discrepancies

- No regenerated numerical result disagreed with `m21c_trends_breakpoints_report.md`.
- The task summary says P9 begins in November 2021, but the authoritative registry starts P9 on 2021-12-01 after P8 ends 2021-11-30. The figures use the registry date.

PNGs are exported at 300 DPI. PDFs retain vector text, coastlines, markers, and intervals; dense map fields and stippling are rasterized within the vector PDF.
