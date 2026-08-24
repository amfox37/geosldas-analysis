# Trends And Observing-System Manuscript Figure Report

These figures are plotting-only products generated from the accepted Phase 1 trend outputs. No analysis was recomputed.

## Outputs And Captions

### Figure 16: long-term state trends

- PNG: `projects/M21C_ls/output/paper_figures/fig16_longterm_rzmc_scf_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/fig16_longterm_rzmc_scf_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/fig16_longterm_rzmc_scf_trends.png`
- Sources: `RZMC_{ol,da,delta}_valid_land_trend_statistics.nc`, `FRLANDSNO_{ol,da,delta}_seasonal_snow_trend_statistics.nc`
- Provenance: Exact production slope; mapped significance is significant_fdr; Robinson projection; 60 S cutoff.

**Draft caption.** Figure 16. Long-term June 2000-May 2024 trends in (a-c) root-zone soil moisture (RZMC) and (d-f) snow-cover fraction (SCF) for the open-loop (OL), data-assimilation (DA), and paired DA-OL series. Trends are exact Theil-Sen slopes after trend-preserving removal of the calendar-month climatology. Black stippling denotes trends significant after Benjamini-Hochberg false-discovery-rate control at 0.05. RZMC uses the valid-land domain and SCF the Northern Hemisphere seasonal-snow domain. The DA-OL panels show the trend of the paired DA-OL series rather than the difference between independently estimated OL and DA trends. Regional RZMC trends are modified substantially by DA, whereas OL and DA exhibit nearly identical long-term SCF trends.

### Supplemental Figure S5: precipitation trends

- PNG: `projects/M21C_ls/output/paper_figures/figS05_precipitation_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figS05_precipitation_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figS05_precipitation_trends.png`
- Sources: `PRECTOTCORRLAND_{ol,da,delta}_valid_land_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on valid land.

**Draft caption.** Figure S5. Long-term precipitation trends for OL, DA, and the paired DA-OL series on valid land. Black stippling denotes production FDR significance. The common OL/DA pattern and null DA-OL result provide a forcing-control check.

### Supplemental Figure S6: SFMC trends

- PNG: `projects/M21C_ls/output/paper_figures/figS06_sfmc_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figS06_sfmc_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figS06_sfmc_trends.png`
- Sources: `SFMC_{ol,da,delta}_valid_land_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on valid land.

**Draft caption.** Figure S6. Long-term surface soil-moisture trends for OL, DA, and the paired DA-OL series on valid land. Black stippling denotes production FDR significance.

### Supplemental Figure S7: snow mass and depth trends

- PNG: `projects/M21C_ls/output/paper_figures/figS07_snow_mass_depth_trends.png`
- PDF: `projects/M21C_ls/output/paper_figures/figS07_snow_mass_depth_trends.pdf`
- Tracked review PNG: `projects/M21C_ls/docs/paper_figures/figS07_snow_mass_depth_trends.png`
- Sources: `SNOMASLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc`, `SNODPLAND_{ol,da,delta}_seasonal_snow_trend_statistics.nc`
- Provenance: Exact production slope and significant_fdr on the Northern Hemisphere seasonal-snow mask.

**Draft caption.** Figure S7. Long-term (a-c) snow-mass and (d-f) snow-depth trends for OL, DA, and the paired DA-OL series on the production Northern Hemisphere seasonal-snow mask. Black stippling denotes production FDR significance.

## Validation

- RZMC FDR-significant counts reproduce 4,909 OL, 10,329 DA, and 7,892 paired DA-OL tiles; significant DA-OL signs reproduce 7,267 positive and 625 negative.
- SCF area-weighted slopes reproduce -0.000554 yr-1 (OL), -0.000549 yr-1 (DA), and -0.000006 yr-1 (DA-OL).
- Precipitation reproduces 3,719 OL and 3,726 DA significant tiles, 3,603 same-sign overlaps, slope correlation 0.9998, and zero DA-OL FDR tiles.
- SFMC reproduces 6,992 OL, 8,966 DA, and 1,412 DA-OL FDR-significant tiles.
- Snow mass reproduces 12 OL, 5 DA, and zero DA-OL FDR tiles; snow depth reproduces 7 OL, 6 DA, and zero DA-OL FDR tiles.
- All maps use exact production Theil-Sen slopes and `significant_fdr`; pointwise confidence-interval exclusion is not used for mapped inference.

## Plotting Choices

- Maps follow the existing report convention: Robinson projection, 60 S cutoff, grey land, thin coastlines, segmented RdBu_r scales centered on a white zero bin, and black stippling.
- OL and DA share a symmetric color scale within each row. Figure 16 retains a separate RZMC DA-OL scale but places SCF DA-OL on the OL/DA scale; snow mass/depth likewise retain the OL/DA scale so negligible differences are not visually exaggerated.

## Discrepancies

- No regenerated numerical result disagreed with `m21c_trends_breakpoints_report.md`.
- The task summary says P9 begins in November 2021, but the authoritative registry starts P9 on 2021-12-01 after P8 ends 2021-11-30. The figures use the registry date.

PNGs are exported at 300 DPI. PDFs retain vector text, coastlines, markers, and intervals; dense map fields and stippling are rasterized within the vector PDF.
