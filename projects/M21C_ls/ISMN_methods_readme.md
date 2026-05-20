# ISMN in-situ soil-moisture validation methods

This note documents the methods behind the ISMN-based in-situ soil-moisture
skill figure produced by:

- `projects/M21C_ls/notebooks/insitu_skill_cached_batch_figures.ipynb`

The figure notebook is cache-only. It reads:

- `outputs_ismn_network_skill/*_raw_timeseries.nc`

Those raw caches were built by:

- `projects/M21C_ls/notebooks/insitu_skill_ismn_network_ol_da.ipynb`

## Draft methods text

In situ soil moisture validation used station observations from the
International Soil Moisture Network (ISMN) local archive. The ISMN archive was
read with the Python `ismn` package, version 1.5.2; the local metadata cache
used by the workflow was `python_metadata/ISMN_data.csv`. The archive itself
does not appear to encode a formal ISMN release/version number in the local
files, so this should be reported as an ISMN archive snapshot/access date if the
download date can be recovered.

We selected ISMN networks with sufficient temporal coverage for the analysis
periods and with sensor depths suitable for surface and, where available,
root-zone comparisons. The networks used in the final hybrid in-situ figure were
SNOTEL, SCAN, USCRN, SMOSMANIA, OZNET, and ARM. SNOTEL, SCAN, and ARM were
evaluated for three periods: pre-ASCAT, 1 June 2000 to 31 May 2007; pre-SMAP,
1 June 2007 to 31 March 2015; and SMAP-era, 1 April 2015 to 1 June 2024. USCRN,
SMOSMANIA, and OZNET were evaluated for the pre-SMAP and SMAP-era periods.

For each station, all available ISMN soil-moisture sensors were inspected.
Surface soil moisture was taken from the available sensor depth nearest 0.05 m.
Root-zone soil moisture was constructed using the current `matlab_strict`
network-specific layer rules used in the cache-building workflow. These strict
composites require the expected set of layer observations at each timestamp;
missing required layers yield missing root-zone values. The strict root-zone
tags used in the current caches include `snotel_n3_c123smv`,
`scan_n4_c1234smv`, `uscrn_n4_c1234smv`, `smosm_n3_c3smv`, and
`oznet_n3_c3smv`.

Quality control retained only ISMN soil-moisture records flagged `G` ("good").
Records with other ISMN flags, including geophysical quality flags associated
with frozen conditions such as in-situ or GLDAS soil/air temperature below
0 degrees C, were excluded through this `G`-flag requirement. No additional
independent model-based snow-cover, snow-depth, or frozen-soil mask was applied
in the plotting workflow.

ISMN observations were averaged to daily means after duplicate timestamps were
averaged, and daily observation timestamps were shifted by 12 hours to align
with the GEOS-LDAS daily files. Model values were daily GEOS-LDAS `SFMC` and
`RZMC` from the OL and DA simulations on the `SMAP_EASEv2_M36_GLOBAL` grid.
Each station was matched to the nearest M36 tile using squared
latitude/longitude distance with a threshold of 0.1 degree squared; the cache
also stores the corresponding great-circle distance.

Statistics were computed from paired, non-missing observation-model samples.
Missing observations, missing model values, and root-zone times with incomplete
required layer availability were omitted pairwise. A station/domain/window was
retained only if it had at least 1000 valid observation days in each required
analysis window. Network summaries used common station membership across windows
for each metric and domain. The final plotted common-site counts were: SNOTEL
34 surface and 14 root-zone sites, SCAN 70 and 39, USCRN 88 and 55, SMOSMANIA
20 and 20, OZNET 16 and 10, and ARM 4 and 0.

For each station, we computed Pearson correlation, anomaly correlation, and
unbiased RMSE. Anomalies were defined by removing the day-of-year climatology
computed from the full paired observation-model record for that
station/domain/experiment, with at least 30 samples required for a day-of-year
climatology estimate. Final plotted values are paired network-mean differences
between DA and OL: `DA - OL` for correlation and anomaly correlation, and
`OL - DA` for unbiased RMSE, so positive values indicate improved DA performance
for all three plotted metrics.

## Short wording note

Do not describe this as simply "nearest depth and nearest time." More accurate
wording is:

> nearest available surface sensor depth to 0.05 m, daily averaging with
> 12-hour timestamp alignment, and exact timestamp pairing after daily
> collocation.

## Items still worth verifying

- Recover the ISMN archive download/access date, if possible.
- Confirm whether the manuscript figure is the hybrid output:
  `outputs_ismn_network_skill/batch_figures/all_networks_hybrid_OL_DA_delta_surface_rz_R_anomR_ubRMSE.png`.
- If reviewers ask specifically about snow-covered days, note that this workflow
  relies on ISMN `G` quality flags and does not apply a separate snow mask.
