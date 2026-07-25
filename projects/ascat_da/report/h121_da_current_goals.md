# H121 ASCAT DA Project: Current Goals and Status

Last updated: 2026-07-24

## Core Question

GEOSldas currently assimilates soil moisture from the legacy EUMETSAT ASCAT
BUFR product, a 25 km per-swath product from Metop-A/B/C. H SAF now provides
H121, a 12.5 km Fibonacci-grid NetCDF Climate Data Record, with H139 as the
near-real-time counterpart.

The central project question is:

**Should GEOSldas switch its ASCAT soil-moisture assimilation source from the
legacy BUFR product to H121 CDR, and does H121 improve the land analysis?**

The work is split between:

- GEOSldas experiment runs:
  `/discover/nobackup/projects/land_da/hsaf_cdr_test/`
- Python analysis, validation, scaling, and evaluation tooling:
  `/discover/nobackup/projects/land_da/geosldas-analysis/`

## Current Experiment Design

The core production comparison is a pair of matched 6-year DA runs:

| Experiment | Period | Assimilated ASCAT source |
| --- | --- | --- |
| `DAv7_M36_ASCAT_type_13_H121` | 2015-04-01 to 2021-04-01 | H SAF H121 CDR |
| `DAv7_M36_ASCAT_type_13_legacy` | 2015-04-01 to 2021-04-01 | Legacy EUMETSAT BUFR |

The two DA runs are intended to be identical except for the ASCAT observation
stream being assimilated. They share forcing, boundary conditions, restart
lineage, and matched obs-error scaling methodology.

Both production runs have completed:

- H121: job `57274390`, finished 2026-07-22
- Legacy: job `57274391`, finished 2026-07-23

A companion open-loop run, `OLv7_M36_MULTI_type_13_H121`, validated the H121
reader end-to-end and provides an OL reference background for O-F diagnostics.
That run is intentionally unscaled (`scale=.false.`), because its role is raw
observation infrastructure validation and background comparison, not DA.

## Work Threads

### 1. H121 Reader and QC in GEOSldas

A new Fortran reader, `read_obs_sm_ASCAT_HSAF`, was added in
`clsm_ensupd_read_obs.F90` on the nested `GEOSldas_GridComp` branch
`feature/amfox/ascat-hsaf-v8`.

The H121 QC sequence was built from the H SAF literature and empirical checks:

- reject open water
- reject bad model/sigma0 processing flags
- reject noisy backscatter using `backscatter40_flag` bit 4 only
- reject wetland/topographic complexity `>= 10%`
- reject subsurface scattering `>= 5%`; fill accepted
- reject soil-moisture sensitivity `<= 1 dB`

Important finding: tightening the subsurface-scattering threshold from 10% to
5% was the dominant data-removal change, and it was net-negative for fit quality
in a 10-day test. It disproportionately removed well-fitting observations
concentrated between 10N and 50N. The sensitivity and backscatter-noise changes
were comparatively minor and better targeted.

### 2. Independent Python Validation

The H121-vs-legacy comparison was reproduced independently from raw files in
`projects/ascat_da`, so QC and product-difference conclusions were not based
only on GEOSldas output.

Key findings:

- H121 provides about 2.8 times denser observation coverage than legacy.
- H121 fits the model background better in aggregate, with lower O-F bias/RMSE.
- H121 and legacy correlate at about 0.80 at matched place/time, so they are
  related but not interchangeable.
- The apparent result that H121 is about 20% saturation drier than legacy above
  50N was reframed. Legacy's per-pass `smpf` QC rejects about 60% of
  high-latitude observations on retrieval-noise flags that H121's
  climatology-based retrieval does not have. Comparing each product against the
  model background showed legacy fit degrading above 50N, while H121 changes
  much less.
- An old 2024-era GEOSldas build double-counted ASCAT obs across consecutive
  assimilation windows, inflating counts by about 35-37% in one archived
  reference run. That bug is absent from the current codebase.

Detailed writeups are in `projects/ascat_da/report/`.

### 3. Obs-Error Scaling

A Python port of the MATLAB z-score climatology generator was built:

```text
obs_scaled = m_mean + (obs_raw - o_mean) * (m_std / o_std)
```

The Python implementation was parity-checked against the MATLAB reference and
used to produce the H121-specific and legacy-specific scaling climatologies for
the production DAv7 runs.

A Fortran filename-truncation issue around 80-character paths was worked around
using short symlinks.

### 4. O-F Diagnostic Evidence

The current O-F diagnostic compares O-F standard deviation, not O-F mean,
because scaling forces the mean toward zero by construction.

Using the common 2015-04 to 2018-03 window available at the time of analysis:

- For independent SMAP Tb diagnostics, both DA runs reduce O-F standard
  deviation versus OL in nearly every month/channel, and H121 consistently beats
  legacy by a small but real margin.
- For each run's own assimilated ASCAT species, H121 produces the larger
  absolute O-F standard-deviation reduction.

These two lines of evidence currently indicate that H121 assimilation is a real
improvement over legacy assimilation, rather than sampling noise.

### 5. Independent Validation and Triple Collocation

The active next phase is a from-scratch Python port of the MATLAB IV/TC
evaluation pipeline over the full 6-year DA output.

The Python port has been verified against live MATLAB runs at each stage:

- daily pairs
- pentad climatology
- independent validation
- triple collocation

A known non-blocking gap remains for the older H119/H120 ASCAT product:
MATLAB's natural-neighbor interpolation differs from Python's linear
`griddata`, causing persistent per-cell coverage differences. This has been
quantified and is not considered a logic bug for the H121-vs-legacy decision.

The current full-span IV/TC jobs are running for DA-H121, DA-legacy, H119/H120,
and SMAP, with a TC job queued for SMOS-IC, ASCAT H121 DA, and SMAP. This will
provide the independent SMOS-IC/SMAP/CYGNSS verdict on H121 versus legacy.

## Decision Evidence So Far

Current evidence favors H121:

- The H121 reader and QC path work end-to-end.
- Independent Python processing reproduces and explains the main product
  differences.
- H121 gives substantially denser coverage.
- H121 has better aggregate fit to the model background.
- O-F standard-deviation diagnostics show H121 DA beating legacy DA against
  both independent SMAP Tb observations and the assimilated ASCAT stream.

The remaining major confirmation step is the full 6-year IV/TC evaluation
against independent references.

## Naming Notes

| Term | Meaning |
| --- | --- |
| OL | Open loop, no assimilation |
| DA | Data assimilation |
| v7 | Current H121 CDR production experiment generation |
| v8 | Older or separate campaign generation in nearby work |
| `type_13` | Grid/config identifier used in OLv7/DAv7 H121 experiment names |
| legacy | EUMETSAT ASCAT BUFR, 25 km, currently operational |
| H121 | H SAF ASCAT CDR, 12.5 km, NetCDF candidate replacement |
| H139 | Near-real-time counterpart to H121 |
| O-F | Observation minus forecast innovation |
| IV/TC | Independent Validation / Triple Collocation |

