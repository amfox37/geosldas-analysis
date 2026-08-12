# H121 ASCAT DA Project: Current Goals and Status

Last updated: 2026-08-12

Detailed evidence, methods, and caveats through this date are recorded in
[`ascat_da_work_summary_2026-07-25_to_2026-08-12.md`](ascat_da_work_summary_2026-07-25_to_2026-08-12.md).

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
Separate DA-matched OL summary products use scaled observations where a
scale-comparable O-F comparison is required; the raw and scaled products must
not be mixed.

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

Using the completed common-period diagnostics:

- For independent SMAP Tb diagnostics, both DA runs reduce O-F standard
  deviation versus OL, and H121 has a small advantage over legacy.
- For each run's own assimilated ASCAT species, H121 produces the larger
  absolute O-F standard-deviation reduction, although this is not independent
  evidence because each run is evaluated against its assimilated observation
  family.

The robust O-F conclusion is that both ASCAT DA experiments improve on OL.
The H121-versus-legacy difference is much smaller than either experiment's
improvement over OL.

### 5. Instrument-Variable and Triple-Collocation Evaluation

The Python port of the MATLAB IV/TC evaluation pipeline has been completed over
the full six-year DA output. The implementation was checked against live MATLAB
runs at each stage:

- daily pairs
- pentad climatology
- instrument-variable estimates
- triple collocation

A known non-blocking gap remains for the older H119/H120 ASCAT product:
MATLAB's natural-neighbor interpolation differs from Python's linear
`griddata`, causing persistent per-cell coverage differences. This has been
quantified and is not considered a logic bug for the H121-vs-legacy decision.

The completed estimates show that both ASCAT DA runs generally improve on OL,
but H121-versus-legacy differences are small and mixed across independent
reference combinations. The IV metric reported in figures is correlation `R`,
not `R2`.

### 6. ISMN In-Situ Evaluation

The common-sample ISMN evaluation is complete for surface and generic 0-1 m
root-zone soil moisture. All DA runs improve on OL. Legacy is slightly stronger
for surface `R` and anomaly `R`, while H121 is slightly stronger for root-zone
`R`, anomaly `R`, and ubRMSE. This supports ASCAT assimilation but does not give
a single all-depth verdict that H121 universally outperforms legacy.

### 7. FOV12.5 and Peat-QC Experiments

The first six months of the original H121, wider-footprint FOV12.5, and
FOV12.5 + peat-QC experiments have been compared. The peat QC causes complete
observation loss in core affected tiles and narrow bands of low-but-nonzero
counts near peat boundaries.

The low-retention edge population is defined from a clear distribution gap:
tiles must have at least 20 FOV12.5 observations and retain more than zero but
no more than 10% after peat QC. This identifies 179 Canada/Alaska tiles and 576
globally. A raw-OFA reconstruction for six representative boundary tiles
reproduces all 2,720 retain/reject decisions using the implemented
footprint-weighted 10% peat threshold. The boundary pattern is expected
per-observation footprint behavior, not a plotting or indexing bug.

## Decision Evidence So Far

Current evidence supports H121 as a credible replacement candidate, but the
comparison is more nuanced than the initial O-F-only assessment:

- The H121 reader and QC path work end-to-end.
- Independent Python processing reproduces and explains the main product
  differences.
- H121 gives substantially denser coverage.
- Both ASCAT DA runs improve OL in independent SMAP O-F, IV/TC, and ISMN
  diagnostics.
- H121 has a small advantage in independent SMAP O-F spread and ISMN root-zone
  correlation metrics.
- Legacy has a small advantage in ISMN surface correlations and some IV/TC
  combinations.
- Diagnostics against each run's assimilated ASCAT family are useful for fit
  assessment but cannot independently select the better observation stream.

The remaining decision work is to choose the smallest clear evidence figure
set, decide whether the FOV12.5/peat-QC configuration has the desired
coverage-quality balance, and formulate a recommendation that reflects the
small and mixed H121-versus-legacy skill differences.

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
| IV/TC | Instrument Variable / Triple Collocation |
