# Replacing Legacy ASCAT with the H SAF CDR in a Global Soil Moisture Data Assimilation System

*Reader development, quality control, and a peatland footprint correction for
H121/H139 ASCAT soil moisture, evaluated against six years of global ensemble
land data assimilation.*

Internal technical note · GEOSldas M36 global configuration · **25 km EUMETSAT
ASCAT** vs. **12.5 km H SAF H121 CDR / H139 ICDR ASCAT** · 2015–2021 study
period

## Summary

We added support for the H SAF H121 Climate Data Record (CDR) and H139
Interim CDR ASCAT soil moisture retrievals to a global land data assimilation
system, replacing the coarser legacy EUMETSAT ASCAT product used previously.
Across six years of parallel open-loop and data-assimilation experiments, the
reprocessed product is denser, fits the model background better, and — once
assimilated — produces a small but consistent improvement over the legacy
product on two independent lines of evidence. A targeted quality-control
extension that screens observations by footprint peatland fraction removes
the least reliable ~6% of observations, concentrated in high-latitude organic
soils, with no detectable cost elsewhere. One regional exception, over the
Hudson Bay Lowlands, remains only partially resolved and is identified here
as a target for further work.

---

## 1. Motivation

The assimilation system's existing ASCAT observation stream derives from the
legacy EUMETSAT surface soil moisture product: a 25 km (SOMO25) retrieval
distributed as BUFR messages, screened against an external, unversioned
land-cover mask. The H SAF programme has since released a reprocessed Climate
Data Record (H121, 2007–2021) and a continuing Interim CDR (H139), both
retrieved on a native 12.5 km Fibonacci swath grid and distributed as NetCDF
with a much richer set of per-observation quality flags — open-water
fraction, wetland and topographic-complexity fractions, subsurface scattering
probability, backscatter noise flags, and soil-moisture retrieval sensitivity
— embedded directly in the file, rather than depending on a separate static
mask.

This motivated three linked pieces of work: (i) building a reader for the new
product family into the assimilation system, (ii) designing quality control
appropriate to its richer flag set, and (iii) running a controlled,
like-for-like comparison against the legacy product to determine whether
switching observation types is actually beneficial before committing to it
for production reprocessing.

## 2. Observation processing and quality control

The new reader ingests H121/H139 swath granules for MetOp-A, -B and -C,
matching each retrieval to its host model tile and 3-hourly assimilation
window. Quality control is applied as a sequential chain, each check
targeting a specific known ASCAT error source:

1. **Open water** and **poor model/backscatter agreement** — standard
   radar-retrieval screening.
2. **Noisy backscatter** — observations flagged for out-of-limits noise in
   the underlying backscatter measurement are rejected outright; a separate,
   milder "slightly degraded" flag is deliberately *not* used to reject,
   since the assimilation applies a fixed observation-error variance with no
   mechanism to down-weight a noisier-but-usable observation.
3. **Wetland and topographic complexity** — observations over tiles with
   high static wetland or topographic-complexity fraction are rejected (both
   known radar-retrieval confounders); missing coverage in the static
   database is also treated as a reject, since safety cannot be confirmed.
4. **Subsurface scattering** — observations with high modelled probability of
   subsurface volume scattering (a dry/sandy/arid-soil artifact) are
   rejected; a missing value here is treated as a safe accept, since it means
   no scattering was detected at that location, not that the field is
   unavailable.
5. **Retrieval sensitivity** — observations with low soil-moisture
   sensitivity (below 1 dB) are rejected, following the retrieval-uncertainty
   characterization in the H SAF product literature (Hahn et al., 2026).

Snow and frozen-ground conditions are deliberately *not* re-screened in the
reader itself; that check already exists as a shared model-state-based
quality check applied to all satellite soil moisture products in the system,
and duplicating it in the reader was found to be redundant.

> **A tuning lesson.** An early, more aggressive version of the
> subsurface-scattering threshold removed roughly 15% of observations but
> slightly *increased* aggregate model-fit error over a short evaluation
> window — it was preferentially discarding observations that already agreed
> well with the model, concentrated in the 10–50°N latitude band. The final,
> more conservative threshold used throughout this report avoids that
> regression.

## 3. Peatland footprint quality control

ASCAT backscatter retrievals are known to degrade over organic, high-porosity
peat soils, where the dielectric response departs from the mineral-soil
assumptions built into the retrieval algorithm. Because a single ASCAT
footprint (radius on the order of tens of kilometres) can straddle a mixture
of mineral and peat land cover, a single tile-level land-cover flag is not
sufficient: a retrieval can be corrupted by nearby peat even when its nominal
center tile is mineral soil.

We instead compute a **footprint peat fraction**: every model tile within the
sensor's effective footprint is weighted by a Gaussian function of its
distance from the footprint center and its own area, and classified as peat
or mineral by its porosity (peat soils have distinctively higher porosity
than the surrounding mineral catchments). An observation is rejected if the
resulting area- and distance-weighted peat fraction exceeds 10%.

> **Mechanism independently verified.** To confirm this footprint-averaging
> logic behaves as intended rather than merely correlating with rejections,
> the algorithm was reimplemented independently from the model's own porosity
> field and tile geometry and tested against every candidate observation at
> six representative peatland-boundary locations. It reproduced the actual
> accept/reject outcome for **100% of 2,720 test observations**, with zero
> disagreements — direct confirmation that footprint peat fraction crossing
> the 10% threshold, not a coincidental correlate, is the mechanism driving
> rejection near peatland boundaries.

## 4. Experimental design

Four experiment configurations, all sharing the same 36 km global tile mesh,
ensemble size (24 members), forcing, and restart lineage, isolate the effect
of each change:

| Configuration | ASCAT product | Assimilated? | Period |
|---|---|---|---|
| **Open loop (OL)** | — | No (model only) | 2015–2021 |
| **DA — legacy** | EUMETSAT, 25 km | Yes | 2015–2021 |
| **DA — H121** | H SAF H121 CDR, 12.5 km | Yes | 2015–2021 |
| **DA — H121 + peatland QC** | H SAF H121 CDR, 12.5 km | Yes | 6-month test, 2015 |

Two complementary diagnostics are used throughout, both expressed as the
standard deviation of observation-minus-forecast (O−F) residuals rather than
the mean, since scaling (CDF/z-score matching) forces the mean toward zero by
construction and is not informative about fit quality:

- **Self-fit** — O−F spread of each experiment's own assimilated ASCAT
  observations. Reflects how well the assimilation is fitting the data it is
  given, but is not independent of the assimilation itself.
- **Independent check** — O−F spread of SMAP brightness-temperature
  observations, which are monitored but never assimilated in any of these
  experiments. Because SMAP is untouched by the ASCAT choice, any change here
  reflects a genuine change in the underlying land-surface state estimate,
  not an assimilation bookkeeping artifact.

Spatial aggregation weights tiles by land area; combining across satellites
or channels within a diagnostic weights by observation count. Where the open
loop (which uses unscaled observations) is compared directly against a DA
experiment (which uses scaled observations), the open loop's forecasts are
re-matched against that experiment's own scaled observation values, so the
comparison is not confounded by the scaling step itself.

## 5. Results

### 5.1 The reprocessed product is denser and fits differently, not just better

| | |
|---|---|
| **2.8×** | more observations per tile-month than legacy, same satellites |
| **0.80** | correlation with legacy at matched place & time |
| **−0.09** | mean offset, H121 reads systematically drier (volumetric units) |
| **18%** | lower background-fit RMSE than legacy (June 2020 sample) |

Direct collocation of the two products at the same satellite, tile, and
3-hour window shows they are related (r ≈ 0.80) but not interchangeable: H121
reads systematically drier than legacy with lower variance, reflecting
genuine retrieval differences from the finer native grid and updated
algorithm, not merely a processing artifact. On its own terms it also agrees
better with the model background than legacy does over the same evaluation
month.

### 5.2 Assimilation improves on the open loop, and H121 improves on legacy

**Independent check (SMAP Tb O−F stdev, 2015-04 to 2018-03 mean):**

| Configuration | O−F stdev (K) |
|---|---|
| Open loop | 3.888 |
| DA · legacy | 3.794 |
| DA · H121 | 3.783 |

Both DA configurations reduce spread relative to the open loop in nearly
every one of 36 months; H121 is consistently, if modestly, ahead of legacy.

**Self-fit (each experiment's own assimilated ASCAT species, m³/m³, 2015-04
to 2018-03):**

| Configuration | O−F stdev |
|---|---|
| Legacy · OL (re-matched) | 0.0217 |
| Legacy · DA | 0.0206 |
| H121 · OL (re-matched) | 0.0245 |
| H121 · DA | 0.0230 |

H121 shows the larger absolute reduction on its own directly-assimilated
data — a second, independent line of evidence pointing the same direction as
the independent check above.

### 5.3 Peatland footprint QC: a targeted correction, not a blunt one

| | |
|---|---|
| **−6.4%** | observations removed globally |
| **40–50%** | removed locally, in peatland tiles |
| **−4.4%** | tighter ASCAT self-fit in affected regions |
| **≈0** | change to the independent SMAP check |

Applying the footprint peatland QC over a 6-month evaluation window removes
1.65 million of 26.0 million global ASCAT observations (−6.4%), but this is
not spread evenly: rejections are concentrated almost entirely over the
Hudson Bay Lowlands, Scandinavia, and West Siberia, where individual tiles
lose 40–50% of their observations while most of the globe is essentially
untouched. Where it does act, ASCAT self-fit improves by 4.4%; the
independent SMAP check is unchanged to three significant figures. This is the
signature of a well-targeted correction: it improves fit exactly where
organic soils are known to degrade the retrieval, and does no detectable harm
anywhere else.

### 5.4 A regional exception: Hudson Bay

**Independent-check degradation, global vs. Hudson Bay Lowlands box (50–65°N,
100–75°W), DA − OL (K):**

| Region / config | DA − OL (K) |
|---|---|
| Global · baseline | −0.147 |
| Global · peat QC | −0.146 |
| Hudson Bay · baseline | **+0.124** |
| Hudson Bay · peat QC | **+0.073** |

*(Negative = assimilation improves the independent check, the global norm;
positive = assimilation degrades it.)*

Averaged over the Hudson Bay Lowlands, assimilating ASCAT *degrades* the
independent SMAP check relative to the open loop — the opposite sign from the
−0.15 K global improvement above. This holds under both the pre-QC baseline
(+0.124 K) and the peatland-QC configuration (+0.073 K): the new QC reduces
the degradation by roughly 41% but does not eliminate it. The residual
pattern sits on the southwest shore of the bay, the same organic-soil belt
the peatland QC targets, suggesting an additional, still-unresolved error
source in that environment — plausibly a residual water-fraction
misclassification or a seasonal freeze/thaw timing effect not captured by the
current QC — that the porosity-based peat screen only partially captures.

## 6. Discussion

Three findings, taken together, make the case that H SAF H121/H139 CDR is a
suitable replacement for the legacy EUMETSAT ASCAT product in this
assimilation system, and should be preferred for future reprocessing:

1. It is a genuinely different, denser retrieval, not a repackaging of the
   same information — a real gain in observational content, not merely
   resolution.
2. Once assimilated, it improves both an independent monitoring diagnostic
   and its own self-consistency fit, and does so slightly better than the
   legacy product it replaces, on both lines of evidence simultaneously.
3. The one known failure mode of ASCAT retrievals over organic soils can be
   addressed with a physically motivated, footprint-scale QC criterion whose
   mechanism has been independently confirmed at the pixel level, and which
   acts in a narrowly targeted way — large local impact, negligible global
   side effect.

The Hudson Bay result is the clearest open question raised by this work. That
assimilating ASCAT can locally *hurt* an independent diagnostic, in a
specific and geographically coherent region, while helping everywhere else,
is not something the current QC chain fully explains. Because the peatland QC
only partially resolves it, the remaining mechanism is evidently not simply
"more peat correction" — it warrants its own investigation, most plausibly
starting from water-fraction and freeze/thaw seasonality in that specific
landscape.

## 7. Status and next steps

- Extend the peatland footprint QC validation from the six representative
  boundary tiles used here to the full Canada/Alaska domain, to confirm the
  pixel-level agreement holds at scale.
- The QC chain has since been extended to cover an additional soil-moisture
  observation type beyond what is evaluated in this report; that extension
  has not yet been tested in a production experiment.
- Investigate the residual Hudson Bay degradation directly, rather than as a
  byproduct of the peatland QC evaluation — candidate next steps include
  isolating water-fraction QC and freeze/thaw-timing effects specifically in
  that region.
- A per-satellite (rather than pooled) observation-error scaling variant for
  H121 has been scoped but not yet run.

---

*Diagnostics use land-area-weighted spatial aggregation and
observation-count-weighted combination across satellites/channels; see
companion engineering log (`ascat_hsaf_dev_log.md`) for implementation and
provenance detail.*
