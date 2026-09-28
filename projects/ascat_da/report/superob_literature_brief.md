# Brief for a literature search: superobbing vs thinning dense ASCAT soil moisture in a land EnKF

## 1. System

- **Assimilation system:** NASA GMAO GEOS-LDAS (GEOS Land Data Assimilation System). The model is the Catchment land surface model, run on land tiles of a ~25 km cubed-sphere grid (C360, "CF0360"). The analysis is a 24-member ensemble Kalman filter with 3-hourly updates.
- **What the analysis updates:** each tile's soil moisture prognostics (catchment deficit, root-zone excess, surface excess), from observations within a local neighbourhood.
- **Observation errors:** assumed spatially correlated, with a Gaussian-type correlation length of 0.3125° (~35 km) in x and y. The per-observation error std is fixed and is ~0.02 m3/m3 after rescaling.
- **Observations:** EUMETSAT H SAF ASCAT surface soil moisture climate data record (H121), from Metop-A, -B and -C. These are swath data sampled on a ~12.5 km grid (a Fibonacci grid), so neighbouring footprints overlap strongly. The spatial resolution is coarser than the sampling (~25 km or more).
- **QC:** open water, bad model or sigma0 flags, noisy backscatter, wetland/topography fraction, sub-surface scattering, and low SSM sensitivity.
- **Bias correction:** before assimilation, obs are rescaled to the model climatology by z-score (mean/std) matching. The rescaling statistics are stored per 0.25° lat/lon cell and per pentad, computed from an open-loop (no-DA) run's obs-minus-forecast archive.
- **Independent validation:** forecast O-F statistics of SMAP and SMOS L-band brightness temperatures and CYGNSS soil moisture, all monitored but not assimilated, compared with an open-loop reference run. Filter diagnostics are also used: O-A vs O-F, normalized innovations, and Desroziers et al. (2005) ratios.

## 2. The problem

Assimilating the raw ASCAT obs at native density was pathological:
- Each obs was mapped to the model tile it falls in, and multiple obs in a tile were averaged.
- The analysis moved *away* from the obs (O-A std > O-F std).
- Increments ran away, producing negative soil moisture and clamping.
- Skill against independent SMAP/SMOS Tb and CYGNSS got progressively worse than the open loop.

Our interpretation: the obs are much denser than the assumed error correlation length, and the real obs errors are correlated more strongly or over longer distances than modelled. Neighbouring obs then carry largely redundant information. The EnKF overweights them, and the local analysis becomes badly conditioned. This is a hypothesis consistent with the diagnostics, not a proof.

## 3. What we tested

All tests used the same restarts, perturbations, error settings and 5-day period (May 2019). Only one thing changed per experiment.

| Method | Grid spacing | Result |
|---|---|---|
| Raw (tile-averaged) | native ~12.5 km | pathological |
| Super-obs | 0.125° (~14 km) | pathological, though less severe |
| Super-obs | 0.25° (~28 km) | healthy |
| Super-obs | 0.5° (~55 km) | healthy |
| Thinning: one raw obs per 0.25° cell, a fixed location chosen before QC | 0.25° | healthy |

- **Where it breaks:** with a 35 km error correlation length, stability breaks between 0.125° and 0.25° obs spacing. That's roughly obs spacing of 0.4–0.8 of the correlation length. This is our inference from three grid sizes.
- **Skill of the healthy runs:** all slightly beat the open loop on independent SMAP/SMOS Tb (by ~0.05–0.1 K in O-F std) and on CYGNSS. The differences between super-obs and thinning are at noise level.
- **Desroziers obs-error ratio:** ~1.1 for super-obs and ~0.95 for thinning.
- **Desroziers background ratio (HBH):** ~0.4, which suggests the ensemble's forecast error variance is too large.

**How our super-obs are built** (in the obs reader, before assimilation):
- **Averaging:** all QC-passed raw obs in a fixed 0.25° lat/lon cell within the 3-hour window are averaged, in value and in time.
- **Error:** the super-ob error std is NOT reduced by the number of obs averaged, because the errors are assumed correlated.
- **Model value:** an area-weighted average of the model tiles within a circular footprint of radius 0.141°, the equal-area radius of a 0.25° cell.
- **Tile assignment:** each super-ob goes to one tile. If two super-obs land on the same tile, the closer one is kept.
- **Location:** we tested placing each super-ob at (a) the cell centre and (b) the mean lat/lon of its raw obs.

## 4. Current findings and preference

- **Super-obs over thinning.** Skill is essentially tied. Super-obs use all the data rather than discarding ~60% of it, and they're done inside the reader with no preprocessed files. They also suit an operational reanalysis and the near-real-time extension of the record.
- **Mean location over cell centre.**
  - It's more physically consistent: the super-ob represents where its obs actually are, which matters in partly covered cells at swath edges, coastlines, or after QC.
  - It's consistent with using the mean time.
  - In our first test it was neutral on skill.
- **An unexplained ~7% penalty.** At the very first analysis (identical model state in all runs), super-ob O-F std is ~7% higher than raw or thinned obs on the same tiles.
  - The rescaling climatology is not the cause. Rebuilding it from super-ob-aggregated statistics recovered only ~1/5 of the gap.
  - Placement (cell centre vs mean location) is not the cause.
  - Remaining candidates: the super-ob's model footprint (area-averaged model vs averaged obs), and a representativeness mismatch between the area the averaged obs cover and the model tile it's compared to.
- **Scaling climatology.** The z-score climatology should be built from an open loop that monitors the obs in the *same* form as assimilated (super-obs at their mean positions). Otherwise cells are empty or the statistics are inconsistent.
- **Grid geometry concern.** A regular lat/lon super-ob grid gets denser toward the poles in terms of obs per unit area. At 60°N a 0.25° cell is ~0.3 of an EASE-v2 36 km (M36) cell. An equal-area super-ob grid, such as EASE-v2 36 km to match SMAP, may be more appropriate.

## 5. Questions for the literature search

1. **Superobbing vs thinning:** examples and comparisons in land-surface or soil moisture DA, especially ASCAT, SMOS, SMAP or AMSR-E in EnKF or EKF systems (e.g., ECMWF's land surface analysis, GMAO's SMAP L4, NASA LIS, CMEMS/Copernicus land systems). What grid or spacing did they choose, and why?
2. **Spacing vs error correlation:** theory or rules of thumb for the obs spacing, or super-ob size, relative to the observation error correlation length. Includes work on optimal thinning distances, e.g., for satellite radiances and atmospheric motion vectors in NWP. Is spacing ~0.5–1 correlation length a known stability or efficiency threshold?
3. **Super-ob location and time:** precedent for mean location vs grid-cell centre, and for mean time.
4. **Super-ob error:** how super-ob error variance is set when the component errors are correlated (not 1/N), and whether representativeness error should be added.
5. **Model side of a super-ob:** how the model equivalent is computed, i.e. the footprint or averaging area of the operator, and whether it should match the actual spread of the averaged obs.
6. **Error correlations:** estimated horizontal observation-error correlation lengths for ASCAT (and other scatterometer or radiometer) soil moisture, and their effect on DA.
7. **Grid choice:** equal-area vs regular lat/lon grids for superobbing or thinning.
8. **Bias correction at super-ob scale:** CDF or z-score matching done at the super-ob scale vs the native obs scale, and building the matching climatology from an open loop that monitors the super-obs.
9. **Background error:** reports of ensemble forecast error variance being overestimated relative to Desroziers diagnostics in land EnKFs.

Suggested seed topics: "superobservations", "observation thinning", "correlated observation errors", "Desroziers diagnostics", and "ASCAT soil moisture assimilation EnKF". Verify any specific references the search returns; none are asserted here.
