# CERISE Land DA Notes

Created: 2026-06-17

These are lightweight working notes from the public CERISE deliverables archived
under `projects/cerise/data/public_deliverables/`. The goal is to identify the
documents most useful for land data assimilation, reanalysis, and seasonal
initialization work. These deliverables are grey-literature technical reports,
not journal papers.

## Archive

- Source page: https://www.cerise-project.eu/deliverables
- Local manifest: `projects/cerise/data/public_deliverables/00_manifest/cerise_deliverables_manifest.csv`
- Local BibTeX: `projects/cerise/data/public_deliverables/00_manifest/cerise_deliverables.bib`
- Download date in manifest: 2026-06-17
- Public PDFs downloaded: 17
- Public deliverables without posted files: 6
- Confidential deliverables retained in manifest but not downloaded: 19

## GMAO Framing

This note should not be read as "go do SEKF." CERISE is ECMWF-centered, so many
implementation details are expressed in SEKF and IFS/ecLand language. For GMAO,
the useful question is broader:

What do the CERISE reports suggest we should try in GEOS-LDAS, land-only
reanalysis, and weakly coupled global reanalysis?

The two GMAO background papers for interpreting this are:

- Reichle et al. 2023 SMAP radiance weakly coupled land-atmosphere analysis:
  `/Users/amfox/Zotero/storage/CFNN9WYZ/Reichle et al. - 2023 - A weakly coupled land surface analysis with SMAP r.pdf`
- CYGNSS soil-moisture DA paper:
  `/Users/amfox/Documents/Publications/CYGNSS DA paper/Revisions/As Submitted/JHM-D-26-0035_R1.pdf`

Current GMAO baseline from those papers:

- GEOS-LDAS already has a mature ensemble land analysis around the Catchment
  land surface model.
- The system can assimilate SMAP brightness temperatures directly through a
  microwave radiative transfer operator.
- It can also assimilate multiple soil-moisture retrieval products, including
  ASCAT, SMAP/SMOS-style products, and CYGNSS retrievals.
- The CYGNSS work shows that even a moderate-quality additional sensor can help
  when it adds temporal/spatial complementarity and is handled carefully.
- The weakly coupled SMAP radiance work shows that land analysis increments can
  improve near-surface atmospheric analyses and forecasts, not just offline
  land states.
- Known GMAO limitations/opportunities include observation-error tuning,
  representativeness, active-sensor behavior in dry/arid regions, direct use of
  GNSS-R/scatterometer signals, and stronger land-atmosphere ensemble coupling.

## Land Reanalysis Questions

For land-only reanalysis, CERISE is useful because it treats land analysis as a
reanalysis product in its own right, not only as an atmospheric forecast support
tool. The key questions for GMAO are:

- Which additional observations should be used in a GEOS-LDAS land reanalysis,
  and should they be assimilated as retrievals, radiances/brightness
  temperatures, backscatter, reflectivity, or gridded pseudo-observations?
- How should observation errors and quality control vary with vegetation,
  aridity, snow, roughness, drydown stage, and known sensor/model mismatch?
- How should land-only reanalysis be evaluated: against in situ networks,
  independent satellite products, triple collocation, ERA5-Land, operational
  analyses, and downstream forecasts?
- How should land-only reanalysis connect to global atmospheric reanalysis:
  as independent land initial conditions, weakly coupled cycling, or a more
  integrated land-atmosphere analysis?
- How should time-varying vegetation, land cover, snow, soil temperature, and
  other boundary/state variables be represented so that the reanalysis is
  physically consistent over decades?

For global reanalysis, the most relevant CERISE idea is not any one filter. It
is the idea that land can be treated as an active part of the global analysis
problem, with land-only and coupled evaluation happening together.

## Things GMAO Should Consider Trying

1. Dynamic observation-error and QC experiments for ASCAT/CYGNSS/SMAP/SMOS.

   The CYGNSS paper already points to regime-dependent behavior, especially in
   arid and semi-arid regions. CERISE reinforces the need to handle
   representativeness and observation errors explicitly. A practical GMAO test
   would vary observation errors or thinning/blacklisting by aridity,
   vegetation, surface roughness, snow/freeze state, drydown phase, and sensor.

2. Direct active-microwave observation operators.

   CYGNSS is currently assimilated as a soil-moisture retrieval, and ASCAT is
   handled through retrieval/scaling pathways. A clear research direction is
   direct assimilation of GNSS-R reflectivity and scatterometer backscatter once
   robust operators exist. CERISE D1.4 is relevant here because it treats
   observation operators as first-class research objects, including ML
   approaches.

3. Benchmark ML microwave operators against existing GEOS operators.

   D1.4 is about ML-based microwave emissivity/brightness-temperature operators
   for low-frequency observations. For GMAO, this suggests a bounded experiment:
   compare an ML operator with the current SMAP radiance/tau-omega-style
   operator using GEOS-LDAS states as predictors, especially over snow,
   complex vegetation, and transitional wet/dry regimes.

4. Multi-sensor weakly coupled atmospheric-impact experiments.

   Reichle et al. show that SMAP Tb assimilation can improve T2m/q2m analyses
   and forecasts. The CYGNSS paper suggests complementary sensors can improve
   land states. A natural next experiment is whether multi-sensor land DA
   combinations produce larger or more robust atmospheric impacts in weakly
   coupled GEOS than SMAP Tb alone.

5. Member-wise land-atmosphere ensemble forcing.

   The SMAP radiance weakly coupled work notes that the land ensemble was
   constrained by ensemble-average atmospheric forcing in that implementation.
   A high-value GMAO question is whether member-wise atmospheric forcing of the
   land ensemble increases useful spread, improves increments, or strengthens
   land-atmosphere coupling in Hybrid-4DEnVar.

6. Footprint-aware representativeness tests.

   CERISE regional work is useful because it takes coarse satellite footprints
   versus finer model grids seriously. GMAO could test whether footprint-aware
   aggregation, static footprint models, or simple ML/GNN-style operators
   improve comparisons and assimilation for M36/EASE-grid land states and
   coarse microwave observations.

7. Vegetation and land-cover variability in reanalysis.

   D6.2 highlights time-varying LAI/LULC as a reanalysis and seasonal-forecast
   issue. For GMAO, this points to tests with time-varying vegetation forcing,
   vegetation-state initialization, or at least sensitivity experiments that
   separate soil-moisture DA gains from vegetation-boundary changes.

8. Reanalysis evaluation protocol.

   D6.1 is worth mining for a formal evaluation template. GMAO already has in
   situ validation, OmF diagnostics, and triple collocation pieces. CERISE adds
   a useful framing around consistency of land-only reanalysis, operational
   analysis, ERA5/ERA5-Land, seasonal initial conditions, and trends.

9. Soil/snow temperature and screen-level information.

   This is where the CERISE SEKF material matters, but the GMAO question is not
   "copy SEKF." The question is whether analysed near-surface atmospheric
   information, LST, or snow-sensitive observations can improve GEOS soil/snow
   thermal states in land-only or coupled reanalysis without degrading humidity
   and energy-budget consistency.

## Translation Table

- CERISE/ECMWF `SEKF`: not a target method for GMAO by itself; translate to
  GEOS-LDAS EnKF/control-vector/update design questions.
- Screen-level `T2m/RH2m` pseudo-observations: possible analogue is using
  analysed atmospheric/screen-level information or LST-like products as
  indirect constraints on land states.
- CMEM / ML emissivity operators: analogue is the GEOS SMAP radiance operator
  and possible future operators for SMOS, AMSR2, GNSS-R, and scatterometers.
- CARRA/SURFEX footprint work: analogue is representativeness between satellite
  footprints, EASE/M36 grids, Catchment tiles, and in situ validation scales.
- EDA-derived flow-dependent information: analogue is how GEOS atmospheric and
  land ensembles communicate uncertainty in weakly coupled or stronger coupled
  cycling.

## CERISE Documents to Mine

### D1.2: Unified global land DA system

Local PDF:
`projects/cerise/data/public_deliverables/WP1_Land_DA_Methodology/D1.2_Unified_ensemble-based_global_land_data_assimilation_system_and_documentati/CERISE_D1.2_2024_Dec-2024_Unified_ensemble-based_global_land_data_assimilation_system_and_documentation_v1.1.pdf`

Why it matters:

- This is the key ECMWF global LDAS implementation document.
- It describes the move toward a unified SEKF-based land DA system.
- It documents how screen-level T2m/RH2m analyses are used as
  pseudo-observations in the SEKF.
- It extends the SEKF control vector beyond soil moisture to include soil
  temperature layers and first-layer snow temperature.
- It discusses replacing 2D-OI screen-level pseudo-observations with 4D-Var
  generated pseudo-observations, but the reported results are mixed.
- It discusses using EDA information to build flow-dependent background error
  information for the SEKF.

Notable implementation details:

- The SEKF control vector is described as including `swvl1`, `swvl2`, `swvl3`,
  `stl1`, `stl2`, `stl3`, and `tsn`.
- The observation vector includes screen-level `T2m` and `RH2m` plus ASCAT and
  SMOS surface soil moisture terms.
- Soil/snow temperature SEKF development is framed as a replacement for the
  previous 1D-OI approach.
- The report says the new SEKF soil/snow temperature analysis improves short
  lead-time T2m scores relative to the 1D-OI control, especially over snow-free
  land.
- The 4D-Var replacement for 2D-OI pseudo-observations is not presented as
  operationally ready yet, partly because of weighting and resolution issues.

Follow-up questions:

- Which of these changes are in ERA6-Land prototypes versus only NWP cycles?
- How exactly are ASCAT/SMOS observation errors and QC handled in the current
  ECMWF SEKF branch?
- Does D4.1, currently confidential, contain the missing ERA6-Land prototype
  implementation details referenced by D6.1/D6.2?

### D2.1: Coupled global land-atmosphere assimilation infrastructure

Local PDF:
`projects/cerise/data/public_deliverables/WP2_Coupled_Surface_Atmosphere_DA/D2.1_Documentation_of_coupled_assimilation_infrastructure_and_methodology_and_pr/CERISE_D2.1_2024_Dec-2024_Documentation_of_coupled_assimilation_infrastructure_and_methodology_and_preliminary_assessment_v2.0.pdf`

Why it matters:

- This is the bridge between the WP1 unified LDAS and coupled global
  reanalysis/forecast DA.
- It explains why running the SEKF as a separate trajectory is expensive and
  how CERISE reduced this cost.
- It describes outer-loop coupling between atmospheric 4D-Var and land DA.
- It documents gridded 2D-OI T2m/RH2m analysis fields as SEKF inputs.

Notable implementation details:

- The report distinguishes weak coupling, one-way outer-loop coupling, and
  two-way land-atmosphere feedback within the same 12-hour DA window.
- It says the SEKF can be activated flexibly in outer loops.
- It describes land feedback to subsequent outer loops through updated first
  guess land fields.
- It notes concerns about robustness and time-critical-path impacts when
  merging SEKF too directly with the atmospheric 4D-Var task.

Follow-up questions:

- Which outer-loop configuration was actually selected for later prototypes?
- How does the land feedback alter atmospheric increments versus only forecast
  initialization?
- Which fields are written for SEKF first-guess departure calculations at each
  time step?

### D2.3: Coupled skin/soil temperature assimilation for global reanalysis

Local PDF:
`projects/cerise/data/public_deliverables/WP2_Coupled_Surface_Atmosphere_DA/D2.3_Documentation_on_coupled_skin_temperature_assimilation_for_coupled_reanalys/CERISE_D2.3_2025_Dec-2025_Documentation_on_coupled_skin_temperature_assimilation_for_coupled_reanalysis_v1.2.pdf`

Why it matters:

- This follows D2.1 and focuses on assimilating surface-sensitive temperature
  information through the coupled outer-loop framework.
- It documents an expanded SEKF observation vector.
- It describes using 4D-Var extended-control-variable analysis fields as
  gridded pseudo-observations for the SEKF.

Notable implementation details:

- The 4D-Var extended control variable framework is used to analyse first-layer
  soil temperature (`STL1`) from surface-sensitive observations.
- The resulting gridded `XCV-STL1` fields can be written at synoptic times and
  assimilated into the SEKF.
- Desroziers-style diagnostics are used to estimate observation error for
  these pseudo-observation fields.
- Initial experiments with deeper soil-temperature updating did not show clear
  benefit, so the tested setup focuses on updating `STL1`.
- The report mentions large speedups from avoiding repeated SEKF nonlinear
  trajectories, which is important for operational feasibility.

Follow-up questions:

- Is `XCV-STL1` planned as a pathway for skin temperature/LST observations?
- How are the pseudo-observation errors varied by outer loop or observation
  source?
- Does the RH2m degradation noted in some experiments constrain the usable
  coupling strength?

### D1.4: ML observation operators for microwave observations

Local PDF:
`projects/cerise/data/public_deliverables/WP1_Land_DA_Methodology/D1.4_Report_on_observation_operator_methodology_ready_for_implementation_in_coup/CERISE_D1.4_2025_Dec-2025_Report_on_observation_operator_methodology_ready_for_implementation_in_coupled_global_and_regio_v1.0.pdf`

Why it matters:

- This is the observation-operator deliverable for WP1.
- It is about ML-based forward operators for low-frequency microwave
  observations over land, rather than about SEKF mechanics directly.
- It is relevant to direct assimilation of microwave radiances/emissivities
  from sensors such as SMOS, SMAP, and AMSR2.
- It is potentially relevant to replacing, augmenting, or benchmarking older
  physics-based microwave operators such as CMEM.
- It is also relevant to the footprint mismatch problem: coarse satellite
  microwave footprints versus finer regional land-model grids.

Main structure:

- Section 3: ML-based observation operator for a global coupled assimilation
  system.
- Section 4: ML-based observation operator for a regional coupled assimilation
  system.
- Section 5: ML-based observation operator for hydrological applications.
- Annex I: code pointers for the global coupled-assimilation observation
  operator.

Global operator notes:

- The target is low-frequency passive microwave information from about
  1.4-36 GHz.
- The report uses satellite-derived microwave emissivities from AMSR2, SMAP,
  and SMOS.
- SMAP provides L-band brightness temperatures at about 1.4 GHz; SMOS provides
  multi-angular/full-polarization L-band observations; AMSR2 contributes higher
  microwave frequencies.
- The workflow derives emissivities from observations, combines them with
  geophysical predictors, and trains neural-network emissivity forward
  operators.
- The global component is framed as a way to link globally available
  geophysical properties to robust satellite-derived emissivities.
- Initial results reportedly improve performance over snow-covered regions for
  the global emissivity operator.

Regional operator notes:

- The regional component is motivated by high-resolution limited-area models
  and coarse microwave observation footprints.
- SURFEX simulations and observations over a Scandinavian/CARRA-like domain are
  used to create training data.
- CARRA forcing is mentioned as input to the SURFEX/pysurfex workflow.
- Candidate ML methods include multilayer perceptrons, XGBoost, footprint CNNs,
  and graph neural networks.
- The graph neural network approach is meant to account for sub-footprint
  heterogeneity in the land state.
- The dynamic graph neural network appears to perform well but is expensive to
  train; the static graph neural network is presented as a more practical
  promising option.
- The report explicitly discusses use in both offline and coupled regional
  reanalysis systems, including CARRA3-Pv1-style contexts.

Hydrological operator notes:

- The hydrological section explores an ML-based microwave observation operator
  for HYPE-style hydrological applications.
- AMSR2 brightness temperatures are mentioned as a target observation type.
- This may be less immediately relevant to GEOS-LDAS, but useful as a reminder
  that the observation-operator work is broader than NWP/reanalysis systems.

Why this is different from D1.2/D2.1:

- D1.2 and D2.1 are mostly about how land DA and coupled DA are structured.
- D1.4 is about the forward operator needed to use satellite microwave
  observations more directly.
- For ASCAT H121 soil-moisture retrieval assimilation this is only indirectly
  relevant, because ASCAT retrievals are already geophysical soil moisture
  products.
- For SMAP/SMOS/AMSR2 brightness-temperature or emissivity assimilation, D1.4
  is directly relevant.

Follow-up questions:

- Does the Annex contain code URLs or implementation details that should be
  archived separately?
- How does the global ML emissivity operator compare against CMEM by frequency,
  snow condition, and land-cover class?
- Which predictors are used in the final global NN operator, and are they
  available from GEOS-LDAS-style outputs?
- Can the regional static-GNN idea inform footprint-aware comparisons between
  satellite products and M36/EASE-grid model fields?
- Is there a trained model artifact available, or only methodology in the
  deliverable?

### D1.3: Unified regional land DA system

Local PDF:
`projects/cerise/data/public_deliverables/WP1_Land_DA_Methodology/D1.3_Unified_ensemble-based_regional_land_data_assimilation_system_and_documenta/CERISE_D1.3_2024_Dec-2024_Unified_ensemble-based_regional_land_data_assimilation_system_and_documentation_v1.1.pdf`

Why it matters:

- This is the regional counterpart to D1.2.
- It covers ensemble-based regional LDAS work in HARMONIE-AROME/SURFEX-like
  settings.
- It has useful details on EnSRKF/LETKF approaches, soil/snow control vectors,
  screen-level observations, and snow analysis.

Notable implementation details:

- EnSRKF regional LDAS experiments compare against sEKF and show different
  behavior in propagating screen-level observation information into deeper soil.
- The LETKF prototype is described as flexible with respect to observation and
  control variables and suitable for future satellite-observation footprint
  operators.
- The LETKF snow work uses localization and explicit handling for bounded
  variables and absent-snow states.
- Synthetic experiments are used to separate filter behavior from observation
  representativeness problems.
- Snow DA performs best in flatter terrain; mountainous regions remain
  challenging because ensemble perturbation structures and representativeness
  are harder.

Follow-up questions:

- Is the referenced `sfcpert` package useful to inspect directly?
- Which regional control variables map cleanly onto GEOS-LDAS snow/soil
  variables?
- How do regional systems handle the "no snow in all members" case in practice?

### D2.4: Coupled skin/soil temperature assimilation for regional reanalyses

Local PDF:
`projects/cerise/data/public_deliverables/WP2_Coupled_Surface_Atmosphere_DA/D2.4_Documentation_on_coupled_skin_temperature_assimilation_for_regional_reanaly/CERISE_D2.4_2025_Dec-2025_Documentation_on_coupled_skin_temperature_assimilation_for_regional_reanalyses_v1.0.0.pdf`

Why it matters:

- This document is useful for regional coupled DA architecture.
- It describes a Harmonie-Arome outer-loop 4D-Var setup with a land surface
  analysis block inside each outer loop.
- It gives a regional counterpart to the global outer-loop LDAS work.

Notable implementation details:

- The regional SEKF assimilates screen-level temperature, relative humidity,
  and snow-depth fields obtained through optimal interpolation from
  quality-controlled observations.
- The SEKF propagates screen-level increments into soil levels using Jacobians
  computed by perturbed offline SURFEX/ISBA-DIF runs.
- Coupling from atmosphere back to soil can occur through rediagnosed
  screen-level fields and updated observation-operator Jacobians in later
  outer loops.

Follow-up questions:

- Which parts are CARRA2/CARRA3-specific versus generally reusable?
- How much tuning is system-specific for Jacobian bounds and observation
  errors?

### D6.1: Evaluation protocol for demonstrator quality

Local PDF:
`projects/cerise/data/public_deliverables/WP6_Demonstrator_Evaluation/D6.1_Report_providing_a_protocol_for_assessing_the_improvement_in_the_quality_in/CERISE_D6.1_2025_Jun-2025_Report_providing_a_protocol_for_assessing_the_improvement_in_the_quality_in_the_demonstrators_v1.0.pdf`

Why it matters:

- This is not a DA-method document, but it is useful for evaluation framing.
- It gives methods for comparing LDAS/reanalysis prototypes against ERA5,
  ERA5-Land, open-loop land reanalysis, and operational analysis.
- It connects land initialization consistency to seasonal forecast behavior.

Notable implementation details:

- One proposed diagnostic compares prototype offline LDAS reanalysis and ERA5
  or ERA5-Land against ECMWF operational analysis as a reference.
- It explicitly discusses consistency between hindcast and forecast land initial
  conditions for soil moisture and snow depth.
- It includes trend diagnostics for land-surface variables, including ERA5-Land
  soil-water trends.

Follow-up questions:

- Can the ERA5-Land comparison workflow in `projects/era5_land` reuse any of
  these diagnostic definitions?
- Which operational-analysis fields are needed to reproduce the consistency
  diagnostics?

### D6.2: Time-varying vegetation in reanalysis and seasonal forecasts

Local PDF:
`projects/cerise/data/public_deliverables/WP6_Demonstrator_Evaluation/D6.2_Report_providing_feedback_on_the_impact_of_time_varying_vegetation_in_reana/CERISE_D6.2_2025_Dec-2025_Report_providing_feedback_on_the_impact_of_time_varying_vegetation_in_reanalysis_on_seasonal_fo_v3.0.pdf`

Why it matters:

- This report links offline LDAS initial conditions, vegetation forcing, and
  seasonal forecast sensitivity.
- It documents ECMWF offline LDAS experiments with ecLand forced by ERA5
  meteorology and land-only observation assimilation.

Notable implementation details:

- `LDAS_CTL` uses a 2010-2019 seasonal LAI climatology and static LULC.
- `LDAS_VAR` uses seasonally and interannually varying LAI plus interannually
  varying LULC.
- ECMWF seasonal reforecasts test CONTROL, LAI initial-condition changes, and
  persistence of LAI anomalies into the forecast.

Follow-up questions:

- Are the CONFESS LAI/LULC boundary datasets publicly accessible and useful for
  GEOS-LDAS comparisons?
- Is vegetation variability treated as a DA state, a boundary forcing, or both
  in later CERISE work?

## Cross-Cutting Themes

### Screen-level observations as land DA information

The ECMWF global and regional documents repeatedly use screen-level 2 m
temperature and humidity as indirect information for soil and snow states.
The current/legacy structure is often:

- analyse screen-level fields with 2D-OI or atmospheric DA,
- use gridded analysis fields or increments as pseudo-observations,
- map these to soil moisture/temperature/snow control variables through the
  SEKF or regional ensemble filters.

This matters for interpreting papers that simply say "screen-level observations
are assimilated into the land analysis." The implementation detail is that the
land DA may not assimilate raw station observations directly; it may assimilate
gridded analysed fields produced by a separate screen-level analysis.

For GMAO, this is mostly interesting as a possible land-reanalysis information
pathway: can analysed near-surface atmospheric states help constrain land
states without creating circular evaluation problems or atmospheric-analysis
dependency issues?

### ECMWF implementation path, not a GMAO prescription

The CERISE reports show a progression:

- operational/static SEKF for soil moisture,
- extended SEKF control vectors for soil and snow temperature,
- use of EDA-derived Jacobians and potentially flow-dependent background
  errors,
- integration with 4D-Var outer loops to avoid repeated nonlinear trajectories,
- possible replacement or augmentation of 2D-OI pseudo-observations with 4D-Var
  generated fields.

The GMAO translation is to ask which state variables, observations, ensemble
couplings, and error models should be added to the existing GEOS-LDAS/EnKF and
weakly coupled GEOS framework.

### Snow remains awkward

Across global and regional documents, snow creates special cases:

- snow/no-snow discontinuities and bounded variables,
- snowpack insulation that weakens screen-level sensitivity to soil states,
- representativeness problems for point snow-depth observations,
- weaker performance in mountainous terrain,
- need for special initialization/post-processing when ensemble members lack
  snow.

### Regional ensemble DA is more explicit about flow-dependent structure

D1.3 is useful because it discusses ensemble-space filters, localization,
synthetic experiments, spread-skill behavior, and observation
representativeness. This is a good complement to the ECMWF global SEKF
documents, which are more tightly tied to IFS infrastructure.

### Land-only reanalysis and coupled reanalysis need different scorecards

CERISE separates land-only LDAS demonstrators, coupled reanalysis
demonstrators, and seasonal forecast impacts. That is a useful discipline for
GMAO. A land-only reanalysis can improve soil moisture or snow metrics without
necessarily improving atmospheric forecasts, and a coupled global reanalysis
can improve near-surface weather while still needing independent land-state
validation.

## Search Terms That Paid Off

Useful terms for PDF/text searches:

- `SEKF`
- `2D-OI`
- `screen-level`
- `pseudo-observations`
- `T2m`
- `RH2m`
- `ASCAT`
- `SMOS`
- `XCV-STL1`
- `outer loop`
- `EDA`
- `flow-dependent`
- `LDAS`
- `ecLand`
- `SURFEX`
- `LETKF`
- `EnSRKF`
- `snow depth`
- `ERA5-Land`
- `land reanalysis`
- `global reanalysis`
- `ERA6`
- `observation operator`
- `emissivity`
- `brightness temperature`
- `reflectivity`
- `backscatter`
- `vegetation`

Example local search workflow:

```bash
mkdir -p /private/tmp/cerise_pdf_text
find projects/cerise/data/public_deliverables -name '*.pdf' -print0 |
  while IFS= read -r -d '' f; do
    base=$(basename "$f" .pdf)
    pdftotext -layout "$f" "/private/tmp/cerise_pdf_text/${base}.txt"
  done

rg -n -i "SEKF|2D.?OI|screen.?level|pseudo.?observ|ASCAT|SMOS|ERA5-Land" \
  /private/tmp/cerise_pdf_text
```

## Next Reading Pass

1. Build a one-page GMAO opportunity matrix: observation/source, current GEOS
   capability, CERISE idea, possible experiment, and evaluation metric.
2. Mine D1.4 first for observation-operator ideas that could matter for SMAP,
   SMOS, CYGNSS, ASCAT, and future direct microwave DA.
3. Mine D6.1 for reusable land-only and global-reanalysis evaluation metrics,
   especially ERA5-Land, operational-analysis, trend, and seasonal-initial-state
   consistency diagnostics.
4. Read D6.2 for vegetation and land-cover variability as a reanalysis design
   issue, not just a seasonal forecast issue.
5. Read D2.1/D2.3 only after that, focusing on how land-only analysis feeds
   coupled global reanalysis and atmospheric forecast impact.
6. Use D1.2/D1.3 as implementation background, translating SEKF/SURFEX/ecLand
   details into GEOS-LDAS/EnKF/Catchment questions.
7. Check whether public future reports D2.5, D3.5, D5.6, D6.4, D6.5, and D8.7
   become available later and rerun the downloader.
