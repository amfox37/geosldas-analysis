# Manuscript Rebuild Plan

> This document records the structural and editorial decisions for the August 2026 rebuild of the M21C land-reanalysis manuscript.
>
> The scientific conclusions themselves are maintained in `current_scientific_state.md`. This file governs **how those conclusions should be turned into the manuscript**.

## 1. Rebuild principle

The manuscript is being **rebuilt rather than incrementally revised**.

The May 2026 draft remains useful source material for:

- background
- references
- established Methods wording
- dataset descriptions
- previously accepted factual material

but its structure and emphasis are no longer authoritative.

New analyses should replace weaker or redundant material rather than simply being appended to the old paper.

The final manuscript should ideally contain **more science but less prose** than the May draft.

## 2. Central manuscript argument

The paper should answer:

> **How does an evolving land observing system change the strength, spatial structure, physical consequences, and temporal consistency of a long-term multi-sensor land reanalysis?**

The manuscript should develop three linked ideas:

1. **Evolving observational influence**
2. **Physical propagation of assimilation corrections**
3. **Temporal consistency of the resulting land record**

The narrative progression should be:

> **observing system evolves → observational influence changes → assimilation corrections propagate through land processes → the correction regime becomes nonstationary → broad long-term state behavior is nevertheless generally preserved**

The paper should not read as a catalogue of sensors or validation scores.

## 3. Relationship to previous papers

The manuscript must avoid unnecessary overlap with the existing soil-moisture assimilation papers.

Use:

- **Fox et al. (2025)** for the established ASCAT + SMAP GEOS LDAS methodology
- **Fox et al. (2026)** for the current multi-sensor soil-moisture configuration, including SMOS and CYGNSS

These papers should carry much of the detailed soil-moisture DA Methods burden.

Do not re-explain in detail:

- the complete EnKF formulation
- all localization mechanics
- all ASCAT rescaling details
- the full SMAP/SMOS radiative-transfer implementation
- complete CYGNSS QC and transformation chains
- generic monitored-only OmF mechanics
- generic multi-sensor interaction principles

The present paper still needs enough information to understand:

- the June 2000–May 2024 experiment
- OL versus DA
- model/grid/ensemble
- forcing
- which observations were assimilated and when
- P1–P9
- MODIS SCF assimilation
- evaluation datasets
- snow-process analyses
- trend and transition analyses

A useful Methods bridge is:

> GEOS LDAS and its multi-sensor soil-moisture assimilation configuration are described in detail by Fox et al. (2025, 2026). Here we summarize only the elements needed to interpret the 2000–2024 experiment and describe aspects specific to the evolving observing system, MODIS snow-cover assimilation, and the process and temporal-consistency analyses introduced here.

## 4. Target manuscript structure

### 1. Introduction

The Introduction should be rewritten around the manuscript problem rather than inherited from the May draft.

It should establish:

1. land reanalyses rely increasingly on heterogeneous satellite observing systems;
2. those observing systems change substantially over multi-decadal records;
3. data assimilation therefore introduces time-varying observational constraint;
4. this raises two distinct concerns:
   - how observational influence evolves;
   - whether changes in the observing system compromise temporal consistency;
5. assimilation corrections can also propagate beyond directly observed quantities through land-model dynamics;
6. the current experiment provides a 24-year test of these issues.

The Introduction should end with **three explicit study questions**, corresponding to:

- evolution of observational influence;
- propagation into land hydrology;
- long-term trend and transition behavior.

Avoid a long sensor-by-sensor literature review.

### 2. Methods and Data

#### 2.1 GEOS LDAS and Catchment model

Keep concise.

Include:

- CLSM
- EASEv2 M36 grid
- 24-member ensemble
- relevant land-state framework
- experiment forcing

Offload established implementation details to prior references.

#### 2.2 Assimilated observing system

Describe:

- MODIS SCF
- ASCAT soil moisture
- SMOS Tb
- SMAP Tb
- CYGNSS soil moisture

The purpose is to explain the evolving observation record, not to reproduce individual instrument-method papers.

Figure 1 / P1–P9 is the authoritative chronology.

#### 2.3 Land data assimilation

Summarize the EnKF and observation treatment only to the level needed for this experiment.

Cite Fox et al. (2025, 2026) heavily.

MODIS snow assimilation requires enough additional detail because it is central to the new hydrologic-process result.

#### 2.4 Experiments and observing-system periods

Define:

- OL
- DA
- June 2000–May 2024
- same forcing/model setup
- DA as the principal difference
- P1–P9

Do not reintroduce the old competing period taxonomy.

Broader validation periods may be described only as explicit aggregations where necessary.

#### 2.5 Evaluation datasets and metrics

Cover only what is required for:

- ISMN soil moisture
- IMS snow cover
- SNOTEL SWE
- GHCN snow depth
- ERA5-Land comparisons

Avoid excessive metric tutorial material if already standard or documented in Fox et al. (2025, 2026).

State clearly that ERA5-Land is a consistency comparison rather than independent truth.

#### 2.6 Snow-DA hydrologic pathway and water-budget analysis

This is a new and important Methods section.

It should explain:

- why 2001–2006 is used for clean snow attribution;
- construction of signed snow-DA water input;
- water-year/common-window accounting;
- runoff partition:
  - surface runoff
  - baseflow
- ET
- end-of-year storage change
- residual
- use of actual TWLAND endpoint change rather than WCHANGELAND
- treatment of snowmelt as a pathway rather than a terminal sink
- spatial-block uncertainty
- controlled-regression corroboration
- strict non-overlapping-window sensitivity

Keep detailed robustness variants in the Supplement.

#### 2.7 Long-term trends and observing-system transitions

Explain separately:

### Long-term trends
- trend-preserving removal of calendar-month climatology
- Theil–Sen slopes
- paired OL, DA, and DA−OL treatment
- FDR control
- relevant domains/masks

### Structural transitions
- P1–P9 known boundaries
- interrupted time-series level changes
- fitted-AR(1) bootstrap intervals
- boundary-family FDR
- independent changepoint analysis as corroboration

Make clear that:

> a long-term trend and a discrete observing-system transition are different estimands.

## 5. Results architecture

### 3.1 Evolution of the assimilated observing system

Figures 1–3.

Purpose:

> establish how the observing system itself changes through time.

Keep descriptive prose short.

Do not turn this into an instrument catalogue.

### 3.2 Evolution of observational influence

Figures 4–6.

Purpose:

> show that the strength and spatial structure of observational constraint evolve with the observing system.

The principal transition should emerge naturally as SMAP is introduced.

Avoid repeating detailed sensor-comparison conclusions from Fox et al. (2025, 2026).

### 3.3 Independent soil-moisture evaluation

Figure 7.

Purpose:

> demonstrate that the evolving observational influence has measurable consequences for analyzed soil moisture.

Emphasize:

- stronger surface response
- weaker but meaningful root-zone response
- changes across observing-system eras

Keep station-by-station or metric-by-metric description out of the main narrative.

### 3.4 Snow-cover, SWE, and snow-depth evaluation

Figures 8–10.

Purpose:

> distinguish strong direct constraint on snow occurrence/SCF from more model-mediated SWE and depth responses.

This section should set up the process question:

> modest direct SWE/depth validation changes do not imply weak hydrological consequences.

That transition is important for Section 3.6.

### 3.5 Comparison with ERA5-Land

Figures 11–13.

Keep this section concise.

Purpose:

> provide an additional large-scale consistency view of soil-moisture and snow behavior.

Do not present ERA5-Land as truth.

Do not let this section interrupt the main narrative.

### 3.6 Hydrologic consequences of snow-cover assimilation

Figures 14–15.

This should be one of the manuscript's strongest sections.

The order should be:

1. **Common-window water budget — headline result**
2. runoff/ET/storage partition
3. independent controlled-regression corroboration
4. monthly pathway / timing
5. RZMC response
6. concise secondary note on later SM-DA counteraction
7. optional short hydro-energy consistency result

The main quantitative message is:

> approximately 64.3% of snow water added by assimilation becomes total runoff.

Do not lead with regression slopes.

Do not headline snow-removal, residence-time, or delayed-infiltration analyses.

### 3.7 Long-term trends and observing-system transitions

Figures 16–17.

This section should make a clear distinction between:

#### Long-term behavior
- broad trends are generally preserved
- SCF decline is essentially unchanged
- snow mass/depth show little DA-induced secular trend
- RZMC trends are regionally modified

and:

#### Structural observing-system transitions
- April 2015 is the dominant transition
- RZMC DA−OL shifts positively
- SM correction activity increases sharply
- multiple independent diagnostics break at the same time
- OL/DA control-state series do not show the corresponding break

The section should end with the central temporal-consistency interpretation:

> **The evolving observing system produces clear structural changes in analysis influence without producing widespread spurious long-term drift in the underlying land-state record.**

## 6. Discussion architecture

The Discussion should be largely rewritten from scratch.

Do not simply restate the Results.

### 4.1 Evolving observational constraint in long-term land analysis

Discuss:

- why changing observing systems matter for reanalysis
- how observational influence grows/changes through the record
- why SMAP is the dominant transition here
- why a nominally fixed DA framework does not imply stationary observational influence
- limits on attributing changes uniquely to intrinsic sensor information content

### 4.2 From direct observational constraint to model-mediated hydrologic response

Use the snow result to revise the older “directness-of-constraint” framing.

Key point:

> the strongest immediate effects remain near observed quantities, but mass/state corrections can propagate strongly through modeled hydrology.

Discuss:

- snowmelt as transit
- runoff-dominated partition
- ET and storage
- RZMC response
- hydro-energy propagation
- implications for interpreting DA effects on unobserved quantities

### 4.3 Observing-system transitions and long-term trend fidelity

This is the conceptual center of the temporal-consistency discussion.

Make the distinction:

> a reanalysis can contain observing-system-induced level/activity transitions while still preserving broad long-term state trends.

Discuss:

- April 2015 transition
- regional RZMC trend modification
- near-null DA-induced SCF trend
- null precipitation control
- why trend fidelity cannot be assessed from a single global mean alone

Avoid claiming the record is perfectly homogeneous.

### 4.4 Limitations and implications for future multi-sensor reanalysis

Include concise limitations:

- evolving observing systems are inherently nonstationary
- DA settings were not optimized independently for each era
- some validation datasets have temporal/spatial limitations
- ERA5-Land is not independent truth
- snow-process attribution is intentionally restricted to the pre-SM-DA period
- residual water-budget error reflects available temporal sampling
- causality at known observing-system boundaries should be expressed carefully

Potential forward-looking implication:

> future long-term land reanalyses should document and diagnose changes in observational influence explicitly rather than assuming a fixed assimilation algorithm guarantees temporal homogeneity.

A brief connection to forthcoming MERRA-21C-Land may be included in the final paragraph, but it should not become a promotional section.

## 7. Conclusions

The Conclusions should be short and synthetic.

They should not reproduce every validation metric.

The principal conclusions should be approximately:

1. The observational constraint in the 24-year GEOS LDAS record changes substantially as the satellite observing system evolves.
2. The introduction of SMAP in April 2015 marks the clearest system-level transition in soil-water analysis behavior.
3. Assimilation impacts are not confined to directly observed variables: MODIS SCF assimilation adds/removes water that propagates through modeled hydrology.
4. For snow-addition cases, most added snow water is subsequently discharged as runoff, with smaller ET and storage components.
5. Despite clear observing-system transitions and regional modification of RZMC trends, broad long-term behavior—including the SCF decline—is generally preserved.
6. Long-term land-reanalysis evaluation therefore requires consideration of both **observational influence** and **temporal consistency**.

Do not oversell perfect homogeneity or unique causal attribution.

## 8. Main figure sequence

The main-text figure sequence is fixed as:

1. Observing-system timeline
2. Mean assimilated observations per day
3. Monthly assimilated-observation counts
4. Full-period OmF improvement
5. Monthly evolution of OmF improvement
6. Period/sensor OmF maps
7. ISMN soil-moisture validation
8. IMS snow-cover evaluation
9. SNOTEL SWE evaluation
10. GHCN snow-depth evaluation
11. ERA5-Land soil-moisture bars
12. ERA5-Land soil-moisture maps
13. ERA5-Land snow comparison
14. Snow-DA water budget
15. Monthly snow-DA hydrologic pathway
16. Long-term RZMC and SCF trends
17. Observing-system transitions

Do not change this sequence without an explicit manuscript-level decision.

## 9. Main text versus Supplement

The main text should contain the result, not every robustness test.

Move or keep in the Supplement:

- detailed M0–M4 regression hierarchy
- alternative trimming
- alternative spatial-block scales
- non-overlap variants
- delayed-infiltration failure
- snow-removal budget
- individual annual budget details
- water-budget technical audits
- RZMC persistence / residence analysis
- detailed signed snow–SM interaction
- precipitation trend maps
- SFMC trend maps
- snow-mass/depth trend maps
- breakpoint-boundary agreement matrix
- full changepoint sensitivity results

A robustness result should enter the main text only if it changes interpretation or is needed to defend the main claim.

## 10. Writing style

Use a restrained GEOS/GMAO scientific style.

Prefer:

- direct statements
- short paragraphs
- quantitative results followed by interpretation
- clear distinction between evidence and inference
- cautious causal wording
- explicit limitations where relevant

Avoid:

- hype
- repeated claims that results are “important” or “novel”
- long signposting paragraphs
- excessive method recapitulation in Results
- figure-by-figure narration
- repeating the same numerical result in Results and Discussion
- describing internal physical consistency as independent validation

Prefer formulations such as:

- “coincides with”
- “is consistent with”
- “supports the interpretation”
- “indicates”
- “suggests”

where unique causation is not established.

## 11. Citation strategy

For now, use normal human-readable author–year citations.

Examples:

- `Fox et al. (2025)`
- `Fox et al. (2026)`
- `(Reichle et al., 2017)`

Do not yet introduce:

- Pandoc citation keys
- CSL processing
- BibTeX requirements
- Zotero field placeholders

The Word/Zotero conversion strategy will be decided after the scientific manuscript stabilizes.

## 12. Drafting workflow

Work section by section.

Recommended sequence:

1. finalize scientific-state note
2. finalize rebuild plan
3. finalize coauthor-feedback ledger
4. draft Methods
5. draft Results
6. draft Discussion
7. rewrite Introduction
8. write Conclusions
9. write Abstract
10. choose/finalize title
11. perform citation audit
12. perform figure/caption audit
13. perform coauthor-comment closure audit
14. perform aggressive redundancy edit
15. convert to final Word/submission format

The Introduction and Abstract should be written late because the manuscript argument will sharpen during the Methods/Results/Discussion rebuild.

## 13. Revision discipline

When new analysis changes interpretation:

1. update `current_scientific_state.md`;
2. update the relevant manuscript section;
3. record the change in `revision_log.md` if substantive.

Do not leave contradictory interpretations in different files.

When an older draft conflicts with `current_scientific_state.md`, the current-state note wins unless a later explicit decision supersedes it.

## 14. Completion criterion

The rebuild is complete when the manuscript:

- has one coherent central argument;
- uses P1–P9 consistently;
- treats the snow-water budget as the primary process result;
- clearly separates long-term trends from observing-system transitions;
- incorporates the current Fig. 1–17 sequence;
- addresses substantive coauthor feedback;
- avoids unnecessary duplication of Fox et al. (2025, 2026);
- contains no superseded analysis claims;
- is no longer than necessary to support the science.
