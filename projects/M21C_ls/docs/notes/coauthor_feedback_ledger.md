# Coauthor Feedback Ledger

> This document tracks substantive feedback from the May–June 2026 coauthor review of the earlier manuscript and records its disposition in the August 2026 rebuild.
>
> The purpose is not to preserve the old manuscript structure. It is to make sure that the concerns raised by Rolf Reichle and Lauren Andrews remain visible while the manuscript is rebuilt around the current science.
>
> Status reflects the **current August 2026 state**, including analyses and figures produced after the original reviews.

## Status categories

- **Resolved** — concern has been addressed by a subsequent manuscript, figure, analysis, or structural decision.
- **Resolved by new analysis** — the original concern led to, or is substantially answered by, analysis completed after the review.
- **Superseded** — subsequent analysis or restructuring changed the scientific framing enough that the original issue no longer applies in its original form.
- **Incorporate during rewrite** — solution is known but must still be implemented in the new manuscript text.
- **Factual check required** — requires confirmation of experiment/configuration details before final wording.
- **Optional / final cleanup** — useful but not central to the scientific rebuild.
- **Not adopted** — considered but intentionally not pursued.

---

# 1. Rolf Reichle — General comments

## RR-G1. Keep the observing-system evolution as the central focus

**Original concern**

Rolf strongly supported focusing the paper on changes in the observing system and noted that this distinguishes the manuscript from the CYGNSS-DA paper.

**Status: Resolved / governing manuscript principle**

This is now the central manuscript question.

The paper is being rebuilt around:

1. evolving observational influence;
2. physical propagation of assimilation corrections;
3. temporal consistency.

The paper should not revert to a generic sensor-comparison or DA-validation structure.

**Current implementation**

- Fig. 1 defines the observing-system chronology.
- Figs. 2–6 characterize changing observation density and OmF behavior.
- Fig. 17 directly quantifies observing-system transitions.
- The Discussion will explicitly address the nonstationarity of observational influence.

---

## RR-G2. Reduce overlap with the CYGNSS-DA manuscript

**Original concern**

There was substantial overlap in system description, data description, Fig. 4-type sensor comparisons, and associated discussion. Rolf recommended cutting these sections substantially and referring to the CYGNSS-DA paper.

**Status: Incorporate during rewrite**

This concern is even easier to address now because the revised CYGNSS manuscript can be cited as **Fox et al. (2026)**.

Use:

- Fox et al. (2025) for established ASCAT + SMAP methods;
- Fox et al. (2026) for the current multi-sensor configuration including SMOS and CYGNSS.

The current manuscript should retain only the methodological information necessary to understand:

- the 2000–2024 experiment;
- the evolving observing system;
- MODIS snow assimilation;
- the new process analyses;
- the trend/transition analyses.

Generic sensor-comparison discussion should be shortened and interpreted mainly through **temporal evolution**.

---

## RR-G3. MODIS OmF treatment needs to be explicit

**Original concern**

MODIS SCF OmFs appeared in the old full-period figure but not consistently in the temporal/spatial OmF figures. If intentional, the manuscript should explain the difference.

**Status: Incorporate during rewrite**

The final Figs. 4–6 and their captions should be reviewed together during the figure/caption audit.

The text should state explicitly:

- which observation types are represented in each OmF figure;
- why MODIS is or is not shown in each diagnostic;
- the applicable period for each OmF statistic.

Do not leave omission/inclusion implicit.

---

## RR-G4. Examine the evolution of the MODIS observing system, especially Terra-only → Terra+Aqua

**Original concern**

The earlier observing-system periods largely ignored the MODIS evolution. Rolf specifically suggested examining the Terra-only to Terra+Aqua transition.

**Status: Resolved**

The authoritative Fig. 1 now begins with the MODIS Terra-only period and explicitly represents the transition to Terra+Aqua.

Additional Terra/Aqua-specific snow evaluation has been produced for the Supplement:

- `fig08_supp_ims_terra_aqua_scope_maps`
- `fig08_supp_ims_terra_aqua_scope_bars`

The observing-system chronology now treats snow observations as part of the evolving observing system rather than as background to the microwave record.

No further major Terra/Aqua analysis is required unless the final manuscript reveals a specific gap.

---

## RR-G5. Tie the soil-moisture and snow analyses together

**Original concern**

The snow section appeared disconnected from the overall objective of the paper.

**Status: Resolved by new analysis**

This concern led to one of the strongest developments in the paper.

The new snow-process work demonstrates:

> snow DA increment → snow storage/melt → runoff + soil water → ET + storage

The common-window differential water budget shows that for snow-addition tile-years approximately:

- 43.1% becomes surface runoff;
- 21.2% becomes baseflow;
- 35.9% becomes ET;
- 3.9% remains as end-of-year storage;
- total runoff is approximately **64.3%**.

The monthly pathway analysis additionally shows a coherent RZMC response.

Later microwave SM assimilation does not simply become more active where snow DA was stronger, but its **signed corrections tend on average to oppose prior snow-added water**.

The manuscript should therefore treat snow and soil moisture as connected parts of the land-water analysis rather than independent validation stories.

---

## RR-G6. Unify the observing-system periods

**Original concern**

Period definitions varied across the manuscript and figures. Rolf recommended:

- defining them in Fig. 1;
- starting with Terra-only and Terra+Aqua;
- relating broader validation periods explicitly to the numbered periods;
- eliminating redundant Table 1.

**Status: Resolved**

Fig. 1 is now the authoritative chronology and defines **P1–P9**.

The registry is authoritative, including:

- P6 begins **1 April 2015**;
- P9 begins **1 December 2021**.

All detailed observing-system references should use P1–P9.

Broader validation periods may remain where needed for statistical sample size but must be explicitly described as aggregations of P1–P9.

The old Table 1 is redundant and should be removed from the main manuscript. It should not be recreated unless a supplemental table later provides information not contained in Fig. 1.

---

## RR-G7. Reduce redundant text

**Original concern**

The manuscript contained substantial repeated explanation.

**Status: Incorporate during rewrite**

The August manuscript is being rebuilt destructively rather than incrementally revised.

Targets:

- shorten Methods through Fox et al. (2025, 2026);
- reduce sensor-by-sensor narrative;
- remove repeated surface-versus-root-zone interpretation;
- avoid repeating figure captions in Results;
- do not repeat numerical Results in Discussion.

The final manuscript should ideally contain more science but less prose than the May draft.

---

## RR-G8. Should the experiment be described as a “sweeper” for M21C-Land?

**Original concern**

Rolf questioned whether to mention that the experiment contributes directly to development of the forthcoming M21C-Land product, particularly given the experiment resolution and data availability.

**Status: Optional / final cleanup**

Do not make this part of the manuscript motivation.

Preferred treatment is one restrained future-work sentence near the end of the Discussion or Conclusions stating that the analysis informs development of MERRA-21C-Land.

Avoid terminology that makes the present experiment appear merely preliminary.

---

## RR-G9. Has evolving land-observation influence been examined similarly in ERA5 or other reanalyses?

**Original concern**

The Introduction should establish whether comparable analyses of changing land observing systems have been conducted, particularly for ERA5.

**Status: Incorporate during rewrite / literature check required**

The rewritten Introduction should contain a concise literature statement addressing previous work on observing-system changes and temporal consistency in land or atmospheric reanalysis.

Do not claim novelty merely because no identical experiment is known.

The literature search should distinguish:

- general reanalysis observing-system discontinuities;
- land-analysis observing-system studies;
- sensor-impact studies;
- explicit multi-decadal land DA transition analyses.

This is a literature task, not a new numerical-analysis task.

---

## RR-G10. Explain ERA5-Land rather than ERA5 as the comparison dataset

**Original concern**

ERA5-Land does not contain a land analysis in the same sense as ERA5, whereas ERA5 assimilates land-relevant observations. The rationale for choosing ERA5-Land therefore needs to be explicit.

**Status: Incorporate during rewrite**

The Methods should include one concise rationale.

The intended framing is:

- ERA5-Land offers a high-resolution, long-term, internally consistent land-model integration driven by ERA5 forcing;
- it provides a useful large-scale reference trajectory without directly assimilating the same land observations in the way ERA5 does;
- comparison is therefore used as a **consistency comparison**, not independent truth.

Final wording should be checked against the exact ERA5/ERA5-Land configurations.

---

# 2. Rolf Reichle — Inline comments and technical checks

## RR-I1. Avoid “replay” terminology

**Original concern**

“Replay” has a specific meaning in the atmospheric DA system and is not necessarily appropriate here.

**Status: Incorporate during rewrite**

Prefer:

- `land reanalysis experiment`
- `land analysis`
- `OL simulation`
- `DA experiment`

Use “replay” only if technically required and explicitly defined.

---

## RR-I2. Put sensor-specific references next to sensor-specific statements

**Status: Incorporate during rewrite**

The rewritten Introduction and Methods should attach mission/product references directly to the corresponding statements rather than relying on distant grouped citations.

---

## RR-I3. Define observing-system periods in Fig. 1

**Status: Resolved**

Implemented through P1–P9 in the new Fig. 1.

---

## RR-I4. Corrected precipitation description may not match M21C-Land PRECTOTCORR

**Original concern**

The draft risked conflating the precipitation treatment used in this experiment with the eventual M21C-Land precipitation forcing.

**Status: Factual check required**

Keep the current experiment description strictly specific to the actual experiment.

Do not imply that its corrected precipitation is identical to the final M21C-Land forcing.

Final wording should be checked with Qing.

---

## RR-I5. CPC Unified may not be part of the actual corrected precipitation

**Status: Factual check required**

Verify the exact precipitation products and correction procedure with Qing before finalizing Methods.

Do not carry the old CPCU/GPCP/IMERG description forward without verification.

---

## RR-I6. Localization and observation errors are not OL settings

**Original concern**

The old manuscript implied that OL and DA shared localization and observation-error settings.

**Status: Incorporate during rewrite**

State that OL and DA share:

- forcing;
- model physics;
- grid;
- parameterization;
- ensemble/perturbation setup where applicable.

Localization and observation-error parameters apply to DA, not to OL.

The analysis interval may matter to monitored OmF sampling, but it should not be described as an OL assimilation setting.

---

## RR-I7. Do not overclaim that fixed configuration isolates intrinsic observation information content

**Original concern**

The DA system was not independently optimized for each observing-system era, so fixing the assimilation configuration does not guarantee that temporal changes represent intrinsic sensor information content alone.

**Status: Resolved conceptually / incorporate during rewrite**

This caveat is now a governing manuscript rule.

Use language such as:

> changes are consistent with the evolving observing system under a common DA configuration

rather than:

> changes uniquely quantify the intrinsic information content of each observing system.

---

## RR-I8. Warm-mask temperature threshold may be 277.15 K

**Status: Factual check required**

Verify against the actual analysis code and Fox et al. (2026).

Do not retain 275.15 K or another value from an old draft without checking.

---

## RR-I9. “Snow-possible” terminology

**Status: Resolved editorially**

Use a precise domain description such as:

- `Northern Hemisphere seasonal-snow domain`, or
- the exact mask definition used by the production analysis.

Avoid informal “snow-possible” wording.

---

## RR-I10. SMAP/SMOS observation counts: H-only or H+V?

**Status: Factual check required**

Verify what Figs. 2–3 count and make the definition consistent with the actual assimilation configuration.

The final caption/method should say explicitly whether counts represent:

- H polarization only;
- H + V separately;
- assimilation windows;
- or another observational unit.

---

## RR-I11. SMOS wording: not “retrieval screening”; include alias-free-zone restriction

**Status: Factual check required / incorporate during rewrite**

Because SMOS Tb is assimilated, use Tb/observation screening terminology rather than soil-moisture retrieval screening.

Verify and document the alias-free-zone restriction from the current configuration / Fox et al. (2026).

---

## RR-I12. Define what constitutes a MODIS observation in Figs. 2–3

**Original concern**

MODIS enters initially on the 0.05° CMG and is aggregated to tile space; the plotted observation count needs an explicit definition.

**Status: Factual check required**

Check the plotting inputs and describe precisely what is counted.

Do not assume the plotted counts refer to raw CMG pixels.

---

## RR-I13. State the period used for full-period OmF panels

**Status: Incorporate during caption audit**

Each OmF statistic should have an explicit temporal scope.

For sensors with differing availability, “full period” means the valid observation period for that sensor unless the analysis code defines otherwise.

Verify before final captions.

---

## RR-I14. Standardize improvement sign/color conventions

**Original concern**

Rolf recommended a consistent convention in which improvement is visually emphasized, preferably red.

**Status: Largely resolved / final figure audit**

Current figure generation uses a consistent positive/red improvement convention where appropriate.

Before submission perform a final figure audit covering:

- DA−OL versus OL−DA definitions;
- correlation versus error metrics;
- color-bar direction;
- caption terminology;
- red = improvement wherever a common convention is meaningful.

Trend maps are different because red/blue represent positive/negative trends rather than skill improvement; captions must make this distinction obvious.

---

## RR-I15. Reduce Fig. 6 white space / improve layout

**Status: Final figure audit**

The current regenerated figure should be evaluated rather than applying old layout edits blindly.

Only make further changes if the current version remains difficult to read.

---

## RR-I16. Fig. 6 legend / unused panels

**Status: Final figure audit**

Review current version. Do not preserve empty `n/a` panels merely because they existed in the May figure.

---

## RR-I17. Terra versus Aqua comparison

**Status: Resolved**

Explicit Terra/Aqua supplementary analyses now exist.

See RR-G4.

---

## RR-I18. ERA5-Land SM map readability / possible splitting

**Status: Final figure audit**

The current Fig. 12 has been rebuilt since the original comment.

Assess the current figure rather than automatically splitting it.

The main requirement is that:

- surface versus root-zone behavior is understandable;
- positive/red means improvement consistently across displayed metrics.

---

## RR-I19. Experiment labels should simply be OL and DA

**Status: Resolved**

Use **OL** and **DA** consistently.

Do not carry old internal experiment names into figure captions or Results unless required for reproducibility.

---

## RR-I20. Use “OmF” consistently

**Status: Incorporate during rewrite**

Prefer `OmF` throughout after first definition.

Use `OmF standard deviation` or `OmF std. dev.` consistently.

Avoid switching unnecessarily among:

- O−F
- forecast–observation mismatch
- innovation
- residual

unless a conceptual distinction is intended.

---

## RR-I21. Distinguish statements specifically about the soil-moisture analysis

**Status: Incorporate during rewrite**

When interpreting ASCAT/SMOS/SMAP/CYGNSS effects, identify whether the statement refers to:

- soil-moisture analysis behavior;
- the broader land analysis;
- snow analysis;
- or all DA corrections.

Avoid generalizing a soil-moisture diagnostic to the entire land system.

---

## RR-I22. Funding / IMVI acknowledgment

**Status: Optional / final cleanup**

Check acknowledgments and funding immediately before submission.

This is not part of the scientific rebuild.

---

# 3. Lauren Andrews — General and inline comments

## LA-1. Abstract snow and soil-moisture results should have consistent quantitative treatment

**Original concern**

The abstract quantified soil-moisture effects but described snow results only qualitatively.

**Status: Superseded / address when Abstract is rewritten**

The final Abstract will be rewritten only after the manuscript stabilizes.

The new snow-water-budget result now provides a much stronger quantitative snow/process statement than the old validation-only description.

Likely abstract-level number:

> approximately **64% of assimilation-added snow water subsequently leaves as runoff**.

Do not force artificial symmetry between SM and snow metrics if the principal scientific results differ.

---

## LA-2. Instrument-detail material belongs in Methods rather than the Introduction

**Status: Incorporate during rewrite**

The Introduction should motivate the evolving-observing-system problem.

Detailed mission characteristics belong in Methods and should be shortened substantially using Fox et al. (2025, 2026).

---

## LA-3. MODIS Terra start date was missing

**Status: Resolved**

The revised Fig. 1 defines the full P1–P9 chronology beginning with the Terra-only SCF period.

Final caption should retain the exact date.

---

## LA-4. Snow-albedo climatology description may need more precision

**Original concern**

The old statement may have oversimplified how the prescribed snow albedo is constructed; Lauren suggested additional detail, possibly related to the upper percentile of the MODIS climatology.

**Status: Factual check required**

Verify the current GEOS LDAS snow-albedo implementation with Lauren/source code before finalizing Section 2.1.

Do not repeat the old “MODIS-derived climatology” description without checking the precise definition.

---

## LA-5. MODIS surface-temperature screening should be described

**Original concern**

Surface-temperature screening is used to reduce poorly screened cloud contamination and should be documented.

**Status: Factual check required / incorporate during rewrite**

Verify the actual threshold and screening sequence in the current configuration.

Include this either:

- in the MODIS observations subsection, or
- in the SCF assimilation subsection,

but not redundantly in both.

---

## LA-6. Precipitation forcing description was repetitive

**Status: Incorporate during rewrite**

Describe the forcing once in the experiment/model section.

Refer back to it rather than repeating the full precipitation-correction description in multiple subsections.

---

## LA-7. MODIS percentages in the old OmF figure looked wrong

**Status: Factual/figure check required**

Do not assume the old normalized percentages are correct.

The current regenerated Fig. 4 should be checked against its production inputs.

If the current Fig. 4 no longer uses the suspect calculation, mark this resolved during the final figure audit.

---

## LA-8. Gray shading in old IMS figure was difficult to see

**Status: Final figure audit**

The current Fig. 8 has been regenerated since the comment.

Inspect the present rendering and adjust only if the relevant mask/domain remains visually ambiguous.

---

## LA-9. Root-zone interpretation was repetitive

**Status: Incorporate during rewrite**

Discuss the weaker/more heterogeneous RZMC response once clearly in Results and develop its implications in Discussion.

Do not repeat essentially the same surface/root-zone paragraph across multiple Results sections.

---

## LA-10. Find a variable showing the combined effect of all assimilation on the land trajectory

**Original concern**

Lauren suggested that a broader state variable, possibly surface temperature, might reveal where the full assimilation system causes DA to diverge from OL and might connect snow and soil-moisture assimilation.

**Status: Resolved by new analysis, using stronger variables than the original suggestion**

Rather than using temperature alone, subsequent analyses examined physically integrated hydrologic and energy responses.

The new analyses show that DA-induced differences propagate into:

- total land water;
- runoff;
- ET;
- RZMC;
- latent heat;
- sensible heat;
- evaporative fraction.

This provides a clearer physical synthesis than a single temperature variable.

No dedicated TS analysis is currently needed for the main manuscript.

---

## LA-11. Investigate interactions between snow and soil-moisture assimilation

**Original concern**

Lauren specifically asked whether soil-moisture assimilation works harder, perhaps in the opposite direction, where SCF assimilation is more active.

**Status: Resolved by new analysis**

The investigation produced two distinct findings.

### Magnitude/activity hypothesis

Not supported robustly:

> stronger snow DA does **not** reliably lead to greater subsequent SM-DA activity.

This should not become a manuscript claim.

### Signed interaction

A weaker but coherent result is present:

> where snow DA adds water, later microwave SM analysis corrections tend on average to be more negative.

Interpret this as partial population-level opposition/counteraction, not deterministic cancellation.

The detailed diagnostics belong in the Supplement, with at most a short synthesis in the main text.

This directly addresses Lauren's original question while preserving the negative result for the activity-magnitude hypothesis.

---

# 4. Major concerns raised during early structural review that are now answered by new August analysis

These items were not all literal coauthor comments, but they emerged directly from the review process and should remain in the ledger because the August analyses were designed to close them.

## A1. The paper argued temporal nonstationarity without directly demonstrating its effect on the analysis

**Status: Resolved by new analysis**

The new trend and transition analyses now directly test temporal consistency.

### Fig. 16

Shows June 2000–May 2024 OL, DA, and paired DA−OL trends for:

- RZMC
- SCF

Key result:

- regional RZMC trends are modified by DA;
- the long-term SCF decline is essentially identical in OL and DA.

### Fig. 17

Shows observing-system transitions in:

- RZMC DA−OL;
- soil-moisture correction magnitude;
- prognostic soil-water correction activity.

The dominant system-level transition occurs at:

> **April 2015 / P5→P6 / introduction of SMAP Tb assimilation**

This is quantified through interrupted time-series estimates and independently supported by changepoint detection.

The manuscript can now distinguish:

> **structural nonstationarity of observational influence**

from

> **widespread secular distortion of the land-state record**.

That distinction should become central to the Discussion.

---

## A2. Need evidence that observing-system transitions are not merely visible by eye

**Status: Resolved by new analysis**

Independent changepoint analysis now provides corroboration.

Accepted results include:

- 37 accepted independent breaks;
- 20 paired DA−OL series;
- 17 correction diagnostics;
- 10 primary series breaking exactly in April 2015;
- 9 of those also significant in the known-date P6 analysis;
- no corresponding accepted April-2015 break in paired OL or DA state controls.

The full breakpoint-boundary matrix remains supplemental.

---

## A3. Concern that the paper might overstate risks to long-term trends

**Status: Resolved by new analysis**

The new trend analysis substantially sharpens and moderates this interpretation.

The paper should **not** argue that changing observing systems render long-term trends unusable.

Instead:

- precipitation provides a null forcing-control case;
- snow mass/depth have essentially no DA-induced trends;
- SCF decline is almost identical in OL and DA;
- RZMC trends are regionally modified;
- no significant global land-mean RZMC trend is introduced.

Preferred conclusion:

> the evolving observing system clearly changes analysis influence and can introduce structural transitions, while broad long-term state behavior is generally preserved.

---

## A4. Old “directness-of-constraint” conclusion was too restrictive

**Status: Superseded by new analysis**

The old draft argued that impacts become progressively weaker for quantities further removed from the observations.

The snow-water-budget result demonstrates that this is incomplete.

New interpretation:

> observational constraints remain strongest and most immediate near the observed quantities, but assimilation-induced state and mass corrections can propagate strongly through modeled hydrology into runoff, evapotranspiration, storage, soil moisture, and surface-energy partitioning.

The old broad statement that integrated hydrologic reservoirs remain largely unaffected should not return.

---

# 5. Current figure responses to coauthor feedback

The manuscript now has the following intended main-text sequence:

1. observing-system timeline / P1–P9
2. mean assimilated observations per day
3. monthly assimilated-observation counts
4. full-period OmF standard-deviation improvement
5. monthly OmF improvement evolution
6. period/sensor OmF improvement maps
7. ISMN soil-moisture validation
8. IMS snow-cover skill
9. SNOTEL SWE skill
10. GHCN snow-depth skill
11. ERA5-Land soil-moisture bars
12. ERA5-Land soil-moisture maps
13. ERA5-Land snow comparison
14. snow-DA water budget
15. monthly snow-DA hydrologic pathway
16. long-term RZMC and SCF trends
17. observing-system transitions

Important feedback closure provided by the revised figure set:

- **period inconsistency** → resolved by Fig. 1 / P1–P9;
- **Terra-only versus Terra+Aqua** → represented in Fig. 1 and supplemental IMS analysis;
- **snow disconnected from SM/overall objective** → resolved by Figs. 14–15 and associated hydrologic analysis;
- **need combined DA trajectory/process response** → addressed through hydrologic and energy diagnostics;
- **need direct temporal-consistency evidence** → resolved by Figs. 16–17;
- **need evidence beyond known transition dates** → supplemental breakpoint-boundary analysis.

The old figure-specific comments should always be checked against the **current regenerated figures**, not automatically applied to superseded May versions.

---

# 6. Remaining open factual checks before final Methods

These items should be resolved during the Methods rebuild:

1. exact corrected-precipitation formulation used in this experiment;
2. distinction between this precipitation and eventual M21C-Land PRECTOTCORR;
3. exact warm-mask/frozen-soil temperature threshold;
4. exact SMOS alias-free-zone and QC terminology;
5. H-only versus H+V treatment in SMOS/SMAP observation counts;
6. exact definition of a counted MODIS observation in Figs. 2–3;
7. MODIS cloud and surface-temperature screening;
8. exact snow-albedo climatology formulation;
9. exact experiment/spin-up dates and final naming;
10. ERA5-Land versus ERA5 rationale;
11. final funding/acknowledgment language.

These are **not unresolved scientific analyses**. They are factual/configuration checks.

---

# 7. Issues explicitly not requiring additional primary analysis

Unless the manuscript rebuild exposes a new problem, do not reopen:

- the P1–P9 period definition;
- Terra versus Aqua as a primary analysis question;
- the basic snow→hydrology pathway;
- whether stronger snow DA causes greater SM-DA activity — it does not robustly;
- the signed snow→SM opposition result;
- the common-window water-budget construction;
- RZMC residence-time analysis as a main result;
- delayed infiltration under strict non-overlap;
- whether snowmelt is a terminal budget sink;
- whether April 2015 is the dominant detected observing-system transition;
- whether DA creates a significant global RZMC trend;
- whether SCF decline is created by DA.

These questions have been investigated sufficiently for the current manuscript.

---

# 8. Coauthor-feedback closure criteria

Before recirculating the rebuilt manuscript to Rolf and Lauren, confirm that:

- [ ] The Introduction is organized around observing-system evolution rather than sensor history.
- [ ] P1–P9 is used consistently.
- [ ] Fig. 1 includes Terra-only and Terra+Aqua.
- [ ] The old redundant observing-system table is removed from the main text.
- [ ] Methods rely heavily on Fox et al. (2025, 2026) rather than repeating prior implementation detail.
- [ ] “Replay” terminology has been removed or justified.
- [ ] Precipitation forcing has been factually verified.
- [ ] SMOS/SMAP polarization/count definitions have been verified.
- [ ] MODIS observation-count and screening definitions have been verified.
- [ ] Snow-albedo description has been verified.
- [ ] ERA5-Land choice is explained.
- [ ] Snow and soil moisture are connected through the hydrologic-process analysis.
- [ ] The negative snow-DA→SM-DA activity result is not overstated.
- [ ] The signed opposing SM-DA tendency is described cautiously.
- [ ] The directness-of-constraint conclusion reflects the new water-budget result.
- [ ] Fig. 16 is used to discuss trend preservation versus regional modification.
- [ ] Fig. 17 is used to discuss the April 2015 observing-system transition.
- [ ] Structural transitions are distinguished explicitly from secular trends.
- [ ] The manuscript does not imply that reanalysis trends are generally invalid.
- [ ] ERA5-Land is described as a consistency comparison rather than truth.
- [ ] All figure sign/color conventions have been audited.
- [ ] Repetitive Results/Discussion prose has been removed.
- [ ] Any M21C-Land connection is brief and placed near the end.
