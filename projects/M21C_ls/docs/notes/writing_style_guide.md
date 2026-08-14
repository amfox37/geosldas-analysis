# Manuscript Writing Style Guide

This note defines the preferred writing voice for the M21C land-reanalysis manuscript.

The target is **Andrew Fox's scientific voice, with the statistical and causal restraint characteristic of Rolf Reichle's writing**.

It is not intended to imitate another author's prose mechanically. It records the style that has emerged through previous Fox papers, coauthor editing, and the current manuscript-development process.

## 1. Overall voice

The manuscript should sound:

- technically confident without being emphatic;
- quantitative;
- physically grounded;
- figure- and diagnostic-driven;
- synthesis-oriented;
- cautious about causality;
- concise despite technical complexity;
- calm and low-drama.

The paper should read as:

> we measured this carefully, here is the quantitative result, here is the physical interpretation, and here are the limits of what it establishes.

Avoid promotional or manifesto-like prose.

## 2. Andrew's base style

Preserve the features characteristic of the author's existing scientific writing:

- explain complicated results rather than merely reporting metrics;
- connect diagnostic behavior to physical mechanisms;
- use figures to organize the evidentiary argument;
- compare experiments explicitly;
- give enough quantitative information that interpretation can be checked;
- end important result sequences with a synthesis statement;
- use technical terminology consistently;
- allow moderately complex sentences where needed to preserve scientific relationships.

The manuscript should not become artificially terse or generic.

There should be an identifiable argument running through each section rather than a sequence of disconnected observations.

## 3. Add Rolf's restraint

Apply a strong filter against claims that exceed the evidence.

Prefer:

- `generally`
- `on average`
- `tends to`
- `is consistent with`
- `supports the interpretation that`
- `suggests`
- `coincides with`
- `is associated with`
- `larger`
- `smaller`
- `more spatially coherent`
- `regional`
- `localized`

Use stronger words such as:

- `demonstrates`
- `robust`
- `strong`
- `clear`
- `systematic`

only when the analysis genuinely supports them.

Avoid:

- `proves`
- `unequivocally`
- `universally`
- `dramatically`
- `fundamentally`
- `without degradation`
- `unique information`

unless literally justified.

## 4. Claims must be traceable

A substantive quantitative statement should normally be traceable to:

- a figure;
- a table;
- a named diagnostic;
- an analysis result;
- or a cited source.

Interpretation should follow the evidence rather than precede it.

A good Results paragraph usually follows:

1. state the result;
2. quantify it;
3. interpret it physically;
4. add the important qualification if one exists.

Do not bury the numerical evidence beneath interpretation.

## 5. Statistical precision

Be unusually careful about distinctions such as:

- mean versus median;
- domain mean versus tile-wise behavior;
- pooled versus within-tile relationships;
- magnitude versus signed relationships;
- absolute versus relative differences;
- pointwise versus FDR-controlled significance;
- OL−DA versus DA−OL;
- state differences versus analysis-correction activity;
- correlation versus causality;
- validation versus internal physical consistency.

If a sentence would become false after replacing `generally` with `everywhere`, retain the qualifier.

If an effect is regional, call it regional.

If a relationship is noisy at the tile-year level but evident in aggregate, say so.

## 6. Causal language

This manuscript contains several analyses where causal language requires care.

For observing-system transitions, prefer:

> The April 2015 transition coincides with the introduction of SMAP brightness-temperature assimilation and is consistent with a substantial change in soil-moisture analysis behavior.

Do not write:

> SMAP caused the April 2015 discontinuity

unless the design uniquely establishes that.

For snow-process analyses, stronger physical attribution is justified because OL and DA share forcing and the clean 2001–2006 period isolates snow DA.

Even there, distinguish:

- direct mass accounting;
- controlled statistical corroboration;
- model-internal physical consistency;
- independent validation.

## 7. Validation language

Use `independent validation` only when the comparison is genuinely independent.

For ERA5-Land use terms such as:

- `comparison`
- `large-scale consistency`
- `agreement`

rather than treating ERA5-Land as truth.

For modeled runoff, ET, storage, and energy responses produced by the snow/process analyses, describe them as:

- `physical consistency`
- `model-internal process response`
- `propagation of assimilation-induced differences`

not as independent validation of those fluxes.

## 8. Results writing

Results should be figure-driven but **not figure-caption narration**.

Do not mechanically begin every paragraph with:

> Figure X shows...

Instead, lead with the scientific result where that reads better, and use the figure citation as support.

For example:

> The observational constraint strengthened markedly after the introduction of SMAP in April 2015 (Fig. 17).

is usually better than:

> Figure 17 shows the observational constraint after the introduction of SMAP.

Important Results sections should feel cumulative:

> result → quantitative evidence → mechanism → synthesis.

Do not explain the same physical interpretation separately for every metric.

## 9. Discussion writing

The Discussion should interpret rather than replay the Results.

Do not repeat long strings of numbers unless the number itself is central to the argument.

Use the Discussion to answer:

- what does this tell us about an evolving land observing system?
- what does it tell us about how DA corrections propagate?
- what does it tell us about temporal consistency?
- what does it not establish?

The Discussion should remain restrained even when the result is strong.

For example, the snow-water budget is strong enough to state directly that the response is runoff dominated.

The trend analysis is more nuanced and should retain the distinction between:

- clear structural changes in analysis influence;
- regional modification of RZMC trends;
- broad preservation of the underlying long-term state behavior.

## 10. Physical interpretation

Whenever possible, explain *why* a diagnostic behaves as observed.

Examples:

- surface SM responds more directly than RZMC because the microwave observations primarily constrain the near-surface state;
- snowmelt is expected to respond strongly to added snow because it is the transit pathway by which added snow water leaves the pack;
- wetter RZMC is associated with greater ET and latent heat and lower sensible heat, consistent with expected land-surface energy partitioning.

Do not invent mechanisms merely to make a result sound complete.

If the mechanism is uncertain, say so.

## 11. Sentence and paragraph style

Use complete, technically precise prose.

Moderately long sentences are acceptable where they preserve a meaningful chain of reasoning, but avoid sentences carrying multiple unrelated claims.

Paragraphs should normally have one principal job.

Useful paragraph structure:

- result/topic sentence;
- quantitative support;
- physical/statistical interpretation;
- qualification or synthesis.

Avoid excessive headings, bullets, framework boxes, or pedagogical devices in the final manuscript.

The prose should look like a journal paper, not a technical report or AI-generated outline.

## 12. Avoid AI-style prose

Do not use:

- `Importantly,`
- `Notably,`
- `Interestingly,`
- `It is worth noting that`
- `Taken together` repeatedly
- rhetorical questions in Results
- inflated transitions
- repeated summary sentences saying the same thing in different words

Use `Taken together` occasionally where a genuine synthesis is being made, not as a generic paragraph ending.

Do not describe every result as `robust`, `compelling`, `striking`, or `strong`.

Let the evidence carry the emphasis.

## 13. Specific manuscript examples

### Preferred

> The paired water-year accounting indicates that 64.3% of snow water added by assimilation subsequently leaves the land column as runoff, including 43.1% as surface runoff and 21.2% as baseflow (Fig. 14). An independent controlled estimator likewise yields a runoff-dominated response, supporting the interpretation that the partitioning is not simply a consequence of spatial climatological covariation.

### Avoid

> Figure 14 clearly demonstrates the dramatic and robust effect of snow assimilation on the hydrological cycle, with most of the added water ultimately becoming runoff.

---

### Preferred

> The introduction of SMAP brightness-temperature assimilation in April 2015 coincides with the largest observing-system transition in the soil-water diagnostics, including a positive shift in RZMC DA−OL and increased analysis-correction activity (Fig. 17).

### Avoid

> SMAP fundamentally changed the land reanalysis in April 2015.

---

### Preferred

> OL and DA exhibit nearly identical long-term SCF declines, whereas paired DA−OL RZMC trends are significant over 7.0% of the valid-land domain (Fig. 16). Thus, changes in observational constraint can modify regional state evolution without necessarily imposing a comparable trend on the global or domain-mean record.

### Avoid

> The changing observing system does not affect climate trends.

## 14. Final test for each important sentence

Before keeping a substantive interpretation, ask:

1. What diagnostic supports this?
2. Is the statistic defined precisely?
3. Is the spatial/temporal scope stated correctly?
4. Am I implying causality beyond the experiment?
5. Is this independent validation or model-internal consistency?
6. Would Rolf object to one adjective or qualifier?
7. Can the sentence be made simpler without losing scientific meaning?

If the answer to question 6 is probably yes, revise it before circulation.
