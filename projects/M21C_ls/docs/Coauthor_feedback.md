Yes. I extracted the Word comments directly from the `.docx` files.

* `draft_052426.docx`: **0 comments**
* `draft_052426_lca.docx`: **13 Lauren comments** 
* `draft_052426_rr.docx`: **24 Rolf comments**, plus a general-comments block at the top 

## Big picture

Lauren’s comments are mostly **snow-methods / clarity / possible extra diagnostic** comments.

Rolf’s comments are more structural and strategic. He likes the observing-system-evolution framing, but wants the paper to be **less overlapping with the CYGNSS-DA paper**, more focused on **temporal evolution of impact**, and cleaner about observing-system periods, terminology, figures, and method details.

## Highest-priority actions

1. **Replace “replay” with “reanalysis” or another safer term.**
   Rolf says “replay” has a specific ADAS meaning and is probably not right here.

2. **Refocus away from catch-all sensor comparison and toward observing-system evolution.**
   Rolf thinks Fig. 4 and related discussion overlap too much with the CYGNSS-DA paper. He suggests deemphasizing Fig. 4 and leaning harder on Fig. 5 / time evolution.

3. **Unify observing-system periods everywhere.**
   Define them in Fig. 1 with shading/vertical lines, include MODIS Terra-only and Terra+Aqua periods, and make Figs. 5, 7, 11, and 12 use consistent naming.

4. **Tie the snow analysis more tightly to the paper’s main objective.**
   Both Lauren and Rolf note that the snow results feel somewhat disconnected. Rolf suggests looking at Terra-only vs Terra+Aqua or Terra vs Aqua metrics; Lauren suggests discussing snow–soil moisture assimilation interactions.

5. **Check technical details before revising methods.**
   Rolf flags several: corrected precipitation details, whether CPCU is really used, 277.15 K threshold, H-pol vs H+V counts, SMOS alias-free-zone restriction, MODIS count aggregation, and use of “retrieval screening” for Tb assimilation.

6. **Improve figures.**
   Rolf wants improvement consistently shown in red, reduced whitespace in Figs. 6 and 12, possible split of Fig. 12, and clearer Fig. 8 shading. Lauren also says Fig. 8 gray shading is hard to see.

## Lauren’s comments distilled

Lauren flags:

* Abstract sentence ending “toward dense multi-sensor microwave” is incomplete.
* Snow statement in abstract should maybe include numbers if soil moisture has numbers.
* Figure 1 / observing-system description may belong more naturally in Methods.
* Terra start date is missing from Fig. 1.
* Snow albedo description may need more detail, maybe “highest 95%” climatology.
* The statement that precipitation correction “substantially improves” snow/soil moisture is maybe too strong.
* Precipitation-forcing description is repeated in Methods/Experiments.
* Add surface-temperature screening for MODIS SCF to address poorly screened clouds.
* MODIS percentages in Fig. 4 bottom-left look suspicious.
* Fig. 8 gray shading is hard to see.
* Root-zone ERA5-Land paragraph repeats earlier text.
* Consider a variable showing DA divergence from OL, maybe surface temperature, because it responds to both snow and soil moisture assimilation.
* Add a short discussion of snow–soil moisture assimilation interactions: does soil moisture DA behave differently where SCF assimilation is active?

## Rolf’s comments distilled

Rolf’s general comments are the most important. He says:

* The observing-system-evolution focus is good and sets this apart from the CYGNSS-DA paper.
* There is still too much overlap with the CYGNSS-DA paper, especially Fig. 4 and sensor-impact discussion.
* Shorten data/system descriptions and refer to CYGNSS-DA paper where possible.
* MODIS SCF is in Fig. 4 but not Figs. 5–6; justify or make consistent.
* Consider evolution of SCF observing system, especially Terra-only to Terra+Aqua.
* Tie soil moisture and snow analyses together better.
* Define observing-system periods consistently.
* Table 1 is redundant with Fig. 1; move to supplement if retained.
* Define numbered periods in Fig. 1.
* Rebrand “pre-ASCAT / pre-SMAP / SMAP” to something based on what observations are present.
* Possibly mention connection to forthcoming M21C-Land only in future work/conclusion.
* Clarify whether anyone has done a similar ERA5 land observing-system evolution analysis.
* Explain why ERA5-Land, not ERA5, is used as reference.

Then his inline comments cover terminology, citations, precipitation details, thresholds, observation counts, figure periods, color conventions, and several caption/wording fixes.

## My read

This is actually very useful feedback. It says the paper is not broken; it says the paper needs a sharper identity.

The revised paper should become:

> not “multi-sensor DA improves land states,”
> but “the impact of land DA evolves as the observing system evolves, and this matters for long-term land reanalysis interpretation.”

That means the next revision should probably focus on **structure and emphasis**, not adding lots of new science. The one exception is that a **small MODIS Terra/Aqua or soil-moisture–snow interaction diagnostic** might be worth doing if it is easy.