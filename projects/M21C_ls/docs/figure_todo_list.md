Yes — here’s the **figure-specific to-do list** from Lauren + Rolf’s feedback.

## Global figure conventions

1. **Unify “improvement” colors across all maps.**
   Rolf suggests using **red = improvement** everywhere, because red is easier to see and most maps are improvement-dominated. This means changing the sign or colorbar labels for Figs. 4, 6, 8, 12, maybe 13.

2. **Define observing-system periods once, preferably in Fig. 1.**
   Add gray shading and vertical lines in Fig. 1, then use the same period names/numbers in Figs. 5, 6, 7, 11, and 12.

3. **Rename periods based on what observations are present, not what is absent.**
   Instead of “pre-ASCAT / pre-SMAP / SMAP era,” consider something like:

   * SCF only
   * SCF + active microwave SM
   * SCF + active/passive microwave SM
   * SCF + SMAP-era microwave
   * SCF + SMAP + CYGNSS
     Rolf specifically disliked the inconsistency across current period definitions. 

4. **Decide whether Table 1 survives.**
   Rolf thinks Table 1 is redundant if Fig. 1 defines the observing-system timeline. If kept, move it to supplement.

---

## Figure 1 — observing-system timeline

**To-do:**

* Add **Terra MODIS start date**. Lauren noticed it is missing.
* Add **Terra-only vs Terra+Aqua** transition. Rolf specifically wants the observing-system periods to start with:

  * MODIS Terra-only SCF
  * MODIS Terra + Aqua SCF
  * then ASCAT, SMOS, SMAP, CYGNSS additions.
* Add shaded/numbered observing-system periods directly on the graphic.
* Use the same period labels later in Figs. 5, 6, 7, 11, and 12.
* Possibly include the sensor names inside each shaded period, but avoid crowding.
* If Fig. 1 becomes complete, delete Table 1 or move it to supplement.

---

## Figure 2 — mean assimilated observations per day

**To-do:**

* Clarify whether SMOS/SMAP counts are **H-pol only** or **H+V**. Rolf thinks probably H-pol only, as in the CYGNSS-DA paper.
* Clarify how MODIS observations are counted:

  * native 0.05° CMG observations?
  * or tile-aggregated observations actually passed to GEOS LDAS?
* For SMOS, remove language about “retrieval screening” because you assimilate **Tb**, not retrievals.
* Add/confirm the **SMOS alias-free-zone restriction** if that is part of the processing.
* Check whether the plotted sampling density for SMAP is inflated if both polarizations/views are counted.
* Make the caption explicit: “counts are for observations entering the analysis after QC and aggregation.”

---

## Figure 3 — temporal evolution of assimilated observations

**To-do:**

* Again clarify **H-pol only vs H+V** for L-band Tb counts.
* Clarify whether MODIS counts are tile-space counts after aggregation.
* Add vertical lines and/or shaded observing-system periods matching Fig. 1.
* Include the Terra-only to Terra+Aqua transition if it is visible.
* If Terra/Aqua transition is not visible or uninteresting, say so briefly in text.
* Clean up any legend or label crowding.
* This figure should become one of the core figures, because Rolf wants the paper centered on temporal evolution.

---

## Figure 4 — full-period OmF maps by sensor

**To-do:**

* Consider **deemphasizing or moving to supplement**, because Rolf thinks this overlaps strongly with the CYGNSS-DA paper.
* If kept in main text, make the role clear: broad sensor-specific context, not the main result.
* Specify the **exact period** used for each subplot, since each sensor has a different availability period.
* Fix/check the **MODIS SCF percentages**; Lauren thinks the bottom-left values look off.
* Decide whether MODIS belongs here if it is excluded from Figs. 5–6.
* Change colors so **improvement is red**, consistent with later maps.
* Use “OmF std-dev” consistently.
* Make sign convention unmistakable in caption and colorbar.

---

## Figure 5 — monthly OmF evolution

**To-do:**

* Treat this as a central figure. Rolf explicitly says this is closer to the real paper than Fig. 4.
* Use the observing-system periods defined in Fig. 1.
* Clean up shaded-period labels; Rolf noted text overlaps with lines in current periods III–VII.
* Decide whether to include MODIS SCF OmF time series.
  If not, explicitly justify: e.g., Fig. 5 focuses on microwave soil-moisture observing systems; MODIS SCF is evaluated separately.
* Make sign convention consistent with Fig. 4/6.
* If positive means improvement here but negative means improvement in maps, fix that inconsistency.

---

## Figure 6 — OmF maps by period and sensor

**To-do:**

* Use the same numbered/shaded periods defined in Fig. 1.
* Include MODIS SCF periods or explicitly explain why MODIS is omitted.
* Reduce whitespace:

  * crop longitude range if acceptable,
  * reduce vertical gaps between rows,
  * move map-average numbers inside axes,
  * remove unused “n/a” panels if possible.
* Consider putting colorbar in unused top-right space.
* Use **red = improvement**.
* Make sure each panel’s period and sensor availability are obvious.
* This figure is important, but it may be visually too heavy; tighten aggressively.

---

## Figure 7 — ISMN skill by period

**To-do:**

* Replace “pre-ASCAT / pre-SMAP / SMAP era” with period names tied to Fig. 1.
* Or explicitly map these broader validation periods onto the numbered Fig. 1 periods.
* Make clear that these are broader validation aggregations, not the full period breakdown used in Fig. 6.
* Consider renaming to something like:

  * SCF-only period
  * microwave pre-SMAP period
  * SMAP-era microwave period
* Confirm sign convention: positive = DA better for all metrics.

---

## Table 2 / Figure 8 — IMS snow-cover skill

**To-do:**

* Improve gray shading visibility; Lauren says it is hard to see.
* Consider whether Terra-only vs Terra+Aqua snow-cover evaluation should be added here or as a small supplementary figure/table.
* If Terra/Aqua transition is not worth showing, consider a simpler Terra vs Aqua OmF/metric comparison, coordinated with Lauren’s paper.
* Make color convention consistent: if red = improvement, then miss-rate and false-alarm-ratio panels need sign flipping or clear relabeling.
* Check caption typo: “comparing OL DA against IMS” should be “comparing OL and DA against IMS.”

---

## Figure 9 — SNOTEL SWE

**To-do:**

* No major direct figure-format comment, but this figure needs to be better connected to the main story.
* Emphasize visually/caption-wise that SWE changes are weak compared with SCF changes.
* Consider adding a note in caption: negative DA−OL differences indicate improvement for RMSE/ubRMSE/|bias|.
* Coordinate with Fig. 8 and Fig. 10 so the snow story is visually coherent: SCF improves; SWE/SD weak.

---

## Figure 10 — GHCN snow depth

**To-do:**

* Same as Fig. 9: make the weak/heterogeneous snow-depth response visually clear.
* Ensure sign convention matches Fig. 9.
* Consider whether Fig. 9 and Fig. 10 could be shortened or one moved to supplement if the main story gets too snow-heavy.
* Keep if needed to support the “SCF ≠ snow mass/depth” point.

---

## Figure 11 — ERA5-Land soil moisture bar plot

**To-do:**

* Rename periods consistently with Fig. 1.
* Add a sentence/caption note explaining why ERA5-Land is used rather than ERA5.
* Make sure period labels match Fig. 7 and Fig. 12.
* Confirm the warm/snow-free masking threshold, especially Rolf’s note that Qing may use **277.15 K**, not the current value.

---

## Figure 12 — ERA5-Land soil moisture maps

**To-do:**

* Reduce whitespace, as for Fig. 6.
* Consider splitting the figure:

  * surface and root-zone as separate figures, or
  * correlation/anomR and ubRMSE as separate figures.
* Use **red = improvement** consistently. This may be easier if maps are grouped by metric.
* Remove repeated text in the results; Lauren flagged the root-zone interpretation as repetitive.
* Period names must match Fig. 1 / Fig. 11.

---

## Figure 13 — ERA5-Land snow comparison

**To-do:**

* Change `LS_OL` and `LS_DA` to simply `OL` and `DA`, unless there are actually other experiment names.
* Fix caption wording: current “Rows Columns show…” typo.
* Consider replacing “ever-snow” in the method/caption with Rolf’s suggestion: **“snow-possible” domain**.
* Make sure the figure reinforces the key snow result: SCF improves clearly; SWE and snow depth weakly.
* Check whether it belongs in main text or supplement depending on length.

---

## Possible new figure / diagnostic

There are two possible additions from the comments, but I would only do them if easy:

1. **MODIS Terra-only → Terra+Aqua transition diagnostic.**
   Could be OmF, observation counts, or IMS skill split by period. This would connect snow more directly to observing-system evolution.

2. **DA − OL time series/map for a state variable affected by both SM and snow DA.**
   Lauren suggested something like **surface temperature**, because it may show the combined influence of soil moisture and snow assimilation. This could help tie snow and soil moisture together.

My priority order would be: **fix Fig. 1 + periods first, then Fig. 5, then Fig. 6/Fig. 12 layout/color conventions, then decide what to do with Fig. 4 and MODIS.**