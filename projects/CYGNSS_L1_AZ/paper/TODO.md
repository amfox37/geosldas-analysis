# CYGNSS L1 paper TODO

*2026-09-29. Companion to [paper_outline.md](paper_outline.md). These are outstanding tasks, not completed analyses. Preserve the reported results until checks justify changes.*

## Priority 1 — Evidence supporting the main claims

- [ ] **Root-zone update pathway and mechanism (§5.3).** Verify the actual EnKF state updates. Diagnose srfexc, rzexc and catdef analysis increments and corresponding SFMC/RZMC changes. Separate direct cross-covariance updates from subsequent model propagation. Compare L1 and L3 timing, observation counts, wetting events and dry-down persistence. Deliver one focused figure or table; do not assume sampling explains the difference.
- [ ] **Root-zone robustness and practical size (Table 4).** Add OL baseline R/anomR/ubRMSE and percentage ubRMSE reductions using consistent samples. Compare station weighting with equal tile weighting and leave-one-tile-out results across the 19 root-zone tiles. Report nominal CI interpretation and the 24-comparison context; keep marginal surface results secondary.
- [ ] **ASCAT penalty and validation coverage (§5.1/§6).** Overlay surface and root-zone stations with affected Four Corners tiles. Compare L1, L3 and SMAP on common monitor samples; summarize by region, terrain slope, land cover, season/frozen-snow QC and distance to the preprocessing boundary. Distinguish actual evidence of an edge effect from proximity alone. State when station coverage cannot adjudicate the discrepancy.
- [ ] **Configuration development versus evaluation (§3).** Document Jan–Jun 2020 screen tuning and the 2020–2022 scaling climatology. Evaluate the frozen screen and headline skill on 2021–2022. Label residual climatology overlap; consider training-period-only scaling if an independent evaluation is feasible.

## Priority 2 — Observation/error and computational audit

- [ ] **Selection and correlation (§§3–4, Fig. 5).** Trace preprocessing versus reader selection, define “best,” and confirm actual coordinates used by R. Label raw/screened/selected pair populations. Calculate retained-pair separation distributions and off-diagonal R magnitudes; assess whether xcorr 0.15° materially changes updates/conditioning or is chiefly a safeguard. Adjacent tiles do not guarantee 36 km observation separation.
- [ ] **Error assumptions (§4).** Distinguish measured innovation correlation from inferred observation-error correlation and list forecast-error separation assumptions. State Desroziers gain/covariance caveats. Recover the documented rationale for errstd 2.75 dB rather than constructing one retrospectively.
- [ ] **Weak L1 own-fit reduction (§5.1).** Examine forecast observation-space variance, normalized innovations, gain and increment distributions alongside the −0.76% own-fit change and L3's −0.54% L1-monitor change. Do not compare percentage own-fit changes across species as information measures.
- [ ] **Accuracy and cost (Table 2).** Verify timing environments and units (Python comparison versus production per member). Report observation counts, total core-hours, hardware, parallel throughput, failed/retried processing and storage where available. Distinguish measured numerical discrepancies from physical/model error; assess larger-domain feasibility from actual preprocessing cost.

## Priority 3 — Complete publication figures and tables

- [ ] Regenerate Figs. 4–5 from per-observation and pair data on Discover; include sample definitions, counts and uncertainty. Check the dashed curve's actual meaning.
- [ ] Regenerate the other figures from the documented post-fix data; check labels, experiment IDs, units and consistency with the revised wording. Finalize Fig. 2 layout.
- [ ] Complete Table 3's missing SMAP L1-area row from source data. Make L1-area comparisons primary; retain domain-wide effects in text or supplement.
- [ ] Add appropriate paired uncertainty to Table 3 O−F changes, preserving spatial/temporal dependence. With only three years, do not rely on a year-only bootstrap without assessing its limitations.
- [ ] Check spring in-situ skill on matched samples and the SMAP Tb monitor's seasonal QC/coverage. Avoid assuming monitor disagreement identifies which system is wrong.
- [ ] Verify Table 4 baselines/differences and all headline numbers against source reports. Label nominal CIs; document any multiplicity analysis.
- [ ] Keep Fig. 8 and §5.2 terminology consistent: accumulated assimilation-induced forecast change, not instantaneous analysis increment.

## Optional experiments — Decide scope before running

- [ ] **Seasonal/adaptive R sensitivity.** Decide whether this answers a central question for this paper. Define choices using a development period and freeze them before evaluation where feasible. α_opt does not directly prescribe R. Report negative or partial outcomes if run; do not select inclusion by whether spring degradation disappears.
- [ ] **Simpler-operator comparison.** Compare with a clearly specified specular-point/single-tile operator under comparable QC, sampling, scaling and error treatment. Test what footprint integration adds after z-score scaling; distinguish operator-only diagnostics from full DA experiments.
- [ ] **Correlation sensitivity.** If retained-pair diagnostics show meaningful off-diagonal terms, test diagonal R or an alternative correlation model on comparable samples.
- [ ] **Larger domain.** Decide whether expansion is essential to the current claims or better deferred after documenting cost and boundary behavior.

## Manuscript decisions and final checks

- [ ] Select journal based on the intended scientific emphasis and audience; set main-text versus appendix operator detail accordingly.
- [ ] Discuss contributions, co-authorship/acknowledgment and permitted formulation/code disclosure with IGOT collaborators.
- [ ] Explain the scalar 3×5 observable and added footprint information early; use “DDM-derived reflectivity” consistently.
- [ ] Keep exactness conditional on the adopted factorization, tile-uniform state and fixed factors; separate vegetation/implementation approximations.
- [ ] Explain why root-zone skill matters and why physical structure can survive z-score scaling.
- [ ] Keep the L1–L3 result scoped to complete configurations in this domain. Retain ASCAT and spring limitations and the 3.9 dB sensitivity.
- [ ] Reconcile source manifest, figure-generation inputs, configuration records and final numerical tables before submission.
