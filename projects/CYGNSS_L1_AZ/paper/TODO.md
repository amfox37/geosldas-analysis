# CYGNSS L1 paper TODO

*2026-09-29. Companion to [paper_outline.md](paper_outline.md). These are outstanding tasks, not completed analyses. Preserve the reported results until checks justify changes.*

## Working decisions and Discover handoff

- **Working target: Remote Sensing of Environment (RSE).** Lead with the terrain-aware physical operator and its practical ensemble implementation, then observation information and independent DA evaluation.
- **Collaboration: USC MiXIL/IGOT team confirmed.** Scientific contributions and formulation/code disclosure still need agreement.
- Start by reading `paper_outline.md`, `../docs/paper_source_manifest.md`, relevant operator/QC/error reports and the actual run configurations. Locate current code and Discover data paths from those records; do not assume uploaded README snapshots describe the current implementation.
- Use the primary experiment `DA_L1_full_xc015_coh040216_fixedop` and its documented OL/L3/SMAP counterparts. Confirm dates, retained observation IDs, coefficient product version, scaling and masks before analysis.
- For each task, record source files, variables/units, sample selection, commands and configuration/commit identifiers. Save reproducible scripts, machine-readable summaries, figures and a short interpretation in the project's established analysis/report locations; link outputs here.
- Begin with existing-data diagnostics. A proposed full DA rerun is a separate scope decision, not implicit authorization from this checklist. If necessary inputs are unavailable, document exactly what is missing and the extraction required; do not mark the task complete.
- Add Melebari et al. (2025), DOI `10.1109/TGRS.2025.3532591`, to the source manifest/reference records. Locate the supplied PDF or obtain an accessible copy; do not infer implementation choices from the paper.

## Priority 1 — Evidence supporting the main claims

- [ ] **Simpler-operator comparison and residual geometry dependence (RSE priority).** Agree a physically meaningful specular-point/single-tile or flat-terrain baseline with MiXIL; these alternatives test different assumptions and should not be conflated. First evaluate the baseline and coefficient operator on the same retained observations and OL states, documenting incidence convention, normalization, vegetation and dielectric model. Apply comparable scaling fitted on the same development sample. Compare held-out innovation mean/spread versus incidence and azimuth, stratified or controlled for tile, season, soil state and spacecraft, with sample counts and uncertainty. Deliver a baseline specification, reproducible comparison script and geometry-diagnostic figure/table. Reduction of residual geometry dependence supports the physical operator; it does not by itself demonstrate better DA skill. Propose a matched full DA comparison separately if warranted.
- [ ] **Sensitivity of the actual 3×5 multi-tile dB observable.** Trace the production formula and notation: tile power reflectivity corresponds to the paper's squared amplitude coefficient, not its complex amplitude. With fixed coefficients, evaluate $S=\sum_t C_t R_t$, $H=10\log_{10}S$ and $\partial H/\partial m_t=(10/\ln10)C_tR'_t/S$. Check numerical derivatives against production evaluations using several finite-difference steps and valid soil-moisture bounds; examine ensemble observation-space spread where member states are available. Stratify by OL soil moisture relative to each tile's actual $m_{vt}=0.02863+0.30673\,\mathrm{clay}$, opacity and geometry, specifying how multi-tile values are weighted. Report native dB sensitivity and any scaling transformation used by DA separately. A common multiplicative attenuation cancels from the dB moisture derivative; tile-dependent attenuation can alter relative weights and attenuation can affect measurement noise. Do not infer weak dB sensitivity directly from declining linear-reflectivity sensitivity in Fig. 15. Deliver sensitivity curves/maps and a short report relating them to the spring diagnostics; associations are not proof of causality.
- [ ] **Fixed-parameter audit with MiXIL.** Recover the small-scale roughness σ_S and intermediate-scale slope parameter σ_L, units, spatial dependence and defaults actually used to construct σ_p. Document Copernicus DEM resolution/processing and whether σ_L calibration from other DEMs transfers. Evaluate the geometry-dependent condition $(q_z\sigma_S)^2\ll1$ over representative retained observations; the paper's 1.25 cm statement is specific to its study geometries. Verify the production dielectric formulation and treatment of the Mironov transition, including continuity and ensembles crossing m_vt. Distinguish approximations from code discrepancies. Deliver a parameter/provenance table and questions for co-authors; holding roughness fixed is an assumption to assess, not justified by geometry-dependent sensitivity.

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
- [ ] **Correlation sensitivity.** If retained-pair diagnostics show meaningful off-diagonal terms, test diagonal R or an alternative correlation model on comparable samples.
- [ ] **Larger domain.** Decide whether expansion is essential to the current claims or better deferred after documenting cost and boundary behavior.

## Manuscript decisions and final checks

- [ ] Align main-text and appendix detail with the working RSE target, retaining the formulation, assumptions and physical added-value evidence in the main text.
- [ ] Agree scientific contributions, author roles and permitted formulation/code disclosure with the confirmed MiXIL/IGOT co-authors.
- [ ] Explain the scalar 3×5 observable and added footprint information early; use “DDM-derived reflectivity” consistently.
- [ ] Keep exactness conditional on the adopted factorization, tile-uniform state and fixed factors; separate vegetation/implementation approximations.
- [ ] Explain why root-zone skill matters and why physical structure can survive z-score scaling.
- [ ] Keep the L1–L3 result scoped to complete configurations in this domain. Retain ASCAT and spring limitations and the 3.9 dB sensitivity.
- [ ] Reconcile source manifest, figure-generation inputs, configuration records and final numerical tables before submission.
