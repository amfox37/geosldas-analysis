# Source manifest for a paper on CYGNSS L1 assimilation in GEOSldas

*Compiled 2026-09-28 as the starting point for drafting a scientific paper on the whole CYGNSS L1 operator project: forward
operator, preprocessing, GEOSldas implementation, and the AZ data assimilation experiments and validation.* It lists every
relevant document in the two git repositories, what each one covers, whether its numbers can be used, and which paper section it
feeds.

## How to use this list

- **Two repositories** (clone both; paths below are relative to each repo root):
  - `CYGNSS_operator` (`git@github.com:amfox37/CYGNSS_operator.git`): forward model, tile-coefficient operator, preprocessing
    pipeline, QC and thinning. Cited below as **[OP]**.
  - `geosldas-analysis` (`git@github.com:amfox37/geosldas-analysis.git`), project `projects/CYGNSS_L1_AZ/` unless stated: DA
    experiments, scoring, in-situ validation. Cited below as **[GA]**.
- **Status labels.** These matter because of the operator bug:
  - **CURRENT**: post-fix (on or after 2026-09-25). The numbers can be used in the paper.
  - **METHOD**: describes methods, design or derivation that are still valid. There are no DA skill numbers to reuse (or they are
    not DA numbers).
  - **PRE-FIX**: produced with the defective GEOSldas L1 operator (fixed 2026-09-25). **Every L1 DA/O-F number in these documents
    is invalid.** Use them only for history, method ideas or motivation. Do not quote their results.
  - **BACKGROUND**: context, prior work or reference material.
- **The operator bug in one line** (needed to understand PRE-FIX): the GEOSldas operator took each tile's footprint coefficients and
  incidence angle from the obs nearest the tile in the whole day file, while the reader took the observed value from the nearest
  obs inside the 3-h window. As a result, O and F came from different obs on about 30% of tile-cycles. It was fixed by GEOSldas_GridComp
  commits `cdfaaf4` + `11dfdb1`; see [GA] `docs/cygl1_operator_test_project_README.md`.
- **IGOT code is private collaborator code** (USC/MiXIL). Do not reproduce IGOT source excerpts in the paper. See [OP]
  `README.IGOTCodeProvenance.md`.

## Suggested reading order

1. [GA] `docs/cygl1_operator_test_project_README.md`: the whole project on one page, with the experiment table, the bug and fix,
   and pointers.
2. [OP] `CYGNSS_IGOT_GEOSldas_Report.md` + `docs/cygnss_operator_forward_model_explainer.md`: what the operator is.
3. [GA] `docs/cygl1_fixedop_2020_arm_comparison_report.md`: the main DA results (O-F), 2020–2022.
4. [GA] `runs/cygl1_fixedop_ismn_validation_2020_2022.md`: in-situ validation, and how the O-F and in-situ results fit together.
5. [GA] `docs/cygl1_obs_quality_vs_qc_report.md` and `docs/cygl1_obs_error_correlation_report.md`: why the final configuration
   (coherency filter, xcorr 0.15°) is what it is.

## Possible paper structure and sources

| Section | Main sources |
|---|---|
| Introduction / motivation | [OP] `README.CYGNSSOperatorDesign.md`, `README.CYGNSSPrototypeScope.md`; prior CYGNSS L3 DA paper (Fox, Reichle & Liu 2026, JHM, submitted; [GA] `projects/cygnss_da/notebooks/Readme.txt`); refs in [OP] `docs/*.pdf` |
| Forward model and tile-coefficient operator | [OP] `CYGNSS_IGOT_GEOSldas_Report.md`, `docs/cygnss_operator_forward_model_explainer.md`, `README.CYGNSSTileCoefficientOperator.md`, `README.IGOTForwardOperator.md`, `README.CYGNSSVegetation.md`, `README.CYGNSSCompletePythonOperator.md`; schematic figure: [GA] `notebooks/cygl1_operator_story_figure.ipynb` |
| Observations, QC, preprocessing | [OP] `README.CYGNSSObservationInputs.md`, `README.CYGNSSObservationQC.md`, `README.CYGNSSCoefficientPreprocessor.md`, `docs/cygnss_m36_window_thinning_notes.md`, `docs/discover_cygnss_preprocessor_quickstart.md` |
| GEOSldas implementation | [OP] `docs/geosldas_cygnss_preprocessed_operator_handoff.md`, `docs/cygnss_preprocessed_operator_validation_20191101.md`, `README.CYGNSSFortranPortPlan.md`, `fortran/README.md`; GEOSldas source (see "Outside these repos") |
| Experiment design (domain, OL, arms, scaling, scoring) | [GA] `docs/cygl1_operator_test_project_README.md`, `runs/configs/`, `docs/cygl1_fixedop_2020_arm_comparison_report.md` (Arms, scoring rule), `projects/obs_scaling_params/README.md` |
| Observation errors and quality control | [GA] `docs/cygl1_obs_quality_vs_qc_report.md`, `docs/cygl1_obs_error_correlation_report.md` |
| Results: O-F (independent monitors) | [GA] `docs/cygl1_fixedop_2020_arm_comparison_report.md` |
| Results: in situ | [GA] `runs/cygl1_fixedop_ismn_validation_2020_2022.md`, `notebooks/cygl1_ismn_skill_figures.ipynb`, `projects/ascat_da/report/ismn_insitu_validation_methods.md` |
| Discussion (L1 vs L3, spring failure, coverage, R) | the two results documents above, plus the observation-error documents |

## [OP] CYGNSS_operator

All of these are tracked in git. Several contain absolute paths on the author's laptop (`/Users/amfox/Desktop/...`); the content is
unaffected.

### Overview and reports

| Path | Date | Content | Status |
|---|---|---|---|
| `README.md` | 2026-06-19 | Repo overview and "start here" pointers. | METHOD |
| `CYGNSS_IGOT_GEOSldas_Report.md` (also `.docx`; PDF in `docs/`) | 2026-05-12 | 1,300-line prototype report: IGOT → tile-coefficient operator `H(x)=Σ_t C_t R_t(x)`, validation against the full-pixel reference, performance, vegetation. The most complete write-up of the operator. It predates GEOSldas integration and the DA runs. | METHOD |
| `CYGNSS_tile_coefficient_operator_status.pptx` | 2026-05 | 18-slide status deck: what C_t is, workflow, QC, validation against the slow reference, M09 vs M36, vegetation folded into C_t, where an obs can update. | METHOD |
| `CYGNSS_tile_coefficient_operator_status_ensemble_update.pptx` | 2026-09-15 | **Byte-identical copy** of the deck above. | duplicate |
| `docs/cygnss_operator_forward_model_explainer.md` | 2026-09-04 | From the IGOT DDM equation to `H(x)=Σ_t C_t R_t`, step by step, with a worked example (sample 100799). Good basis for the methods section and a schematic figure. | METHOD |

### Forward model and operator design

| Path | Content | Status |
|---|---|---|
| `README.IGOTForwardOperator.md` | The IGOT forward model framed as an observation operator (DDM, BRCS, reflectivity). | METHOD |
| `README.CYGNSSTileCoefficientOperator.md` | Theory of the tile-coefficient reduction, the coefficient formula, and comparison with the per-pixel reference. | METHOD |
| `README.CYGNSSCompletePythonOperator.md` | The complete Python reference operator vs the fast coefficient path. | METHOD |
| `README.CYGNSSOperatorDesign.md` | Design: a CYGNSS operator analogous to the GEOSldas mwRTM Tb operator. | METHOD |
| `README.CYGNSSOperatorAPI.md` | Proposed operator API (`cygIGOT_get_obs(cache, state, params)`). | METHOD |
| `README.CYGNSSPrototypeScope.md` | Scope and open design decisions. | METHOD / motivation |
| `README.CYGNSSGeometryComplexity.md` | Geometry cost of CYGNSS/IGOT vs mwRTM; what to precompute per obs. | METHOD |
| `README.CYGNSSVegetation.md` | Vegetation treatment (opacity folded into C_t). | METHOD |
| `README.TbForwardOperator.md` | Reference description of the GEOSldas L-band Tb (tau-omega) operator, which is the architectural template. | BACKGROUND |
| `README.CYGNSSFortranPortPlan.md`, `fortran/README.md` | Fortran port plan and the first Fortran scalar operator (`cygIGOT_get_obs`). | METHOD |
| `README.CYGNSSCodeGuide.md`, `cygnss_operator/README.md` | Code guide for the Python package. | METHOD (code reference) |
| `README.IGOTCodeProvenance.md` | IGOT code is private USC/MiXIL code: sharing restrictions. | BACKGROUND (**constraint**) |
| `README.GEOSldasRestartSFMC.md` | Utility to derive SFMC from Catchment restarts for operator development. | METHOD (minor) |
| `scripts/plot_cygnss_good_obs_operator_story.py` | 3×3 single-obs "operator anatomy" figure at tile level. Its data loaders are reused by the paper schematic notebook ([GA] `notebooks/cygl1_operator_story_figure.ipynb`, see the [GA] section). Its original 2019-11-01 12z output (`artifacts/out_images/cygnss_good_obs_operator_story_20191101_1200z/`) is PRE-FIX. | METHOD (code) |
| `artifacts/cygl1_story_figure_bundle/` (untracked) | Discover bundle for the schematic (13 MB), with `MANIFEST.txt` listing sources, md5s and the obs selection. Built from `docs/discover_operator_story_figure_handoff.md` ([GA]). | CURRENT (inputs) |

### Observations, QC, preprocessing, thinning

| Path | Date | Content | Status |
|---|---|---|---|
| `README.CYGNSSObservationInputs.md` | 2026-05-12 | Which CYGNSS L1 v3.2 variables the operator uses. | METHOD |
| `README.CYGNSSObservationQC.md` | 2026-07-02 | L1 v3.2 land QC strategy (quality flags, land confidence, coherency, SRTM, water fraction). | METHOD |
| `README.CYGNSSCoefficientPreprocessor.md` | 2026-06-24 | The preprocessor that writes the sparse coefficient NetCDF product, and the GEOSldas-side simulator. | METHOD |
| `README.CYGNSSCoefficientBuildPerformance.md` | 2026-05-13 | Build and evaluation cost (about 22 s/obs to build, about 0.04–0.09 ms/obs/member to evaluate). | METHOD |
| `docs/cygnss_m36_window_thinning_notes.md` | 2026-07-02 | Best-obs-per-M36-tile-per-window selection over SW CONUS. | METHOD |
| `docs/discover_cygnss_preprocessor_quickstart.md` | 2026-06-21 | How the Discover production pipeline is run. The pipeline region is **Arizona + 200 km**, which sets the L1 coverage. | METHOD |
| `docs/cygnss_copernicus_dem_cache_merge_notes.md` | 2026-06-19 | Copernicus DEM gradient cache (operational detail). | METHOD (minor) |

### GEOSldas integration and validation

| Path | Date | Content | Status |
|---|---|---|---|
| `docs/geosldas_cygnss_preprocessed_operator_handoff.md` | 2026-06-24 | How the coefficient product is read and H(x) is computed inside GEOSldas `get_obs_pred()`. | METHOD (the implementation was later bug-fixed, see the operator bug) |
| `docs/cygnss_preprocessed_operator_validation_20191101.md` | 2026-05-22 | First check of the GEOSldas Fortran operator against the Python simulator, `H(x)=10 log10(Σ C_t R_t)`. | METHOD (single-case validation) |
| `docs/obs_error_variance_diagnosis.md` | 2026-08-20 | R vs forecast-error variance (fcstvar) decomposition from ObsFcstAna, on the old `OLv8_M36_all_sensors_AZ`. | PRE-FIX (the method is valid: `fcstvar` = P^f_yy) |
| `artifacts/data/geosldas/tile_quality_da_maps/README.md` | 2026-08-26 | Per-tile L1 quality vs DA performance data (SNR predicts DA performance). | PRE-FIX |

### Reference papers (`docs/`)

| Path | Content |
|---|---|
| `docs/Campbell et al. - 2020 - Modeling the Effects of Topography on Delay-Doppler Maps.pdf` | Topography effects on DDMs. |
| `docs/Melebari et al. - 2023 - Improved Geometric Optics with Topography (IGOT) Model ....pdf` | The IGOT model paper. |
| `docs/forward_model_inputs.pdf` | Forward model inputs summary. |

## [GA] geosldas-analysis, `projects/CYGNSS_L1_AZ/`

### Post-fix results and methods (CURRENT)

| Path | Date | Content |
|---|---|---|
| `docs/cygl1_operator_test_project_README.md` | 2026-09-28 | **Project overview:** bug and fix, build, table of all 9 fixed-operator experiments, layout, templates, scoring rule, gotchas. A copy of the project directory's README. |
| `docs/cygl1_fixedop_2020_arm_comparison_report.md` | 2026-09-27/28 | **Main O-F results.** 7 arms in 2020 and 4 arms extended through 2022, vs a single unscaled OL with cross-masked scaled obs. Monthly, yearly and pooled tables; noise-vs-gain split (corr, α_opt). Best L1 arm: coherency filter + errstd 3.9. The spring degradation recurs every year, with a moving window and two mechanisms (signal collapse vs overshoot). Revised next steps. |
| `runs/cygl1_fixedop_ismn_validation_2020_2022.md` | 2026-09-28 | **In-situ validation** over the L1 area (tiles with ≥100 L1 obs: 42 surface / 33 root-zone stations on 25 / 19 tiles), tile-cluster bootstrap. L1 improves root-zone R by +0.03 to +0.04 and beats L3; L3 is neutral. Bias/ubRMSE/RMSE. O-F restricted to the same tiles (L1 and L3 tie on Tb; L3's lead is ASCAT plus its own fit). SOILSCAPE per-tile section. L1 coverage limit (Arizona + 200 km). Appendix: all-station numbers (misleading). |
| `docs/cygl1_obs_quality_vs_qc_report.md` | 2026-09-26 | L1 obs quality vs raw-granule QC fields (Jan–Jun 2020, fixed build). A U-shaped dependence on coherency_ratio motivates the two-sided 0.40–2.16 filter. |
| `docs/cygl1_obs_error_correlation_report.md` | 2026-09-26 | L1 obs-error spatial correlation from the fixed-operator OL (2020–2022). The correlation is mostly along-track and falls to zero beyond about 0.15°. xcorr 0.625° makes R near-singular (condition number about 1e7). This is why xcorr is 0.15°. |
| `docs/cygl1_fixedop_2020_2022_stats_README.md` | 2026-09-28 | README for the O-F stats bundle (pkl/nc4 formats, cross-mask convention, headline table). |
| `docs/cygl1_fixedop_2020_7arm_stats_README.md` | 2026-09-27 | README for the earlier 2020-only bundle (superseded by the one above). |
| `notebooks/cygl1_ismn_skill_figures.ipynb` | 2026-09-28 | CYGNSS-paper-style in-situ figures (raw skill, Δ vs OL, Δ vs L3; paired-t CIs). |
| `notebooks/cygl1_operator_story_figure.ipynb` | 2026-09-29 | **Methods schematic** (4×3). One assimilated L1 obs (cyg02 sample 163383 ch 3, 2020-11-16 00z window, SW Arizona) traced from native ~30 m IGOT pixels (σ factor, footprint, pixel contribution to H(x)) to M36 tile coefficients and the increment, in both candidate L1 arms. It checks the recomputed C_t against the product. Needs the `cygnss` env and the bundle in [OP] `artifacts/cygl1_story_figure_bundle/`. Output: `output/operator_story_20201116_0000z/`. |
| `docs/discover_operator_story_figure_handoff.md` | 2026-09-29 | Discover task brief for the schematic: obs selection criteria and bundle contents. |
| `paper/paper_outline.md` | 2026-09-29 | Paper outline: storyline, section key points with numbers and sources, figure/table list with status, open decisions. |
| `runs/configs/` | 2026-09-28 | exeinp, bat_inp and namelist templates for every fixed-operator experiment (reproducibility / supplement). |

### Pre-fix documents (numbers invalid; history, methods and ideas only)

| Path | Date | Content |
|---|---|---|
| `README.md` | 2026-09-15 | Old project README (monitoring run `OLv8_M36_all_sensors_AZ`). Superseded by `docs/cygl1_operator_test_project_README.md`. |
| `runs/OLv8_M36_all_sensors_AZ.md` | 2026-07-25 | The original 13-species monitoring OL (reader validation). |
| `runs/cygl1_assim_R_sweep.md` | 2026-08-20 | R sweep (full / half / quarter R). |
| `runs/cygl1_operator_diagnosis.md` | 2026-08-20 | L1 vs L3: intrinsic noise vs crude operator. |
| `runs/cygl1_coherency_stratification.md` | 2026-08-25 | coherency_ratio as an L1 fit predictor (the idea was later confirmed post-fix). |
| `runs/cygl1_thinning_and_localization_summary.md` | 2026-08-25 | Thinning vs localization for dense L1. |
| `runs/weekly_update_2026-08-27.md` | 2026-08-27 | Weekly update: paired thinning, hard gate, coherency screening. |
| `runs/cygl1_coh05_density_spectrum_and_24mo_result.md` | 2026-09-08 | Coherency multi-seed and the 24-month dense075_coh05 "best" result. Shown to be an operator artifact. |
| `dense075_coh05_provenance_and_methods.md` | 2026-09-15 | End-to-end provenance of the dense075_coh05 figures. Useful as a template for a provenance/methods trace across both repos. |
| `runs/cygl1_dense075_coh05_ismn_validation.md` | 2026-09-15 | First ISMN check (pre-fix). Superseded by the fixed-operator version. |
| `notebooks/final_dense075_coh05_omf_figures.ipynb`, `coherency_screening_omf_figures.ipynb`, `cygnss_l1_quality_and_error.ipynb` | 2026-08/09 | Pre-fix figure notebooks. Their plotting code is reusable. |

### Related methods elsewhere in geosldas-analysis

| Path | Content | Status |
|---|---|---|
| `projects/ascat_da/report/ismn_insitu_validation_methods.md` | Methods and provenance of the ISMN skill driver (`run_ismn_ol_da_skill.py`) used for the in-situ validation. | METHOD |
| `projects/obs_scaling_params/README.md` (+ `docs/`) | z-score observation scaling climatologies (used for all species; the L1 clim is built with `run_cygnss_l1_scaling_params.py`). | METHOD |
| `projects/cygnss_da/notebooks/Readme.txt` | Data README for the prior CYGNSS **L3** soil-moisture DA paper: Fox, Reichle & Liu (2026), JHM, submitted. | BACKGROUND (prior work to cite and build on) |
| `projects/cygnss_da/notebooks/CYG_insitu_plotter_100325.ipynb` | In-situ figure code of that paper (paired-t CI convention reused here). | METHOD |

## Headline numbers to use (all post-fix)

These are from the CURRENT documents above; check them against the source tables before use.

- **O-F, whole period 2020–2022, % change of O-F stdv vs OL** (negative = better). Filter + errstd 3.9: SMOS −0.97, SMAP −1.20, ASCAT
  +0.55, CYGNSS L3 −2.08, L1 own −0.68. L3-only: −1.14 / −1.48 / −3.61 / −6.01 (own fit). SMAP-only: −14.1 / −11.5 / +12.6 / +2.9.
- **Over the L1 area only,** O-F Tb is a tie: SMAP −1.55 / −1.66 (L1 arms) vs −1.76 (L3).
- **In situ over the L1 area** (vs OL, 95% CI excludes 0 unless marked n.s.): L1 root-zone R +0.030 to +0.044, root-zone anomR +0.019 to
  +0.027, root-zone ubRMSE −0.4 to −0.5 × 10⁻³ m³/m³. L3 is n.s. everywhere. L1 − L3: surface anomR +0.012, root-zone R +0.032 (errstd 2.75).
  SMAP-only: surface anomR +0.088, root-zone R +0.066.
- **The spring degradation recurs every year** (May–Jun 2020, Apr–Jun 2021, Feb–May 2022). Signal collapse in 2020/2022 (corr < 0.1,
  α_opt ≈ 0); overshoot in 2021 (α_opt ≈ 0.45–0.75).

## Outside these two repos (for completeness)

- **GEOSldas operator source (Fortran):** GEOSldas_GridComp, branch `feature/amfox/cygnss-ascat-hsaf-v8`,
  `GEOSlandassim_GridComp/cygnss_preprocessed_obs.F90` (operator) and `clsm_ensupd_read_obs.F90` (reader, `read_obs_cygnss_l1_scalar`).
  The bug-fix commits are `cdfaaf4`, `11dfdb1`. The Discover checkout is `/gpfsm/dnb34/amfox/GEOSldas_cygnss_operator/GEOSldas/`.
- **Data bundles on Discover** (`/discover/nobackup/projects/land_da/cygl1_operator_test/`):
  - `cygl1_fixedop_2020_2022_stats.tar.gz`: O-F pkl/nc4, 2020–2022.
  - `cygl1_ismn_fixedop_2020_2022_inputs.tar.gz`: ISMN station skill CSV + per-tile L1 counts, for the figures notebook.
- **Pre-fix docs archived as tarballs** in the project directory's `archive/`. Not needed.
