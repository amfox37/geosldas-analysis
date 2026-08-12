# ASCAT DA recent work summary: 2026-07-25 to 2026-08-12

## Purpose

This note records the main ASCAT data-assimilation work completed since the
July 24 project-status note. It explains what was built, why each diagnostic
was needed, what interpretation problems were uncovered, what we currently
believe, and which caveats must travel with the figures.

The central project question remains:

**Should GEOSldas replace the legacy EUMETSAT ASCAT BUFR observation stream
with the H SAF H121 CDR stream, and what changes to footprint handling and QC
are justified?**

The work in this period moved the project from completed production runs and
first-pass O-F evidence into a broader evaluation using:

- full-period O-F diagnostics;
- instrument-variable (IV) and triple-collocation (TC) estimates;
- independent ISMN surface and root-zone measurements;
- controlled H121 footprint and peat-QC experiments; and
- a raw-OFA case study of individual observations at peat boundaries.

The main lesson has been methodological as much as scientific: observation
sampling, scaling, and units have to be made explicit before differences in
O, F, or O-F can be interpreted.

## Executive summary

1. **Both ASCAT DA experiments improve the model over OL.** This is supported
   by independent SMAP O-F diagnostics, IV/TC results, and paired ISMN skill.
2. **H121 is not uniformly superior to legacy in every independent metric.**
   H121 has a small advantage in independent SMAP O-F spread and in ISMN
   root-zone correlations, while legacy has a small advantage in ISMN surface
   correlations and in some IV/TC combinations. Many H121-versus-legacy
   differences are much smaller than either run's improvement over OL.
3. **Diagnostics against an assimilated observation family are useful but not
   independent.** H121 DA fits H121 observations especially well, but that is
   expected and cannot by itself decide between observation streams.
4. **The original OL legacy/H121 O and F comparison mixed different temporal
   samples.** The two OL summary products used the same OL trajectory but were
   independently cross-masked against different DA experiments. The resulting
   forecast means could differ substantially even though the underlying OL run
   was identical.
5. **A new jointly matched OL dataset resolved that sampling ambiguity.** On
   identical tile/cycle support, H121-minus-legacy F is essentially zero, while
   the raw observation difference remains large and strongly latitude
   dependent.
6. **The jointly matched dataset has different units from the DA-matched OL
   products.** Its O values are raw ASCAT wetness, not scaled soil moisture.
   Consequently, O-F in that product is a mixed-unit diagnostic and must not be
   described as a physical innovation bias.
7. **The peat-QC boundary count pattern is real, not a plotting artifact.** A
   reconstruction from raw OFA records and the exact GEOSldas porosity/footprint
   logic predicts all 2,720 sampled retain/reject decisions correctly.
8. **The low-but-nonzero boundary population has a measurable, non-arbitrary
   breakpoint.** Among partially retained tiles with at least 20 baseline
   observations, the low-retention cluster ends near 7% of the FOV12.5 count.
   A 10% retained-count threshold lies in a clear gap and identifies 179
   Canada/Alaska edge tiles and 576 globally.

## Experiment and data inventory

### Full six-year comparison

| Analysis label | Experiment | Role |
| --- | --- | --- |
| `OL` | `OLv7_M36_MULTI_type_13_H121` | No-assimilation baseline |
| `DA_H121` | `DAv7_M36_ASCAT_type_13_H121` | Assimilates H SAF H121 ASCAT |
| `DA_legacy` | `DAv7_M36_ASCAT_type_13_legacy` | Assimilates legacy BUFR ASCAT |
| `DA_SMAP` | `DAv7_M36_SMAP_type_13_comb_fp_scaled` | Independent comparison experiment assimilating SMAP |

The common evaluation period is 2015-04-01 through 2021-03-31. CYGNSS-based
verification begins in 2018-08 because that is when CYGNSS observations become
available.

The main local bundles are:

- [`data/step4_h121_cdr_test_20260725`](../../../data/step4_h121_cdr_test_20260725/README.md):
  final IV and TC output;
- [`data/ismn_ol_da_skill_bundle`](../../../data/ismn_ol_da_skill_bundle/README.md):
  station, network, and inventory tables from the ISMN evaluation; and
- `data/omf_compare_sums`: monthly and full-period ObsFcstAna summary files.

### First-six-month footprint/QC comparison

The footprint/QC analysis uses 2015-04-01 through 2015-09-30:

| Label | Experiment | Difference being tested |
| --- | --- | --- |
| `DA_H121` | original H121 experiment | Existing H121 configuration |
| `DA_baseline_FOV12p5` | H121 FOV12.5 experiment | Wider 12.5 km Gaussian footprint scale |
| `DA_peatlandqc` | H121 FOV12.5 + peat QC | FOV12.5 plus new peat/wetland QC |

The two new experiments share the same model configuration, observation
parameters, scaling climatology, and period. Their intended controlled
difference is the new QC logic. The local provenance is in
[`data/omf_compare_sums/stats_first6mo/README.md`](../../../data/omf_compare_sums/stats_first6mo/README.md).

Here and below, an **observed tile** means a tile with at least one pooled
H SAF Metop-A/B observation in the FOV12.5 experiment. An **impacted tile**
is an observed tile whose pooled count is lower after peat QC.

## Common figure and aggregation conventions

The three full-period figure notebooks were made visually and statistically
consistent:

- global maps stop at 60 degrees south;
- global maps use Robinson projection when Cartopy is available;
- the FOV12.5/peat-QC regional and global detail maps use native M36 EASEv2
  grid-cell `pcolormesh` rendering rather than point markers;
- land is grey where no plotted value is present;
- mapped and global means use GEOSldas tile area, not a plain tile count;
- diverging color bars are segmented and have a white neutral interval around
  zero, normally 2% of the displayed half-range;
- panels carry `(a)`, `(b)`, and subsequent labels;
- figure markdown states what the figure contains and which color/direction is
  favorable;
- titles do not advertise `AW`, although the documented means remain area
  weighted;
- figures are displayed in the notebook and written as PNG/CSV products; PDFs
  are not produced; and
- run labels use `DA SMAP`, not `SMAP + ASCAT`.

Where observations are combined across platforms, the notebooks call the
collection a **species group**, not a family. The observation-support figure
was also changed to a vertical 3-by-1 layout so each global map is large enough
to inspect.

These choices make the decision figures easier to compare and prevent a map
with many small tiles from dominating a global mean merely through tile count.

## 1. Full-period O-F diagnostic notebook

Primary notebook:
[`h121_legacy_omf_figures.ipynb`](../notebooks/h121_legacy_omf_figures.ipynb)

### Why this work was needed

The initial O-F summaries established that assimilation reduces innovation
spread, but they did not make observation support, individual species, or the
legacy/H121 observation differences easy to inspect. The notebook was expanded
into a main figure sequence plus supporting diagnostics so that a headline
result can always be traced back to sampling, platform, O, and F.

### Main comparison logic

Each observation group is compared with the OL/background summary containing
the corresponding observations:

- SMAP uses `OL_vs_SMAPobs`;
- legacy ASCAT uses `OL_vs_legacyobs`; and
- H121 ASCAT uses `OL_vs_H121obs`.

This matters because the OL summary is not a universal interchangeable
baseline. Its observations and retained tile/cycle population are determined
by the cross-mask used to construct it.

O-F standard deviation is the primary fit diagnostic. O-F mean is retained as
context, but scaling tends to force mean innovations toward zero and makes it
less discriminating.

### Main O-F result

Area-weighted full-period summaries show:

| Verification group | DA H121 improvement vs OL | DA legacy improvement vs OL |
| --- | ---: | ---: |
| Independent SMAP Tb O-F stddev | 5.68% | 4.49% |
| Legacy ASCAT O-F stddev | 17.31% | 16.96% |
| H121 ASCAT O-F stddev | 22.21% | 13.30% |

The independent SMAP comparison gives H121 a small advantage. The large H121
advantage against H121 observations is useful confirmation that the H121 DA
path behaves as intended, but it is not independent evidence because H121 DA
assimilates that observation family. Likewise, legacy-observation diagnostics
must be interpreted with their own circularity and sampling context.

### Added diagnostic coverage

The notebook now includes:

- monthly and mapped O-F stddev changes relative to OL;
- absolute monthly O-F stddev and O-F mean;
- observation support maps and monthly counts by observation species group;
- OL O and F time series and maps;
- common-finite-tile maps;
- H121-minus-legacy O and F maps with aligned 5-degree latitude profiles;
- seasonal O differences;
- tile scatter, histogram/CDF, O-F mean, and per-platform diagnostics;
- species-level monthly and full-period O-F improvement; and
- long-form CSV exports for every available mean and standard-deviation metric.

The support diagnostics are deliberately more extensive than the likely final
paper figure set. Their role is to reveal whether a result originates in the
observation product, model trajectory, sampling, or a single Metop platform.

## 2. The OL sampling problem and jointly matched dataset

### What looked wrong

In the original OL support figures, H121-minus-legacy O and F differences were
surprisingly similar, while the corresponding O-F difference was very small.
Because both products were described as coming from the same OL run, it was
tempting to conclude that F should be identical and that the calculation was
wrong.

The fact that a small O-F difference algebraically implies similar O and F
differences was not accepted as proof. That argument is circular and does not
establish that O and F were extracted or sampled correctly.

### Root cause

`OL_vs_legacyobs` and `OL_vs_H121obs` both use the single OL trajectory, but
each is independently cross-masked against its corresponding DA experiment.
They therefore retain different days/cycles at a given tile. Different samples
from the same time-varying forecast trajectory can have different F means.

Restricting maps to the same finite tiles (70,584 cells in the earlier paired
map) only aligns spatial support. It does not reconstruct a common cycle-level
sample, and monthly group summaries cannot be retrospectively cycle-matched
once they have already been aggregated.

Another important clarification followed: O in these DA-matched OL products
comes from the scaled DA observation stream. These are not the original raw
ASCAT wetness observations.

### Corrective dataset

The `OL_legacy_h121_xmask` product was generated from the single OL run alone,
without using a DA run. Legacy/H121 platform pairs contribute only when both
species are present for the same owner tile and assimilation cycle:

- legacy species 5 paired with H121 species 8 (Metop-A);
- legacy species 6 paired with H121 species 9 (Metop-B); and
- legacy species 7 paired with H121 species 10 (Metop-C).

The preprocessor/runtime population is at most one selected observation per
owner tile, species, and window. Matching is therefore performed on that
tile/cycle population rather than by accumulating all raw specular or swath
points.

The pair check is exact: monthly and full-period `N_data` differences are zero
for all three platform pairs.

The generator was developed in the Discover-only `GEOSldas_GridComp` worktree
on branch `feature/amfox/obsfcstana-nc4-postproc`, commit `2fad400`. It adds a
self-cross-masking postprocessor while retaining the existing sums/stats file
format. The generated local files are:

- `data/omf_compare_sums/OL_legacy_h121_xmask/OL_legacy_h121_xmask_monthly_stats.nc4`
- `data/omf_compare_sums/OL_legacy_h121_xmask/OL_legacy_h121_xmask_temporal_stats.nc4`

### What the corrected product shows

On 69,689 common valid tiles:

| H121 minus legacy quantity | Area-weighted mean |
| --- | ---: |
| Raw O mean | -0.05942 raw ASCAT wetness units |
| OL F mean | +0.000294 m3 m-3 |
| Mixed-unit O-F mean | -0.05971 |

The near-zero F difference is the expected result when both observation
families sample the same OL forecast trajectory on identical tile/cycle
support. The raw observation difference remains large. It is especially
negative at northern latitudes, reaching an area-weighted mean near -0.18 in
the 50-55 degree north band.

### Units warning

In `OL_legacy_h121_xmask`:

- O is raw ASCAT wetness;
- F is model soil moisture in `m3 m-3`.

Therefore O-F is not a physical innovation in one common unit. It is kept only
as a bookkeeping diagnostic and is explicitly labelled as mixed-unit. This
dataset is appropriate for comparing raw legacy/H121 O and for verifying that
matched F is nearly identical. It is not a replacement for the scaled,
DA-matched OL baselines used in the primary DA-versus-OL O-F figures.

## 3. Instrument-variable and triple-collocation figures

Primary notebook:
[`h121_iv_tc_skill_figures.ipynb`](../notebooks/h121_iv_tc_skill_figures.ipynb)

The notebook calls the first method **instrument variable (IV)**. The source
bundle's older documentation sometimes calls the stage "independent
validation," but the implemented method is the lagged instrument-variable
framework. The figure text now uses the intended terminology.

### IV method and figures

The key displayed metric is model R, obtained from the reported model IV R2.
R was chosen instead of R2 for readability and consistency with the other
correlation diagnostics. The bundle uses a two-day lag and requires at least
100 paired samples per cell.

The IV figures provide:

- area-weighted model R by verification sensor and run;
- DA-minus-OL R maps; and
- H121-minus-legacy R maps.

The mapped IV comparisons intentionally use only SMAP L3, SMOS-IC, and CYGNSS
L3. The earlier fourth H121-observation panel was removed from the direct
H121-minus-legacy map because it is circular for DA H121.

Headline area-weighted model R values include:

| IV sensor | OL | DA legacy | DA H121 |
| --- | ---: | ---: | ---: |
| SMAP L3 | 0.528 | 0.576 | 0.571 |
| SMOS-IC | 0.593 | 0.636 | 0.636 |
| CYGNSS L3 | 0.641 | 0.654 | 0.641 |
| ASCAT H121 | 0.513 | 0.609 | 0.666 |
| ASCAT H119/H120 | 0.578 | 0.711 | 0.728 |

The independent SMAP/SMOS/CYGNSS rows indicate that both ASCAT DA runs improve
on OL, but they do not show a universal H121 advantage. ASCAT-family rows are
valuable diagnostics but are not fully independent when a run assimilates the
same family. The H119/H120 IV results also retain a documented coverage/
geometry caveat from an older pair-generation path.

### TC method and figures

The TC figures use model fractional mean squared error:

```text
fMSE = sigma2_error / variance
```

Lower is better. Values are filtered to the physically interpretable range
`0 <= fMSE <= 1` before averaging because TC can be unstable when a covariance
denominator approaches zero.

The figures show:

- area-weighted model fMSE by valid triplet and run;
- OL-minus-DA fMSE maps, where positive means assimilation reduces error; and
- legacy-minus-H121 maps where a direct non-circular comparison exists.

Area-weighted means for the two triplets containing both ASCAT DA runs are:

| TC triplet | OL | DA legacy | DA H121 |
| --- | ---: | ---: | ---: |
| CYGNSS L3 + model + SMOS-IC | 0.592 | 0.549 | 0.552 |
| SMAP L3 + model + SMOS-IC | 0.708 | 0.647 | 0.643 |

Both runs improve over OL. Legacy is marginally lower for the CYGNSS triplet;
H121 is marginally lower for the SMAP triplet. Circular combinations are
excluded by design, so the triplet containing ASCAT H121 does not report the
DA H121 run.

## 4. ISMN in-situ evaluation

Primary method note:
[`ismn_insitu_validation_methods.md`](ismn_insitu_validation_methods.md)

Primary figure notebook:
[`h121_ismn_skill_figures.ipynb`](../notebooks/h121_ismn_skill_figures.ipynb)

### Why this work was needed

ISMN provides an independent station-based check of whether the DA experiments
improve surface and root-zone soil moisture, rather than merely fitting the
satellite streams used in assimilation or validation.

A Discover batch workflow was added because the complete problem is too large
for an interactive notebook. It parses the ISMN archive directly, builds
cached daily observation/model series, and scores all four runs over the full
experiment window.

### Station selection and common-sample design

- All available networks are considered rather than using the six-network
  M21C list.
- A station/domain must have at least 1,000 paired days in every run. This
  common-run requirement prevents a DA comparison from changing its station
  population.
- Surface uses the sensor nearest 0.05 m among sensors no deeper than 0.10 m.
- Root zone is a profile-weighted 0-1 m composite. It requires at least three
  finite layers at a timestamp, including a shallow sensor with depth <= 0.20 m
  and a deep sensor with depth >= 0.50 m.
- Only ISMN `G`-flagged observations are retained.
- Observations are daily averaged and shifted 12 hours to align with the model
  day.
- Model `sm_surface` and `sm_rootzone` come from `SMAP_L4_SM_gph`, averaging
  all eight three-hourly values per day for every run.

The generic root-zone method is intentionally different from the older M21C
network-specific strict layer rules. Those rules were hand-tuned for six
networks and do not generalize. The tradeoff is that absolute root-zone
definitions vary somewhat with station sensor layout; paired OL-versus-DA
changes at a fixed station remain much more defensible.

### Metrics and error bars

The station metrics are R, anomaly R, bias, RMSE, and ubRMSE. Anomaly R removes
a day-of-year climatology using a 31-day circular window.

The main figures emphasize R, anomaly R, and ubRMSE. Bias remains available in
the station/network tables but was not added to the mirrored Figure 3/Figure 6
layout.

All plotted deltas are signed so positive means better:

- R and anomaly R: `run - OL`;
- ubRMSE: `OL - run`; and
- H121 advantage: the corresponding H121 improvement minus legacy improvement.

Error bars are approximate 95% confidence intervals around the tile-area-
weighted paired station mean. The variance is calculated from station deltas,
and the standard error uses Kish effective sample size. Pairing means each
station serves as its own OL control.

Figure 6 was revised to mirror Figure 3 exactly: the same networks, surface on
top, root zone on the bottom, and R/anomaly R/ubRMSE in the same column order.
This removed the confusing appearance of a different network population and
restored R to the ranking plot.

### ISMN result

The final area-weighted common-station summaries use 936 surface stations and
578 root-zone stations (927 and 572 respectively for anomaly R):

| Domain and metric | OL | DA legacy | DA H121 | DA SMAP |
| --- | ---: | ---: | ---: | ---: |
| Surface R | 0.6220 | 0.6337 | 0.6299 | 0.6708 |
| Surface anomaly R | 0.4668 | 0.4775 | 0.4713 | 0.5548 |
| Surface ubRMSE | 0.06408 | 0.06290 | 0.06291 | 0.06105 |
| Root-zone R | 0.6437 | 0.6663 | 0.6731 | 0.6766 |
| Root-zone anomaly R | 0.4684 | 0.4993 | 0.5087 | 0.5422 |
| Root-zone ubRMSE | 0.04688 | 0.04598 | 0.04582 | 0.04585 |

All DA runs improve on OL. Legacy is slightly stronger at the surface for R
and anomaly R, while H121 is stronger in the root zone for both correlation
metrics and has a small root-zone ubRMSE advantage. DA SMAP is the strongest
surface experiment and is comparable to H121 in root-zone ubRMSE.

The result is therefore nuanced: ISMN supports the value of ASCAT assimilation
but does not provide a single all-depth verdict that H121 dominates legacy.

## 5. H121 FOV12.5 and peat-QC comparison

Primary notebook:
[`h121_fov_peatlandqc_omf_figures.ipynb`](../notebooks/h121_fov_peatlandqc_omf_figures.ipynb)

### Why this work was needed

The original H121 setup, a wider 12.5 km footprint scale, and the same wider
footprint plus peat QC needed to be compared on a common map set. The immediate
questions were:

- how much does FOV12.5 change observation support and mean sampled O?
- what additional observations does peat QC remove?
- do those changes alter O-F spread for assimilated H SAF observations?
- do independent monitor-only SMAP diagnostics change?
- what is happening in Canada and Alaska, where peat coverage is extensive?

### Figure set

The notebook pools H SAF Metop-A/B using component observation counts and
produces:

- observation-count maps and pairwise differences;
- scaled O-mean maps and pairwise differences;
- O-F standard-deviation maps and pairwise differences;
- experiment-minus-OL relative O-F standard-deviation maps for H SAF and
  monitor-only SMAP; and
- a Canada/Alaska zoom repeating counts, O means, count/O differences, and
  relative O-F spread for FOV12.5 and FOV12.5 + peat QC;
- Canada/Alaska and global maps of complete-loss tiles plus low-but-nonzero
  edge-tile counts; and
- a final 1-by-2 Canada/Alaska SMAP-only relative O-F-spread comparison so the
  independent monitor observations are not visually compressed by the larger
  H SAF differences.

The O-mean and O-F-spread pairwise maps use common finite tiles, but they are
not matched observation by observation. They include both changed sampling
and changed DA trajectories.

For Figure 7 and the regional repeat, `OL_vs_H121obs` is used so H SAF O is on
the scaled basis used by the DA experiment. This reference is exact for the
original H121 experiment and for monitor-only SMAP. It is only approximate for
the FOV12.5 and peat-QC H SAF comparisons because those observation samples are
not cycle-matched to the OL reference.

### Canada/Alaska boundary pattern

The Canada/Alaska count maps showed narrow bands of low-but-nonzero counts at
the edge of peat-QC regions. This initially looked suspicious because a simple
tile mask might be expected to produce an abrupt all-or-nothing boundary.

The QC is not a simple owner-tile mask. It is evaluated for each observation
at its actual footprint center. Observations assigned to the same M36 owner
tile can be shifted within that tile and therefore overlap neighboring peat
tiles with different weights. Most observations can be rejected while a few
remain below the threshold, producing exactly the observed low-count edge.

The final edge-tile definition was derived from the count distribution rather
than from an arbitrary color-bar limit. It requires:

- at least 20 pooled H SAF observations in the FOV12.5 baseline;
- at least one observation retained after peat QC; and
- no more than 10% of the baseline count retained after peat QC.

In Canada/Alaska, the low-retention cluster ends at 6.97% and the next
partially retained tile is at 39.68%. Globally, the cluster ends at 7.04% and
the next tile is at 15.79%. The 10% threshold therefore lies in an empty gap
in both populations. The selected edge tiles retain at most 28 observations,
so their maps use a segmented 0-30 count scale.

The resulting six-month count summary is:

| Domain | Tiles with FOV12.5 observations | Impacted tiles | Impacted fraction | Zero retained among impacted | Edge tiles among impacted | Area-weighted mean count reduction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Canada/Alaska | 9,687 | 1,918 | 19.8% | 84.2% | 179 (9.3%) | 93.3% |
| Global, 60S-85N | 82,045 | 4,590 | 5.6% | 76.9% | 576 (12.5%) | 90.0% |

Grey cells on these maps are a **complete-loss peat-QC impact proxy**: the
FOV12.5 experiment had observations and peat QC retained none. They are not an
independent map of `POROS`, because the temporal-statistics files do not carry
the production porosity field. Colored cells show retained counts only for the
edge population defined above.

## 6. Individual-observation peat-QC case study

Case-study bundle:
[`footprint_case_study/README.md`](../footprint_case_study/README.md)

### Exact weighting logic

For an observation center, candidate model tiles are those whose centers fall
within a 25 km ellipse. For each candidate tile:

```text
unnormalized_weight_i = tile_area_i * exp(-0.5 * (distance_i / 12.5 km)^2)
weight_i = unnormalized_weight_i / sum(unnormalized_weight)
```

The 12.5 km value is the Gaussian distance scale, not a hard boundary. At
12.5 km the distance-only factor is about 0.607. The 25 km value is the hard
candidate cutoff (`2 * FOV`); at that distance the Gaussian factor is about
0.135 before area weighting. Distances are converted to longitude/latitude
semi-axes, which is why the boundaries appear elliptical on a lon/lat plot.
The labelled distances are radii, not diameters.

A tile is classified as PEATCLSM when `POROS >= 0.90`. The footprint peat
fraction is the normalized weighted mean of this binary peat indicator. The
observation is rejected when the footprint peat fraction is at least 0.10.

### Raw-OFA audit

The Discover extraction directly compared:

- `DAv7_M36_ASCAT_type_13_H121_FOV12p5`; and
- `DAv7_M36_ASCAT_type_13_H121_FOV12p5_peatlandqc`.

It opened 1,463 OFA files for each experiment over April-September 2015 and
examined H SAF Metop-A/B observations for six peat-boundary tiles:
`3823`, `4338`, `4813`, `5230`, `6548`, and `9449`.

The baseline contained 2,720 valid candidate observations. Only 30 remained in
the peat-QC experiment; 2,690 were lost. The extraction verified that row
existence, finite O, finite O and F, and `assim_flag == 1` all reproduced the
production `N_data` counts for every tile/species/experiment combination.

The portable Python reconstruction then hard-asserted that:

- every candidate had a reconstructable footprint;
- every reconstructed QC decision matched the actual OFA retain/reject
  outcome; and
- no retained or rejected observation appeared on the wrong side of the 0.10
  threshold.

All 2,720 decisions agreed (100%). This establishes a causal mechanism for the
sampled boundary tiles, not merely a visual correlation.

### Visual products

The committed bundle includes:

- six-panel porosity/observation-center maps;
- the reconstructed threshold separation;
- an actual observation time series for tile 3823 (446 candidates, 5 retained);
  and
- one paired map for each target tile, showing the nearest actual retained and
  rejected observation locations over a POROS pcolormesh.

On the paired maps:

- the solid ellipse marks the 12.5 km Gaussian scale;
- the dashed ellipse marks the 25 km contribution cutoff;
- magenta rings mark contributing tile centers and label normalized weights;
- red outlines identify PEATCLSM display cells; and
- the observation center and owner-tile center are shown separately.

The pcolormesh rectangles are inferred around tile centers for visualization;
they are not exact catchment polygons.

### Scope

The raw audit proves the mechanism for six tiles and six months. It is not a
full-domain reconstruction. The summary-statistics notebook now quantifies
complete-loss and low-retention tiles across Canada/Alaska and globally, but it
does not recompute a footprint peat fraction for every observation in those
domains. Extending the raw-OFA reconstruction would therefore strengthen the
domain-wide attribution, although the six-tile audit is sufficient to show
that the observed boundary mechanism is genuine and follows the implemented
per-observation footprint peat weighting.

## What changed in our interpretation

Several corrections made during this period should be preserved explicitly:

1. **IV terminology:** figure text now says instrument variable, not
   independent validation.
2. **IV display metric:** maps and scorecards use R, not R2.
3. **O in DA-matched OL summaries:** these O values are scaled observations,
   not necessarily the original retrieval values.
4. **O in the jointly matched OL summary:** these O values are raw wetness and
   are not commensurate with F in `m3 m-3`.
5. **Common map cells are not common observations:** tile-level cross-masking
   cannot fix different day/cycle samples.
6. **A small O-F difference is not proof that O and F were computed correctly:**
   raw OFA and exact support matching were needed to establish the result.
7. **Peat QC is observation-footprint based:** owner tiles near a peat boundary
   can contain both retained and rejected observations.
8. **ISMN Figure 3 and Figure 6 must use the same network population and panel
   layout:** otherwise apparent differences can be a plotting-population
   artifact.
9. **Root-zone ISMN is generic rather than network-specific:** this broadens
   coverage but requires caution when comparing absolute skill across networks.
10. **Peat-QC edge tiles are defined by retained fraction, not an absolute
    count:** the final criterion is supported baseline count of at least 20 and
    retention of more than zero but no more than 10%.

## Current scientific assessment

The evidence still supports moving forward with H121 as a credible replacement
for legacy ASCAT, but the recent independent evaluations argue for a measured
rather than absolute conclusion.

Strong evidence in favor of H121 includes:

- a working end-to-end reader/QC/scaling/assimilation path;
- denser observation support than legacy;
- a small improvement over legacy in independent SMAP O-F spread;
- stronger H121-observation fit, as expected for the H121 DA experiment;
- better ISMN root-zone R and anomaly R than legacy; and
- a raw jointly matched comparison that cleanly separates product differences
  from forecast-sampling differences.

Evidence that prevents claiming universal H121 dominance includes:

- slightly better legacy ISMN surface R and anomaly R;
- mixed, small H121-versus-legacy differences across independent IV/TC
  combinations;
- circularity in diagnostics that use the assimilated ASCAT family; and
- sensitivity of O/F summaries to scaling and temporal support.

The robust conclusion is that both ASCAT DA configurations improve OL, and the
choice between H121 and legacy should consider coverage, product continuity,
reader/QC quality, and regional behavior alongside small differences in global
skill means.

## Remaining work

1. Decide which O-F, IV/TC, and ISMN panels form the smallest clear decision
   figure set; retain the larger support suite for audit/provenance.
2. Evaluate whether the FOV12.5 and peat-QC H SAF O-F comparisons require a
   dedicated cycle-matched OL product, analogous to the legacy/H121 matching
   correction.
3. Quantify the FOV12.5 and peat-QC impact beyond the first six months if these
   configurations are candidates for production.
4. Consider extending the peat-boundary raw-OFA reconstruction from six tiles
   to the full edge population if observation-level domain-wide attribution is
   needed; the retained-count population itself is now quantified.
5. Decide whether the 10% footprint peat threshold and 12.5 km Gaussian scale
   give the desired balance between excluding peat-contaminated retrievals and
   preserving useful near-boundary observations.
6. Keep scaled and raw observation products visibly separated in filenames,
   legends, units, and prose.
7. Update `h121_da_current_goals.md` once the FOV/peat-QC decision and final
   H121-versus-legacy recommendation are settled.

## Repository artifacts and commit trail

Key commits in this period:

| Date | Commit | Purpose |
| --- | --- | --- |
| 2026-07-25 | `923c428` | Add the ISMN OL/DA validation workflow and methods |
| 2026-07-25 | `19e7029` | Document H121 comparison status and initial O-F evidence |
| 2026-07-26 | `ea228da` | Add the IV/TC and ISMN figure notebooks |
| 2026-07-26 | `9f68e1f` | Expand the legacy/H121 O-F diagnostics |
| 2026-07-27 | `75be7ce` | Refine O/F support and sampling diagnostics |
| 2026-08-10 | `8841876` | Standardize and correct the H121 skill figure notebooks |
| 2026-08-10 | `adc4d65` | Add jointly matched OL and expanded O/F diagnostics |
| 2026-08-10 | `d3727a3` | Add FOV12.5 and peatland-QC comparison maps |
| 2026-08-10 | `11a55a5` | Add the reproducible raw-OFA peat-boundary case study |

The generated notebook output directories under `projects/ascat_da/output/`
are intentionally ignored by Git. The footprint case study is different: its
small portable inputs and final evidence figures are committed so the QC
reconstruction can be rerun without Discover access. Its regenerated
per-observation result table remains ignored under
`projects/ascat_da/footprint_case_study/output/`.
