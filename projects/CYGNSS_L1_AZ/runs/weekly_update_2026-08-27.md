# CYGNSS L1 AZ weekly update notes

Prepared 2026-08-27. Scope: this week's update excluding the CYGNSS L1 R-sweep
material, which was already presented last week.

## Rule-of-three summary

1. Extended and evaluated the paired thinning-density CYGNSS L1 experiment.
2. Diagnosed the full-density failure mode and validated the hard-gate fix.
3. Tested coherency-based observation-quality screening against a random
   count-matched control.

## 1. Paired thinning-density experiment

The paired thinning-density experiment compares an open-loop control against
CYGNSS L1 DA arms with increasing local observation density. The key distinction
is whether the L1 observations are kept isolated/sparse, allowed limited local
overlap/intermediate, or assimilated at full/unthinned density.

The useful report result this week is the 22-month extension for the
intermediate-density arm. The exact-population cross-mask comparison gives a
clean CygL1 own-fit improvement without broad degradation of independent
verification species. In the 2020-2021 xmask summary, CygL1 OmF StDev improves
by about 9.2 percent while SMOS, SMAP, ASCAT, and CYGNSS L3 are essentially
flat.

Recommended slide bullets:

- Intermediate-density CYGNSS L1 remains stable when extended to the full
  Jan 2020-Dec 2021 record.
- The exact matched-population comparison shows about 9 percent CygL1 own-fit
  improvement.
- Independent verification species remain near neutral, so the result is not a
  short-window artifact.

Supporting files:

- `projects/CYGNSS_L1_AZ/output/thinning_expts/figures/intermediate_vs_ol_xmask_202001_202112_omf_stdv_percent_maps.png`
- `projects/CYGNSS_L1_AZ/output/thinning_expts/figures/intermediate_xmask_202001_202112_omf_oma_monthly_timeseries.png`
- `projects/CYGNSS_L1_AZ/output/thinning_expts/figures/intermediate_vs_ol_xmask_202001_202112_omf_stdv_summary.csv`
- `projects/CYGNSS_L1_AZ/scripts/plot_thinning_omf_stdv_comparison.py`
- `projects/CYGNSS_L1_AZ/scripts/plot_thinning_oma_timeseries.py`

## 2. Hard-gate fix for dense-density harm

The original full-density dense run showed cross-species harm: independent
SMOS, SMAP, ASCAT, and CYGNSS L3 OmF StDev degraded when the unthinned CYGNSS L1
stream was assimilated. That was the failure mode the hard gate was designed to
address.

The hard-gated dense run keeps the full observation stream but applies a much
tighter CYGNSS L1 admission rule around each tile update. In the currently
staged gated-dense summaries, the verification species are near neutral rather
than broadly degraded. CygL1 own-fit is only weakly improved in the full-period
mean, so the honest phrasing is that the hard gate removes active harm; it does
not by itself establish strong CygL1 skill.

Useful numbers:

- Original dense vs OL, CYGNSS L1: about -0.29 percent reduction by the
  thinning summary convention, i.e. slightly worse.
- Original dense vs OL degradation: SMOS -5.4 percent, SMAP -6.1 percent,
  ASCAT -9.8 percent, CYGNSS L3 -3.8 percent reduction by the same convention,
  i.e. all worse.
- Hard-gated dense, common CYGNSS-covered domain: CygL1 about -0.65 percent
  using `(DA - OL) / OL`, while SMOS and SMAP are near zero and ASCAT is modestly
  positive/worse in the lat-masked area-weighted map summary.

Recommended slide bullets:

- Full-density assimilation originally damaged independent verification species,
  showing a real dense-localization failure mode.
- The hard gate keeps full density but limits which L1 observations jointly
  affect each tile update.
- The gated dense run is best described as no longer actively harmful; CygL1
  own-fit is roughly neutral to slightly improved.

Supporting files:

- `projects/CYGNSS_L1_AZ/output/thinning_expts/figures/dense_vs_ol_omf_stdv_percent_maps.png`
- `projects/CYGNSS_L1_AZ/output/thinning_expts/figures/dense_vs_ol_omf_stdv_summary.csv`
- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/gated_dense_omf_stdv_maps_5x3.png`
- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/omf_stdv_percent_maps_area_weighted_means_lat_lt_37p5.csv`

Sign-convention note: the thinning figures use `100 * (OL - DA) / OL`, where
positive means a reduction in OmF StDev. The coherency-screening figures use
`100 * (DA - OL) / OL`, where negative means a reduction in OmF StDev. Do not
compare colors or signs across those two figure families without checking the
axis label.

## 3. Coherency-based quality screening

The quality/error analysis used open-loop innovations to estimate implied
observation error by quality-field bin. Coherency ratio emerged as the strongest
available per-observation quality signal: the lowest coherency bin has roughly
3.5x the implied error variance of the best middle bins in the open loop, which
is about 1.9x in error standard deviation. Coherency also carries a systematic
bias structure, so it is not just a scatter metric.

The DA screening experiment then compared a coherency-screened arm against a
random count-matched control. In the full-period lat-masked summaries,
coherency screening improves CygL1 own-fit more than the random matched control
(-8.0 percent vs -4.2 percent OmF StDev change from matching OL). The
independent Tb and soil-moisture verification groups remain close to the noise
floor, so the result is real but narrow.

Recommended slide bullets:

- Open-loop innovation diagnostics identify low coherency as a high-error
  CYGNSS L1 subset.
- Coherency screening beats one random count-matched draw on CygL1 own-fit.
- The verification species do not show a clear advantage, so the result is
  promising but not yet broad forecast-skill evidence.

Supporting files:

- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/report_coherency_implied_error_and_bias.png`
- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/report_coherency_screening_full_period_outcome.png`
- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/coherency_screened_omf_stdv_maps_5x3.png`
- `projects/CYGNSS_L1_AZ/output/coherency_screening_figures/coherency_randmatch_omf_stdv_maps_5x3.png`
- `projects/CYGNSS_L1_AZ/output/cygnss_da_quality_bundle/cygl1_coherency_screening_experiment_spec.md`
- `projects/CYGNSS_L1_AZ/notebooks/cygnss_l1_quality_and_error.ipynb`
- `projects/CYGNSS_L1_AZ/notebooks/coherency_screening_omf_figures.ipynb`

## Suggested short deck structure

Intro slide:

- Paired thinning-density: intermediate stays stable and useful over 22 months.
- Hard gate: full-density harm is removed, but strong CygL1 skill is not yet
  demonstrated.
- Coherency screening: quality selection improves CygL1 own-fit relative to one
  random matched draw, but verification impact remains near neutral.

Then one slide per topic, each with one primary figure and the three bullets
above.
