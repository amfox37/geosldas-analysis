# CYGNSS L1 assimilation: observation-error (R) sweep

Four-experiment comparison testing how the analysis responds to the assumed
CYGNSS L1 observation-error variance, over calendar 2020 on the AZ box.

Analysed 2026-08-20. **The runs themselves are undocumented** — see
[Provenance gaps](#provenance-gaps) before building on this.

## Experiments

Stats products under `output/stats_output/` (git-ignored):

| label | stats tag | CYGNSS L1 `errstd` | R (variance) | R / R_full |
|---|---|---:|---:|---:|
| OL baseline | `OL_scaled_baseline` | — | — | — |
| full R | `DA_cygl1assim_fullR` | 4.4 dB | 19.36 | 1 |
| half R | `DA_cygl1assim_halfR` | 2.2 dB | 4.84 | 1/4 |
| quarter R | `DA_cygl1assim_quarterR` | 1.1 dB | 1.21 | 1/16 |

**The labels describe the standard deviation, not the variance.** `errstd` in
`obs_param_nml(56)` is the observation error standard deviation ("default obs
error std (before scaling)"), so halving and quartering it scales the
Kalman-gain-relevant variance by 1/4 and 1/16. The step from half R to quarter
R is a much larger change in gain than the naming suggests.

Each has a tile-space `temporal_stats_*_20200101_20201231.nc4` (909 tiles x 13
species) and a domain-aggregated `spatial_stats_*_202001_202012.pkl`
(12 months x 13 species).

CYGNSS L1 (species 12) is the only assimilated species. The other twelve are
monitored, so their O-F and O-A are independent checks on the analysis.

## Species and domain

909 tiles, roughly 117.4-106.2 W and 29.3-39.6 N, same tile set as the
`OLv8_M36_all_sensors_AZ` monitoring run.

| index0 | species | units | errstd (monitoring config) |
|---:|---|---|---:|
| 0-3 | SMOS Tbh/Tbv asc/desc | K | 4.0 |
| 4-7 | SMAP L1C Tbh/Tbv asc/desc | K | 4.0 |
| 8-10 | ASCAT HSAF MetOp-A/B/C | % | 9.0 |
| 11 | `CYGNSS_SM_6hr` (L3) | m3 m-3 | 0.04 |
| 12 | `CYGNSS_L1_DDM3X5_CROP_SCALAR` | dB | 3.0 |

The DA arms use `errstd` 4.4 / 2.2 / 1.1 dB for species 12, not the 3.0 dB of
the monitoring configuration.

Observation counts agree across the four experiments to within **0.05%**, so
all comparisons below are like-for-like. Mean observations per day and the
number of tiles with data: SMOS 0.60 / 794, SMAP 1.47 / 841, ASCAT 1.23 / 581,
CYGNSS L3 0.80 / 677, CYGNSS L1 0.82 / 573.

## Headline result

**Assimilating CYGNSS L1 degrades the fit to every independent observing
system, and the degradation scales with how much weight the observation is
given.** Change in O-F standard deviation relative to the OL baseline:

| group | full R | half R | quarter R |
|---|---:|---:|---:|
| SMOS | +3.9% | +6.8% | +10.5% |
| SMAP | +4.5% | +8.0% | +12.6% |
| ASCAT | +5.8% | +13.3% | +23.1% |
| CYGNSS L3 | +2.1% | +5.9% | +11.2% |
| CYGNSS L1 (assimilated) | +0.2% | +0.3% | +0.6% |

Nothing improves at any R. The relationship is a clean power law
(`D = a * r^-c`) in the true variance ratio, with exponents 0.35 (SMOS), 0.37
(SMAP), 0.50 (ASCAT), 0.60 (CYGNSS L3). Extrapolated:

| group | variance x1.5 | errstd x1.5 (variance x2.25) |
|---|---:|---:|
| SMOS | 3.5% | 3.0% |
| SMAP | 3.9% | 3.4% |
| ASCAT | 5.0% | 4.1% |
| CYGNSS L3 | 1.8% | 1.4% |

Still degrading under either reading of "increase R by 50%".

The damage is seasonal: SMAP peaks at **+31% in May** at quarter R against
about +5% in January-March; ASCAT peaks near +45% in October. Annual means
understate the worst months substantially.

## The assimilation is working; the observation is uninformative

Three independent lines of evidence.

**The analysis does what it is asked to.** O-A relative to O-F for CYGNSS L1 is
-2.9% / -5.9% / -9.7% across the sweep — monotonic, correctly signed, scaling
with the Kalman gain. The pull is spatially uniform across the whole box. In
the OL baseline O-A equals O-F exactly, confirming the fields are read
correctly.

**The same analysis step pushes every other sensor away.** At quarter R, O-A
exceeds O-F by +1.9% (SMOS), +2.4% (SMAP), +4.6% (ASCAT), +6.7% (CYGNSS L3).
This is instantaneous, at the update, not a forecast-model drift. The two soil
moisture retrievals suffer more than the brightness temperatures, consistent
with the increment acting on soil moisture directly while Tb is one RTM step
removed.

**The increments add variance, not bias correction.** The forecast *mean* is
almost identical across all four experiments, while O-F *stdv* degrades by up
to 30% in the same months. Meanwhile the biases that exist go uncorrected:
SMOS and SMAP forecasts run 2-3 K cold from April to October, and CYGNSS L1's
own forecast is high by about 0.6 dB in January-April, a seasonal offset that
changes sign later in the year.

## Error budget

Normalised O-F standard deviation (1.0 = innovation variance consistent with
the assumed error budget):

| group | OL | full R | half R | quarter R |
|---|---:|---:|---:|---:|
| SMOS | 0.721 | 0.761 | 0.814 | 0.874 |
| SMAP | 0.650 | 0.689 | 0.742 | 0.803 |
| ASCAT | 0.614 | 0.650 | 0.732 | 0.836 |
| CYGNSS L3 | 0.400 | 0.411 | 0.435 | 0.465 |
| CYGNSS L1 | 0.736 | 0.743 | **1.355** | **2.218** |

Every species sits below 1.0 in the baseline, so the error specification is
loose across the board, worst for CYGNSS L3.

At full R the assumed observation error variance alone (19.36) **exceeds the
entire measured innovation variance** (17.40, which already includes background
error). That is impossible if R is correct, so R is over-specified at full R
without needing any further assumption. The `errstd` giving a consistent budget
is about **3.0 dB** — which is what the monitoring run uses. The sweep therefore
brackets the consistent value at 4.4 (too large), 2.2 and 1.1 (too small), and
the *least* damaging arm is the most over-specified one.

**This is the opposite direction from what skill wants.** When the
innovation-consistency and independent-skill criteria disagree, R is not the
parameter at fault: it indicates a systematic component in the innovations that
the EnKF is treating as state error.

### Ensemble spread did not collapse

These stats files carry **no direct forecast-error information**. `F_stdv` and
`A_stdv` are temporal standard deviations over 2020, not ensemble spread —
`O_stdv` is identical to four decimals across all four experiments, and
`F_stdv` is the same order as `O_stdv`. The only handle on background error is
the normalised innovation, and it is indirect.

Implied `HPH` from `(OmF_stdv / OmF_norm_stdv)^2 - R`: -1.96 (full R), +0.42
(half R), +0.76 (quarter R). Small at every setting and slightly *rising* as R
falls, so there is no evidence of spread collapse under stronger assimilation.
The mildly negative full-R value is about 10% of R and is within what this
statistic can produce, since `OmF_norm_stdv` is the standard deviation of a
ratio rather than a ratio of standard deviations.

An earlier version of this note claimed spread collapse. That was an artefact
of assuming the sweep scaled variance rather than standard deviation, which
forced `HPH` negative at two of the three points.

## Why it fails: well scaled, but weakly informative

CYGNSS L1 is **not** poorly scaled. Checked directly from the open-loop first
and second moments — z-score scaling forces both `O_mean -> F_mean` and
`O_stdv -> F_stdv`, so a scaled observation should match its forecast in both:

| group | O-F mean / O_stdv | O_stdv / F_stdv |
|---|---:|---:|
| SMOS (4 species) | 8.7-14.2% | 1.08-1.21 |
| SMAP (4 species) | 6.4-13.4% | 1.06-1.19 |
| ASCAT (3 species) | 1.7-19.8% | 0.87-0.88 |
| CYGNSS L3 | 8.2% | 0.911 |
| **CYGNSS L1** | **5.0%** | **0.986** |

CYGNSS L1 has the smallest normalised bias of all thirteen species and a
variance ratio of 0.986, the closest to unity in the set. It is the
best-matched observation in the system.

**The problem is information content.** Recovering the observation-forecast
correlation from `var(O-F) = var(O) + var(F) - 2 cov`:

| group | implied r |
|---|---:|
| SMAP | 0.933 |
| SMOS | 0.922 |
| CYGNSS L3 | 0.658 |
| ASCAT | 0.639 |
| **CYGNSS L1** | **0.352** |

CYGNSS L1's innovation spread (3.095 dB) **exceeds both** its own variability
(2.701) and the forecast's (2.738), which only occurs when observation and
forecast are close to independent.

So the observation is correctly scaled but carries little information about
model error. Its innovations are correctly-sized noise, and the analysis
injects that noise into soil moisture. This accounts for every symptom
together: the forecast mean barely moves because the noise is zero-mean, the
variance degrades because noise adds variance, every independent sensor
registers it at the analysis step, and the damage scales with the Kalman gain.

**No value of R fixes this.** R controls how hard an increment is applied, not
its signal-to-noise ratio. The four other species correlate at 0.64-0.93 with
the same model over the same tiles and period, so this is specific to the
observation rather than a model representativeness problem.

Caveat: correlation against the forecast is not correlation against truth, and
these are aggregate temporal statistics, not paired time series. A tile-level
or event-level analysis would be needed to characterise where and when the
observation does carry signal.

Note on the earlier hypothesis: species 12 has empty scaling file paths in the
monitoring run's `obsparam` while other species point at z-score products. That
difference is real but does not explain the result, and the moment-matching
above indicates the observation reaching the analysis is scaled.

## Provenance gaps

Nothing about these runs is recorded anywhere, and the stats files carry no
global attributes. To close:

1. **The DA runs' own `ldas_obsparam.txt`.** Needed to confirm `scale`, the
   scaling file paths, and the actual `errstd` used in each arm. This is the
   single most valuable missing item.
2. **Discover experiment paths and run dates** for all four experiments.
3. Whether the OL baseline is the same model configuration as the DA arms
   apart from assimilation being switched off.

Resolved 2026-08-20: the `errstd` values and the fact that the sweep scales
standard deviation rather than variance were confirmed from the DA runs' nml
(`obs_param_nml(56)`).

Adding the DA obsparam files to `OLv8_M36_all_sensors_AZ_describe/`, or a
sibling describe directory, would match the provenance convention already used
for the monitoring run.

## Reproducing the figures

```bash
python projects/CYGNSS_L1_AZ/scripts/plot_cygl1_R_sweep.py [--nmin 20]
```

Reads `output/stats_output/` and the monitoring run's tilecoord, writes to
`output/R_sweep_figures/`:

- `R_sweep_fig01_omf_stdv_maps.png` — O-F stdv, open loop plus each DA run minus OL
- `R_sweep_fig02_omf_stdv_monthly.png` — monthly O-F stdv, absolute and percent change
- `R_sweep_fig03_observation_counts.png` — observations per day, maps and monthly totals
- `R_sweep_fig04_mean_obs_and_forecast.png` — mean observation against forecast mean
- `R_sweep_fig05_oma_vs_omf.png` — O-A against O-F, the analysis-impact check

Conventions follow `plot_az_omf_summary.py`: species grouped into observing
systems, `N_data`-weighted statistics, tiles below `--nmin` excluded from
weighted means, maps drawn with `pcolormesh` on the native EASEv2 M36 grid.
