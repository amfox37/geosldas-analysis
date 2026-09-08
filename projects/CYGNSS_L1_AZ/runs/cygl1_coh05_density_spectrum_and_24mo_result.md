# CYGNSS L1 AZ: coherency multi-seed close-out, the coh05 density spectrum, and the 24-month dense075_coh05 result

Covers 2026-08-28 through 2026-09-08 — picks up exactly where `runs/weekly_update_2026-08-27.md`
left off. That note's three items (paired thinning-density 22mo, the hard-gate fix, and coherency
screening) were all still open threads; this note closes two of them out and adds a new,
higher-confidence result on top.

## Rule-of-three summary

1. Coherency-screening's Arm-A-vs-Arm-B result (last week's item 3) is now confirmed real and
   seed-robust, not a single-seed artifact — but the verdict itself (narrow, not decisive) is
   unchanged.
2. A new experiment family combines spatial thinning (from item 1) with coherency screening (from
   item 3) across three obs-density tiers. The densest of these, `dense075_coh05`, was extended to
   a full 24-month record and is now the project's most robust positive CYGNSS L1 skill result to
   date.
3. Coherency filtering *alone*, without thinning or the hard gate, still actively harms
   independent verification species at full density — confirming that density reduction
   specifically (not just per-obs quality screening) is what avoids cross-species harm.

## 1. Coherency screening: multi-seed close-out (2026-08-28)

Last week's report flagged Arm A (coherency-screened, `coherency_ratio >= 0.523`) beating Arm B
(random obs-count-matched control) roughly 2x on CygNSS L1's own O-F skill (-8.0% vs -4.2%
`(DA-OL)/OL` stdv) as a single-seed result needing confirmation. Two more independent random-draw
seeds for Arm B were run (full year 2020, same gated binary, same `errstd=3.0dB`) and evaluated.

**Result: the ~2x gap replicates cleanly across all 3 B seeds — a real, seed-robust effect.**

| group | Arm A | B seed1 | B seed2 | B seed3 | B mean (range) |
|---|---|---|---|---|---|
| CygL1 own skill (dB, `(DA-OL)/OL`) | **-8.03%** | -4.20% | -4.66% | -4.51% | -4.46% (-4.20 to -4.66%) |
| Tb monitor (K) | +0.89% | +0.87% | +0.88% | +0.84% | +0.86% (+0.84 to +0.88%) |
| SM monitor (m3/m3) | ~0.00% | ~0.00% | ~0.00% | ~0.00% | ~0.00% |

Event-level diagnostics (gain_proxy/%improved/%toward_obs) are statistically indistinguishable
between all 4 arms — the screen changes *which* obs get admitted, not *how* the update
mechanistically behaves once obs are in.

**Verdict unchanged in kind, strengthened in confidence**: coherency screening gives a real ~2x
CygL1-own-skill edge over random matching, but the co-primary withheld Tb/SM metrics — which the
R-sweep result identified as the more decisive metrics historically — show no practically
meaningful advantage for either arm over OL. This closes the "was it a fluke" question without
changing the "narrow, real but not decisive" verdict.

## 2. Context for what came next: per-tile obs quality predicts DA performance (2026-08-26)

Built the same week as the coherency work above, and part of what motivated it: per-tile (557
AZ-domain tiles) cross-correlation of two candidate CYGNSS L1 obs-quality metrics against three
DA-performance metrics on the gated-dense run (Spearman rho):

| quality metric | vs frac_improved | vs frac_rightdir | vs gain_proxy |
|---|---|---|---|
| obs-fcst correlation `r` | -0.04 | -0.01 | +0.22 |
| SNR (`fcstvar/obsvar`) | **+0.37** | **+0.41** | +0.83 (partly mechanical) |

**SNR predicts tile-level DA performance; obs-fcst correlation does not.** This is a different
metric from `coherency_ratio` (the per-obs quality field used in the screening experiment above),
but the same underlying question. Not yet followed up as its own screening experiment — noted here
for completeness since it was previously undocumented.

## 3. New experiment family: the coherency-filtered density spectrum ("coh05" arms)

The original paired thinning-density experiment (in last week's report) asked "how much local obs
density survives before DA behavior degrades," using pure spatial thinning (sparse/intermediate/
dense). This new family asks the same question with the coherency filter (`coherency_ratio >=
0.5`) from Section 1 applied *on top of* three density tiers, on the **ungated** binary (a
separate pre-`406206a` build, surgically relinked from the gated binary's own object files — see
[[cygl1_paired_thinning_density_experiment]] memory for the full technique).

**Three arms, Jan-Jun 2020 initially:**

| arm | source thinning | N kept (Jan-Jun 2020) | vs full stream |
|---|---|---|---|
| `intermediate_coh05` | DA-intermediate (min_sep=2.40°) + coherency filter | 3,940 | — |
| `dense075_coh05` | new tier, min_sep=0.75° + coherency filter | 19,802 | 5.0x intermediate |
| `dense_coh05` | full/unthinned stream + coherency filter only | 84,567 | — |

The `dense075` tier itself was newly calibrated this period, filling a gap between intermediate
and full density (a prior obs-spacing check had found intermediate's own median neighbor spacing
looser than assumed):

| min_sep_deg | N (pre-coherency) | vs intermediate | median NN dist | local-interaction median/p90 |
|---|---|---|---|---|
| 1.25 | 12,756 | 2.4x | 1.336° | 4/7 |
| 1.00 | 17,842 | 3.4x | 1.091° | 7/11 |
| **0.75 (chosen)** | **26,340** | **5.0x** | **0.824° (~1.8-2.3 grid cells)** | **12/18** |
| 0.60 | 34,284 | 6.5x | 0.681° | 16/25 |
| 0.50 | 42,044 | 8.0x | 0.578° | 21/32 |

**Config**: `errstd=2.75`, `xcorr=ycorr=0.625` for all three arms (identical otherwise —
one-intervention-per-experiment). One real gotcha worth recording: the first submission at
`xcorr=ycorr=1.25` (the value originally intended, matching the template default) crashed both
jobs immediately with `LDAS ERROR (3000) from check_compact` — `check_compact()` requires
`xcompact/ycompact >= 2 x` that *species'* own `xcorr/ycorr`, checked **per species**, and
`xcompact/ycompact` stayed at the project-standard 1.25° (shared with every other arm), which only
tolerates up to `xcorr/ycorr=0.625°` for this one species. Corrected to 0.625° (`xcompact/2`,
exactly at the limit) and resubmitted successfully.

### 6-month result (own-arm OmF vs OL cross-mask, N-weighted stdv)

| arm | Tb `(DA-OL)/OL` | SM `(DA-OL)/OL` | CygL1 own `(DA-OL)/OL` |
|---|---|---|---|
| `intermediate_coh05` | +0.06% | -0.05% | -8.06% |
| `dense075_coh05` | +0.30% | +0.35% | -8.71% |
| `dense_coh05` (full stream) | **+16.23%** | **+23.26%** | -3.30% |

The two thinned+filtered arms are clean (monitor species at noise floor); the full-stream
coherency-only arm is **not** — real, large cross-species harm survives the coherency filter when
there's no thinning and no hard gate underneath it.

## 4. Headline: dense075_coh05's 24-month result

`intermediate_coh05` and `dense075_coh05` were extended to the full Jan-Dec 2020 record
(2026-09-07), then `dense075_coh05` alone was extended a further 12 months (a from-BEG_DATE
restart to add `catch_progn_incr`/`inst3_1d_lndfcstana_Nt` output collections, completed cleanly
2026-09-08) to the full **Jan 2020 - Dec 2021, 24 months**:

| period | Tb `(DA-OL)/OL` | SM `(DA-OL)/OL` | CygL1 own `(DA-OL)/OL` |
|---|---|---|---|
| 6mo (Jan-Jun 2020) | +0.30% | +0.35% | -8.71% |
| 12mo (Jan-Dec 2020) | +0.08% | +0.32% | -11.80% |
| **24mo (Jan 2020-Dec 2021)** | **+0.45%** | **+0.63%** | **-12.88%** |

**The CygL1 skill signal (roughly -12 to -13% OmF stdv reduction vs OL) is durable across 4x the
record length**, not a short-window artifact, while monitor Tb/SM stay pinned at noise floor
(≤1.2% in either direction) throughout. `intermediate_coh05` was deliberately not extended past
12 months (stays at -13.18%/Tb -0.10%/SM -0.23%, essentially matching dense075_coh05).

### How this compares to every other intervention tried in this project

| approach | own CygL1 skill vs OL | monitor Tb/SM cross-species impact |
|---|---|---|
| thinning alone, no coherency filter (`intermediate`, 22mo) | -9.2% | neutral |
| **thinning + coherency filter (`dense075_coh05`, 24mo)** | **-12.9%** | **neutral** |
| hard gate, full density, no thinning (`dense_gated`, 22mo) | -0.4% (~flat) | neutral (fixed from active harm — see `runs/weekly_update_2026-08-27.md` §2) |
| coherency filter alone, full density, no gate (`dense_coh05`, 6mo) | -3.3% | **actively harmful** (Tb +16%, SM +23%) |

**Thinning — with or without the coherency filter layered on top — is the only approach so far
that both avoids cross-species harm and shows a meaningfully large CygL1-own skill improvement.**
The hard gate alone stops the harm but doesn't yet show real skill; coherency filtering alone
doesn't stop the harm at all. `dense075_coh05` combines both ideas and is currently this project's
best evidence that CYGNSS L1 assimilation can add real, durable skill without collateral damage.

Recommended slide bullets:

- The coh05 density-spectrum experiment confirms: obs density reduction (thinning), not per-obs
  quality screening alone, is what prevents cross-species harm at full observation density.
- `dense075_coh05` (thinning + coherency screening) shows a ~13% CygL1 own-fit skill improvement,
  stable from 6 to 24 months, with independent Tb/SM verification species staying neutral
  throughout.
- This is the strongest and most durable positive CYGNSS L1 result produced by this project to
  date.

## Not yet done / open items

- `dense_coh05` (full-stream, coherency-only) was not extended past 6 months — not planned, since
  it's already shown to underperform every thinned alternative.
- The obs-spacing anomaly flagged 2026-09-06 (some `intermediate`-tier obs pairs closer than the
  documented 2.40° minimum-separation admission floor) is still unexplained. It bears on the
  "clean isolation" characterization used throughout the broader thinning-density thread and
  should be revisited before leaning on it further.
- Two divergent copies of the project's shared `postproc_ObsFcstAna` toolkit exist in this
  environment (an older one in the `dnb34` GEOSldas checkout, a newer NC4-capable one sitting in a
  git-unregistered `hsaf_cdr_test` worktree) — not yet reconciled.

## Supporting files

Stats/data (gitignored, not in this repo — paths for reference):

- `output/postproc_paired_density/stats_output/spatial_stats_DA_paired_dense075_coh05_202001_202112.pkl`
  / `temporal_stats_DA_paired_dense075_coh05_20200101_20211231.nc4` (own-arm, 24mo)
- `output/postproc_paired_density/stats_output/spatial_stats_OL_paired_monitor_xmask_dense075_coh05_202001_202112.pkl`
  / matching `.nc4` (OL cross-masked to dense075_coh05's obs population, 24mo)
- Equivalent `_202001_202012` (12mo, intermediate_coh05 + dense075_coh05) and `_202001_202006`
  (6mo, all 3 arms) files for the shorter periods.

Scripts (committed this period, `geosldas-analysis` commits `9c840e6`/`56dceeb`/`54bbedf`, not yet
pushed to `origin/main`):

- `scripts/build_cygl1_dense075_thinning.py` — builds the dense075 tier for a given date range.
- `scripts/filter_cygl1_by_coherency.py` — coherency filter, any tier/date range.
- `scripts/compare_cygl1_coh05_omf.py` — the comparison table above, `--period {6mo,12mo,24mo,all}`.
- `scripts/postproc_drivers/run_cygl1_paired_density_coh05_24mo.py` /
  `run_cygl1_paired_OL_xmask_coh05_24mo.py` — the 24mo postproc driver scripts (relocated out of
  the broken `hsaf_cdr_test` worktree this period).
